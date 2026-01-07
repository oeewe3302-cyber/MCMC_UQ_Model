function [t_all, y_all] = simulate_multi_fedbatch(batch_data, constants, params)
% SIMULATE_MULTI_FEDBATCH Simulate fed-batch with multiple feeding system
%
% Version: CLEAN & FAST - Optimized for quick test simulation
%
% Syntax:
%   [t_all, y_all] = simulate_multi_fedbatch(batch_data, constants, params)
%
% Inputs:
%   batch_data - Batch structure with experimental data
%   constants  - Structure with kinetic constants
%   params     - Structure with parameters to estimate
%
% Outputs:
%   t_all - Time vector (h)
%   y_all - State matrix [time x 4] = [X, G, Q, A]

    t_all = [];
    y_all = [];
    
    % Initial state (Day 0, prefeed) - 4 states only
    x = [batch_data.x(1); 
         batch_data.glc(1); 
         batch_data.gln(1); 
         batch_data.nh4(1)];
    
    % *** EXTREMELY FAST solver for MCMC ***
    % Very relaxed settings - bad parameters will be rejected anyway
    RelTol = 1e-3;      % 1% tolerance (very relaxed!)
    AbsTol = 1e-4;      % Very relaxed
    MaxStep = 1;       % Very large steps
    InitialStep = 0.001;  % Start with very big steps  
    
    % Loop through each day
    for d = 0:(length(batch_data.t)-2)
        t_start = batch_data.t(d+1);
        t_end = batch_data.t(d+2);
        
        % Current volume
        V_current = batch_data.vol(d+1);
        
        % ============================================
        % Apply multiple feeding at t_start
        % ============================================
        Vf_glc = batch_data.feed_glc_vol(d+1);
        Vf_feeds = batch_data.feed_vols(d+1, :);
        Vf_total = Vf_glc + sum(Vf_feeds);
        
        if Vf_total > 0
            V_after_feed = V_current + Vf_total;
            
            glc_from_pure = Vf_glc * batch_data.glc_stock;
            glc_from_feeds = 0;
            gln_from_feeds = 0;
            
            for f = 1:batch_data.n_feeds
                glc_from_feeds = glc_from_feeds + Vf_feeds(f) * batch_data.feed_glc(f);
                gln_from_feeds = gln_from_feeds + Vf_feeds(f) * batch_data.feed_gln(f);
            end
            
            glc_total_added = glc_from_pure + glc_from_feeds;
            gln_total_added = gln_from_feeds;
            
            % Update state (feeding effect)
            x(1) = x(1) * V_current / V_after_feed;                    % X
            x(2) = (x(2)*V_current + glc_total_added) / V_after_feed;  % G
            x(3) = (x(3)*V_current + gln_total_added) / V_after_feed;  % Q
            x(4) = x(4) * V_current / V_after_feed;                    % A
        end
        
        % ============================================
        % ODE Integration with MinStep protection (prevents hanging)
        % ============================================
        
        % Turn warnings into errors temporarily (so we can catch ODE failures)
        warning_state = warning('query', 'MATLAB:ode15s:IntegrationTolNotMet');
        warning('error', 'MATLAB:ode15s:IntegrationTolNotMet');
        
        opts = odeset('RelTol', RelTol, ...
                      'AbsTol', AbsTol, ...
                      'MaxStep', MaxStep, ...
                      'InitialStep', InitialStep, ...
                      'MinStep', 0.0001, ...
                      'NonNegative', [1 2 3 4]);

        % Integration with error handling
        try
            [t_temp, y_temp] = ode15s(@(t,y) ode_fedbatch(t, y, constants, params), ...
                                       [t_start, t_end], x, opts);
            
            % Restore warning state
            warning(warning_state.state, 'MATLAB:ode15s:IntegrationTolNotMet');
            
        catch ME
            % Restore warning state even on error
            warning(warning_state.state, 'MATLAB:ode15s:IntegrationTolNotMet');
            
            % Re-throw to be caught by likelihood function
            rethrow(ME);
        end
        
        % Store results
        t_all = [t_all; t_temp];
        y_all = [y_all; y_temp];
        
        % Update to next prefeed state
        x = y_temp(end, :)';
    end
end
