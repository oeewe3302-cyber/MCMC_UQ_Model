function loglik = likelihood_fedbatch(p, batches, constants)
% LIKELIHOOD_FEDBATCH Calculate log-likelihood for MCMC
%
% Version: STATISTICALLY CORRECT - No artificial weights
%
% Syntax:
%   loglik = likelihood_fedbatch(p, batches, constants)
%
% Inputs:
%   p         - Parameter vector [8 x 1]:
%               [Kglc, Kgln, KIamm, mglc, a1, a2, dgln, ramm]
%   batches   - Cell array of batch data (concentrations in mM)
%   constants - Fitted constants structure
%
% Output:
%   loglik - Log-likelihood (higher is better)
%
% Method:
%   - Standard Gaussian likelihood for each measurement
%   - State-specific adaptive sigma (measurement noise model)
%   - No artificial weights (statistically correct)

    loglik = 0;
    
    % Convert parameter vector to structure
    params = vec2par(p);
    
    % Loop through all batches
    for i = 1:numel(batches)
        b = batches{i};
        
        try
            % Simulate with current parameters
            [t_sim, y_sim] = simulate_multi_fedbatch(b, constants, params);
            
            % Compare at measurement times
            nDays = length(b.t);
            for d = 1:nDays
                % Find closest simulation time
                [~, idx] = min(abs(t_sim - b.t(d)));
                
                % Predicted values (4 states only)
                Yp = y_sim(idx, :)';
                
                % Observed values (all in mM, Xv in 10^6 cells/mL)
                Yo = [b.x(d); b.glc(d); b.gln(d); b.nh4(d)];
                
                % Calculate log-likelihood for each state variable
                for j = 1:4
                    if isnan(Yo(j))
                        continue; 
                    end
                    
                    % *** State-specific measurement noise model ***
                    % Reflects typical experimental uncertainty for each variable
                    % Smaller sigma = stricter fitting for that variable
                    if j == 1  % Xv (10^6 cells/mL)
                        % Cell count: ~5% relative + 0.3 baseline
                        sigma_j = 0.05 * abs(Yo(j)) + 0.3;
                    elseif j == 2  % GLC (mM)
                        % Glucose assay: ~5% relative + 0.5 baseline (less critical)
                        sigma_j = 0.05 * abs(Yo(j)) + 0.5;
                    elseif j == 3  % GLN (mM)
                        % Glutamine: ~5% relative + 0.2 baseline (more critical)
                        sigma_j = 0.05 * abs(Yo(j)) + 0.2;
                    else  % j == 4, NH4 (mM)
                        % Ammonia: ~5% relative + 0.3 baseline
                        sigma_j = 0.05 * abs(Yo(j)) + 0.3;
                    end
                    
                    % Prevent division by zero
                    sigma_j = max(sigma_j, 1e-6);
                    
                    % Standardized residual
                    r = (Yo(j) - Yp(j)) / sigma_j;
                    
                    % Standard Gaussian log-likelihood
                    % log p(y|theta) = -0.5 * (y-y_pred)^2/sigma^2 - log(sigma*sqrt(2*pi))
                    loglik = loglik - 0.5 * r^2 - log(sigma_j * sqrt(2*pi));
                end
            end
        catch ME
            % If simulation fails, return very low likelihood
            % This causes MCMC to reject these parameters
            loglik = loglik - 1e10;
        end
    end
    
    % Prevent NaN or Inf in final result
    if isnan(loglik) || isinf(loglik)
        loglik = -1e10;
    end
end


function params = vec2par(p)
% VEC2PAR Convert parameter vector to structure
% 8 parameters (no KIlac)
    params = struct();
    params.Kglc = p(1);
    params.Kgln = p(2);
    params.KIamm = p(3);
    params.mglc = p(4);
    params.a1 = p(5);
    params.a2 = p(6);
    params.dgln = p(7);
    params.ramm = p(8);
end
