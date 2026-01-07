function [p_est, posterior, diagnostics] = run_mcmc_with_psrf(batches, constants, status_handle)
% RUN_MCMC_WITH_PSRF MCMC with convergence diagnostics (PSRF)
%
% Version: Command-line compatible (status_handle optional)
% WITHOUT LACTATE (8 parameters, KIlac removed)
%
% Runs 3 independent MCMC chains and calculates PSRF for convergence
%
% Syntax:
%   [p_est, posterior, diagnostics] = run_mcmc_with_psrf(batches, constants)
%   [p_est, posterior, diagnostics] = run_mcmc_with_psrf(batches, constants, status_handle)
%
% Inputs:
%   batches - Cell array of batch data structures
%   constants - Structure with kinetic constants
%   status_handle - (optional) GUI status text handle for updates
%
% Outputs:
%   p_est - Median parameter estimates from all chains combined
%   posterior - Combined posterior samples from all chains
%   diagnostics - Structure with PSRF values and chain data
%       .chains - Cell array of 3 chains
%       .PSRF - PSRF values at each checkpoint
%       .iterations - Iteration numbers for PSRF
%       .final_PSRF - Final PSRF values

    % Handle optional status_handle
    if nargin < 3
        status_handle = [];
    end

    % MCMC settings
    nIter = 12000;
    burn = 6000;
    n_chains = 3;
    

    % Order: [Kglc, Kgln, KIamm, mglc, a1, a2, dgln, ramm] - KIlac REMOVED
    % Bounds - Based on working previous code (more conservative)
    LB = [0.1,  0.01,   1,   0.001, 0.001, 0.01,  0.001, 0.001];
    UB = [1,    0.3,     80,   0.1,     0.01,     0.5,     0.05,   0.05];

    % Initial values for 3 chains - CENTERED in bounds to avoid violations
    % Based on MCMC diagnostics: moved problem parameters away from bounds
    % Order: [Kglc, Kgln, KIamm, mglc, a1, a2, dgln, ramm]
    p_init = zeros(n_chains, 8);
    p_init(1, :) = [0.2,  0.01,  5,   0.008,   0.002,  0.02,  0.005, 0.005];   % Conservative
    p_init(2, :) = [0.5,  0.02,  20,  0.01,   0.005,  0.1,  0.01,  0.01];   % Middle
    p_init(3, :) = [0.8,  0.1,  50,  0.08,   0.007,   0.3,  0.02,  0.02];   % Exploratory
    
    % Proposal step sizes - DRASTICALLY REDUCED to prevent bounds violations
    % Based on MCMC diagnostics: Kgln, a1, mglc were hitting bounds too often
    % Rule: Small enough that ±3*prop stays well within bounds
    % Order: [Kglc, Kgln, KIamm, mglc, a1, a2, dgln, ramm]
    prop = [0.02, 0.01, 2, 0.002, 0.0002, 0.002, 0.003, 0.003];
    
    % Storage
    chains = cell(n_chains, 1);
    acceptance_rates = zeros(n_chains, 1);
    true_acceptance_rates = zeros(n_chains, 1);
    bounds_rejections = zeros(n_chains, 1);
    bounds_per_param_all = zeros(n_chains, 8);
    
    % PSRF checkpoints
    psrf_checkpoints = 1000:500:nIter;
    n_checkpoints = length(psrf_checkpoints);
    PSRF_history = zeros(n_checkpoints, 8);
    
    fprintf('\n=== Running %d MCMC chains ===\n', n_chains);
    
    % Run chains
    for k = 1:n_chains
        fprintf('Chain %d/%d starting...\n', k, n_chains);
        fprintf('  [Progress will update every 50 iterations]\n');
        
        [chains{k}, acceptance_rates(k), true_acceptance_rates(k), bounds_rejections(k), bounds_per_param_all(k,:)] = ...
            run_single_chain(p_init(k, :), batches, constants, ...
                            LB, UB, prop, nIter, status_handle, k);
        
        fprintf('Chain %d/%d completed!\n', k, n_chains);
        fprintf('  Overall acceptance: %.1f%%\n', acceptance_rates(k) * 100);
        fprintf('  True acceptance: %.1f%% (excluding bounds)\n', true_acceptance_rates(k) * 100);
        fprintf('  Bounds rejections: %d (%.1f%%)\n', bounds_rejections(k), ...
                100*bounds_rejections(k)/nIter);
        
        if bounds_rejections(k) > 0
            fprintf('  Bounds violations by parameter:\n');
            param_names = {'Kglc', 'Kgln', 'KIamm', 'mglc', 'a1', 'a2', 'dgln', 'ramm'};
            for p = 1:8
                if bounds_per_param_all(k,p) > 0
                    pct = 100 * bounds_per_param_all(k,p) / bounds_rejections(k);
                    fprintf('    %s: %d (%.1f%%)\n', param_names{p}, bounds_per_param_all(k,p), pct);
                end
            end
        end
        fprintf('\n');
    end
    
    % Calculate PSRF at checkpoints
    fprintf('\n=== Calculating PSRF convergence diagnostics ===\n');
    checkpoint_idx = 1;
    
    for iter = psrf_checkpoints
        fprintf('PSRF at iteration %d:\n', iter);
        
        samples = cell(n_chains, 1);
        for k = 1:n_chains
            samples{k} = chains{k}(burn+1:iter, :);
        end
        
        % Calculate PSRF for all parameters at once
        psrf_values = calculate_psrf(samples);
        
        PSRF_history(checkpoint_idx, :) = psrf_values;
        checkpoint_idx = checkpoint_idx + 1;
        
        param_names = {'Kglc', 'Kgln', 'KIamm', 'mglc', 'a1', 'a2', 'dgln', 'ramm'};
        for i = 1:8
            fprintf('  %8s: %.4f', param_names{i}, psrf_values(i));
            if psrf_values(i) < 1.05
                fprintf(' ✓\n');
            elseif psrf_values(i) < 1.1
                fprintf(' ✓\n');
            elseif psrf_values(i) < 1.2
                fprintf(' ⚠\n');
            else
                fprintf(' ✗\n');
            end
        end
        fprintf('\n');
    end
    
    % Combine chains
    fprintf('=== Combining chains ===\n');
    posterior = [];
    for k = 1:n_chains
        posterior = [posterior; chains{k}(burn+1:end, :)];
    end
    fprintf('Combined posterior: %d samples\n', size(posterior, 1));
    
    % Parameter estimates
    p_est = median(posterior, 1);
    
    % Final PSRF
    samples_final = cell(n_chains, 1);
    for k = 1:n_chains
        samples_final{k} = chains{k}(burn+1:end, :);
    end
    
    % Calculate final PSRF for all parameters at once
    final_PSRF = calculate_psrf(samples_final);
    
    % Store diagnostics
    diagnostics = struct();
    diagnostics.chains = chains;
    diagnostics.PSRF_history = PSRF_history;  % Fixed: was "PSRF"
    diagnostics.iterations = psrf_checkpoints;
    diagnostics.final_PSRF = final_PSRF;
    diagnostics.acceptance_rates = acceptance_rates;
    diagnostics.true_acceptance_rates = true_acceptance_rates;
    diagnostics.bounds_rejections = bounds_rejections;
    diagnostics.burn = burn;  % Added
    diagnostics.nIter = nIter;  % Added
    diagnostics.n_chains = n_chains;  % Added for completeness
    
    param_names = {'Kglc', 'Kgln', 'KIamm', 'mglc', 'a1', 'a2', 'dgln', 'ramm'};
    
    fprintf('\n=== Acceptance Rates Summary ===\n');
    for k = 1:n_chains
        fprintf('Chain %d: %.1f%% overall, %.1f%% true acceptance\n', ...
                k, 100*acceptance_rates(k), 100*true_acceptance_rates(k));
    end
    
    total_bounds = sum(bounds_per_param_all, 1);
    if any(total_bounds > 0)
        fprintf('\n=== Parameters frequently hitting bounds ===\n');
        for p = 1:8
            if total_bounds(p) > 100
                fprintf('  %s: %d times | LB=%.4f, UB=%.4f\n', param_names{p}, total_bounds(p), LB(p), UB(p));
            end
        end
        fprintf('\n');
    end
    
    fprintf('\nFinal PSRF values (target < 1.1):\n');
    for i = 1:8
        psrf_val = diagnostics.final_PSRF(i);
        status = '';
        if psrf_val < 1.05
            status = '✓ Excellent';
        elseif psrf_val < 1.1
            status = '✓ Good';
        elseif psrf_val < 1.2
            status = '⚠ Marginal';
        else
            status = '✗ Poor';
        end
        fprintf('  %8s: %.4f %s\n', param_names{i}, psrf_val, status);
    end
    
    fprintf('\nEstimates:\n');
    for i = 1:8
        fprintf('  %8s: %.4f\n', param_names{i}, p_est(i));
    end
    fprintf('\n');
end


function [chain, acceptance_rate, true_acceptance_rate, rejected_bounds, bounds_per_param] = run_single_chain(p_init, batches, constants, ...
                                                      LB, UB, prop, nIter, status_handle, chain_num)
% RUN_SINGLE_CHAIN Run a single MCMC chain with progress reporting
    
    p_curr = p_init;
    chain = zeros(nIter, 8);
    acc = 0;
    rejected_bounds = 0;
    bounds_per_param = zeros(1, 8);
    
    fprintf('  Chain %d: Calculating initial likelihood...\n', chain_num);
    log_curr = likelihood_fedbatch(p_curr, batches, constants);
    fprintf('  Chain %d: Starting iterations...\n', chain_num);
    
    for i = 1:nIter
        p_prop = p_curr + prop .* randn(1, 8);
        
        below_LB = p_prop < LB;
        above_UB = p_prop > UB;
        
        if any(below_LB) || any(above_UB)
            chain(i, :) = p_curr;
            rejected_bounds = rejected_bounds + 1;
            bounds_per_param = bounds_per_param + (below_LB | above_UB);
            continue;
        end
        
        log_prop = likelihood_fedbatch(p_prop, batches, constants);
        
        alpha = exp(log_prop - log_curr);
        if rand < alpha
            p_curr = p_prop;
            log_curr = log_prop;
            acc = acc + 1;
        end
        
        chain(i, :) = p_curr;
        
        % Update progress - EVERY 50 iterations (like previous code)
        if mod(i, 50) == 0
            n_evaluated = i - rejected_bounds;
            true_acc_rate = 0;
            if n_evaluated > 0
                true_acc_rate = acc / n_evaluated;
            end
            
            progress_pct = 100 * i / nIter;
            
            fprintf('  Chain %d: %d/%d (%.0f%%) | Accept: %.1f%% | True: %.1f%%\n', ...
                    chain_num, i, nIter, progress_pct, 100*acc/i, 100*true_acc_rate);
            
            if ~isempty(status_handle) && ishandle(status_handle)
                set(status_handle, 'String', sprintf('⏳ Chain %d: %d/%d (%.0f%%) | Accept: %.1f%% (True: %.1f%%)', ...
                    chain_num, i, nIter, progress_pct, 100*acc/i, 100*true_acc_rate));
                drawnow;
            end
        end
    end
    
    acceptance_rate = acc / nIter;
    n_evaluated = nIter - rejected_bounds;
    true_acceptance_rate = acc / n_evaluated;
end
