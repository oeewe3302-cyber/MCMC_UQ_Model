function PSRF = calculate_psrf(chains)
% CALCULATE_PSRF Calculate Potential Scale Reduction Factor (Gelman-Rubin)
%
% Inputs:
%   chains - Cell array of chains {n_chains x 1}
%            Each chain is [n_samples x n_params] matrix
%
% Output:
%   PSRF - PSRF values for each parameter [1 x n_params]
%          Values < 1.1 indicate convergence

    n_chains = length(chains);
    n_samples = size(chains{1}, 1);
    n_params = size(chains{1}, 2);
    
    PSRF = zeros(1, n_params);
    
    for p = 1:n_params
        % Extract parameter values from all chains
        chain_means = zeros(n_chains, 1);
        chain_vars = zeros(n_chains, 1);
        
        for c = 1:n_chains
            values = chains{c}(:, p);
            chain_means(c) = mean(values);
            chain_vars(c) = var(values);
        end
        
        % Grand mean
        grand_mean = mean(chain_means);
        
        % Between-chain variance (B)
        B = n_samples * var(chain_means);
        
        % Within-chain variance (W)
        W = mean(chain_vars);
        
        % Pooled variance estimate
        V_hat = ((n_samples - 1) / n_samples) * W + (B / n_samples);
        
        % PSRF (R-hat)
        if W > 0
            PSRF(p) = sqrt(V_hat / W);
        else
            PSRF(p) = 1.0;  % If no variance, assume converged
        end
        
        % Prevent inf/nan
        if ~isfinite(PSRF(p))
            PSRF(p) = 999;  % Very large value indicating problem
        end
    end
end
