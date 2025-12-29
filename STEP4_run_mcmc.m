%% STEP4_run_mcmc.m
% Run MCMC parameter estimation with PSRF diagnostics
%
% Requires:
%   batches, constants (from STEP2)
%
% Usage:
%   >> STEP4_run_mcmc
%
% Output:
%   p_est - Estimated parameters (median)
%   posterior - All posterior samples
%   diagnostics - PSRF and convergence diagnostics
%
% Note: This takes ~20-40 minutes depending on your computer

clear; clc

fprintf('\n========================================\n');
fprintf('STEP 4: MCMC PARAMETER ESTIMATION\n');
fprintf('========================================\n\n');

% Load from STEP2
if ~exist('workspace_step2.mat', 'file')
    error('❌ workspace_step2.mat not found! Run STEP2_fit_constants first.');
end

load('workspace_step2.mat', 'batches', 'constants');
fprintf('✅ Loaded batches and constants\n\n');

% MCMC settings
fprintf('=== MCMC Settings ===\n');
fprintf('Chains: 3\n');
fprintf('Iterations: 4000 per chain\n');
fprintf('Burn-in: 2000\n');
fprintf('Total samples: 12000 (after burn-in)\n');
fprintf('Estimated time: 20-40 minutes\n\n');

% Confirm to run
fprintf('⚠️  This will take a while. Continue? (Ctrl+C to cancel)\n');
pause(3);

% Run MCMC
fprintf('\n⏳ Starting MCMC...\n');
fprintf('   Progress will be updated every 100 iterations\n\n');

try
    tic;
    % Pass empty status_handle for command-line mode
    [p_est, posterior, diagnostics] = run_mcmc_with_psrf(batches, constants, []);
    elapsed = toc;
    
    fprintf('\n✅ MCMC COMPLETE in %.1f minutes!\n', elapsed/60);
    
    % Display results
    param_names = {'Kglc', 'Kgln', 'KIamm', 'mglc', 'a1', 'a2', 'dgln', 'ramm'};
    
    fprintf('\n=== Estimated Parameters ===\n');
    for i = 1:8
        fprintf('%8s = %.4f\n', param_names{i}, p_est(i));
    end
    
    fprintf('\n=== Convergence (PSRF < 1.1 = Good) ===\n');
    for i = 1:8
        psrf_val = diagnostics.final_PSRF(i);
        status = '';
        if psrf_val < 1.05
            status = '✅ Excellent';
        elseif psrf_val < 1.1
            status = '✅ Good';
        elseif psrf_val < 1.2
            status = '⚠️  Marginal';
        else
            status = '❌ Poor';
        end
        fprintf('%8s: %.4f %s\n', param_names{i}, psrf_val, status);
    end
    
    % Save results
    save('workspace_step4.mat', 'batches', 'constants', ...
         'p_est', 'posterior', 'diagnostics');
    fprintf('\n💾 Saved to: workspace_step4.mat\n');
    
    % Save diagnostics
    if ~exist('mcmc_diagnostics', 'dir')
        mkdir('mcmc_diagnostics');
    end
    save_psrf_diagnostics(diagnostics, 'mcmc_diagnostics');
    fprintf('💾 Diagnostics saved to: mcmc_diagnostics/\n');
    
    fprintf('\n✅ STEP 4 COMPLETE!\n');
    fprintf('Next: Run STEP5_final_simulation\n\n');
    
catch ME
    fprintf('❌ ERROR: %s\n', ME.message);
    rethrow(ME);
end
