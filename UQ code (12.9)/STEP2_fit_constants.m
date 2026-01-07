%% STEP2_fit_constants.m
% Fit kinetic constants from experimental data
%
% Requires:
%   batches (from STEP1)
%
% Usage:
%   >> STEP2_fit_constants
%
% Output:
%   constants - Fitted constants structure
%   MAT - Matrix of fitted values for each batch

clear; clc

fprintf('\n========================================\n');
fprintf('STEP 2: FITTING CONSTANTS\n');
fprintf('========================================\n\n');

% Load batches from STEP1
if ~exist('workspace_step1.mat', 'file')
    error('❌ workspace_step1.mat not found! Run STEP1_load_data first.');
end

load('workspace_step1.mat', 'batches');
fprintf('✅ Loaded %d batches\n\n', length(batches));

% Fit constants
fprintf('⏳ Fitting constants (90-200h exponential phase)...\n');
fprintf('   This may take 5-10 seconds...\n\n');

try
    [constants, MAT] = fit_constants_fedbatch(batches);
    
    % Display results
    fprintf('\n=== Fitted Constants (4-ODE, NO LACTATE) ===\n');
    fprintf('μ_max     = %.4f h⁻¹\n', constants.mu_max);
    fprintf('Y_xv/glc  = %.2e (10^6 cells/mL)/mM\n', constants.Yxv_glc);
    fprintf('Y_xv/gln  = %.2e (10^6 cells/mL)/mM\n', constants.Yxv_gln);
    fprintf('Y_amm/gln = %.3f mM/mM\n\n', constants.Yamm_gln);
    
    % Check for NaN
    if any(isnan([constants.mu_max, constants.Yxv_glc, constants.Yxv_gln, constants.Yamm_gln]))
        warning('⚠️  Some constants are NaN! Check your data.');
    end
    
    % Save to workspace
    save('workspace_step2.mat', 'batches', 'constants', 'MAT');
    fprintf('💾 Saved to: workspace_step2.mat\n');
    
    fprintf('\n✅ STEP 2 COMPLETE!\n');
    fprintf('Next: Run STEP3_test_simulation\n\n');
    
catch ME
    fprintf('❌ ERROR: %s\n', ME.message);
    rethrow(ME);
end
