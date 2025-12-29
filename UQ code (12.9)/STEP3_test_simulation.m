%% STEP3_test_simulation.m
% Test simulation with initial parameter guesses
%
% Requires:
%   batches, constants (from STEP2)
%
% Usage:
%   >> STEP3_test_simulation
%
% Output:
%   Plots comparing model vs data

clear; clc

fprintf('\n========================================\n');
fprintf('STEP 3: TEST SIMULATION\n');
fprintf('========================================\n\n');

% Load from STEP2
if ~exist('workspace_step2.mat', 'file')
    error('❌ workspace_step2.mat not found! Run STEP2_fit_constants first.');
end

load('workspace_step2.mat', 'batches', 'constants');
fprintf('✅ Loaded batches and constants\n\n');

% ============================================================
% USER PARAMETERS - Initial guesses for simulation
% ============================================================
fprintf('=== Initial Parameter Guesses ===\n');
params = struct();
params.Kglc = 0.2;    % Glucose half-saturation (mM) %0.1
params.Kgln = 0.01;    % Glutamine half-saturation (mM) %0.01
params.KIamm = 10;     % Ammonia inhibition constant (mM) %3
params.mglc = 0.02;     % Glucose maintenance (mmol/(10^6 cells·h)) %0.01
params.a1 = 0.001;      % Glutamine maintenance numerator
params.a2 = 0.2;      % Glutamine maintenance denominator (mM)
params.dgln = 0.01;    % Glutamine degradation rate (h⁻¹)
params.ramm = 0.1;    % Ammonia removal rate (mmol/(10^6 cells·h)) %0.01

fprintf('Kglc   = %.3f mM\n', params.Kglc);
fprintf('Kgln   = %.3f mM\n', params.Kgln);
fprintf('KIamm  = %.1f mM\n', params.KIamm);
fprintf('mglc   = %.3f\n', params.mglc);
fprintf('a1     = %.3f\n', params.a1);
fprintf('a2     = %.3f\n', params.a2);
fprintf('dgln   = %.4f\n', params.dgln);
fprintf('ramm   = %.3f\n\n', params.ramm);

% Select batch to simulate (default: first batch)
batch_idx = 1;
fprintf('Simulating: %s\n', batches{batch_idx}.name);

% Run simulation
fprintf('\n⏳ Running simulation...\n');
tic;
try
    [t, y] = simulate_multi_fedbatch(batches{batch_idx}, constants, params);
    elapsed = toc;
    fprintf('✅ Simulation complete in %.2f sec\n', elapsed);
    fprintf('   Generated %d time points\n\n', length(t));
    
    % Check for NaN
    if any(isnan(y(:)))
        warning('⚠️  Simulation contains NaN values!');
    end
    
    % Plot results
    fprintf('📊 Creating plots...\n');
    
    figure('Position', [100 100 1200 800], 'Name', 'Test Simulation Results');
    
    vars = {'x', 'glc', 'gln', 'nh4'};
    titles = {'Viable Cells (Xv)', 'Glucose (GLC)', 'Glutamine (GLN)', 'Ammonia (NH4)'};
    ylabels = {'10^6 cells/mL', 'mM', 'mM', 'mM'};
    
    for i = 1:4
        subplot(2, 2, i);
        
        % Plot model
        plot(t, y(:,i), 'b-', 'LineWidth', 2.5, 'DisplayName', 'Model');
        hold on;
        
        % Plot data
        b = batches{batch_idx};
        scatter(b.t, b.(vars{i}), 80, 'r', 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
                'DisplayName', 'Data');
        
        xlabel('Time (h)', 'FontSize', 11);
        ylabel(ylabels{i}, 'FontSize', 11);
        title(titles{i}, 'FontSize', 12, 'FontWeight', 'bold');
        legend('Location', 'best', 'FontSize', 10);
        grid on;
        box on;
    end
    
    sgtitle(sprintf('Test Simulation: %s', batches{batch_idx}.name), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, 'test_simulation.png');
    fprintf('💾 Figure saved: test_simulation.png\n');
    
    % Save results
    save('workspace_step3.mat', 'batches', 'constants', 'params', 't', 'y');
    fprintf('💾 Saved to: workspace_step3.mat\n');
    
    fprintf('\n✅ STEP 3 COMPLETE!\n');
    fprintf('Next: Run STEP4_run_mcmc (optional, takes ~30 min)\n');
    fprintf('      Or skip to STEP5_final_simulation with these params\n\n');
    
catch ME
    fprintf('❌ ERROR: %s\n', ME.message);
    rethrow(ME);
end
