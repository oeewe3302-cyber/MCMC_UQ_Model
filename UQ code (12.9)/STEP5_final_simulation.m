%% STEP5_final_simulation.m
% Final simulation with estimated parameters and export results
%
% Requires:
%   batches, constants, p_est (from STEP4)
%   OR batches, constants, params (from STEP3)
%
% Usage:
%   >> STEP5_final_simulation
%
% Output:
%   Plots for all batches
%   Excel file with simulation results

clear; clc

fprintf('\n========================================\n');
fprintf('STEP 5: FINAL SIMULATION & EXPORT\n');
fprintf('========================================\n\n');

% Load workspace
if exist('workspace_step4.mat', 'file')
    fprintf('✅ Loading MCMC results from STEP4...\n');
    load('workspace_step4.mat', 'batches', 'constants', 'p_est');
    use_mcmc = true;
else
    if exist('workspace_step3.mat', 'file')
        fprintf('✅ Loading test parameters from STEP3...\n');
        load('workspace_step3.mat', 'batches', 'constants', 'params');
        use_mcmc = false;
    else
        error('❌ No workspace found! Run STEP2-4 first.');
    end
end

% Convert parameters
if use_mcmc
    fprintf('   Using MCMC estimated parameters\n\n');
    params = struct();
    params.Kglc = p_est(1);
    params.Kgln = p_est(2);
    params.KIamm = p_est(3);
    params.mglc = p_est(4);
    params.a1 = p_est(5);
    params.a2 = p_est(6);
    params.dgln = p_est(7);
    params.ramm = p_est(8);
else
    fprintf('   Using test parameters (no MCMC)\n\n');
end

% Display parameters
fprintf('=== Final Parameters ===\n');
fprintf('Kglc   = %.4f mM\n', params.Kglc);
fprintf('Kgln   = %.4f mM\n', params.Kgln);
fprintf('KIamm  = %.4f mM\n', params.KIamm);
fprintf('mglc   = %.4f\n', params.mglc);
fprintf('a1     = %.4f\n', params.a1);
fprintf('a2     = %.4f\n', params.a2);
fprintf('dgln   = %.5f\n', params.dgln);
fprintf('ramm   = %.4f\n\n', params.ramm);

% Simulate all batches
n_batches = length(batches);
fprintf('⏳ Simulating %d batches...\n', n_batches);

results = cell(n_batches, 1);

for i = 1:n_batches
    fprintf('  [%d/%d] %s... ', i, n_batches, batches{i}.name);
    tic;
    
    try
        [t, y] = simulate_multi_fedbatch(batches{i}, constants, params);
        results{i}.t = t;
        results{i}.y = y;
        results{i}.name = batches{i}.name;
        results{i}.success = true;
        
        elapsed = toc;
        fprintf('✅ %.2f sec\n', elapsed);
    catch ME
        results{i}.success = false;
        results{i}.error = ME.message;
        fprintf('❌ FAILED\n');
    end
end

fprintf('\n✅ All simulations complete!\n\n');

% Plot results
fprintf('📊 Creating plots...\n');

vars = {'x', 'glc', 'gln', 'nh4'};
titles = {'Viable Cells (Xv)', 'Glucose (GLC)', 'Glutamine (GLN)', 'Ammonia (NH4)'};
ylabels = {'10^6 cells/mL', 'mM', 'mM', 'mM'};
colors = lines(n_batches);

for v = 1:4
    figure('Position', [100 + (v-1)*50, 100 + (v-1)*50, 1000, 600], ...
           'Name', titles{v});
    
    for i = 1:n_batches
        if ~results{i}.success
            continue;
        end
        
        % Plot model
        plot(results{i}.t, results{i}.y(:,v), '-', ...
             'LineWidth', 2, 'Color', colors(i,:), ...
             'DisplayName', sprintf('%s (model)', results{i}.name));
        hold on;
        
        % Plot data
        b = batches{i};
        scatter(b.t, b.(vars{v}), 60, colors(i,:), 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
                'DisplayName', sprintf('%s (data)', results{i}.name));
    end
    
    xlabel('Time (h)', 'FontSize', 12);
    ylabel(ylabels{v}, 'FontSize', 12);
    title(titles{v}, 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 9);
    grid on;
    box on;
    
    % Save figure
    filename = sprintf('final_simulation_%s.png', vars{v});
    saveas(gcf, filename);
    fprintf('   💾 Saved: %s\n', filename);
end

% Export to Excel
fprintf('\n📊 Exporting results to Excel...\n');

for i = 1:n_batches
    if ~results{i}.success
        continue;
    end
    
    % Create table
    T = table(results{i}.t, ...
              results{i}.y(:,1), ...
              results{i}.y(:,2), ...
              results{i}.y(:,3), ...
              results{i}.y(:,4), ...
              'VariableNames', {'Time_h', 'Xv_1e6cellsmL', 'GLC_mM', 'GLN_mM', 'NH4_mM'});
    
    % Write to Excel
    filename = sprintf('simulation_%s.xlsx', results{i}.name);
    writetable(T, filename);
    fprintf('   💾 Saved: %s\n', filename);
end

% Save final workspace
save('workspace_final.mat', 'batches', 'constants', 'params', 'results');
fprintf('\n💾 Final workspace saved: workspace_final.mat\n');

fprintf('\n========================================\n');
fprintf('✅✅✅ ALL STEPS COMPLETE! ✅✅✅\n');
fprintf('========================================\n\n');

fprintf('Generated files:\n');
fprintf('  - Figures: final_simulation_*.png\n');
fprintf('  - Excel: simulation_*.xlsx\n');
fprintf('  - Workspace: workspace_final.mat\n');
if use_mcmc
    fprintf('  - MCMC diagnostics: mcmc_diagnostics/\n');
end
fprintf('\n');
