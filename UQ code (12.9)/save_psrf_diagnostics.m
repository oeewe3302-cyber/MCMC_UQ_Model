function save_psrf_diagnostics(diagnostics, output_dir)
% SAVE_PSRF_DIAGNOSTICS Save PSRF convergence diagnostics to files
%
% Creates:
%   1. psrf_values.csv - PSRF values at each checkpoint
%   2. psrf_summary.txt - Convergence summary
%   3. psrf_plot_data.mat - Data for plotting
%
% Inputs:
%   diagnostics - Output from run_mcmc_with_psrf
%   output_dir - Directory to save files (default: current)

    if nargin < 2
        output_dir = '.';
    end
    
    % Create output directory if needed
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    
    param_names = {'Kglc', 'Kgln', 'KIamm', 'mglc', 'a1', 'a2', 'dgln', 'ramm'};  % KIlac removed
    
    %% 1. Save PSRF values to CSV
    csv_file = fullfile(output_dir, 'psrf_values.csv');
    
    % Header
    fid = fopen(csv_file, 'w');
    fprintf(fid, 'Iteration');
    for i = 1:8  % Changed from 9 to 8
        fprintf(fid, ',%s', param_names{i});
    end
    fprintf(fid, '\n');
    
    % Data
    for row = 1:size(diagnostics.PSRF_history, 1)
        fprintf(fid, '%d', diagnostics.iterations(row));
        for col = 1:8  % Changed from 9 to 8
            fprintf(fid, ',%.6f', diagnostics.PSRF_history(row, col));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
    
    fprintf('✓ Saved PSRF values to: %s\n', csv_file);
    
    %% 2. Save summary text
    txt_file = fullfile(output_dir, 'psrf_summary.txt');
    fid = fopen(txt_file, 'w');
    
    fprintf(fid, '=== MCMC Convergence Diagnostics ===\n\n');
    fprintf(fid, 'Number of chains: %d\n', length(diagnostics.chains));
    fprintf(fid, 'Burn-in: %d\n', diagnostics.burn);
    fprintf(fid, 'Total samples per chain: %d\n', size(diagnostics.chains{1}, 1));
    fprintf(fid, '\n');
    
    fprintf(fid, 'Acceptance rates:\n');
    for k = 1:length(diagnostics.acceptance_rates)
        fprintf(fid, '  Chain %d:\n', k);
        fprintf(fid, '    Overall: %.2f%% (includes bounds rejections)\n', ...
                diagnostics.acceptance_rates(k) * 100);
        fprintf(fid, '    True: %.2f%% (excludes bounds rejections)\n', ...
                diagnostics.true_acceptance_rates(k) * 100);
        fprintf(fid, '    Bounds rejected: %d (%.2f%%)\n', ...
                diagnostics.bounds_rejections(k), ...
                100*diagnostics.bounds_rejections(k)/(diagnostics.burn*2));
    end
    fprintf(fid, '\n');
    
    fprintf(fid, 'Final PSRF values (iteration %d):\n', diagnostics.iterations(end));
    fprintf(fid, 'Target: < 1.1 for convergence, < 1.05 for excellent convergence\n\n');
    
    for i = 1:8  % Changed from 9 to 8
        psrf_val = diagnostics.final_PSRF(i);
        status = '';
        if psrf_val < 1.05
            status = 'Excellent';
        elseif psrf_val < 1.1
            status = 'Good';
        elseif psrf_val < 1.2
            status = 'Marginal';
        else
            status = 'Poor - needs more iterations';
        end
        fprintf(fid, '  %8s: %.4f - %s\n', param_names{i}, psrf_val, status);
    end
    
    fprintf(fid, '\n');
    fprintf(fid, 'Convergence assessment:\n');
    if all(diagnostics.final_PSRF < 1.1)
        fprintf(fid, '✓ All parameters converged (PSRF < 1.1)\n');
    elseif all(diagnostics.final_PSRF < 1.2)
        fprintf(fid, '⚠ Marginal convergence (PSRF < 1.2)\n');
        fprintf(fid, '  Recommend more iterations or different initial values\n');
    else
        fprintf(fid, '✗ Poor convergence detected\n');
        fprintf(fid, '  Action needed: Increase iterations or check model\n');
    end
    
    fclose(fid);
    fprintf('✓ Saved summary to: %s\n', txt_file);
    
    %% 3. Save MATLAB data for plotting
    mat_file = fullfile(output_dir, 'psrf_plot_data.mat');
    
    plot_data.iterations = diagnostics.iterations;
    plot_data.PSRF_history = diagnostics.PSRF_history;
    plot_data.param_names = param_names;
    plot_data.final_PSRF = diagnostics.final_PSRF;
    plot_data.chains = diagnostics.chains;
    
    save(mat_file, 'plot_data');
    fprintf('✓ Saved plot data to: %s\n', mat_file);
    
    %% 4. Instructions for plotting
    fprintf('\n=== How to plot PSRF ===\n');
    fprintf('Run the following code:\n\n');
    fprintf('  load(''%s'');\n', mat_file);
    fprintf('  figure;\n');
    fprintf('  for i = 1:8\n');  % Changed from 9 to 8
    fprintf('      subplot(3, 3, i);\n');
    fprintf('      plot(plot_data.iterations, plot_data.PSRF_history(:, i), ''LineWidth'', 2);\n');
    fprintf('      hold on;\n');
    fprintf('      yline(1.1, ''--r'', ''Target'');\n');
    fprintf('      yline(1.0, ''-k'');\n');
    fprintf('      xlabel(''Iteration'');\n');
    fprintf('      ylabel(''PSRF'');\n');
    fprintf('      title(plot_data.param_names{i});\n');
    fprintf('      grid on;\n');
    fprintf('  end\n');
    fprintf('  sgtitle(''PSRF Convergence Diagnostics'');\n\n');
    
    fprintf('Or use: plot_psrf_from_file(''%s'');\n\n', mat_file);
end
