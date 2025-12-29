%% EXAMPLE_USAGE.m
% Examples of how to use the command-line workflow

%% Example 1: Full workflow (automatic)
% Runs all 5 steps including MCMC (~30-40 minutes)
RUN_ALL

%% Example 2: Quick workflow (skip MCMC)
% Runs steps 1-3, 5 only (~15 seconds)
RUN_ALL('quick')

%% Example 3: Manual step-by-step
% Run each step separately for more control

% Step 1: Load data
STEP1_load_data

% Step 2: Fit constants
STEP2_fit_constants

% Step 3: Test simulation
STEP3_test_simulation

% Step 4: MCMC (optional, takes ~30 min)
% STEP4_run_mcmc  % Uncomment to run

% Step 5: Final simulation
STEP5_final_simulation

%% Example 4: Load previous results
% Load workspace from previous run
load workspace_final.mat

% Check parameters
params

% Re-simulate a specific batch
batch_idx = 1;
[t, y] = simulate_multi_fedbatch(batches{batch_idx}, constants, params);

% Plot
figure;
plot(t, y(:,1));  % Xv
hold on;
scatter(batches{batch_idx}.t, batches{batch_idx}.x);
xlabel('Time (h)');
ylabel('Xv (10^6 cells/mL)');
title('Manual simulation');

%% Example 5: Modify parameters and re-simulate
load workspace_step2.mat  % Load batches and constants

% Define your own parameters
my_params = struct();
my_params.Kglc = 0.2;   % Your value
my_params.Kgln = 0.08;  % Your value
my_params.KIamm = 12;   % Your value
my_params.mglc = 0.15;  % Your value
my_params.a1 = 0.015;   % Your value
my_params.a2 = 0.025;   % Your value
my_params.dgln = 0.012; % Your value
my_params.ramm = 0.015; % Your value

% Simulate
[t, y] = simulate_multi_fedbatch(batches{1}, constants, my_params);

% Plot
figure;
subplot(2,2,1); plot(t, y(:,1)); title('Xv');
subplot(2,2,2); plot(t, y(:,2)); title('GLC');
subplot(2,2,3); plot(t, y(:,3)); title('GLN');
subplot(2,2,4); plot(t, y(:,4)); title('NH4');

%% Example 6: Check MCMC convergence
% If you ran STEP4, check convergence
load workspace_step4.mat

% View final PSRF
param_names = {'Kglc', 'Kgln', 'KIamm', 'mglc', 'a1', 'a2', 'dgln', 'ramm'};
for i = 1:8
    fprintf('%s: PSRF = %.4f\n', param_names{i}, diagnostics.final_PSRF(i));
end

% Plot PSRF evolution
load mcmc_diagnostics/psrf_plot_data.mat
figure;
for i = 1:8
    subplot(3,3,i);
    plot(plot_data.iterations, plot_data.PSRF_history(:,i), 'LineWidth', 2);
    hold on;
    yline(1.1, '--r', 'Target');
    yline(1.0, '-k');
    title(plot_data.param_names{i});
    ylabel('PSRF');
    xlabel('Iteration');
    grid on;
end
sgtitle('PSRF Convergence');

%% Example 7: Compare different parameter sets
load workspace_step2.mat

% Test parameters
params_test = struct();
params_test.Kglc = 0.15;
params_test.Kgln = 0.05;
params_test.KIamm = 10;
params_test.mglc = 0.1;
params_test.a1 = 0.01;
params_test.a2 = 0.02;
params_test.dgln = 0.01;
params_test.ramm = 0.01;

% MCMC parameters (if available)
if exist('workspace_step4.mat', 'file')
    load workspace_step4.mat p_est
    params_mcmc = struct();
    params_mcmc.Kglc = p_est(1);
    params_mcmc.Kgln = p_est(2);
    params_mcmc.KIamm = p_est(3);
    params_mcmc.mglc = p_est(4);
    params_mcmc.a1 = p_est(5);
    params_mcmc.a2 = p_est(6);
    params_mcmc.dgln = p_est(7);
    params_mcmc.ramm = p_est(8);
    
    % Simulate both
    [t1, y1] = simulate_multi_fedbatch(batches{1}, constants, params_test);
    [t2, y2] = simulate_multi_fedbatch(batches{1}, constants, params_mcmc);
    
    % Compare
    figure('Position', [100 100 1200 800]);
    vars = {'Xv', 'GLC', 'GLN', 'NH4'};
    for i = 1:4
        subplot(2,2,i);
        plot(t1, y1(:,i), 'b-', 'LineWidth', 2, 'DisplayName', 'Test params');
        hold on;
        plot(t2, y2(:,i), 'r-', 'LineWidth', 2, 'DisplayName', 'MCMC params');
        scatter(batches{1}.t, batches{1}.(lower(vars{i})), 60, 'k', 'filled', ...
                'DisplayName', 'Data');
        xlabel('Time (h)');
        ylabel(vars{i});
        title(vars{i});
        legend('Location', 'best');
        grid on;
    end
    sgtitle('Test vs MCMC Parameters');
end
