%% STEP1_load_data.m
% Load experimental data from Excel file
%
% Usage:
%   >> STEP1_load_data
%
% Output:
%   batches - Cell array of batch structures (saved to workspace)

clear all
clc

fprintf('\n========================================\n');
fprintf('STEP 1: LOADING EXCEL DATA\n');
fprintf('========================================\n\n');

% File path
filepath = 'data.xlsx';

% Check if file exists
if ~exist(filepath, 'file')
    error('❌ data.xlsx not found! Please place it in the current directory.');
end

% Load data
fprintf('📂 Loading: %s\n', filepath);
try
    batches = load_excel_data(filepath);
    fprintf('✅ Successfully loaded %d batches!\n\n', length(batches));
    
    % Display batch info
    fprintf('=== Batch Information ===\n');
    for i = 1:length(batches)
        b = batches{i};
        fprintf('  Batch %d: %s\n', i, b.name);
        fprintf('    Time points: %d (%.0f - %.0f h)\n', length(b.t), b.t(1), b.t(end));
        fprintf('    Initial Xv: %.2f (10^6 cells/mL)\n', b.x(1));
        fprintf('    Initial GLC: %.2f mM\n', b.glc(1));
        fprintf('    Initial GLN: %.2f mM\n', b.gln(1));
        fprintf('    Feed media: %d types\n\n', b.n_feeds);
    end
    
    % Save to workspace
    save('workspace_step1.mat', 'batches');
    fprintf('💾 Saved to: workspace_step1.mat\n');
    
    fprintf('\n✅ STEP 1 COMPLETE!\n');
    fprintf('Next: Run STEP2_fit_constants\n\n');
    
catch ME
    fprintf('❌ ERROR: %s\n', ME.message);
    fprintf('   Check if data.xlsx has correct format\n\n');
    rethrow(ME);
end
