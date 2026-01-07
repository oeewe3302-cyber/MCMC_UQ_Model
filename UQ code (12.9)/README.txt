========================================
CHO CELL FED-BATCH MCMC ANALYSIS
Command-Line Version (4-ODE, NO LACTATE)
========================================

📋 OVERVIEW
-----------
This package provides a complete workflow for CHO cell fed-batch 
bioprocess modeling and parameter estimation using MCMC.

✅ Features:
- 4-ODE system (Xv, GLC, GLN, NH4)
- Multiple feeding system support
- MCMC with PSRF convergence diagnostics
- Fast solver (2-3 sec per simulation)
- Command-line interface (no GUI needed!)


📁 FILES
--------
Main workflow scripts:
  STEP1_load_data.m         - Load Excel data
  STEP2_fit_constants.m     - Fit kinetic constants
  STEP3_test_simulation.m   - Test simulation with initial guess
  STEP4_run_mcmc.m          - MCMC parameter estimation (~30 min)
  STEP5_final_simulation.m  - Final simulation & export
  RUN_ALL.m                 - Run all steps automatically

Core functions:
  load_excel_data.m         - Excel data loader
  fit_constants_fedbatch.m  - Constants fitting (mu_max, yields)
  simulate_multi_fedbatch.m - Fed-batch simulator
  ode_fedbatch.m            - 4-ODE system
  likelihood_fedbatch.m     - Log-likelihood for MCMC
  run_mcmc_with_psrf.m      - MCMC with 3 chains
  calculate_psrf.m          - PSRF calculation
  save_psrf_diagnostics.m   - Save diagnostics


🚀 QUICK START
--------------
1. Place data.xlsx in the same folder as these scripts
2. Open MATLAB in this folder
3. Run:
   >> RUN_ALL           % Full workflow (includes MCMC)
   or
   >> RUN_ALL('quick')  % Skip MCMC (faster)


📝 STEP-BY-STEP USAGE
---------------------

Option 1: Run all at once
--------------------------
>> RUN_ALL           % Takes ~30-40 minutes (includes MCMC)
>> RUN_ALL('quick')  % Takes ~1 minute (skips MCMC)


Option 2: Run each step manually
----------------------------------
>> STEP1_load_data          % 1-2 seconds
   → Loads data.xlsx
   → Saves: workspace_step1.mat

>> STEP2_fit_constants      % 3-5 seconds
   → Fits mu_max, yields from 90-200h data
   → Saves: workspace_step2.mat

>> STEP3_test_simulation    % 2-3 seconds
   → Runs simulation with initial parameter guesses
   → Creates: test_simulation.png
   → Saves: workspace_step3.mat

>> STEP4_run_mcmc           % 20-40 minutes (OPTIONAL)
   → Estimates 8 parameters using MCMC
   → Creates: mcmc_diagnostics/
   → Saves: workspace_step4.mat
   
   ⚠️ This step is OPTIONAL!
      You can skip it and use test parameters

>> STEP5_final_simulation   % 5-10 seconds
   → Simulates all batches with final parameters
   → Creates: final_simulation_*.png
   → Creates: simulation_*.xlsx
   → Saves: workspace_final.mat


📊 EXPECTED OUTPUT
------------------
After running all steps:

Figures:
  test_simulation.png           - Test simulation results
  final_simulation_x.png        - Xv for all batches
  final_simulation_glc.png      - GLC for all batches
  final_simulation_gln.png      - GLN for all batches
  final_simulation_nh4.png      - NH4 for all batches

Excel files:
  simulation_Batch1.xlsx        - Simulation data for Batch1
  simulation_Batch2.xlsx        - Simulation data for Batch2
  ...

Workspaces:
  workspace_step1.mat           - After data loading
  workspace_step2.mat           - After constants fitting
  workspace_step3.mat           - After test simulation
  workspace_step4.mat           - After MCMC (if ran)
  workspace_final.mat           - Final results

MCMC diagnostics (if STEP4 ran):
  mcmc_diagnostics/psrf_values.csv
  mcmc_diagnostics/psrf_summary.txt
  mcmc_diagnostics/psrf_plot_data.mat


🔧 CUSTOMIZATION
----------------

1. Change initial parameter guesses:
   Edit STEP3_test_simulation.m, lines 26-33

2. Change MCMC settings:
   Edit run_mcmc_with_psrf.m:
   - nIter (line 17): number of iterations per chain
   - burn (line 18): burn-in period
   - n_chains (line 19): number of chains

3. Change solver settings:
   Edit simulate_multi_fedbatch.m, lines 29-32:
   - RelTol: relative tolerance (smaller = more accurate, slower)
   - AbsTol: absolute tolerance
   - MaxStep: maximum time step

4. Change ODE parameters:
   Edit ode_fedbatch.m, lines 17-19:
   - AMM_THRESHOLD: ammonia inhibition threshold
   - MIN_G: minimum glucose constraint
   - MIN_Q: minimum glutamine constraint


⚡ PERFORMANCE
--------------
Typical timing on modern computer:

STEP1: 1-2 seconds
STEP2: 3-5 seconds  
STEP3: 2-3 seconds
STEP4: 20-40 minutes (OPTIONAL!)
STEP5: 5-10 seconds

Total with MCMC: ~30-45 minutes
Total without MCMC: ~15 seconds (quick mode)


❓ TROUBLESHOOTING
------------------

Problem: "data.xlsx not found"
Solution: Make sure data.xlsx is in the current MATLAB directory
          Check with: >> pwd

Problem: "workspace_stepX.mat not found"
Solution: Run previous steps first
          Or run: >> RUN_ALL

Problem: Simulation very slow (>30 sec)
Solution: Check solver settings in simulate_multi_fedbatch.m
          Default should be fast (RelTol=1e-3, AbsTol=1e-6)

Problem: Constants contain NaN
Solution: Check your data.xlsx
          Make sure it has valid data from 90-200h

Problem: MCMC not converging (PSRF > 1.2)
Solution: Increase iterations in run_mcmc_with_psrf.m
          Or adjust parameter bounds


📧 SUPPORT
----------
Check generated files for debugging:
- workspace_*.mat: Contains all variables at each step
- mcmc_diagnostics/: MCMC convergence information

Load workspace to inspect:
>> load workspace_step2.mat
>> constants  % View fitted constants
>> batches{1} % View first batch data


========================================
GOOD LUCK WITH YOUR ANALYSIS! 🔬
========================================
