function [constants, MAT] = fit_constants_fedbatch(batches)
% FIT_CONSTANTS_FEDBATCH Fit kinetic constants with multiple feeding correction
%
% Version: WITHOUT LACTATE (4 yields only)
%
% Syntax:
%   [constants, MAT] = fit_constants_fedbatch(batches)
%
% Inputs:
%   batches - Cell array of batch data structures (concentrations in mM)
%
% Outputs:
%   constants - Structure with fitted constants:
%       .mu_max    - Maximum specific growth rate (h^-1)
%       .Yxv_glc   - Yield of cells on glucose ((10^6 cells/mL)/mM)
%       .Yxv_gln   - Yield of cells on glutamine ((10^6 cells/mL)/mM)
%       .Yamm_gln  - Yield of ammonia on glutamine (mM/mM)
%   MAT - Matrix of fitted values for each batch [n x 4]
%
% Method:
%   - μ_max: Linear regression on log(Xv) from 90-200h
%   - Yields: Mass balance with multiple feeding correction (90-200h)
%   - All concentrations in mM (converted from g/L in load_excel_data)

    % ============================================================
    % ============== USER PARAMETERS (MODIFY HERE) ===============
    % ============================================================
    
    % Time range for fitting (hours)
    FIT_TIME_START = 0;   % Start time for exponential phase fitting (h)
    FIT_TIME_END = 70;    % End time for exponential phase fitting (h)
    
    % Time range for YIELD calculation (Batch phase: 0h to 144h, following Xing 2009)
    YIELD_TIME_START = 0;  % Start time for Yield Mass Balance (h)
    YIELD_TIME_END = 144;  % End time for Yield Mass Balance (h) - 144h (6 days)

    % Minimum data points required
    MIN_DATA_POINTS = 3;   % Minimum points in time range for valid fit
    
    % ============================================================
    % ==================== END USER PARAMETERS ===================
    % ============================================================

    n = numel(batches);
    MAT = zeros(n, 4);  % Changed from 5 to 4 (no lactate)
    
    % ========================================
    % μ_max: Calculate for each batch individually
    % ========================================
    mu_max_values = [];
    
    for b = 1:n
        bb = batches{b};
        idx_exp = bb.t >= FIT_TIME_START & bb.t <= FIT_TIME_END;
        
        if sum(idx_exp) > MIN_DATA_POINTS
            p = polyfit(bb.t(idx_exp), log(max(bb.x(idx_exp), eps)), 1);
            mu_max_values = [mu_max_values; p(1)];  % Collect all mu_max
        end
    end
    
    % Use 95th percentile (upper bound)
    if ~isempty(mu_max_values)
        mu_max_95th = prctile(mu_max_values, 95);
    else
        mu_max_95th = 0;
    end
    
    for b = 1:n
        bb = batches{b};
        
        % ========================================
        % μ_max: Individual batch value
        % ========================================
        idx_exp = bb.t >= FIT_TIME_START & bb.t <= FIT_TIME_END;
        if sum(idx_exp) > MIN_DATA_POINTS
            p = polyfit(bb.t(idx_exp), log(max(bb.x(idx_exp), eps)), 1);
            MAT(b, 1) = p(1);  % Individual batch μ_max
        end
        
        % ========================================
        % Yields: Mass balance with multiple feeding
        % ========================================

        % *** YIELD 계산을 위한 새로운 인덱스 정의 (0h ~ 144h) ***
    idx_yield = bb.t >= YIELD_TIME_START & bb.t <= YIELD_TIME_END;
    day_idx = find(idx_yield); % 이 인덱스가 Y 계산에 사용됨

if length(day_idx) < 2
        continue; 
    end
    
    t_start_idx = day_idx(1); % t = 0h에 해당하는 인덱스
    t_end_idx = day_idx(end); % t = 144h에 해당하는 인덱스
        
        % ================================================
        % Calculate pre-feeding volumes
        % Measured volumes are POST-FEEDING, but concentrations are PRE-FEEDING
        % ================================================
        V_start_post = bb.vol(t_start_idx);  % mL (post-feeding)
        V_end_post = bb.vol(t_end_idx);      % mL (post-feeding)
        
        % Feeding at each time point
        feed_at_start = bb.feed_glc_vol(t_start_idx) + sum(bb.feed_vols(t_start_idx, :));  % mL
        feed_at_end = bb.feed_glc_vol(t_end_idx) + sum(bb.feed_vols(t_end_idx, :));        % mL
        
        % Pre-feeding volumes (to match with pre-feeding concentrations)
        V_start = V_start_post - feed_at_start;  % mL (pre-feeding)
        V_end = V_end_post - feed_at_end;        % mL (pre-feeding)
        
        % ================================================
        % Calculate total glucose fed (pure + all feeds)
        % Units: mM × mL = μmol
        % ================================================
        % Pure glucose feed
        feed_glc_in_period = sum(bb.feed_glc_vol(t_start_idx:t_end_idx-1));  % mL
        GLC_fed_pure = feed_glc_in_period * bb.glc_stock;  % μmol (mM × mL)
        
        % Feed media (multiple feeds)
        GLC_fed_media = 0;
        GLN_fed_media = 0;
        
        for f = 1:bb.n_feeds
            feed_vol = sum(bb.feed_vols(t_start_idx:t_end_idx-1, f));  % mL
            GLC_fed_media = GLC_fed_media + feed_vol * bb.feed_glc(f);  % μmol
            GLN_fed_media = GLN_fed_media + feed_vol * bb.feed_gln(f);  % μmol
        end
        
        % Total fed (μmol)
        GLC_fed = GLC_fed_pure + GLC_fed_media;
        GLN_fed = GLN_fed_media;  % Only from feed media
        
        % ================================================
        % Mass balance - TOTAL AMOUNT BASED (following literature)
        % ================================================
        % Cell calculation: 
        %   bb.x: 10^6 cells/mL (measured pre-feeding)
        %   V: mL (calculated pre-feeding volume)
        %   bb.x * V = (10^6 cells/mL) × (mL) = 10^6 cells ✓
        dXv = bb.x(t_end_idx)*V_end - bb.x(t_start_idx)*V_start;  % 10^6 cells
        
        % Substrate/metabolite calculations: mM × mL = μmol, convert to mmol
        dGLC_consumed = (bb.glc(t_start_idx)*V_start - bb.glc(t_end_idx)*V_end + GLC_fed) / 1000;  % mmol
        dGLN_consumed = (bb.gln(t_start_idx)*V_start - bb.gln(t_end_idx)*V_end + GLN_fed) / 1000;  % mmol
        % dLAC_produced REMOVED - lactate not modeled
        dNH4_produced = (bb.nh4(t_end_idx)*V_end - bb.nh4(t_start_idx)*V_start) / 1000;  % mmol
        
        % ================================================
        % Yields - CONCENTRATION-BASED for ODE compatibility
        % ================================================
        % Convert from total basis to concentration basis (mL→L)
        % Y_total = (10^6 cells)/(mmol) from mass balance
        % Y_conc = Y_total / 1000 for use with mM (mmol/L) units in ODE
        % This accounts for X being in cells/mL and G being in mmol/L
        if dGLC_consumed > 0
            MAT(b, 2) = dXv / dGLC_consumed / 1000;  % Y_xv/glc ((10^6 cells/mL)/mM)
        end
        if dGLN_consumed > 0
            MAT(b, 3) = dXv / dGLN_consumed / 1000;  % Y_xv/gln ((10^6 cells/mL)/mM)
        end
        % Y_lac/glc REMOVED - lactate not modeled
        if dGLN_consumed > 0
            MAT(b, 4) = abs(dNH4_produced / dGLN_consumed);  % Y_amm/gln (mmol/mmol) - always positive
        end
    end
    
    % Average yields across batches, but use 95th percentile for mu_max
    avg = mean(MAT, 1);
    
    constants.mu_max = mu_max_95th;  % Use 95th percentile, not mean!
    constants.Yxv_glc = avg(2);
    constants.Yxv_gln = avg(3);
    % constants.Ylac_glc REMOVED - lactate not modeled
    constants.Yamm_gln = avg(4);  % Changed from avg(5) to avg(4)
    
    % Display results
    fprintf('\n=== Fitted Constants (WITHOUT LACTATE) ===\n');
    fprintf('μ_max     = %.4f h^-1 (95th percentile)\n', constants.mu_max);
    fprintf('  Individual batch values: min=%.4f, mean=%.4f, max=%.4f\n', ...
            min(MAT(:,1)), mean(MAT(:,1)), max(MAT(:,1)));
    fprintf('Y_xv/glc  = %.2e (10^6 cells/mL)/mM (concentration-based)\n', constants.Yxv_glc);
    fprintf('Y_xv/gln  = %.2e (10^6 cells/mL)/mM (concentration-based)\n', constants.Yxv_gln);
    % fprintf('Y_lac/glc = ...') REMOVED
    fprintf('Y_amm/gln = %.3f mM/mM\n\n', constants.Yamm_gln);
end
