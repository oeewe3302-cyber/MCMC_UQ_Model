function dydt = ode_fedbatch(~, y, constants, params)
% ODE_FEDBATCH 4-ODE system for CHO cell fed-batch culture
%
% Version: SAFE - NaN prevention for MCMC robustness
% WITHOUT LACTATE - Improved parameter estimation
%
% State Variables:
%   y(1) = Xv  - Viable cell density (10^6 cells/mL)
%   y(2) = GLC - Glucose concentration (mM)
%   y(3) = GLN - Glutamine concentration (mM)
%   y(4) = NH4 - Ammonia concentration (mM)

    % ============================================================
    % USER PARAMETERS
    % ============================================================
    AMM_THRESHOLD = 9;    % Ammonia inhibition threshold (mM)
    MIN_G = 1;            % Minimum glucose (mM)
    MIN_Q = 0.2;          % Minimum glutamine (mM)
    
    % ============================================================
    % SAFETY CHECK: Input validation
    % ============================================================
    % Check for NaN/Inf in state variables
    if any(isnan(y)) || any(isinf(y))
        dydt = zeros(4, 1);
        return;
    end
    
    % Check for NaN/Inf in constants
    const_values = [constants.mu_max, constants.Yxv_glc, constants.Yxv_gln, constants.Yamm_gln];
    if any(isnan(const_values)) || any(isinf(const_values)) || any(const_values <= 0)
        dydt = zeros(4, 1);
        return;
    end
    
    % Check for NaN/Inf in parameters
    param_values = [params.Kglc, params.Kgln, params.KIamm, params.mglc, ...
                    params.a1, params.a2, params.dgln, params.ramm];
    if any(isnan(param_values)) || any(isinf(param_values))
        dydt = zeros(4, 1);
        return;
    end
    
    % ============================================================
    % Extract state variables with MINIMUM CONSTRAINTS
    % ============================================================
    X = max(y(1), 0.1);       % Viable cell density (10^6 cells/mL)
    G = max(y(2), MIN_G);   % Glucose concentration (mM)
    Q = max(y(3), MIN_Q);   % Glutamine concentration (mM)
    A = max(y(4), 0.005);       % Ammonia concentration (mM)
    
    % Additional safety: prevent extreme values
    X = min(X, 1e3);  % Cap at 1000 x 10^6 cells/mL
    G = min(G, 1e3);  % Cap at 1000 mM
    Q = min(Q, 1e2);  % Cap at 100 mM
    A = min(A, 1e2);  % Cap at 100 mM
    
    % ============================================================
    % Growth rate (Monod + ammonia inhibition only)
    % ============================================================
    % Prevent division by very small numbers
    Kglc_safe = max(params.Kglc, 1e-6);
    Kgln_safe = max(params.Kgln, 1e-6);
    KIamm_safe = max(params.KIamm, 1e-6);
    
    mu = constants.mu_max * ...
        (G / (Kglc_safe + G)) * ...
        (Q / (Kgln_safe + Q));
    
    if A >= AMM_THRESHOLD
        mu = mu * (KIamm_safe / (KIamm_safe + A));
    end
    
    mu = max(mu, 0);
    mu = min(mu, 1.0);  % Cap growth rate at 1.0 h^-1 (very high already)
    
    % ============================================================
    % Glutamine maintenance
    % ============================================================
    a2_safe = max(params.a2, 1e-6);  % Prevent division by zero
    mgln = (params.a1 * Q) / (a2_safe + Q);
    mgln = max(mgln, 0);  % Non-negative
    mgln = min(mgln, 1.0);  % Cap at reasonable value
    
    % ============================================================
    % Prevent division by zero in yields
    % ============================================================
    Yxv_glc_safe = max(constants.Yxv_glc, 1e-10);
    Yxv_gln_safe = max(constants.Yxv_gln, 1e-10);
    
    % ============================================================
    % ODE System (4 equations only)
    % ============================================================
    dX = mu * X;
    dG = -(mu / Yxv_glc_safe + params.mglc) * X;
    dQ = -(mu / Yxv_gln_safe + mgln) * X - params.dgln * Q;
    dA = constants.Yamm_gln * (mu / Yxv_gln_safe) * X - ...
         params.ramm * X + params.dgln * Q;
    
    % ★★★ ENFORCE MINIMUM CONSTRAINTS ★★★
    if Q <= MIN_Q && dQ < 0
        dQ = 0;  % Q cannot decrease below MIN_Q
        dA = constants.Yamm_gln * (mu / Yxv_gln_safe) * X - params.ramm * X;
    end
    
    if G <= MIN_G && dG < 0
        dG = 0;  % G cannot decrease below MIN_G
    end
    
    % Stack derivatives (4 equations)
    dydt = [dX; dG; dQ; dA];
    
    % ============================================================
    % FINAL SAFETY CHECK: Prevent NaN/Inf in output
    % ============================================================
    if any(isnan(dydt)) || any(isinf(dydt))
        dydt = zeros(4, 1);
        return;
    end
    
    % Clip to reasonable range
    dydt = max(min(dydt, 1e3), -1e3);
end
