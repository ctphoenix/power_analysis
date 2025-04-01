#!/usr/bin/env python3

"""
multi_arm_biomarker_sim.py

Simulates a multi-arm RCT with K experimental arms. Each arm k has a
corresponding biomarker X_k. Generates data under a linear model with
optional correlations among biomarkers, fits the model, and applies
multiple-testing correction for the K interaction tests (Tk:Xk).

Features:
  - K experimental arms (no SoC).
  - Different effect sizes for each biomarker (beta_int array).
  - Holm or Benjamini–Hochberg corrections for multiple comparisons.
  - Success criterion can be 'any' or 'all' significant interactions.

Dependencies:
    numpy, pandas, statsmodels, (optional) matplotlib
"""

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# --------------------------------------------------------------------------
# 1) DATA GENERATION
# --------------------------------------------------------------------------

def generate_data(n_per_arm, K, beta_main, beta_int, mu_X, gamma_X, sigma_X, sigma_y, p=None):
    """
    Generate one dataset for a multi-arm trial:
      - K experimental arms (indexed 1 to K).
      - Each arm k has a paired biomarker X_k with interaction effect beta_int[k-1].
    
    Args:
        n_per_arm : int
            Sample size PER ARM. Total N = K * n_per_arm.
        K : int
            Number of experimental arms.
        beta_main : array-like of length K
            Main effect of each treatment k (relative to implicit baseline).
        beta_int : array-like of length K
            Interaction effect for each arm-biomarker pair (Tk * Xk).
        mu_X : array-like of length K
            Mean of each biomarker X_k.
        gamma_X : array-like of length K
            Loading for latent factor U -> X_k for each biomarker, controlling correlation.
            If all 0 => X_k are independent.
        sigma_X : array-like of length K
            Residual std dev for each biomarker's random component.
        sigma_y : float
            Residual std dev for the outcome.
        p : array-like of length K, optional
            Randomization probabilities for arms 1..K. If None, default is equal allocation.

    Returns:
        pd.DataFrame with columns: ["T","Y","X1",...,"XK"].
            T=k => experimental arm k (1..K).
    """

    if p is None:
        # Equal allocation across K arms
        p = np.ones(K) / K
    assert len(p) == K, f"Probabilities array p should have length K={K}"
    assert abs(sum(p) - 1) < 1e-8, "Probabilities must sum to 1."

    N = K * n_per_arm

    # Treatment assignment (arms 1 to K)
    T = np.random.choice(a=range(1, K + 1), size=N, p=p)

    # Latent factor for correlation among biomarkers
    U = np.random.normal(loc=0, scale=1, size=N)

    # Generate biomarkers
    X = np.zeros((N, K))
    for k in range(K):
        # X_{ik} = mu_X[k] + gamma_X[k]*U_i + sigma_X[k]*epsilon
        eps = np.random.normal(0, 1, N)
        X[:,k] = mu_X[k] + gamma_X[k]*U + sigma_X[k]*eps

    # Generate outcome
    beta0 = 0.0
    # (Optionally you could have a "main effect" of each biomarker X_k,
    #  but let's set them = 0 here for simplicity.)
    beta_2 = np.zeros(K)

    Y = np.zeros(N)
    noise = np.random.normal(0, sigma_y, N)

    for i in range(N):
        t_i = T[i] # Arm index from 1 to K
        y_val = beta0
        # Add main effect for the assigned arm k (indexed 0 to K-1 in beta_main)
        y_val += beta_main[t_i - 1]
        # Add biomarker main effects (beta_2) and interaction effects (beta_int)
        for k in range(K): # Iterate through biomarker indices 0 to K-1
            # Main effect of X_k (X[i, k]) for everyone
            y_val += beta_2[k] * X[i, k]
            # Interaction effect if assigned to arm k+1 (T_i = k+1)
            # Interaction beta_int[k] corresponds to biomarker X[i, k]
            if t_i == (k + 1):
                y_val += beta_int[k] * X[i, k]
        y_val += noise[i]
        Y[i] = y_val

    # Build DataFrame
    data = pd.DataFrame({
        "T": T,
        "Y": Y,
    })
    for k in range(K):
        data[f"X{k+1}"] = X[:,k]

    return data

# --------------------------------------------------------------------------
# 2) FIT MODEL & EXTRACT P-VALUES
# --------------------------------------------------------------------------

def fit_linear_model_and_get_pvalues(data, K):
    """
    Fit the linear model using manually created interaction terms:
        Y = beta0 + sum_{k=2..K} main_k*I(T=k)
           + sum_{j=1..K} beta2j*Xj
           + sum_{k=1..K} beta_int_k * I(T=k)*Xk
        + error
    where T=1 is the reference arm by default using dummy coding.

    Returns:
        (coefs, pvals) for the K interaction terms (Tk*Xk).
        coefs, pvals each length K in the order k=1..K.
    """
    # Create dummy variables for T (arms 1 to K)
    # drop_first=True makes T1 the reference category
    T_dummies = pd.get_dummies(data['T'], prefix='T', drop_first=True)
    df_fit = pd.concat([data, T_dummies], axis=1)

    # --- Formula Construction ---
    # Start with main effects of treatments T2...TK (relative to T1)
    treatment_terms = list(T_dummies.columns) # ['T2', 'T3', ..., 'TK']

    # Add main effects of all biomarkers X1...XK
    biomarker_terms = [f"X{k_i}" for k_i in range(1, K + 1)]

    # Create specific interaction terms Tk*Xk for k=1...K
    interaction_terms = []
    interaction_term_names_in_formula = []
    for k_i in range(1, K + 1):
        term_name = f"T{k_i}X{k_i}"
        interaction_terms.append(term_name)
        # Create the actual interaction column in the dataframe
        # I(T=k_i) * Xk_i
        df_fit[term_name] = (df_fit['T'] == k_i).astype(float) * df_fit[f"X{k_i}"]
        interaction_term_names_in_formula.append(term_name) # We need coefficients for these

    # Combine all terms for the formula string
    # Model: Y ~ T2+...+TK + X1+...+XK + T1X1+...+TKXK
    all_terms = treatment_terms + biomarker_terms + interaction_term_names_in_formula
    formula = f"Y ~ {' + '.join(all_terms)}"

    # --- Fit Model ---
    try:
        model = smf.ols(formula=formula, data=df_fit).fit()
    except Exception as e:
        print(f"Error fitting model with formula: {formula}")
        print(f"Data columns: {df_fit.columns}")
        print(f"Exception: {e}")
        # Return NaNs or raise error if fitting fails
        nan_array = np.full(K, np.nan)
        return nan_array, nan_array


    # --- Extract Interaction P-values ---
    coefs = []
    pvals = []
    # Get results for the interaction terms T1X1, T2X2, ..., TKXK
    for term_name in interaction_term_names_in_formula:
        coefs.append(model.params.get(term_name, np.nan))
        pvals.append(model.pvalues.get(term_name, np.nan))

    # Check if any p-values couldn't be extracted
    if any(np.isnan(pvals)):
         print(f"Warning: Could not extract all interaction p-values. Model summary:")
         # print(model.summary()) # Optional: print summary for debugging
         # Find which terms were missing
         missing_terms = [term for term, p in zip(interaction_term_names_in_formula, pvals) if np.isnan(p)]
         print(f"Missing p-values for terms: {missing_terms}")
         print(f"Available model coefficients: {list(model.params.index)}")


    return np.array(coefs), np.array(pvals)

# --------------------------------------------------------------------------
# 3) MULTIPLE-TESTING CORRECTIONS
# --------------------------------------------------------------------------

def holm_step_down(pvals, alpha=0.05):
    """
    Holm's step-down procedure for controlling the FWER.

    Returns an array of booleans of the same length as pvals, indicating which
    hypotheses are rejected at level alpha. 
    """
    m = len(pvals)
    # sort p-values ascending
    order = np.argsort(pvals)
    sorted_p = pvals[order]
    reject = np.zeros(m, dtype=bool)

    for i in range(m):
        # Holm threshold = alpha/(m-i)
        if sorted_p[i] <= alpha/(m - i):
            # reject
            reject[order[i]] = True
        else:
            # once we fail to reject at position i, higher pvals won't be rejected in Holm step-down
            # but we must still check if a smaller pval after re-sorting can pass. 
            # Actually in Holm's step-down, we do "continue" to check the others, 
            # but the logic is that once we find a fail, all subsequent are fails as well.
            break
    # all subsequent tests in the sorted list are also not rejected
    return reject

def benjamini_hochberg(pvals, alpha=0.05):
    """
    Benjamini–Hochberg procedure for controlling FDR.
    
    Returns a boolean array (same length as pvals) indicating which 
    hypotheses are rejected.
    """
    m = len(pvals)
    order = np.argsort(pvals)
    sorted_p = pvals[order]

    # BH critical values: alpha * (i/m)
    threshold = np.array([(i+1)/m * alpha for i in range(m)])
    below = sorted_p <= threshold

    # The largest i for which p_i <= alpha*(i/m)
    if not np.any(below):
        # no rejections
        return np.zeros(m, dtype=bool)

    max_index = np.where(below)[0][-1]  # last (highest) index where condition holds
    reject = np.zeros(m, dtype=bool)
    reject[order[:max_index+1]] = True
    return reject

# --------------------------------------------------------------------------
# 4) APPLY CORRECTION & DEFINE "SUCCESS" CRITERIA
# --------------------------------------------------------------------------

def test_biomarkers(pvals, alpha=0.05, method="holm", success_criteria="any"):
    """
    Given a list of p-values for each biomarker's interaction,
    apply a multiple-testing method and then define "success."

    Args:
        pvals : array-like of length K
        alpha : float
        method: "holm" or "BH"
        success_criteria: "any" or "all"
            - "any" => success if at least one biomarker is significant
            - "all" => success if all biomarkers are significant
    
    Returns:
        (reject_array, success_boolean)
        reject_array is bool array of length K (which biomarkers are significant)
        success_boolean is True/False based on the success_criteria
    """
    if method.lower() == "holm":
        reject_array = holm_step_down(pvals, alpha=alpha)
    elif method.lower() in ["bh", "benjamini-hochberg"]:
        reject_array = benjamini_hochberg(pvals, alpha=alpha)
    else:
        # No correction
        reject_array = (pvals < alpha)

    if success_criteria == "any":
        success = np.any(reject_array)
    elif success_criteria == "all":
        success = np.all(reject_array)
    else:
        raise ValueError("Unrecognized success_criteria. Use 'any' or 'all'.")
    
    return reject_array, success

# --------------------------------------------------------------------------
# 5) POWER ESTIMATION BY SIMULATION
# --------------------------------------------------------------------------

def estimate_power(K, n_per_arm, beta_main, beta_int, mu_X, gamma_X, sigma_X, sigma_y,
                   alpha=0.05, method="holm", success_criteria="any",
                   n_sims=1000):
    """
    Estimate power by running multiple simulations.
    Uses the updated generate_data and fit_linear_model functions.
    """
    success_count = 0
    all_reject_arrays = [] # Store reject arrays from each sim

    for sim in range(n_sims):
        if (sim + 1) % 100 == 0:
            print(f"Running simulation {sim + 1} / {n_sims}")

        # 1. Generate data (K arms, N = K*n_per_arm total)
        data = generate_data(n_per_arm, K, beta_main, beta_int, mu_X, gamma_X, sigma_X, sigma_y)

        # 2. Fit model and get interaction p-values
        coefs, pvals = fit_linear_model_and_get_pvalues(data, K)

        # Handle potential NaN p-values from model fitting issues
        if np.any(np.isnan(pvals)):
            print(f"Warning: NaN p-values in sim {sim+1}. Skipping this simulation for power estimate.")
            # Optionally, count failures or handle differently
            continue # Skip this sim from power calculation

        # 3. Test biomarkers using specified method and criteria
        reject_array, success = test_biomarkers(pvals, alpha, method, success_criteria)

        all_reject_arrays.append(reject_array)
        if success:
            success_count += 1

    power = success_count / n_sims # Be careful if sims were skipped

    # Calculate per-biomarker power (proportion of times each biomarker was rejected)
    if all_reject_arrays:
      per_biomarker_power = np.mean(np.array(all_reject_arrays), axis=0)
    else:
      per_biomarker_power = np.zeros(K) # Or handle as appropriate

    print(f"\nPower Estimation Complete ({n_sims} simulations):")
    print(f"  Overall Power ({success_criteria} criteria, {method} correction): {power:.3f}")
    print(f"  Per-Biomarker Power ({method} correction):")
    for k in range(K):
        print(f"    Biomarker {k+1}: {per_biomarker_power[k]:.3f}")

    return power, per_biomarker_power

# --------------------------------------------------------------------------
# 6) DEMO / MAIN
# --------------------------------------------------------------------------

def main():
    np.random.seed(42)

    # --- Simulation Setup ---
    K = 4 # Number of EXPERIMENTAL arms
    n_per_arm_list = [25, 50, 100, 150, 200, 250] # Sample size PER ARM
    N_total_list = [K * n for n in n_per_arm_list] # Total N

    # Effect sizes (length K)
    beta_main = np.array([0.1, 0.15, 0.2, 0.25]) # Main effects vs implicit baseline
    beta_int = np.array([0.2, 0.25, 0.3, 0.35]) # Interaction effects Tk*Xk

    # Biomarker properties (length K)
    mu_X = np.zeros(K) # Mean of each biomarker
    gamma_X = np.full(K, 0.3) # Correlation factor (e.g., 0.3 => moderate corr)
    sigma_X = np.ones(K) # Std dev of biomarker noise

    # Outcome noise
    sigma_y = 1.0

    # Simulation parameters
    alpha = 0.05
    n_sims = 500 # Reduced for quicker testing; increase for accuracy (e.g., 1000+)

    # --- Run Power Simulations ---
    results_holm_any = []
    results_bh_any = []
    per_biomarker_holm = {n: [] for n in n_per_arm_list}
    per_biomarker_bh = {n: [] for n in n_per_arm_list}

    print(f"Starting power simulation: K={K} arms, N per arm={n_per_arm_list}")

    for n_per_arm in n_per_arm_list:
        print(f"\n--- N per arm = {n_per_arm} (Total N = {K * n_per_arm}) ---")
        # Holm correction, 'any' success
        power_holm, pbp_holm = estimate_power(
            K, n_per_arm, beta_main, beta_int, mu_X, gamma_X, sigma_X, sigma_y,
            alpha=alpha, method="holm", success_criteria="any", n_sims=n_sims
        )
        results_holm_any.append(power_holm)
        per_biomarker_holm[n_per_arm] = pbp_holm

        # BH correction, 'any' success
        power_bh, pbp_bh = estimate_power(
            K, n_per_arm, beta_main, beta_int, mu_X, gamma_X, sigma_X, sigma_y,
            alpha=alpha, method="bh", success_criteria="any", n_sims=n_sims
        )
        results_bh_any.append(power_bh)
        per_biomarker_bh[n_per_arm] = pbp_bh

    # --- Plotting (Optional) ---
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nMatplotlib not found. Skipping plots.")
        plt = None

    if plt:
        plt.figure(figsize=(12, 10))

        # Plot overall power
        plt.subplot(2, 1, 1)
        plt.plot(N_total_list, results_holm_any, marker='o', linestyle='-', label="Overall Power (Holm, Any Success)")
        plt.plot(N_total_list, results_bh_any, marker='s', linestyle='--', label="Overall Power (BH, Any Success)")
        plt.title(f'Overall Power vs. Total Sample Size (K={K})')
        plt.xlabel('Total Sample Size (N)')
        plt.ylabel('Power (Any Significant Biomarker)')
        plt.grid(True, linestyle=':')
        plt.ylim(0, 1.05)
        plt.legend()

        # Plot per-biomarker power (using BH correction as example)
        plt.subplot(2, 1, 2)
        pbp_bh_array = np.array([per_biomarker_bh[n] for n in n_per_arm_list]) # Shape (n_sizes, K)
        for k in range(K):
            plt.plot(N_total_list, pbp_bh_array[:, k], marker='^', linestyle=':', label=f'Biomarker {k+1} Power (BH)')

        plt.title(f'Per-Biomarker Power vs. Total Sample Size (BH Correction)')
        plt.xlabel('Total Sample Size (N)')
        plt.ylabel('Power')
        plt.grid(True, linestyle=':')
        plt.ylim(0, 1.05)
        plt.legend()


        plt.tight_layout()
        # Save the figure
        figname = f"figs/power_curves_K{K}_sims{n_sims}.png"
        plt.savefig(figname)
        print(f"\nPlots saved to {figname}")
        # plt.show() # Optionally display plot interactively

if __name__ == "__main__":
    main()
