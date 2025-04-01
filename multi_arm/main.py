#!/usr/bin/env python3

"""
run_power_curves.py

An example 'execution script' that uses the multi-arm biomarker simulation
code to:
  - Sweep over various sample sizes
  - Estimate power (per biomarker + overall success) for each sample size
  - Plot power curves

Prerequisites:
  - multi_arm_biomarker_sim.py in the same directory
  - python packages: numpy, matplotlib, statsmodels, pandas
"""

import numpy as np
import matplotlib.pyplot as plt
from multi_arm_biomarker_sim import estimate_power

def main():
    # Random seed for reproducibility
    np.random.seed(42)
    
    # Suppose we have K=4 experimental arms, each with its own biomarker:
    #  - e.g. TMS, Ketamine, ECT, Voice-biomarker
    #  - hypothetical interaction effects based on prior literature
    K = 4
    beta_main = [0.2, 0.4, 0.3, 0.2]     # main effects vs. implicit baseline (e.g., Arm 1)
    beta_int  = [0.3, 0.6, 0.2, 0.3]     # distinct biomarker interactions
    
    # Biomarker generation parameters:
    # - mu_X: mean for each biomarker
    # - gamma_X: correlation strength with latent factor (can be different per biomarker)
    # - sigma_X: residual standard deviation (can be different per biomarker)
    mu_X = [0.0, 0.0, 0.0, 0.0]          # means for each biomarker
    gamma_X = [0.5, 0.7, 0.3, 0.4]        # correlation strengths with latent factor
    sigma_X = [1.0, 0.8, 1.2, 1.0]        # residual standard deviations
    sigma_y = 1.0                          # outcome noise

    # Extended sample sizes per arm, including up to 600:
    sample_sizes = [50, 100, 150, 200, 300, 400, 500, 600]

    alpha = 0.05
    n_sims = 2000
    
    # We'll compare both "any" and "all" success criteria
    success_criteria_list = ["any", "all"]
    
    # We can compare two multiple-testing methods: "holm" vs "BH"
    methods = ["holm", "BH"]
    
    for success_criteria in success_criteria_list:
        for method in methods:
            # Prepare arrays to store results
            #  shape: (len(sample_sizes), K) for biomarker-level power
            power_biomarkers = np.zeros((len(sample_sizes), K))
            # shape: (len(sample_sizes),) for overall success
            success_rates = np.zeros(len(sample_sizes))
            
            for i, n_per_arm in enumerate(sample_sizes):
                power_per_biom, overall_success = estimate_power(
                    K=K,
                    n_per_arm=n_per_arm,
                    beta_main=beta_main,
                    beta_int=beta_int,
                    mu_X=mu_X,
                    gamma_X=gamma_X,
                    sigma_X=sigma_X,
                    sigma_y=sigma_y,
                    alpha=alpha,
                    method=method,
                    success_criteria=success_criteria,
                    n_sims=n_sims
                )
                power_biomarkers[i, :] = power_per_biom
                success_rates[i] = overall_success
                print(f"[{method.upper()}] criteria='{success_criteria}', n={n_per_arm} => "
                      f"per-biomarker={power_per_biom}, overall={overall_success:.3f}")

            # --- Make a separate figure for each method and success criteria ---
            # 1) Plot per-biomarker power
            plt.figure(figsize=(10, 6))
            for k_i in range(K):
                plt.plot(sample_sizes, power_biomarkers[:, k_i],
                        marker='o', linestyle='-',
                        label=f"Biomarker {k_i+1} (β={beta_int[k_i]})")
                        
            plt.xlabel("Sample Size per Arm")
            plt.ylabel("Power (per biomarker)")
            plt.title(f"Per-Biomarker Power (method={method})")  # Remove success criteria
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.savefig(f"figs/biomarker_power_{method}.png", dpi=300)  # Remove success criteria from filename
            plt.show()

            # 2) Plot overall success rate
            plt.figure(figsize=(10, 6))
            total_sample_sizes = [K * n for n in sample_sizes]
            plt.plot(total_sample_sizes, success_rates, marker='o', linestyle='-', color='darkblue')
            plt.xlabel("Total Sample Size (N)")
            plt.ylabel("Overall Success Rate")
            plt.title(f"Overall Success Rate ('{success_criteria}') (method={method})")
            plt.grid(True, alpha=0.3)
            plt.savefig(f"figs/overall_success_{method}_{success_criteria}.png", dpi=300)
            plt.show()

if __name__ == "__main__":
    main()
