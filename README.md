# Power Analysis

This repository contains code, simulations, and documentation related to statistical power analysis for clinical trial designs.

---

## Overview

This repository currently houses two main projects:

1.  **`multi_arm/`**: Python simulations and LaTeX documents for designing and analyzing multi-arm clinical trials aimed at validating treatment-predictive biomarkers.
2.  **`power_calc_app/`**: An R Shiny application potentially related to power calculations (details inferred from file names).

---

## Multi-Arm Biomarker Trial Simulations (`multi_arm/`)

### Purpose

This project provides tools to simulate multi-arm clinical trials where the goal is to confirm if specific, pre-defined patient biomarkers predict who responds best to different experimental treatments. It helps estimate the statistical power and required sample size for such complex trial designs.

### Contents

*   **Python Scripts (`.py`)**: Code for generating simulated trial data and analyzing the results to estimate power.
    *   `multi_arm_biomarker_sim.py`: Core simulation functions.
    *   `main.py`: Example script to run simulations and generate plots.
*   **LaTeX Documents (`.tex`)**: Detailed manuscript and executive summary describing the methodology, rationale, and simulation results.
*   **Bibliography (`.bib`)**: References for the LaTeX documents.
*   **Figures (`figs/`)**: Output directory for power curve plots generated by the Python scripts.
*   **Slides (`slides.md`, `slides.pdf`)**: Presentation summarizing the project.
*   **Plan (`plan.md`)**: Technical specification for a future, extended simulation framework.

### Running the Simulations

1.  Navigate to the `multi_arm` directory.
2.  Ensure you have Python 3 and necessary libraries (likely `numpy`, `pandas`, `statsmodels`, `matplotlib`) installed.
3.  Run the example simulation script:
    ```bash
    python main.py 
    ```
4.  Check the `figs/` directory for output plots.
5.  Compile the LaTeX documents (`main.tex`, `executive_summary.tex`) using a LaTeX distribution (like TeX Live, MacTeX, MiKTeX) or an online editor like Overleaf.

---

## Power Calculation Shiny App (`power_calc_app/`)

### Purpose

This directory contains an R Shiny application, likely providing an interactive interface for performing certain types of power calculations.

### Contents

*   **`app.R`**: The R script containing the Shiny application's UI (User Interface) and server logic.
*   **`run.R`**: An R script used to deploy the application to shinyapps.io (using the `rsconnect` package).
*   **`power.tex`**: A LaTeX document, possibly related to the theoretical basis of the calculations in the app.
*   **`rsconnect/`**: Directory containing deployment configuration files for shinyapps.io.

### Running/Deploying the App

*   **To run locally (requires R and Shiny package):**
    1.  Open R/RStudio.
    2.  Set the working directory to `power_calc_app/`.
    3.  Install necessary packages if not already present (e.g., `install.packages('shiny')`).
    4.  Run the command: `shiny::runApp()`
*   **To re-deploy to shinyapps.io:**
    1.  Ensure you have the `rsconnect` package installed and configured with appropriate account credentials (token/secret).
    2.  The script `run.R` contains the deployment command. You can execute it from R/RStudio:
        ```R
        # Make sure your working directory is the repository root
        # or adjust the path in run.R if needed.
        source("power_calc_app/run.R") 
        ```

--- 