library(shiny)
library(rhandsontable) # For the editable grid

# --- Core Simulation Logic (Placeholder) ---
# This function will contain the R translation of the Python simulation
run_power_simulation <- function(J, K, n_per_arm, beta_int_list, gamma, n_sims = 1000, alpha = 0.05) {
  # Placeholder implementation
  
  # 1. Define fixed parameters (as discussed)
  beta_txt <- rep(0.3, J) # Fixed main treatment effects (cosmetic)
  beta_bio <- rep(0, K)   # Fixed main biomarker effects
  mu_X <- rep(0, K)
  sigma_X <- rep(1, K)
  sigma_y <- 1
  
  # 2. Get list of interactions to test from beta_int_list
  hypotheses_to_test <- beta_int_list # List of lists: list(j=j, k=k, beta=beta)
  num_hypotheses <- length(hypotheses_to_test)
  
  if (num_hypotheses == 0) {
    return(list(
      power_per_interaction = data.frame(j=integer(), k=integer(), Interaction=character(), Power=numeric(), SE=numeric()),
      power_holm_any = list(Power=NA, SE=NA),
      power_holm_all = list(Power=NA, SE=NA)
    ))
  }

  # Store rejection status: rows=sims, cols=hypotheses
  rejected_status <- matrix(FALSE, nrow = n_sims, ncol = num_hypotheses)
  
  # --- Simulation Loop ---
  for (i in 1:n_sims) {
    # --- A. Generate Data ---
    N_total <- sum(n_per_arm)
    # Assign treatments - ensure it's a factor for lm()
    treatment_assignments <- factor(sample(1:J, size = N_total, replace = TRUE, prob = n_per_arm / N_total))
    
    # Generate latent factor
    U <- rnorm(N_total)
    
    # Generate K biomarkers
    X <- matrix(0, nrow = N_total, ncol = K)
    for (k in 1:K) {
      eps <- rnorm(N_total)
      X[, k] <- mu_X[k] + gamma * U + sigma_X[k] * eps
    }
    colnames(X) <- paste0("X", 1:K)
    
    # Generate Outcome Y
    Y <- numeric(N_total)
    beta0 <- 0.0
    
    # Add main treatment effects (relative to treatment 1 as baseline)
    # We now rely on lm() to handle this via the factor T
    # But Y must still be *generated* with the effects
    Y <- Y + beta0
    for(j in 2:J) {
        Y <- Y + (as.numeric(treatment_assignments) == j) * beta_txt[j] 
    }
    # Or using model.matrix for generation:
    # T_dummies_gen <- model.matrix(~ treatment_assignments) # includes intercept
    # if (J > 1) {
    #    main_eff_contrib_gen <- T_dummies_gen[, 2:J, drop=FALSE] %*% beta_txt[2:J] 
    # } else {
    #    main_eff_contrib_gen <- 0
    # }
    # Y <- Y + beta0 + main_eff_contrib_gen
        
    # Add main biomarker effects
    Y <- Y + X %*% beta_bio
    
    # Add specified interaction effects AND create columns for model
    interaction_term_values <- matrix(0, nrow = N_total, ncol = num_hypotheses)
    interaction_names <- character(num_hypotheses)
    colnames(interaction_term_values) <- interaction_names # Initialize column names

    if (num_hypotheses > 0) {
      for (h_idx in 1:num_hypotheses) {
        hyp <- hypotheses_to_test[[h_idx]]
        j <- hyp$j
        k <- hyp$k
        beta_val <- hyp$beta
        # Use a safe name for the column/term (e.g., no special chars)
        term_col_name <- paste0("Interaction_T", j, "_X", k)
        interaction_names[h_idx] <- term_col_name
        
        # Calculate T_ij * X_ik for this interaction
        term_value <- (as.numeric(treatment_assignments) == j) * X[, k]
        interaction_term_values[, h_idx] <- term_value
        Y <- Y + term_value * beta_val # Add effect to generated Y
      }
      colnames(interaction_term_values) <- interaction_names # Assign final names
    }

    # Add noise
    Y <- Y + rnorm(N_total, 0, sigma_y)
    
    # Combine into a data frame for fitting
    sim_data <- data.frame(Y = Y, T = treatment_assignments, X)
    # Add the *manually calculated* interaction terms needed for the model
    if (num_hypotheses > 0) {
       sim_data <- cbind(sim_data, interaction_term_values)
    }

    # --- B. Fit Model --- 
    # Construct formula: Y ~ T + X1 + ... + XK + Interaction_Tj_Xk_1 + ...
    # Let lm handle dummy coding for factor T
    biomarker_terms <- paste0("X", 1:K)
    
    # Only include terms if they exist
    term_parts <- c("T", biomarker_terms)
    if (num_hypotheses > 0) {
      term_parts <- c(term_parts, interaction_names)
    }
    
    formula_str <- paste("Y ~", paste(term_parts, collapse = " + "))
    
    # Fit the model
    fit <- tryCatch({
      lm(as.formula(formula_str), data = sim_data)
    }, error = function(e) {
      # warning("Model fit failed in sim ", i, ": ", e$message) # Reduce warnings
      NULL
    })

    if (is.null(fit)) {
      rejected_status[i, ] <- FALSE 
      next 
    }
    
    # --- C. Extract P-values ---
    summary_fit <- summary(fit)
    coefs <- summary_fit$coefficients
    
    pvals <- numeric(num_hypotheses)
    if (num_hypotheses > 0) {
        for (h_idx in 1:num_hypotheses) {
          term_name <- interaction_names[h_idx] # The name we gave the column
          if (term_name %in% rownames(coefs)) {
            pvals[h_idx] <- coefs[term_name, "Pr(>|t|)"]
          } else {
            # warning("Interaction term ", term_name, " not found in model results for sim ", i) # Reduce warnings
            pvals[h_idx] <- 1.0 
          }
        }
    }
    
    pvals[is.na(pvals)] <- 1.0

    # --- D. Apply Holm Correction ---
    if (num_hypotheses > 0) {
       corrected_pvals <- p.adjust(pvals, method = "holm")
       rejected_status[i, ] <- (corrected_pvals <= alpha)
    } else {
       # No hypotheses, nothing to reject
    }

  } # End simulation loop

  # --- Calculate Results ---
  # Power per interaction
  if (num_hypotheses > 0) {
    power_per_hyp <- colMeans(rejected_status, na.rm = TRUE)
    se_per_hyp <- sqrt(power_per_hyp * (1 - power_per_hyp) / n_sims)
    
    # Add j and k indices to the results data frame for sorting
    j_indices <- sapply(hypotheses_to_test, function(hyp) hyp$j)
    k_indices <- sapply(hypotheses_to_test, function(hyp) hyp$k)
    interaction_labels <- paste0("T", j_indices, " * X", k_indices)
    
    power_df <- data.frame(
      j = j_indices,
      k = k_indices,
      Interaction = interaction_labels,
      Power = power_per_hyp,
      SE = se_per_hyp
    )
    
    sim_success_any <- rowSums(rejected_status) > 0
    sim_success_all <- rowSums(rejected_status) == num_hypotheses
    
    power_holm_any <- mean(sim_success_any, na.rm = TRUE)
    se_holm_any <- sqrt(power_holm_any * (1 - power_holm_any) / n_sims)
    
    power_holm_all <- mean(sim_success_all, na.rm = TRUE)
    se_holm_all <- sqrt(power_holm_all * (1 - power_holm_all) / n_sims)
    
  } else {
    # Handle case with no hypotheses
    # Include j and k columns for consistency
    power_df <- data.frame(j=integer(), k=integer(), Interaction=character(), Power=numeric(), SE=numeric())
    power_holm_any <- list(Power=NA, SE=NA)
    power_holm_all <- list(Power=NA, SE=NA)
  }

  return(list(
    power_per_interaction = power_df,
    power_holm_any = list(Power = power_holm_any, SE = se_holm_any),
    power_holm_all = list(Power = power_holm_all, SE = se_holm_all)
  ))
}


# --- Shiny UI ---
ui <- fluidPage(
  titlePanel("Multi-Arm Biomarker Trial Power Simulation"),

  fluidRow(
    # --- Column 1: Design & Sample Size ---
    column(5, 
      h4("Study Design"),
      fluidRow(
          column(6, numericInput("num_treatments", "Treatments (J):", value = 4, min = 1, step = 1)),
          column(6, numericInput("num_biomarkers", "Biomarkers (K):", value = 6, min = 1, step = 1))
      ),
      
      hr(),
      h4("Sample Sizes per Arm"),
      uiOutput("sample_size_sliders") # Dynamic sliders based on J
    ),

    # --- Column 2: Interactions, Params, Button & Results ---
    column(7, 
      h4("Interaction Effects"),
      helpText("Enter hypothesized interaction coefficients. Rows=Biomarkers, Cols=Treatments. Leave blank or 0 for no interaction."),
      rHandsontableOutput("interaction_grid"),

      hr(),
      h4("Simulation Parameters"),
      sliderInput("gamma", "Biomarker Correlation Factor (gamma):", min = 0, max = 2, value = 0.3, step = 0.1, width="100%"),
      sliderInput("n_sims", "Number of Simulations:", min = 100, max = 5000, value = 500, step = 100, width="100%"), 
      
      hr(),
      actionButton("run_sim", "Run Simulation", icon = icon("play"), width = "100%"), 
      
      hr(),
      h4("Simulation Results (alpha = 0.05)"),
      verbatimTextOutput("results_text")
    )
  )
)

# --- Shiny Server ---
server <- function(input, output, session) {

  # Reactive value to store interaction grid data
  # Initialize with defaults
  vals <- reactiveValues(
    interaction_data = {
      J_init <- 4; K_init <- 6
      df <- as.data.frame(matrix(0.0, nrow = K_init, ncol = J_init))
      colnames(df) <- paste0("T", 1:J_init); rownames(df) <- paste0("X", 1:K_init)
      if (J_init >= 4 && K_init >= 5) {
        df[1, 1] <- 0.3; df[2, 2] <- 0.6; df[3, 3] <- 0.2; df[4, 4] <- 0.3
        df[5, 3] <- -0.3; df[5, 4] <- -0.3
      }
      df
    }
  )

  # Render the editable interaction grid (original logic with resizing inside)
  output$interaction_grid <- renderRHandsontable({
    J <- input$num_treatments
    K <- input$num_biomarkers
    req(J, K, J >= 1, K >= 1)
    
    current_df <- vals$interaction_data
    
    # Adjust dimensions if needed, preserving existing data
    if (is.null(current_df) || ncol(current_df) != J || nrow(current_df) != K) {
        new_df <- as.data.frame(matrix(0.0, nrow = K, ncol = J))
        colnames(new_df) <- paste0("T", 1:J)
        rownames(new_df) <- paste0("X", 1:K)
        if (!is.null(current_df)) { # Preserve if possible
            common_rows <- min(nrow(current_df), K)
            common_cols <- min(ncol(current_df), J)
            if (common_rows > 0 && common_cols > 0) {
                new_df[1:common_rows, 1:common_cols] <- current_df[1:common_rows, 1:common_cols, drop = FALSE]
            }
        }
        # Update the reactive value if resize occurs
        # Note: This might cause re-rendering issues, but matches previous state
        vals$interaction_data <- new_df 
        current_df <- new_df # Use the new df for rendering this time
    } else if (!is.null(input$interaction_grid)) {
        # If dimensions match, reflect user edits if available
        # This was implicit before, let's try making it explicit
         edited_df <- hot_to_r(input$interaction_grid)
         if(!is.null(edited_df) && all(dim(edited_df) == dim(current_df))) {
             vals$interaction_data <- edited_df
             current_df <- edited_df
         }
    }
    
    rhandsontable(current_df, 
                  rowHeaders = rownames(current_df), 
                  colHeaders = colnames(current_df),
                  stretchH = "all", height = 300) %>%
      hot_cols(type = "numeric", format = "0.0")
  })
  
  # Update the reactive grid data when the table is edited by the user
  # This observer might be needed to capture edits properly before button click
  observeEvent(input$interaction_grid, {
      if (!is.null(input$interaction_grid)) {
          # Check dimensions before updating to avoid issues during resize
          current_J <- isolate(input$num_treatments)
          current_K <- isolate(input$num_biomarkers)
          edited_data <- hot_to_r(input$interaction_grid)
          if(!is.null(edited_data) && nrow(edited_data) == current_K && ncol(edited_data) == current_J) {
             vals$interaction_data <- edited_data
          }
      }
  }, ignoreNULL = TRUE)

  # Dynamically generate sample size sliders
  output$sample_size_sliders <- renderUI({
    J <- input$num_treatments
    req(J, J >= 1)
    
    # Define default values for the first 4 sliders
    default_values <- c(150, 50, 300, 450) 
    
    sliders <- lapply(1:J, function(j) {
      input_id <- paste0("n_trt_", j)
      
      # Set default value based on j, fallback to 100 if j > 4
      value_default <- if(j <= length(default_values)) default_values[j] else 100
      
      # Read previous value if it exists, otherwise use the calculated default
      current_val <- input[[input_id]]
      value_final <- if(!is.null(current_val)) current_val else value_default
      
      sliderInput(inputId = input_id, 
                  label = paste("N for Treatment", j, ":"), 
                  min = 50,           # Updated min value
                  max = 1000,       
                  value = value_final, 
                  step = 50)        
    })
    do.call(tagList, sliders)
  })

  # Reactive expression to run simulation (simplified input reading)
  simulation_results <- eventReactive(input$run_sim, {
      
      J <- input$num_treatments
      K <- input$num_biomarkers
      gamma_val <- input$gamma
      n_sims_val <- input$n_sims
      # Read grid data from the reactive value, assume it's updated by observer/render
      grid_data_current <- vals$interaction_data 
      
      # Basic validation (can keep this)
      req(J, K, gamma_val, n_sims_val, grid_data_current, 
          J >= 1, K >= 1, n_sims_val >= 100)
      req(nrow(grid_data_current) == K && ncol(grid_data_current) == J)
          
      # Read sample sizes (simpler)
      n_per_arm <- vapply(1:J, function(j) {
          input_id <- paste0("n_trt_", j)
          # Assume input exists if req(J) passed
          as.integer(input[[input_id]] %||% 100) # Use %||% for safety
      }, FUN.VALUE = integer(1))
      
      # Process Grid (same as before)
      beta_int_list <- list()
      for (row_k in 1:K) {
          for (col_j in 1:J) {
              beta_val <- tryCatch(as.numeric(grid_data_current[row_k, col_j]), warning = function(w) NA)
              if (!is.na(beta_val) && beta_val != 0) {
                  beta_int_list[[length(beta_int_list) + 1]] <- list(j = col_j, k = row_k, beta = beta_val)
              }
          }
      }
      
      # Run Simulation with Progress (simpler version)
      results <- NULL
      withProgress(message = 'Running Simulation...', value = 0, {
           results <- tryCatch({
               # Use simpler progress update
               incProgress(0.3, detail = "Simulating...")
               res <- run_power_simulation(
                  J = J, K = K, n_per_arm = n_per_arm, beta_int_list = beta_int_list, 
                  gamma = gamma_val, n_sims = n_sims_val, alpha = 0.05
               )
               incProgress(0.7, detail = "Done.")
               Sys.sleep(0.5)
               res 
           }, error = function(e) {
              showNotification(paste("Simulation Error:", e$message), type = "error", duration = 15)
              NULL 
           })
      })
      
      return(results)
      
  })

  # Render the results (ensure it handles NULL gracefully)
  output$results_text <- renderPrint({
      res <- simulation_results()
      # Basic check for NULL
      if(is.null(res)) return("Simulation failed or has not run yet.")
      
      cat("Power per Interaction (Holm Corrected):\n")
      # Sort the results before printing
      power_df_sorted <- res$power_per_interaction
      if (nrow(power_df_sorted) > 0) {
          # Order by j, then k
          power_df_sorted <- power_df_sorted[order(power_df_sorted$j, power_df_sorted$k), ]
          
          # Create the formatted power string
          formatted_power_str <- sprintf("%.3f (SE: %.3f)", 
                                     power_df_sorted$Power, 
                                     power_df_sorted$SE)
                                     
          # Create data frame for printing, applying format() for alignment
          power_output_df <- data.frame(
              Interaction = format(power_df_sorted$Interaction, justify = "left"),
              Power = format(formatted_power_str, justify = "right") # Right-align power
          )
          # Print with options to suppress quotes and row names
          print(power_output_df, quote = FALSE, row.names = FALSE)
          
      } else {
          cat("  (No non-zero interactions specified in the grid)\n")
      }
      
      cat("\nOverall Power (Holm Correction):\n")
      pow_any <- if(is.list(res$power_holm_any) && !is.na(res$power_holm_any$Power)) res$power_holm_any$Power else NA
      se_any <- if(is.list(res$power_holm_any) && !is.na(res$power_holm_any$SE)) res$power_holm_any$SE else NA
      pow_all <- if(is.list(res$power_holm_all) && !is.na(res$power_holm_all$Power)) res$power_holm_all$Power else NA
      se_all <- if(is.list(res$power_holm_all) && !is.na(res$power_holm_all$SE)) res$power_holm_all$SE else NA

      cat(sprintf("  Detecting ANY interaction: %.3f (SE: %.3f)\n", pow_any, se_any))
      cat(sprintf("  Detecting ALL interactions: %.3f (SE: %.3f)\n", pow_all, se_all))
  })

} # End server

# Run the application 
shinyApp(ui = ui, server = server) 
# --- TEST CODE (Comment out shinyApp() call above when testing) ---
# message("--- Running Direct Test of run_power_simulation ---")
# test_results <- tryCatch({
#   run_power_simulation(
#     J = 2, 
#     K = 2, 
#     n_per_arm = c(50, 50), 
#     beta_int_list = list(list(j = 1, k = 1, beta = 0.5)), # Test T1*X1
#     gamma = 0.3, 
#     n_sims = 20, # Small n_sims for quick test
#     alpha = 0.05
#   )
# }, error = function(e) {
#   message("Error during test execution: ", e$message) 
#   NULL
# })
# 
# message("--- Test Results ---")
# print(test_results)
# message("--- Test Complete ---")
# --------------------------------------------------------- 