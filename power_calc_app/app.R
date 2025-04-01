library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(shinydashboard)
library(shinyWidgets)

# Helper functions for biomarker-based power analysis

# Calculate probability of selection P(X > X_t)
calc_selection_prob <- function(X_t, sigma_Z = 1) {
  1 - pnorm(X_t)
}

# Calculate E[Z|X > X_t]
calc_expected_Z <- function(X_t, sigma_Z = 1) {
  standardized_threshold <- X_t
  mills_ratio <- dnorm(standardized_threshold) / (1 - pnorm(standardized_threshold))
  mills_ratio
}

# Calculate correlation between X and Y for control group
calc_rho_XY_control <- function(beta_2, sigma_Z, sigma_Y) {
  beta_2 / (sigma_Y * sqrt(1 + sigma_Z^2))
}

# Calculate correlation between X and Y for treatment group
calc_rho_XY_treatment <- function(beta_2, beta_3, sigma_Z, sigma_Y) {
  (beta_2 + beta_3) / (sigma_Y * sqrt(1 + sigma_Z^2))
}

# Convert from Cohen's d and correlations to beta parameters
convert_params_to_betas <- function(d_overall, rho_XY_control, rho_XY_treatment, sigma_Y, sigma_Z) {
  # Calculate beta parameters from Cohen's d and correlation
  beta_1 <- d_overall * sigma_Y
  beta_2 <- rho_XY_control * sigma_Y * sqrt(1 + sigma_Z^2)
  beta_3 <- (rho_XY_treatment - rho_XY_control) * sigma_Y * sqrt(1 + sigma_Z^2)
  
  list(beta_1 = beta_1, beta_2 = beta_2, beta_3 = beta_3)
}

# Calculate effect size delta in the selected subpopulation
calc_effect_size <- function(d_overall, rho_XY_control, rho_XY_treatment, X_t, sigma_Z = 1, sigma_Y = 1) {
  beta_1 <- d_overall
  beta_3 <- rho_XY_treatment - rho_XY_control
  
  # Calculate expected Z in selected subpopulation
  expected_Z <- calc_expected_Z(X_t)
  
  # Calculate effect size
  beta_1 + beta_3 * expected_Z
}

# Calculate Cohen's d
calc_cohens_d <- function(d_overall, rho_XY_control, rho_XY_treatment, X_t, sigma_Z = 1, sigma_Y = 1) {
  calc_effect_size(d_overall, rho_XY_control, rho_XY_treatment, X_t)
}

# Calculate power for two-group comparison
calc_power_two_groups <- function(d_overall, rho_XY_control, rho_XY_treatment, X_t, 
                                  sigma_Z, sigma_Y, n, p_T, alpha = 0.05) {
  # Effect size
  delta <- calc_effect_size(d_overall, rho_XY_control, rho_XY_treatment, X_t, sigma_Z, sigma_Y)
  
  # Selection probability
  p_select <- calc_selection_prob(X_t, sigma_Z)
  
  # Subgroup sample sizes
  n_sub <- n * p_select
  n_T <- n_sub * p_T
  n_C <- n_sub * (1 - p_T)
  
  # Standard error
  SE <- sigma_Y * sqrt((1/n_T) + (1/n_C))
  
  # Critical value
  z_crit <- qnorm(1 - alpha/2)
  
  # Power calculation
  pnorm(abs(delta)/SE - z_crit)
}

# Calculate power with known control response
calc_power_known_control <- function(d_overall, rho_XY_control, rho_XY_treatment, X_t, 
                                   sigma_Z, sigma_Y, n, p_T, alpha = 0.05) {
  
  delta <- calc_effect_size(d_overall, rho_XY_control, rho_XY_treatment, X_t, sigma_Z, sigma_Y)
  d <- delta / sigma_Y
  
  # Selection probability
  p_select <- calc_selection_prob(X_t, sigma_Z)
  
  # Sample size in treatment group
  n_T <- n * p_select * p_T
  
  # Power calculation for two-sided test
  z_crit <- qnorm(1 - alpha/2)
  pnorm(d * sqrt(n_T) - z_crit)
}

# Calculate true responder rate based on effect size
calc_true_responder_rate <- function(d_overall, rho_XY_control, rho_XY_treatment, X_t, sigma_Z = 1, sigma_Y = 1) {
  # Calculate effect size in the selected subpopulation
  delta <- calc_effect_size(d_overall, rho_XY_control, rho_XY_treatment, X_t)
  
  # Assuming a responder is defined as someone with a clinically significant improvement
  # (e.g., 50% reduction in symptoms or 1 SD improvement)
  clinical_threshold <- 1.0  # 1 SD improvement
  
  # Calculate probability of response
  pnorm(delta - clinical_threshold)
}

# Calculate economic metrics
calc_economic_metrics <- function(d_overall, rho_XY_control, rho_XY_treatment, X_t, 
                                 total_market, cost_per_test, value_per_responder, cost_per_treatment) {
  # Calculate selection probability
  p_select <- calc_selection_prob(X_t)
  
  # Calculate responder rate in selected population
  responder_rate <- calc_true_responder_rate(d_overall, rho_XY_control, rho_XY_treatment, X_t)
  
  # Number of patients selected
  n_selected <- total_market * p_select
  
  # Number of true responders
  n_responders <- n_selected * responder_rate
  
  # Testing costs (applied to all patients)
  testing_cost <- total_market * cost_per_test
  
  # Treatment costs (applied only to selected patients)
  treatment_cost <- n_selected * cost_per_treatment
  
  # Value generated from responders
  value_generated <- n_responders * value_per_responder
  
  # Net benefit
  net_benefit <- value_generated - testing_cost - treatment_cost
  
  # Return on investment (ROI)
  total_cost <- testing_cost + treatment_cost
  roi <- (value_generated - total_cost) / total_cost
  
  # Cost per responder
  cost_per_responder <- total_cost / n_responders
  
  list(
    selection_prob = p_select,
    responder_rate = responder_rate,
    n_selected = n_selected,
    n_responders = n_responders,
    testing_cost = testing_cost,
    treatment_cost = treatment_cost,
    value_generated = value_generated,
    net_benefit = net_benefit,
    roi = roi,
    cost_per_responder = cost_per_responder
  )
}

# UI Definition
ui <- dashboardPage(
  dashboardHeader(title = "MADRS Biomarker Power Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Power Analysis", tabName = "power", icon = icon("chart-line")),
      menuItem("Economic Analysis", tabName = "economic", icon = icon("dollar-sign")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Power Analysis Tab
      tabItem(tabName = "power",
              fluidRow(
                box(
                  title = "Statistical Parameters", width = 4, status = "primary",
                  sliderInput("d_overall", "Overall Treatment Effect (Cohen's d):", 
                              min = -1, max = 1, value = 0.3, step = 0.05),
                  
                  sliderInput("rho_XY_treatment", "Biomarker-Outcome Correlation (ρ_XY) in Treatment:", 
                              min = -0.5, max = 0.5, value = 0.25, step = 0.05),
                  helpText("This represents how well the biomarker predicts the outcome in the treatment group."),
                  
                  sliderInput("rho_XY_control", "Biomarker-Outcome Correlation (ρ_XY) in Control:", 
                              min = -0.5, max = 0.5, value = 0.0, step = 0.05),
                  helpText("The difference between treatment and control correlations indicates biomarker's predictive value."),
                  
                  sliderInput("n", "Total Sample Size (n):", 
                            min = 10, max = 1000, 
                            value = 200, 
                            step = 10),
                  sliderInput("p_T", "Proportion in Treatment Group:", min = 0.1, max = 0.9, value = 0.5, step = 0.05),
                  sliderInput("alpha", "Significance Level (α):", min = 0.01, max = 0.1, value = 0.05, step = 0.01)
                ),
                
                tabBox(
                  title = "Power Analysis Results", width = 8,
                  tabPanel("Power vs Threshold", 
                           plotlyOutput("power_threshold_plot"),
                           checkboxInput("show_known_control", "Show Power with Known Control", TRUE)),
                  tabPanel("Power vs Sample Size", 
                           plotlyOutput("power_sample_plot"),
                           sliderInput("X_t", "Biomarker Threshold", min = -3, max = 3, value = 0, step = 0.1)),
                  tabPanel("Effect Size", 
                           plotlyOutput("effect_size_plot"))
                )
              ),
              
              # fluidRow(
              #   box(
              #     title = "Threshold Analysis", width = 12, status = "info",
              #     sliderInput("X_t", "Biomarker Threshold (X_t):", min = -3, max = 3, value = 0, step = 0.1),
              #     fluidRow(
              #       valueBoxOutput("selection_box", width = 3),
              #       valueBoxOutput("effect_size_box", width = 3),
              #       valueBoxOutput("cohens_d_box", width = 3),
              #       valueBoxOutput("power_box", width = 3)
              #     )
              #   )
              # )
      ),
      
      # Economic Analysis Tab
      tabItem(tabName = "economic",
              fluidRow(
                box(
                  title = "Economic Parameters", width = 4, status = "primary",
                  numericInput("total_market", "Total Patient Population:", 100000, min = 1000, max = 10000000, step = 10000),
                  numericInput("cost_per_test", "Cost per Biomarker Test ($):", 2000, min = 10, max = 10000, step = 50),
                  numericInput("value_per_responder", "Value per Treatment Responder ($):", 50000, min = 1000, max = 100000, step = 1000),
                  numericInput("cost_per_treatment", "Cost per Treatment ($):", 10000, min = 100, max = 50000, step = 500),
                  helpText("These parameters define the economic model for biomarker-guided treatment."),
                  helpText("Adjust them to reflect your specific scenario.")
                ),
                box(
                  title = "Cost-Benefit Analysis", width = 8,
                  plotlyOutput("cost_benefit_plot", height = "400px")
                )
              ),
              fluidRow(
                box(
                  title = "Economic Summary", width = 12, status = "info",
                  tableOutput("economic_summary_table")
                )
              )
      ),
      
      # About Tab
      tabItem(tabName = "about",
              fluidRow(
                box(
                  title = "About This App", width = 12,
                  p("This application helps design and analyze biomarker-based clinical trials in depression."),
                  p("Key Features:"),
                  tags$ul(
                    tags$li("All measurements are standardized (mean 0, variance 1) for simplicity"),
                    tags$li("Effect sizes are reported as Cohen's d"),
                    tags$li("Power calculations account for reduced sample size in selected populations")
                  ),
                  
                  h4("Economic Analysis Explanation:"),
                  p("The economic analysis tab helps evaluate the financial implications of using a biomarker to select patients for treatment."),
                  
                  h5("Key Economic Metrics:"),
                  tags$ul(
                    tags$li(strong("Testing Cost"), " = Total Patient Population × Cost per Test"),
                    tags$li("This is the cost of testing all patients in the population, regardless of threshold"),
                    
                    tags$li(strong("Treatment Cost"), " = Number of Selected Patients × Cost per Treatment"),
                    tags$li("This is the cost of treating only those patients who meet the biomarker threshold"),
                    
                    tags$li(strong("Value Generated"), " = Number of Responders × Value per Responder"),
                    tags$li("This represents the economic value created by patients who respond to treatment"),
                    
                    tags$li(strong("Net Benefit"), " = Value Generated − Testing Cost − Treatment Cost"),
                    tags$li("This is the overall economic benefit after accounting for all costs")
                  ),
                  
                  h5("How Threshold Affects Economics:"),
                  tags$ul(
                    tags$li("Higher thresholds select fewer patients, reducing treatment costs but potentially missing responders"),
                    tags$li("Lower thresholds treat more patients, increasing costs but potentially capturing more responders"),
                    tags$li("The optimal threshold maximizes net benefit by balancing these tradeoffs")
                  ),
                  
                  h5("Mathematical Formulas:"),
                  tags$ul(
                    tags$li("Selection Probability: P(X > X_t) = 1 - Φ(X_t) where Φ is the standard normal CDF"),
                    tags$li("Number of Selected Patients = Total Population × Selection Probability"),
                    tags$li("Responder Rate = P(Response | X > X_t) based on effect size in selected population"),
                    tags$li("Number of Responders = Number of Selected Patients × Responder Rate")
                  )
                )
              )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Calculate rho_XY for treatment group based on control correlation and beta_3
  rho_XY_treatment <- reactive({
    # First calculate beta_2 from rho_XY_control
    beta_2 <- input$rho_XY_control * input$sigma_Y * sqrt(1 + input$sigma_Z^2)
    
    # Then calculate rho_XY_treatment
    (beta_2 + input$beta_3) / (input$sigma_Y * sqrt(1 + input$sigma_Z^2))
  })
  
  # Calculate key metrics based on current inputs
  current_metrics <- reactive({
    list(
      selection_prob = calc_selection_prob(input$X_t, input$sigma_Z),
      expected_Z = calc_expected_Z(input$X_t, input$sigma_Z),
      effect_size = calc_effect_size(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, input$X_t, input$sigma_Z, input$sigma_Y),
      cohens_d = calc_cohens_d(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, input$X_t, input$sigma_Z, input$sigma_Y),
      power_two_groups = calc_power_two_groups(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, input$X_t, 
                                            input$sigma_Z, input$sigma_Y, 
                                            input$n, input$p_T, input$alpha),
      power_known_control = calc_power_known_control(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, input$X_t, 
                                                  input$sigma_Z, input$sigma_Y, 
                                                  input$n, input$p_T, input$alpha)
    )
  })
  
  # Economic metrics
  economic_metrics <- reactive({
    calc_economic_metrics(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, input$X_t,
                          input$total_market, input$cost_per_test, input$value_per_responder, input$cost_per_treatment)
  })
  
  # Value boxes for threshold analysis
  output$selection_box <- renderValueBox({
    valueBox(
      paste0(round(current_metrics()$selection_prob * 100, 1), "%"),
      "Population Selected",
      icon = icon("users"),
      color = "blue"
    )
  })
  
  output$effect_size_box <- renderValueBox({
    valueBox(
      paste0(round(current_metrics()$effect_size, 1), " points"),
      "MADRS Reduction",
      icon = icon("chart-bar"),
      color = "green"
    )
  })
  
  output$cohens_d_box <- renderValueBox({
    valueBox(
      round(current_metrics()$cohens_d, 2),
      "Cohen's d",
      icon = icon("balance-scale"),
      color = "purple"
    )
  })
  
  output$power_box <- renderValueBox({
    valueBox(
      paste0(round(current_metrics()$power_two_groups * 100, 1), "%"),
      "Statistical Power",
      icon = icon("bolt"),
      color = "red"
    )
  })
  
  # Cost-Benefit Plot
  output$cost_benefit_plot <- renderPlotly({
    thresholds <- seq(-3, 3, by = 0.1)
    
    economic_data <- data.frame(
      Threshold = thresholds,
      Testing_Cost = sapply(thresholds, function(x) {
        calc_economic_metrics(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                             x, input$total_market, input$cost_per_test, 
                             input$value_per_responder, input$cost_per_treatment)$testing_cost
      }),
      Treatment_Cost = sapply(thresholds, function(x) {
        calc_economic_metrics(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                             x, input$total_market, input$cost_per_test, 
                             input$value_per_responder, input$cost_per_treatment)$treatment_cost
      }),
      Value_Generated = sapply(thresholds, function(x) {
        calc_economic_metrics(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                             x, input$total_market, input$cost_per_test, 
                             input$value_per_responder, input$cost_per_treatment)$value_generated
      }),
      Net_Benefit = sapply(thresholds, function(x) {
        calc_economic_metrics(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                             x, input$total_market, input$cost_per_test, 
                             input$value_per_responder, input$cost_per_treatment)$net_benefit
      })
    )
    
    # Convert values to millions for plotting
    economic_data_millions <- economic_data
    economic_data_millions$Testing_Cost <- economic_data_millions$Testing_Cost / 1000000
    economic_data_millions$Treatment_Cost <- economic_data_millions$Treatment_Cost / 1000000
    economic_data_millions$Value_Generated <- economic_data_millions$Value_Generated / 1000000
    economic_data_millions$Net_Benefit <- economic_data_millions$Net_Benefit / 1000000
    
    plot_data <- economic_data_millions %>%
      pivot_longer(cols = c(Testing_Cost, Treatment_Cost, Value_Generated, Net_Benefit),
                   names_to = "Metric", values_to = "Value")
    
    # Rename metrics for better display
    plot_data$Metric <- gsub("_", " ", plot_data$Metric)
    
    # Create the plot directly with plotly instead of ggplot
    plot_ly(plot_data, x = ~Threshold, y = ~Value, color = ~Metric, type = "scatter", mode = "lines",
            line = list(width = 4),
            hovertemplate = "%{fullData.name}: $%{y:.1f} Million<extra></extra>") %>%
      add_segments(x = input$X_t, xend = input$X_t, y = min(plot_data$Value), yend = max(plot_data$Value),
                   line = list(color = "red", dash = "dash", width = 2.5), showlegend = FALSE, # Thicker dash line
                   hovertemplate = "Threshold: %{x:.2f}<extra></extra>") %>%
      add_segments(x = min(thresholds), xend = max(thresholds), y = 0, yend = 0,
                   line = list(color = "darkgrey", dash = "dash", width = 2.5), showlegend = FALSE, # Thicker zero line
                   hoverinfo = "none") %>%
      layout(title = "Economic Analysis vs Biomarker Threshold",
             xaxis = list(title = "Biomarker Threshold (X_t)"),
             yaxis = list(title = "Value ($ millions)"),
             hovermode = "closest")
  })
  
  # Economic Summary Table
  output$economic_summary_table <- renderTable({
    # Calculate metrics at current threshold
    econ <- calc_economic_metrics(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                                input$X_t, input$total_market, input$cost_per_test, 
                                input$value_per_responder, input$cost_per_treatment)
    
    # Create summary table
    data.frame(
      Metric = c("Biomarker Threshold", 
                "Patients Selected", 
                "Selection Rate", 
                "Responder Rate", 
                "Number of Responders",
                "Testing Cost", 
                "Treatment Cost", 
                "Total Cost",
                "Value Generated", 
                "Net Benefit", 
                "Return on Investment (ROI)",
                "Cost per Responder"),
      Value = c(
        sprintf("%.2f", input$X_t),
        sprintf("%d", round(econ$n_selected)),
        sprintf("%.1f%%", 100 * econ$selection_prob),
        sprintf("%.1f%%", 100 * econ$responder_rate),
        sprintf("%d", round(econ$n_responders)),
        sprintf("$%s", format(round(econ$testing_cost), big.mark=",")),
        sprintf("$%s", format(round(econ$treatment_cost), big.mark=",")),
        sprintf("$%s", format(round(econ$testing_cost + econ$treatment_cost), big.mark=",")),
        sprintf("$%s", format(round(econ$value_generated), big.mark=",")),
        sprintf("$%s", format(round(econ$net_benefit), big.mark=",")),
        sprintf("%.1f%%", 100 * econ$roi),
        sprintf("$%s", format(round(econ$cost_per_responder), big.mark=","))
      )
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = 'l')
  
  # Power vs Threshold Plot
  output$power_threshold_plot <- renderPlotly({
    thresholds <- seq(-3, 3, by = 0.1)
    
    power_data <- data.frame(
      Threshold = thresholds,
      Power_Two_Groups = sapply(thresholds, function(x) {
        calc_power_two_groups(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                            x, 1, 1, # Using standardized values
                            input$n, input$p_T, input$alpha)
      })
    )
    
    if (input$show_known_control) {
      power_data$Power_Known_Control <- sapply(thresholds, function(x) {
        calc_power_known_control(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                               x, 1, 1, # Using standardized values
                               input$n, input$p_T, input$alpha)
      })
      
      plot_data <- power_data %>%
        pivot_longer(cols = c(Power_Two_Groups, Power_Known_Control),
                     names_to = "Method", values_to = "Power")
    } else {
      plot_data <- power_data %>%
        select(Threshold, Power = Power_Two_Groups)
    }
    
    p <- ggplot(plot_data) +
      geom_line(aes(x = Threshold, y = Power, color = if(exists("Method", plot_data)) Method else "Power"), size = 1) +
      geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgrey") +
      geom_vline(xintercept = input$X_t, linetype = "dashed", color = "red") +
      scale_y_continuous(limits = c(0, 1)) +
      labs(title = "Power vs Biomarker Threshold",
           x = "Biomarker Threshold (X_t)",
           y = "Statistical Power") +
      theme_minimal()
    
    ggplotly(p)
  })
  
  # Power vs Sample Size Plot
  output$power_sample_plot <- renderPlotly({
    # Generate custom sample sizes with logarithmic spacing
    sample_sizes <- c(
      # 10-100 by 10s
      seq(10, 100, by = 10),
      # 100-1000 by 100s
      seq(100, 1000, by = 100)[-1],
      # 1000-10000 by 1000s
      seq(1000, 10000, by = 1000)[-1]
    )
    
    power_data <- data.frame(
      Sample_Size = sample_sizes,
      Power_Two_Groups = sapply(sample_sizes, function(n) {
        calc_power_two_groups(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                             input$X_t, 1, 1, # Using standardized values
                             n, input$p_T, input$alpha)
      })
    )
    
    if (input$show_known_control) {
      power_data$Power_Known_Control <- sapply(sample_sizes, function(n) {
        calc_power_known_control(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, 
                                input$X_t, 1, 1, # Using standardized values
                                n, input$p_T, input$alpha)
      })
    }
    
    # Create plot
    p <- plot_ly() %>%
      add_trace(
        data = power_data,
        x = ~Sample_Size, 
        y = ~Power_Two_Groups,
        type = 'scatter', 
        mode = 'lines',
        line = list(color = 'blue', width = 3.5), # Increased line thickness
        name = "Two-Group Comparison"
      ) 
    
    if (input$show_known_control) {
      p <- p %>% add_trace(
        data = power_data,
        x = ~Sample_Size, 
        y = ~Power_Known_Control,
        type = 'scatter', 
        mode = 'lines',
        line = list(color = 'green', width = 3.5, dash = 'dash'), # Increased line thickness
        name = "Known Control"
      )
    }
    
    # Add reference line for target power
    p <- p %>% add_segments(
      x = 10, 
      xend = 10000, 
      y = 0.8, 
      yend = 0.8, 
      line = list(color = 'red', dash = 'dash', width = 2.5), # Increased line thickness
      name = "80% Power"
    )
    
    # Set layout with logarithmic x-axis
    p <- p %>% layout(
      title = "Power vs Sample Size",
      xaxis = list(
        title = "Sample Size (n)", 
        type = "log", 
        tickvals = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000),
        ticktext = c("10", "20", "50", "100", "200", "500", "1K", "2K", "5K", "10K")
      ),
      yaxis = list(title = "Power", range = c(0, 1.1)),
      legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.02), # Lowered legend position
      showlegend = TRUE
    )
    
    p
  })
  
  # Effect Size Plot
  output$effect_size_plot <- renderPlotly({
    thresholds <- seq(-3, 3, by = 0.1)
    
    effect_data <- data.frame(
      Threshold = thresholds,
      Effect_Size = sapply(thresholds, function(x) {
        calc_effect_size(input$d_overall, input$rho_XY_control, input$rho_XY_treatment, x, 1, 1)
      })
    )
    
    # Calculate standard errors for confidence bounds
    effect_data$SE <- sapply(thresholds, function(x) {
      # Selection probability
      p_select <- calc_selection_prob(x)
      
      # Subgroup sample sizes
      n_sub <- input$n * p_select
      n_T <- n_sub * input$p_T
      n_C <- n_sub * (1 - input$p_T)
      
      # Only calculate SE if we have enough subjects in each group
      if(n_T < 2 || n_C < 2) {
        return(NA)
      }
      
      # Standard error of effect size (approximation)
      sqrt((1/n_T + 1/n_C) * (1 + effect_data$Effect_Size[which(thresholds == x)]^2/8))
    })
    
    # Add upper and lower bounds
    effect_data$Upper <- effect_data$Effect_Size + 1.96 * effect_data$SE
    effect_data$Lower <- effect_data$Effect_Size - 1.96 * effect_data$SE
    
    # Create the plot with standard confidence interval bands
    p <- plot_ly() %>%
      # Add CI shading
      add_trace(
        data = effect_data,
        x = c(effect_data$Threshold, rev(effect_data$Threshold)),
        y = c(effect_data$Lower, rev(effect_data$Upper)),
        fill = 'toself',
        fillcolor = 'rgba(0, 0, 255, 0.2)',
        line = list(color = 'transparent'),
        name = "95% CI",
        hoverinfo = "none",
        showlegend = TRUE
      ) %>%
      # Add effect size line on top (with no markers/dots)
      add_trace(
        data = effect_data,
        x = ~Threshold, 
        y = ~Effect_Size,
        type = 'scatter', 
        mode = 'lines',
        line = list(color = 'blue', width = 2),
        name = "Effect Size",
        hoverinfo = "skip" # Removes hover dots
      ) %>%
      # Add threshold line
      add_segments(
        x = input$X_t, 
        xend = input$X_t, 
        y = -0.5, 
        yend = 1.5, 
        line = list(color = 'red', dash = 'dash'),
        name = "Current Threshold"
      ) %>%
      layout(
        title = "Effect Size vs Biomarker Threshold",
        xaxis = list(title = "Biomarker Threshold (X_t)"),
        yaxis = list(title = "Effect Size (Cohen's d)", range = c(-0.5, 1.5)),
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.1),
        showlegend = TRUE,
        hovermode = FALSE # Disable hover mode entirely
      )
    
    p
  })
}

# Run the application
shinyApp(ui = ui, server = server)
