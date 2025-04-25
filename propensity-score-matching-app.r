# app.R

# Load required packages
library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(genieBPC)
library(MatchIt)
library(tableone)
library(scales)

# UI definition
ui <- fluidPage(
  titlePanel("Propensity Score Matching for Treatment Cohorts"),
  
  # Setup and data loading section
  sidebarLayout(
    sidebarPanel(
      h3("Data Selection"),
      
      # GENIE BPC Data Connection
      selectInput("consortium", 
                  "Select GENIE BPC Consortium:", 
                  choices = c("NSCLC", "CRC", "BrCa", "Prostate", "Pancreas")),
      
      # Data version selection
      selectInput("version", 
                  "Select Data Version:", 
                  choices = NULL),
      
      actionButton("loadData", "Load Data", class = "btn-primary"),
      hr(),
      
      # Only show these options after data is loaded
      conditionalPanel(
        condition = "output.dataLoaded",
        h3("Cohort Selection"),
        selectInput("treatmentAgent", "Select Treatment Agent:", choices = NULL),
        hr(),
        
        h3("Propensity Score Settings"),
        selectInput("numericalCovariates", "Select Numerical Covariates:", 
                    choices = NULL, multiple = TRUE),
        selectInput("categoricalCovariates", "Select Categorical Covariates:",
                    choices = NULL, multiple = TRUE),
        sliderInput("caliper", "Caliper Width:", min = 0.1, max = 1.0, value = 0.2, step = 0.1),
        actionButton("runMatching", "Run Propensity Score Matching", class = "btn-success")
      )
    ),
    
    mainPanel(
      # Initial data summary
      conditionalPanel(
        condition = "output.dataLoaded",
        h3("Dataset Summary"),
        fluidRow(
          column(4, valueBoxOutput("totalPatientsBox", width = NULL)),
          column(4, valueBoxOutput("totalSamplesBox", width = NULL)),
          column(4, valueBoxOutput("treatmentAgentsBox", width = NULL))
        ),
        
        # Debug: Check treatment data
        tabsetPanel(
          tabPanel("Treatment Data", 
                  h4("Sample of Treatment Data"),
                  DTOutput("treatmentDataTable"),
                  h4("Top Agents in Treatment Data"),
                  plotlyOutput("topAgentsPlot")
          ),
          
          # Results only shown after matching
          tabPanel("Matching Results",
                  conditionalPanel(
                    condition = "output.matchingComplete",
                    h3("Matching Summary"),
                    fluidRow(
                      column(4, valueBoxOutput("originalCohortBox", width = NULL)),
                      column(4, valueBoxOutput("matchedCohortBox", width = NULL)),
                      column(4, valueBoxOutput("matchRateBox", width = NULL))
                    ),
                    
                    h3("Propensity Score Distributions"),
                    plotlyOutput("psDistPlot"),
                    
                    h3("Covariate Balance Assessment"),
                    plotlyOutput("balancePlot"),
                    h4("Tabular Balance Assessment"),
                    DTOutput("balanceTable"),
                    
                    h3("Cohort Characteristics"),
                    DTOutput("characteristicsTable"),
                    
                    h3("Download Matched Cohorts"),
                    fluidRow(
                      column(6, downloadButton("downloadTreated", "Download Query Cohort CSV")),
                      column(6, downloadButton("downloadControl", "Download Matched Controls CSV"))
                    )
                  )
          )
        )
      ),
      
      # Show this when data is not yet loaded
      conditionalPanel(
        condition = "!output.dataLoaded",
        h3("Instructions"),
        p("Welcome to the GENIE BPC Propensity Score Matching application."),
        p("To begin:"),
        tags$ol(
          tags$li("Select a GENIE BPC cohort from the dropdown menu"),
          tags$li("Select a data version"),
          tags$li("Click 'Load Data' to retrieve the dataset"),
          tags$li("You will then be able to select treatment agents and covariates for matching")
        )
      )
    )
  )
)

# Server definition
server <- function(input, output, session) {
  # Reactive values to store our data and analysis results
  rv <- reactiveValues(
    samples_df = NULL,
    treatment_df = NULL,
    clinical_df = NULL,
    available_agents = NULL,
    processed_data = NULL,
    treated_with_ps = NULL,
    control_with_ps = NULL,
    matched_treated = NULL,
    matched_control = NULL,
    ps_model = NULL,
    dataLoaded = FALSE,
    matchingComplete = FALSE
  )
  
  # Flag for data loading status
  output$dataLoaded <- reactive({
    return(rv$dataLoaded)
  })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  
  # Flag for matching completion status
  output$matchingComplete <- reactive({
    return(rv$matchingComplete)
  })
  outputOptions(output, "matchingComplete", suspendWhenHidden = FALSE)
  
  # Update version choices based on consortium selection
  observe({
    # In a real implementation, we would query available versions
    # For now, we'll provide placeholder versions
    updateSelectInput(session, "version",
                     choices = c("v2.0-public", "v1.1-public", "v1.0-public"))
  })
  
  # Load data action
  observeEvent(input$loadData, {
    # Reset previous data
    rv$dataLoaded <- FALSE
    rv$matchingComplete <- FALSE
    
    # Show loading message
    showNotification("Loading GENIE BPC data...", type = "message", duration = NULL, id = "loading")
    
    # Try to load data - for demonstration we'll use tryCatch to handle potential errors
    tryCatch({
      # In a real implementation, use genieBPC functions to load data
      # For example:
      # consortium <- input$consortium
      # version <- input$version
      # clin_data <- pull_data_clinical(cohort = consortium, version = version)
      # treat_data <- pull_data_treatments(cohort = consortium, version = version)
      
      # For demonstration purposes, we'll create sample data that resembles the expected structure
      # (In your real app, replace this with actual genieBPC data loading)
      
      # Sample clinical data
      rv$clinical_df <- data.frame(
        PATIENT_ID = paste0("P", 1:500),
        GENDER = sample(c("Male", "Female"), 500, replace = TRUE),
        RACE = sample(c("White", "Black", "Asian", "Other"), 500, replace = TRUE),
        ETHNICITY = sample(c("Hispanic", "Non-Hispanic"), 500, replace = TRUE),
        CURRENT_AGE_DEID = round(rnorm(500, mean = 65, sd = 12)),
        STAGE_HIGHEST_RECORDED = sample(c("I", "II", "III", "IV"), 500, replace = TRUE),
        stringsAsFactors = FALSE
      )
      
      # Sample sample data (can be same as clinical for simplicity)
      rv$samples_df <- rv$clinical_df %>%
        select(PATIENT_ID) %>%
        mutate(SAMPLE_ID = paste0("S", 1:500),
               OS_MONTHS = round(rexp(500, 1/30), 1)) # Random survival times
      
      # Sample treatment data
      agents <- c("Carboplatin", "Cisplatin", "Paclitaxel", "Gemcitabine", 
                 "Pembrolizumab", "Nivolumab", "Atezolizumab", "Docetaxel", "Pemetrexed")
      
      rv$treatment_df <- data.frame(
        PATIENT_ID = sample(rv$clinical_df$PATIENT_ID, 800, replace = TRUE),
        AGENT = sample(agents, 800, replace = TRUE),
        EVENT_TYPE = "TREATMENT",
        stringsAsFactors = FALSE
      )
      
      # Process the loaded data
      rv$available_agents <- unique(rv$treatment_df$AGENT)
      
      # Update treatment agent selection
      updateSelectInput(session, "treatmentAgent",
                       choices = sort(rv$available_agents))
      
      # Update covariate selections
      updateSelectInput(session, "numericalCovariates",
                       choices = c("CURRENT_AGE_DEID", "OS_MONTHS"),
                       selected = c("CURRENT_AGE_DEID", "OS_MONTHS"))
      
      updateSelectInput(session, "categoricalCovariates",
                       choices = c("GENDER", "RACE", "ETHNICITY", "STAGE_HIGHEST_RECORDED"),
                       selected = c("GENDER", "RACE", "STAGE_HIGHEST_RECORDED"))
      
      # Flag data as loaded
      rv$dataLoaded <- TRUE
      
      # Process data for analysis
      rv$processed_data <- preprocess_data(rv$samples_df, rv$treatment_df, rv$clinical_df)
      
      # Remove loading notification
      removeNotification(id = "loading")
      showNotification("Data loaded successfully!", type = "message")
      
    }, error = function(e) {
      # Handle errors
      removeNotification(id = "loading")
      showNotification(paste("Error loading data:", e$message), type = "error", duration = NULL)
    })
  })
  
  # Data preprocessing function
  preprocess_data <- function(samples_df, treatment_df, clinical_df) {
    # Get all treatments for each patient
    patient_treatments <- treatment_df %>%
      filter(!is.na(AGENT)) %>%
      group_by(PATIENT_ID) %>%
      summarize(TREATMENTS = paste(unique(AGENT), collapse = ","))
    
    # Merge clinical and sample data
    merged_df <- clinical_df %>%
      inner_join(samples_df, by = "PATIENT_ID")
    
    # Add treatment info
    final_df <- merged_df %>%
      left_join(patient_treatments, by = "PATIENT_ID") %>%
      mutate(TREATMENTS = ifelse(is.na(TREATMENTS), "No Treatment", TREATMENTS))
    
    # Create treatment indicator columns
    for (agent in unique(treatment_df$AGENT)) {
      # Clean agent name for column
      agent_col <- paste0("RECEIVED_", gsub("[^[:alnum:]]", "_", agent))
      
      # Create indicator column
      final_df[[agent_col]] <- ifelse(grepl(agent, final_df$TREATMENTS), 1, 0)
    }
    
    return(final_df)
  }
  
  # Display total patients metric
  output$totalPatientsBox <- renderValueBox({
    req(rv$processed_data)
    valueBox(
      length(unique(rv$processed_data$PATIENT_ID)),
      "Total Patients",
      icon = icon("users"),
      color = "navy"
    )
  })
  
  # Display total samples metric
  output$totalSamplesBox <- renderValueBox({
    req(rv$processed_data)
    valueBox(
      nrow(rv$processed_data),
      "Total Samples", 
      icon = icon("vials"),
      color = "light-blue"
    )
  })
  
  # Display treatment agents metric
  output$treatmentAgentsBox <- renderValueBox({
    req(rv$available_agents)
    valueBox(
      length(rv$available_agents),
      "Treatment Agents",
      icon = icon("prescription-bottle-medical"),
      color = "aqua"
    )
  })
  
  # Display treatment data table
  output$treatmentDataTable <- renderDT({
    req(rv$treatment_df)
    datatable(
      head(rv$treatment_df, 10),
      options = list(scrollX = TRUE, dom = 'frtip')
    )
  })
  
  # Display top agents plot
  output$topAgentsPlot <- renderPlotly({
    req(rv$treatment_df)
    top_agents <- rv$treatment_df %>%
      filter(!is.na(AGENT)) %>%
      count(AGENT) %>%
      arrange(desc(n)) %>%
      head(10)
    
    plot_ly(top_agents, x = ~AGENT, y = ~n, type = "bar", marker = list(color = "#1E88E5")) %>%
      layout(title = "Top 10 Agents in Treatment Data",
             xaxis = list(title = "Agent"),
             yaxis = list(title = "Count"))
  })
  
  # Run propensity score matching
  observeEvent(input$runMatching, {
    req(rv$processed_data, input$treatmentAgent)
    
    # Reset previous matching results
    rv$matchingComplete <- FALSE
    
    # Check if covariates are selected
    selected_covariates <- c(input$numericalCovariates, input$categoricalCovariates)
    if (length(selected_covariates) == 0) {
      showNotification("Please select at least one covariate for matching.", type = "error")
      return()
    }
    
    # Show loading message
    showNotification("Calculating propensity scores and matching cohorts...", 
                    type = "message", duration = NULL, id = "matching")
    
    # Try to run the matching process
    tryCatch({
      # Get selected agent
      selected_agent <- input$treatmentAgent
      cleaned_agent <- gsub("[^[:alnum:]]", "_", selected_agent)
      selected_agent_col <- paste0("RECEIVED_", cleaned_agent)
      
      # Define query and control cohorts
      if (!selected_agent_col %in% colnames(rv$processed_data)) {
        removeNotification(id = "matching")
        showNotification(paste("Column", selected_agent_col, "not found in the data."), 
                        type = "error", duration = NULL)
        return()
      }
      
      query_cohort <- rv$processed_data %>% filter(.data[[selected_agent_col]] == 1)
      control_pool <- rv$processed_data %>% filter(.data[[selected_agent_col]] == 0)
      
      # Check cohort sizes
      if (nrow(query_cohort) == 0) {
        removeNotification(id = "matching")
        showNotification(paste("No patients found in the query cohort for", selected_agent), 
                        type = "warning", duration = NULL)
        return()
      }
      
      if (nrow(control_pool) == 0) {
        removeNotification(id = "matching")
        showNotification("No patients found in the control pool. Cannot perform matching.", 
                        type = "warning", duration = NULL)
        return()
      }
      
      # Calculate propensity scores and perform matching
      results <- calculate_and_match_ps(
        query_cohort, 
        control_pool, 
        selected_covariates, 
        input$numericalCovariates,
        input$caliper
      )
      
      # Store results
      rv$treated_with_ps <- results$treated_with_ps
      rv$control_with_ps <- results$control_with_ps
      rv$matched_treated <- results$matched_treated
      rv$matched_control <- results$matched_control
      rv$ps_model <- results$ps_model
      
      # Check if matching was successful
      if (is.null(rv$matched_treated) || nrow(rv$matched_treated) == 0) {
        removeNotification(id = "matching")
        showNotification("Matching failed. Please try adjusting the caliper or select different covariates.", 
                        type = "error", duration = NULL)
        return()
      }
      
      # Flag matching as complete
      rv$matchingComplete <- TRUE
      
      # Remove loading notification and show success message
      removeNotification(id = "matching")
      showNotification(paste("Successfully matched", nrow(rv$matched_treated), "patients!"), 
                      type = "message")
      
    }, error = function(e) {
      # Handle errors
      removeNotification(id = "matching")
      showNotification(paste("Error during matching:", e$message), type = "error", duration = NULL)
    })
  })
  
  # Function to calculate propensity scores and perform matching
  calculate_and_match_ps <- function(treated_df, control_df, covariates, numerical_covariates, caliper) {
    # Combine datasets for modeling
    treated_df$TREATMENT <- 1
    control_df$TREATMENT <- 0
    combined_df <- rbind(treated_df, control_df)
    
    # Create formula for matching
    ps_formula <- as.formula(paste("TREATMENT ~", paste(covariates, collapse = " + ")))
    
    # Handle missing data
    combined_df <- combined_df %>% drop_na(all_of(covariates))
    
    # Check if there are enough samples after dropping missing values
    if (sum(combined_df$TREATMENT == 1) == 0 || sum(combined_df$TREATMENT == 0) == 0) {
      return(list(
        treated_with_ps = NULL,
        control_with_ps = NULL,
        matched_treated = NULL,
        matched_control = NULL,
        ps_model = NULL
      ))
    }
    
    # Run matching using MatchIt
    m.out <- matchit(
      ps_formula,
      data = combined_df,
      method = "nearest",   # Change to nearest neighbor matching
      caliper = caliper,    # Optional: Caliper to limit the distance
      distance = "glm"
    )
    
    # Get matched data
    matched_data <- match.data(m.out)
    
    # Extract propensity scores (for later use if needed)
    combined_df$PROPENSITY_SCORE <- m.out$distance
    
    # Split back into treated and control groups (before matching)
    treated_with_ps <- combined_df %>% filter(TREATMENT == 1)
    control_with_ps <- combined_df %>% filter(TREATMENT == 0)
    
    # Extract matched treated and control (after matching)
    matched_treated <- matched_data %>% filter(TREATMENT == 1)
    matched_control <- matched_data %>% filter(TREATMENT == 0)
    
    return(list(
      treated_with_ps = treated_with_ps,
      control_with_ps = control_with_ps,
      matched_treated = matched_treated,
      matched_control = matched_control,
      ps_model = m.out
    ))
  }
  
  # Render matching summary boxes
  output$originalCohortBox <- renderValueBox({
    req(rv$treated_with_ps)
    valueBox(
      nrow(rv$treated_with_ps),
      "Original Query Cohort",
      icon = icon("users"),
      color = "navy"
    )
  })
  
  output$matchedCohortBox <- renderValueBox({
    req(rv$matched_treated)
    valueBox(
      nrow(rv$matched_treated),
      "Matched Query Patients",
      icon = icon("user-check"),
      color = "light-blue"
    )
  })
  
  output$matchRateBox <- renderValueBox({
    req(rv$matched_treated, rv$treated_with_ps)
    match_rate <- round(nrow(rv$matched_treated) / nrow(rv$treated_with_ps) * 100, 1)
    valueBox(
      paste0(match_rate, "%"),
      "Match Rate",
      icon = icon("percentage"),
      color = "aqua"
    )
  })
  
  safe_rbind <- function(...) {
    data_frames <- list(...)
    valid_data_frames <- data_frames[sapply(data_frames, nrow) > 0]
    
    if (length(valid_data_frames) > 0) {
      do.call(rbind, valid_data_frames)
    } else {
      # Return an empty data frame with expected columns (adjust as needed)
      data.frame(PS = numeric(0), Group = character(0))
    }
  }
  
  # Render propensity score distribution plot
  output$psDistPlot <- renderPlotly({
    req(rv$treated_with_ps, rv$control_with_ps, rv$matched_treated, rv$matched_control)
    
    rv$matched_treated$PROPENSITY_SCORE <- rv$matched_treated$distance
    rv$matched_control$PROPENSITY_SCORE <- rv$matched_control$distance
    
    # Define colors
    colors <- c(
      "Query Cohort (Before Matching)" = "blue",
      "Control Pool (Before Matching)" = "red",
      "Query Cohort (After Matching)" = "lightblue",
      "Matched Controls" = "lightcoral"
    )
    
    # Before matching plot
    p_before <- plot_ly() %>%
      add_histogram(
        data = rv$treated_with_ps,
        x = ~PROPENSITY_SCORE,
        name = "Query Cohort (Before Matching)",
        marker = list(color = colors["Query Cohort (Before Matching)"]),
        opacity = 0.7
      ) %>%
      add_histogram(
        data = rv$control_with_ps,
        x = ~PROPENSITY_SCORE,
        name = "Control Pool (Before Matching)",
        marker = list(color = colors["Control Pool (Before Matching)"]),
        opacity = 0.7
      ) %>%
      layout(
        title = "",  # Remove the title from individual plots
        xaxis = list(title = "Propensity Score"),
        yaxis = list(title = "Count"),
        barmode = "overlay",
        showlegend = TRUE
      )
    
    # After matching plot
    p_after <- plot_ly() %>%
      add_histogram(
        data = rv$matched_treated,
        x = ~PROPENSITY_SCORE,
        name = "Query Cohort (After Matching)",
        marker = list(color = colors["Query Cohort (After Matching)"]),
        opacity = 0.7
      ) %>%
      add_histogram(
        data = rv$matched_control,
        x = ~PROPENSITY_SCORE,
        name = "Matched Controls",
        marker = list(color = colors["Matched Controls"]),
        opacity = 0.7
      ) %>%
      layout(
        title = "",  # Remove the title from individual plots
        xaxis = list(title = "Propensity Score"),
        yaxis = list(title = "Count"),
        barmode = "overlay",
        showlegend = TRUE
      )
    
    # Combine using subplot and add annotations for titles at the top
    subplot(p_before, p_after, nrows = 1, shareY = TRUE, titleX = TRUE, titleY = TRUE) %>%
      layout(
        annotations = list(
          list(
            x = 0,  # Position of "Before Matching" title
            y = 1.05,  # Position slightly above the plot area
            text = "Before Matching",
            showarrow = FALSE,
            font = list(size = 14, color = "black"),
            align = "center",
            xref = "paper", 
            yref = "paper"   
          ),
          list(
            x = 0.75,  # Position of "After Matching" title
            y = 1.05,  # Position slightly above the plot area
            text = "After Matching",
            showarrow = FALSE,
            font = list(size = 14, color = "black"),
            align = "center",
            xref = "paper",  
            yref = "paper"   
          )
        )
      )
  })
  
  # Calculate standardized mean difference
  # for numeric variables 
  calculate_smd <- function(x1, x2) {
    if (length(x1) == 0 || length(x2) == 0) return(NA)
    
    mean1 <- mean(x1, na.rm = TRUE)
    mean2 <- mean(x2, na.rm = TRUE)
    var1 <- var(x1, na.rm = TRUE)
    var2 <- var(x2, na.rm = TRUE)
    
    pooled_sd <- sqrt((var1 + var2) / 2)
    
    if (pooled_sd == 0) return(0)
    return(abs(mean1 - mean2) / pooled_sd)
  }
  
  # and for categorical variables 
  categorical_smd <- function(x1, x2) {
    x1 <- as.factor(x1)
    x2 <- as.factor(x2)
    
    # Align factor levels
    all_levels <- union(levels(x1), levels(x2))
    x1 <- factor(x1, levels = all_levels)
    x2 <- factor(x2, levels = all_levels)
    
    p1 <- prop.table(table(x1))
    p2 <- prop.table(table(x2))
    
    # Ensure all levels are represented
    p1 <- p1[all_levels]; p1[is.na(p1)] <- 0
    p2 <- p2[all_levels]; p2[is.na(p2)] <- 0
    
    return(max(abs(p1 - p2)))
  }
  
  # Balance assessment data
  balance_data <- reactive({
    req(rv$treated_with_ps, rv$control_with_ps, rv$matched_treated, rv$matched_control)
    req(input$numericalCovariates, input$categoricalCovariates)
    
    covariates <- c(input$numericalCovariates, input$categoricalCovariates)
    
    balance_df <- data.frame(
      Covariate = character(),
      SMD_Before = numeric(),
      SMD_After = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (cov in covariates) {
      tryCatch({
        x1 <- rv$treated_with_ps[[cov]]
        x2 <- rv$control_with_ps[[cov]]
        y1 <- rv$matched_treated[[cov]]
        y2 <- rv$matched_control[[cov]]
        
        # Use numeric SMD if both vectors are numeric
        if (is.numeric(x1) && is.numeric(x2) && is.numeric(y1) && is.numeric(y2)) {
          smd_before <- calculate_smd(x1, x2)
          smd_after <- calculate_smd(y1, y2)
        } else {
          # For categorical: use prop-based SMD approximation (e.g., max abs diff)
          smd_before <- categorical_smd(x1, x2)
          smd_after <- categorical_smd(y1, y2)
        }
        
        balance_df <- rbind(balance_df, data.frame(
          Covariate = cov,
          SMD_Before = smd_before,
          SMD_After = smd_after
        ))
      }, error = function(e) {
        message("❌ Error for covariate: ", cov)
        message("   ", conditionMessage(e))
      })
    }
    
    return(balance_df)
  })
  
  # Render balance plot
  output$balancePlot <- renderPlotly({
    req(balance_data())
    
    balance_df <- balance_data()
    
    plot_ly() %>%
      add_trace(
        type = "bar",
        y = balance_df$Covariate,
        x = balance_df$SMD_Before,
        orientation = "h",
        name = "Before Matching",
        marker = list(color = "red")
      ) %>%
      add_trace(
        type = "bar",
        y = balance_df$Covariate,
        x = balance_df$SMD_After,
        orientation = "h",
        name = "After Matching",
        marker = list(color = "green")
      ) %>%
      layout(
        title = "Standardized Mean Differences Before and After Matching",
        xaxis = list(
          title = "Standardized Mean Difference",
          range = c(0, max(c(balance_df$SMD_Before, balance_df$SMD_After), na.rm = TRUE) * 1.1)
        ),
        yaxis = list(title = "Covariate"),
        barmode = "group",
        shapes = list(
          list(
            type = "line",
            x0 = 0.1, x1 = 0.1,
            y0 = -0.5, y1 = nrow(balance_df) - 0.5,
            line = list(color = "gray", width = 2, dash = "dash")
          )
        )
      )
  })
  
  # Render balance table
  output$balanceTable <- renderDT({
    req(balance_data())
    
    balance_df <- balance_data()
    
    datatable(
      balance_df,
      options = list(
        scrollX = TRUE,
        dom = 'frtip'
      ),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("SMD_Before", "SMD_After"), digits = 3) %>%
      formatStyle(
        "SMD_After",
        backgroundColor = styleInterval(0.1, c("white", "#ffcccc"))
      )
  })
  
  # Cohort characteristics comparison
  output$characteristicsTable <- renderDT({
    req(rv$matched_treated, rv$matched_control)
    req(input$numericalCovariates, input$categoricalCovariates)
    
    covariates <- c(input$numericalCovariates, input$categoricalCovariates)
    
    # Create TableOne object
    vars <- covariates
    factorVars <- input$categoricalCovariates
    
    # Combine data for comparison
    matched_data <- rbind(
      cbind(rv$matched_treated, Group = "Treated"),
      cbind(rv$matched_control, Group = "Control")
    )
    
    # Generate table one
    tab1 <- CreateTableOne(
      vars = vars,
      strata = "Group",
      data = matched_data,
      factorVars = factorVars,
      test = TRUE
    )
    
    # Convert to data frame
    tab_df <- as.data.frame(print(tab1, printToggle = FALSE, showAllLevels = TRUE))
    
    # Clean up and reshape
    tab_df$Variable <- rownames(tab_df)
    rownames(tab_df) <- NULL
    
    # Reorder columns for better presentation
    cols <- c("Variable", "Treated", "Control", "p")
    tab_df <- tab_df[, cols[cols %in% colnames(tab_df)]]
    
    datatable(
      tab_df,
      options = list(
        scrollX = TRUE,
        dom = 'frtip'
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        "p",
        backgroundColor = styleInterval(0.05, c("#ffcccc", "white"))
      )
  })
  
  # Download handlers
  output$downloadTreated <- downloadHandler(
    filename = function() {
      agent <- input$treatmentAgent
      cleaned_agent <- gsub("[^[:alnum:]]", "_", agent)
      paste0(cleaned_agent, "_matched_query_cohort.csv")
    },
    content = function(file) {
      write.csv(rv$matched_treated, file, row.names = FALSE)
    }
  )
  
  output$downloadControl <- downloadHandler(
    filename = function() {
      agent <- input$treatmentAgent
      cleaned_agent <- gsub("[^[:alnum:]]", "_", agent)
      paste0(cleaned_agent, "_matched_controls.csv")
    },
    content = function(file) {
      write.csv(rv$matched_control, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
