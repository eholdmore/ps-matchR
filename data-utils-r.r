# utils.R

#' Load clinical data from genieBPC package
#' 
#' @param cohort The BPC cohort name
#' @param version The BPC data version
#' @return A data frame with clinical patient data
load_clinical_data <- function(cohort, version) {
  # Check if parameters are valid
  if (is.null(cohort) || is.null(version)) {
    return(NULL)
  }
  
  # Load clinical data using genieBPC
  tryCatch({
    clinical_data <- genieBPC::pull_data_clinical(
      cohort = cohort, 
      version = version
    )
    return(clinical_data)
  }, error = function(e) {
    warning(paste("Error loading clinical data:", e$message))
    return(NULL)
  })
}

#' Load treatment data from genieBPC package
#' 
#' @param cohort The BPC cohort name
#' @param version The BPC data version
#' @return A data frame with treatment data
load_treatment_data <- function(cohort, version) {
  # Check if parameters are valid
  if (is.null(cohort) || is.null(version)) {
    return(NULL)
  }
  
  # Load treatment data using genieBPC
  tryCatch({
    treatment_data <- genieBPC::pull_data_medications(
      cohort = cohort, 
      version = version
    )
    
    # In GENIE BPC, medications data is equivalent to the treatment data in your previous app
    # We'll make sure the column names match what your app expects
    if ("medication" %in% colnames(treatment_data)) {
      treatment_data$AGENT <- treatment_data$medication
    }
    
    return(treatment_data)
  }, error = function(e) {
    warning(paste("Error loading treatment data:", e$message))
    return(NULL)
  })
}

#' Load sample data from genieBPC package
#' 
#' @param cohort The BPC cohort name
#' @param version The BPC data version
#' @return A data frame with sample data
load_sample_data <- function(cohort, version) {
  # Check if parameters are valid
  if (is.null(cohort) || is.null(version)) {
    return(NULL)
  }
  
  # Load sample data using genieBPC
  tryCatch({
    # For GENIE BPC, we'll use the biomarker data which contains sample information
    sample_data <- genieBPC::pull_data_genomic(
      cohort = cohort, 
      version = version
    )
    
    # Create a simple sample data frame with PATIENT_ID
    if (!is.null(sample_data) && "record_id" %in% colnames(sample_data)) {
      # Map patient IDs to appropriate column name
      samples <- unique(sample_data[c("record_id", "sample_id")])
      names(samples)[names(samples) == "record_id"] <- "PATIENT_ID"
      names(samples)[names(samples) == "sample_id"] <- "SAMPLE_ID"
      
      return(samples)
    }
    
    # If we don't have genomic data, construct sample data from clinical data
    clinical_data <- load_clinical_data(cohort, version)
    if (!is.null(clinical_data) && "record_id" %in% colnames(clinical_data)) {
      # Use record_id as patient_id and create sample_id
      samples <- data.frame(
        PATIENT_ID = clinical_data$record_id,
        SAMPLE_ID = paste0("S_", clinical_data$record_id)
      )
      return(samples)
    }
    
    return(NULL)
  }, error = function(e) {
    warning(paste("Error loading sample data:", e$message))
    return(NULL)
  })
}

#' Preprocess data for propensity score matching
#' 
#' @param samples_df Sample data frame
#' @param treatment_df Treatment data frame
#' @param clinical_df Clinical data frame
#' @return A processed data frame ready for matching
preprocess_data <- function(samples_df, treatment_df, clinical_df) {
  # Standardize column names to match the expected format
  required_columns <- c("PATIENT_ID", "SAMPLE_ID")
  clinical_required <- c("PATIENT_ID")
  
  # Map common GENIE BPC column names to your expected format
  column_mapping <- list(
    "record_id" = "PATIENT_ID",
    "sample_id" = "SAMPLE_ID",
    "medication" = "AGENT",
    "sex" = "GENDER",
    "age_at_seq_rptd" = "CURRENT_AGE_DEID",
    "stage" = "STAGE_HIGHEST_RECORDED",
    "os_months" = "OS_MONTHS"
  )
  
  # Apply column mapping to each dataframe
  for (old_col in names(column_mapping)) {
    new_col <- column_mapping[[old_col]]
    
    if (old_col %in% colnames(samples_df)) {
      samples_df[[new_col]] <- samples_df[[old_col]]
    }
    
    if (old_col %in% colnames(treatment_df)) {
      treatment_df[[new_col]] <- treatment_df[[old_col]]
    }
    
    if (old_col %in% colnames(clinical_df)) {
      clinical_df[[new_col]] <- clinical_df[[old_col]]
    }
  }
  
  # Get unique agents from the treatment data
  if ('AGENT' %in% colnames(treatment_df)) {
    unique_agents <- unique(treatment_df$AGENT[!is.na(treatment_df$AGENT)])
  } else {
    unique_agents <- c()
  }
  
  # Merge clinical data with sample data
  merged_df <- clinical_df %>%
    dplyr::inner_join(samples_df, by = "PATIENT_ID")
  
  # Create a dictionary to hold patient treatment information
  patient_treatments <- treatment_df %>%
    filter(!is.na(AGENT)) %>%
    group_by(PATIENT_ID) %>%
    summarize(TREATMENTS = paste(unique(AGENT), collapse = ","))
  
  # Add treatment info to the merged data
  final_df <- merged_df %>%
    left_join(patient_treatments, by = "PATIENT_ID") %>%
    mutate(TREATMENTS = ifelse(is.na(TREATMENTS), "No Treatment", TREATMENTS))
  
  # Create treatment indicator columns
  for (agent in unique_agents) {
    agent_col <- paste0("RECEIVED_", gsub("[^[:alnum:]]", "_", agent))
    final_df[[agent_col]] <- ifelse(grepl(agent, final_df$TREATMENTS), 1, 0)
  }
  
  # Convert categorical variables to numeric for modeling
  for (col in c('GENDER', 'RACE', 'ETHNICITY', 'STAGE_HIGHEST_RECORDED')) {
    if (col %in% colnames(final_df)) {
      final_df[[col]] <- as.factor(final_df[[col]])
      final_df[[paste0(col, "_CODE")]] <- as.numeric(final_df[[col]])
    }
  }
  
  return(final_df)
}

#' Calculate standardized mean difference between two groups
#' 
#' @param x1 Vector of values from first group
#' @param x2 Vector of values from second group
#' @return The standardized mean difference
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

#' Calculate propensity scores and perform matching
#' 
#' @param treated_df Data frame for treated (query) group
#' @param control_df Data frame for control group
#' @param covariates Vector of covariate names
#' @param numerical_covariates Vector of numerical covariate names
#' @param caliper Caliper width for matching
#' @return A list containing matched cohorts and model information
calculate_and_match_ps <- function(treated_df, control_df, covariates, numerical_covariates, caliper) {
  # Combine datasets for modeling
  treated_df$TREATMENT <- 1
  control_df$TREATMENT <- 0
  combined_df <- rbind(treated_df, control_df)
  
  # Create formula for matching
  ps_formula <- as.formula(paste("TREATMENT ~", paste(covariates, collapse = " + ")))
  
  # Handle missing data
  combined_df <- combined_df %>% tidyr::drop_na(all_of(covariates))
  
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
  m.out <- MatchIt::matchit(
    ps_formula,
    data = combined_df,
    method = "nearest",
    distance = "glm",
    caliper = caliper,
    std.caliper = TRUE
  )
  
  # Get matched data
  matched_data <- MatchIt::match.data(m.out)
  
  # Extract propensity scores
  combined_df$PROPENSITY_SCORE <- m.out$distance
  
  # Split back into treated and control groups
  treated_with_ps <- combined_df %>% filter(TREATMENT == 1)
  control_with_ps <- combined_df %>% filter(TREATMENT == 0)
  
  # Extract matched treated and control
  if (nrow(matched_data) > 0) {
    matched_treated <- matched_data %>% filter(TREATMENT == 1)
    matched_control <- matched_data %>% filter(TREATMENT == 0)
    
    # Calculate match distance
    matched_treated$MATCH_DISTANCE <- matched_data %>% 
      group_by(subclass) %>% 
      summarize(dist = abs(diff(PROPENSITY_SCORE))) %>%
      pull(dist)
    
    matched_control$MATCH_DISTANCE <- matched_treated$MATCH_DISTANCE
  } else {
    matched_treated <- NULL
    matched_control <- NULL
  }
  
  return(list(
    treated_with_ps = treated_with_ps,
    control_with_ps = control_with_ps,
    matched_treated = matched_treated,
    matched_control = matched_control,
    ps_model = m.out
  ))
}
