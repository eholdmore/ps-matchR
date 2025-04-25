# Propensity Score Matching App for GENIE BPC

This application allows researchers to perform propensity score matching on GENIE BPC data to create balanced cohorts for treatment effect analysis.

## 🌐 Live Demo

You can access the deployed Shiny app here:  
🔗 [ps-matchR on shinyapps.io](https://eholdmore.shinyapps.io/ps-matchR)

## Project Structure

```
propensity_score_matcher/
├── app.R           # Main Shiny application 
├── utils.R         # Utility functions for data processing
├── README.md       # This file
└── .gitignore      # Git ignore file
```

## Key Differences from the Streamlit Version

This R/Shiny version has several important differences from the original Python/Streamlit app:

1. **Data Source**: Uses the genieBPC R package instead of local MSK-CHORD data files
2. **Technology Stack**: R/Shiny instead of Python/Streamlit
3. **Library Usage**: 
   - MatchIt for propensity score matching (instead of scikit-learn)
   - plotly for interactive visualizations
   - tableone for cohort characteristics comparison
   - dplyr for data manipulation

## Setup Instructions

### Prerequisites

- R (version 4.0.0 or higher)
- RStudio (recommended)

### Installation

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/ps_matchR.git
   cd ps_matchR
   ```

2. Install required packages:
   ```R
   install.packages(c("shiny", "dplyr", "ggplot2", "plotly", "DT", "MatchIt", "tableone", "scales"))
   
   # Install genieBPC from GitHub
   if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
   remotes::install_github("genie-bpc/genieBPC")
   ```

### Running the Application

1. Open `app.R` in RStudio
2. Click the "Run App" button or run:
   ```R
   shiny::runApp()
   ```

## Using the Application

1. Select a GENIE BPC cohort and data version
2. Click "Load Data" to retrieve the dataset
3. Select a treatment agent to define your query cohort
4. Choose covariates for propensity score calculation
5. Adjust the caliper width as needed
6. Click "Run Propensity Score Matching" to perform matching
7. Review the matching results and download the matched cohorts

## Implementation Notes

### Data Loading

The app uses functions from the [genieBPC package](https://genie-bpc.github.io/genieBPC/) to load clinical, treatment, and genomic data. The `pull_data_clinical()`, `pull_data_medications()`, and `pull_data_genomic()` functions are used to fetch the required datasets.

### Propensity Score Calculation

Propensity scores are calculated using logistic regression via the [MatchIt package](https://cran.r-project.org/web/packages/MatchIt/vignettes/MatchIt.html). The app fits a model predicting treatment assignment based on selected covariates.

### Matching Process

The app uses nearest-neighbor matching with a user-specified caliper width. Patients from the control group are matched to patients in the treatment group based on their propensity scores.

### Balance Assessment

The app assesses covariate balance by calculating standardized mean differences (SMD) before and after matching. An SMD less than 0.1 indicates good balance for that covariate.

## Adapting to Your Needs

If you need to modify the app for your specific needs:

- **Change available consortia**: Update the `choices` parameter in the `selectInput("consortium", ...)` function
- **Add more covariates**: Modify the preprocessing function to extract additional clinical features
- **Customize visualizations**: Adjust the plotly graph objects in the appropriate render functions

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Original [Python/Streamlit implementation](https://github.com/eholdmore/propensity_score_matcher) by eholdmore
- GENIE BPC data provided by AACR Project GENIE
