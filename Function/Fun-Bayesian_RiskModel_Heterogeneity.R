
############################## Required Functions ##############################

### Stochastic Search Variable Selection
source("~/Thesis/Bayesian_RiskModel_TEH/Function/Helper_Script-Function/Fun-SSVS.R")

### Risk Calculations
source("~/Thesis/Bayesian_RiskModel_TEH/Function/Helper_Script-Function/Fun-Risk_Calculations.R")

### Evaluate Heterogeneity
source("~/Thesis/Bayesian_RiskModel_TEH/Function/Helper_Script-Function/Fun-Risk_Heterogeneity.R")



################## Full Function for Evaluating Heterogeneity ##################

Risk_Heterogeneity <- function(Trial_Data, Burn_in, Iterations, 
                               SSVS_script, Risk_script, Heterogeneity_script,
                               Optimized = TRUE) {
  
  ### Using Stochastic Search Variable Selection (SSVS) to select prognostic variables
  if(Optimized == TRUE) {
    SSVS_results <- SSVS_jags(Trial_Data, Burn_in, Iterations, model_string.SSVS = SSVS_script)
    included <- SSVS_results$included # Variables selected through SSVS
    ### Checking for included to not be empty
    if (any(included == "N/A")){
      summary <- list(Heterogeneity = list(Posterior_Means = "N/A",
                                           Intervals = "N/A",
                                           Heterogeneity = "N/A",
                                           R_hat = "N/A"),
                      Risk = list(Risk = list(Risk = "N/A",
                                              Posterior_Means = "N/A",
                                              R_hat = "N/A"),
                                  Optimized_Risk = if_else(Optimized == TRUE, 
                                                           "Optimized Risk Model",
                                                           "Full Risk Model")),
                      SSVS = "N/A",
                      Selection_Probability = "N/A",
                      Variable_Selections = "N/A",
                      Iterations = "N/A",
                      Samples = "N/A")
      return(summary)
      
    }
  } else {
    included = "0"
  }
  
  
  
  ### Calculating Risk based on Bayesian risk model
  Risk_results <- Risk_jags(Trial_Data, included, Burn_in, Iterations, model_string.Risk = Risk_script, Optimized = Optimized)
  Risk <- Risk_results$Risk # Predicted patient risk
  Heterogeneity.Data <- bind_cols(Trial_Data, Risk) # Adding risk to trial data
  
  ### Evaluating treatment effect heterogeneity from treatment-risk-interaction model
  Heterogeneity <- Eval_Heterogeneity(Heterogeneity.Data, Burn_in, Iterations, model_string.Heterogeneity)
  
  ### Returning Summary 
  
  if(Optimized == TRUE){
    summary <- list(Heterogeneity = Heterogeneity,
                    Risk = list(Risk = Risk_results,
                                Optimized_Risk = if_else(Optimized == TRUE, 
                                                         "Optimized Risk Model",
                                                         "Full Risk Model")),
                    SSVS = SSVS_results$included,
                    Selection_Probability = (SSVS_results$Selection_Probability),
                    Variable_Selections = (SSVS_results$Variable_Selections),
                    Iterations = Iterations,
                    Samples = Iterations * 4)
  } else {
    summary <- list(Heterogeneity = Heterogeneity,
                    Risk = list(Risk = (Risk),
                                Optimized_Risk = if_else(Optimized == TRUE, 
                                                         "Optimized Risk Model",
                                                         "Full Risk Model")),
                    Iterations = Iterations,
                    Samples = Iterations * 4)
  }
                     
                     
  
  
  return(summary)
}
