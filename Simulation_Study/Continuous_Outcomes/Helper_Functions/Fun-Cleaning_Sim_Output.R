
########################## Cleaning Simulation Output ##########################


Clean_Output <- function(Trial_Data, Heterogeneity, Heterogeneity_Noise, Alternative_Methods, N_trial) {
  
  if(any(Heterogeneity$Heterogeneity$Posterior_Means == "N/A")){
    return("N/A")
  }
  
  
  ### Defining Heterogeneity Outcomes ------------------------------------------
  
  ### From the heterogeneity model
  Beta_Trt <- Heterogeneity$Heterogeneity$Posterior_Means$beta_Trt
  Beta_Risk <- Heterogeneity$Heterogeneity$Posterior_Means$beta_Risk
  Beta_intx <- Heterogeneity$Heterogeneity$Posterior_Means$beta_intx
  
  HDI <- Heterogeneity$Heterogeneity$Intervals$beta_intx
  
  ### From Risk model
  Prediction_Error <- (sum((Trial_Data$Outcome - Heterogeneity$Risk$Risk$Risk)^2)) / N_trial
  
  ### From SSVS model
  included <- Heterogeneity$SSVS
  
  
  ### Defining Noisy Heterogeneity Outcomes ------------------------------------
  
  ### From heterogeneity (noise) model
  Beta_intx_Noise <- Heterogeneity_Noise$Heterogeneity$Posterior_Means$beta_intx
  
  HDI_Noise <- Heterogeneity_Noise$Heterogeneity$Intervals$beta_intx
  
  ### From risk (noise) model
  Prediction_Error_Noise <- sum((Trial_Data$Outcome - Heterogeneity_Noise$Risk$Risk$Risk)^2) / N_trial
  
  
  ### Defining Alternative Methods to Heterogeneity Outcomes--------------------
  
  ### From interaction model
  Intx_pval <- Alternative_Results$Intx_pval
  
  ### From SDIR
  SDIR = Alternative_Results$SDIR
  SDIR_pval = Alternative_Results$SDIR_pval
  
  
  ### Creating Data Frame
  Outcomes <- data.frame(### Optimized Risk Model ------------------------------
                         Beta_Trt = Beta_Trt,
                         Beta_Risk = Beta_Risk,
                         Beta_intx = Beta_intx_Noise,
                         Lower.CI = HDI[1],
                         Upper.CI = HDI[2],                     
                         rejection.Risk = if_else((HDI[1] < 0 & 0 < HDI[2]), 
                                                  "Accept", "Reject"),
                         SSVS_results = paste(included, collapse = ", " ),
                         Prediction_Error = Prediction_Error,
                         ### Noisy Risk Model ----------------------------------
                         Beta_intx_Noise = Beta_intx_Noise,
                         Lower.CI_Noise = HDI_Noise[1],
                         Upper.CI_Noise = HDI_Noise[2],
                         rejection.Risk_Noise = if_else(
                           Noise.type == "Noise",
                           if_else(HDI_Noise[1] < 0 & 0 < HDI_Noise[2], 
                                   "Accept", "Reject"), 
                           NA),
                         Prediction_Error_Noise = Prediction_Error_Noise,
                         ### Interaction Model ---------------------------------
                         rejection.interaction = if_else(Intx_pval < 0.05,
                                                         "Reject", "Accept"),
                         ### SDIR ----------------------------------------------
                         rejection.SDIR = if_else(SDIR_pval <= 0.05,                                                        
                                                  "Reject", "Accept"),
                         SDIR = SDIR)
}
