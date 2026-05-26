
########################## Cleaning Simulation Output ##########################


Clean_Output <- function(Trial_Data, Heterogeneity, Heterogeneity_Noise, Alternative_Results, N_trial, Noise.Vars) {
  
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
  Y <- mean(Trial_Data$Outcome)
  Y_hat <- sum(Heterogeneity$Risk$Risk$Risk) / N_trial
  Bias <- sum(Heterogeneity$Risk$Risk$Risk - Trial_Data$Outcome) / N_trial
  abs_Bias <- sum(abs(Heterogeneity$Risk$Risk$Risk - Trial_Data$Outcome)) / N_trial
  MSE  <- sum((Heterogeneity$Risk$Risk$Risk - Trial_Data$Outcome)^2) /N_trial
  waic <- Heterogeneity$Risk$Risk$waic
  
  
  ### From SSVS model
  included <- Heterogeneity$SSVS
  
  
  ### Defining Noisy Heterogeneity Outcomes ------------------------------------
  
  ### From heterogeneity (noise) model
  Beta_intx_Noise <- Heterogeneity_Noise$Heterogeneity$Posterior_Means$beta_intx
  
  HDI_Noise <- Heterogeneity_Noise$Heterogeneity$Intervals$beta_intx
  
  ### From risk (noise) model
  Y_hat_noise <- mean(Heterogeneity_Noise$Risk$Risk$Risk)
  Bias_noise <- sum(Heterogeneity_Noise$Risk$Risk$Risk - Trial_Data$Outcome) / N_trial
  abs_Bias_noise <- sum(abs(Heterogeneity_Noise$Risk$Risk$Risk - Trial_Data$Outcome)) / N_trial
  MSE_noise  <- sum((Heterogeneity_Noise$Risk$Risk$Risk - Trial_Data$Outcome)^2) /N_trial
  waic_noise <- Heterogeneity_Noise$Risk$waic
  
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
                         Y = Y,
                         Y_hat = Y_hat,
                         Bias = Bias,
                         abs_Bias = abs_Bias,
                         MSE = MSE,
                         waic = waic,
                         ### Noisy Risk Model ----------------------------------
                         Beta_intx_Noise = Beta_intx_Noise,
                         Lower.CI_Noise = HDI_Noise[1],
                         Upper.CI_Noise = HDI_Noise[2],
                         rejection.Risk_Noise = if_else(
                           Noise.Vars >= 1,
                           if_else(HDI_Noise[1] < 0 & 0 < HDI_Noise[2], 
                                   "Accept", "Reject"), 
                           NA),
                         Y_hat_noise = Y_hat_noise,
                         Bias_noise = Bias_noise,
                         abs_Bias_noise = abs_Bias_noise,
                         MSE_noise = MSE_noise,
                         waic_noise = waic_noise,
                         ### Multiple Interactions (Optimized) -----------------
                         rejection.Opt_MultInt = if_else(Heterogeneity$Rejection_OptMultInt,
                                                         "Reject", "Accept"),
                         ### Interaction Model ---------------------------------
                         rejection.interaction = if_else(Intx_pval < 0.05,
                                                         "Reject", "Accept"),
                         rejection.interaction_adj = if_else(Intx_pval < (0.05 / (4 + Noise.Vars)),
                                                             "Reject", "Accept"),
                         ### SDIR ----------------------------------------------
                         rejection.SDIR = if_else(SDIR_pval <= 0.05,                                                        
                                                  "Reject", "Accept"),
                         SDIR = SDIR)
}
