
############################# Generating Outcomes ##############################

calculating_Outcomes <- function(Simulated_Data, Heterogeneity.type, ProgVar.type, Noise.type) {
  
  ### Defining Coefficients for outcome generation -----------------------------

  beta.GeneratingOutcome <- case_when(
    ### No Heterogeneity and Categorical Prognostic Variable
    Heterogeneity.type == "No_Heterogeneity" & 
      ProgVar.type == "Categorical" ~  c(1.00,    # Intercept
                                         0.50,    # Treatment
                                         0.75,    # Categorical 
                                         0.00,    # Treatment * Categorical Predictor
                                         -0.75,   # Gender
                                         -0.75),  # Age
    ### Heterogeneity and Categorical Prognostic Variable
    Heterogeneity.type == "Heterogeneity" & 
      ProgVar.type == "Categorical" ~  c(1.00,    # Intercept
                                         0.50,    # Treatment
                                         0.75,    # Categorical 
                                         1.00,    # Treatment * Categorical Predictor
                                         -0.75,   # Gende
                                         -0.75),  # Age
    ### No Heterogeneity and Continous Prognostic Variable
    Heterogeneity.type == "No_Heterogeneity" & 
      ProgVar.type == "Continuous" ~   c(1.00,    # Intercept
                                         1.00,    # Treatment
                                         0.75,    # Continuous 
                                         0.00,    # Treatment * Continuous Predictor
                                         -0.75,   # Gender
                                         -0.75),  # Age
    ### Heterogeneity and Continuous Prognostic Variable
    Heterogeneity.type == "Heterogeneity" & 
      ProgVar.type == "Continuous" ~   c(1.00,    # Intercept
                                         1.00,    # Treatment
                                         0.75,    # Continuous 
                                         1.00,    # Treatment * Continuous Predictor
                                         -0.75,   # Gender
                                         -0.75),  # Age
    TRUE ~ NA)
  
  
  ### Generating Outcomes ------------------------------------------------------
  
  ### Calculating Trial Outcome
  Trial_Data <- Simulated_Data |>
    rowwise() |>
    mutate(mu = Outcome_BL +
             beta.GeneratingOutcome[1] + 
             beta.GeneratingOutcome[2] * Treatment + 
             beta.GeneratingOutcome[3] * ProgVar + 
             beta.GeneratingOutcome[4] * Treatment * ProgVar + 
             beta.GeneratingOutcome[5] * Gender +
             beta.GeneratingOutcome[6] * Age,
           Outcome = rnorm(1, mean = mu, sd = 0.5)) |> 
    select(-mu)
  
  
  ### Generating (weakly) correlated noise variables ---------------------------
  
  if(Noise.type == "Noise") {
    Trial_Data <- Trial_Data |>
      mutate(X1 = Age * 0.8 + Outcome_BL * 0.8 + rnorm(n, sd = 0.5),
             X2 = rbinom(1, size = 1, prob = 0.3 + 0.25 * Gender ))
  }
  
  return(Trial_Data)
}





