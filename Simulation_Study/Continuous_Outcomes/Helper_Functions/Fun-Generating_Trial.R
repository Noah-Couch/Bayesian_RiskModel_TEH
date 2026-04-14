
############################## Required Functions ##############################

### Generating baseline characteristics
source("~/Thesis/Bayesian_RiskModel_TEH/Simulation_Study/Continuous_Outcomes/Helper_Functions/Fun-Generating_Data.R")

### Outcome Generation
source("~/Thesis/Bayesian_RiskModel_TEH/Simulation_Study/Continuous_Outcomes/Helper_Functions/Fun-Calculating_Outcomes.R")

############################ Generating Trial Data #############################

Gen_Trial_Data <- function(N_trial, seed, Heterogeneity.type, ProgVar.type, Noise.type) {
  
  ### Patient Identification and Assignment ------------------------------------
  ### Simulated Study Design: Single treatment group v. control
  ###                         N; n1 = n2 = N/2
  
  ### Variables that remain constant across simulations (PID and Treatment Assignment)
  
  ### Participant ID
  PID <- seq(from = 1, to = N_trial, by = 1)
  
  ### Treatment Assignment
  Treatment <- c( rep( 1, N_trial/2 ),
                  rep( 0, N_trial/2 ))
  
  
  ### Generating Baseline Characteristics --------------------------------------
  
  Simulated_Data <- generating_data(seed, PID, Treatment, ProgVar.type)
  
  
  ### Generating Outcomes ------------------------------------------------------
  
  Trial_Data <- calculating_Outcomes(Simulated_Data, Heterogeneity.type, ProgVar.type, Noise.type)
  
  return(Trial_Data)
}






