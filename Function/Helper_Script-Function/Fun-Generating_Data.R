
############################ Generating Trial Data #############################

generating_data <- function(seed, PID, Treatment, ProgVar.type) {
  
  ### Generating baseline characteristics --------------------------------------
  
  ### Gender
  Gender <- rbinom(200, size = 1, prob = 0.5)
  
  ### Age
  Age <- rnorm(200, mean = 0, sd = 1)
  
  ### Baseline Measurement of Outcome (Y_0)
  Outcome_BL <- rnorm(200, mean = 0, sd = 1)
  
  ### Prognostic Variable (Categorical or Continuous)
  ProgVar <- case_when(ProgVar.type == "Categorical" ~ rbinom(200, size = 1, prob = 0.5),
                       ProgVar.type == "Continuous"  ~ rnorm(200, mean = 0, sd = 1))
  
  ### Combining variables to create simulated trial data -----------------------
  
  Simulated_Data <- tibble(PID, Treatment, Gender, Age, Outcome_BL, ProgVar)
  
  return(Simulated_Data)
}






