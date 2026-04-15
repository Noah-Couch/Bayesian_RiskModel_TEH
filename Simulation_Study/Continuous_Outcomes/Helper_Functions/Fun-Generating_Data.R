
##################### Generating Baseline Characteristics ######################

generating_data <- function(PID, Treatment, Z4.type) {
  
  ### Generating baseline characteristics --------------------------------------
  
  ### Z1
  Z1 <- rbinom(N_trial, size = 1, prob = 0.5)
  
  ### Z2
  Z2 <- rnorm(N_trial, mean = 0, sd = 1)
  
  ### Baseline Measurement of Outcome (Y_0)
  Z3 <- rnorm(N_trial, mean = 0, sd = 1)
  
  ### Prognostic Variable (Categorical or Continuous)
  Z4 <- case_when(Z4.type == "Categorical" ~ rbinom(N_trial, size = 1, prob = 0.5),
                  Z4.type == "Continuous"  ~ rnorm(N_trial, mean = 0, sd = 1))
  
  ### Combining variables to create simulated trial data -----------------------
  
  Simulated_Data <- tibble(PID, Treatment, Z1, Z2, Z3, Z4)
  
  return(Simulated_Data)
}






