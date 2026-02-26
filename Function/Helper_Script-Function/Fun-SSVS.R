
source("/deac/sta/beaversGrp/coucns24/Bayesian_RiskModel_TEH/Function/Helper_Script-Function/JAGS_Scripts.R")

########### Implementing Stochastic Search Variable Selection (SSVS) ###########

require(rjags)

SSVS_jags <- function(Trial_Data, Burn_in, Iterations, model_string.SSVS) {
  
  ### Checking for Columns -----------------------------------------------------
  
  col_check <- c("PID", "Treatment", "Outcome")
  if(all(col_check %in% colnames(Trial_Data))){} else{
    print("Trial Data should include 'PID', 'Treatment', and 'Outcome'")
    break
  }
  
  
  ### Setting initial values ---------------------------------------------------
  
  ### Creating the Design Matrix and Outcome variable for JAGS script (SSVS)
  ### Ordering of variables:
  ###     - Gender
  ###     - Age
  ###     - Outcome_BL
  ###     - Prognostic Variable
  ###     - X1 (Noise)
  ###     - X2 (Noise)
  
  ### Design
  X <- Trial_Data |> 
    select(-c(PID, Treatment, Outcome)) |> 
    as.matrix()
  
  ### Outcome
  delta_Y <- Trial_Data$Outcome - Trial_Data$Outcome_BL
  
  ### Data dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  ### List of initial values for JAGS samples
  initial_values <- list(Y = delta_Y,
                         n = n,
                         X = X,
                         p = p)
  
  
  ### Initializing model -------------------------------------------------------
  
  model <- jags.model(textConnection(model_string.SSVS),
                      data = initial_values,
                      n.chains=3)
  
  
  ### Sampling from posterior --------------------------------------------------
  
  update(model, Burn_in)
  JAGS_samps  <- coda.samples(model, 
                              variable.names=c("delta"), 
                              n.iter= Iterations)
  
  
  ### Returning posterior selection numbers ------------------------------------
  
  ### Combining the samples of all three chains
  Samps <- JAGS_samps[[1]] |> rbind(JAGS_samps[[2]], JAGS_samps[[3]])
  
  ### Pulling deltas
  delta <- Samps
  
  ### Defining the column names for which covariates each beta and delta represent
  colnames(delta) <- colnames(X)
  
  return(delta)
}















### Sampling from the posterior
