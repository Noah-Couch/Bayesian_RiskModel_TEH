
############################### Calculating Risk ###############################

require(rjags)

Risk_jags <- function(Trial_Data, included, Burn_in, Iterations,
                      model_string.Risk, Optimized = TRUE) {
  
  ### Checking for optimized risk or noisy risk --------------------------------
  
  if(Optimized == FALSE) {
    included <- Trial_Data |> 
      select(-c(PID, Treatment, Outcome)) |>
      colnames()
  }
  
  
  ### Setting initial values ---------------------------------------------------
  
  ### Design
  X <- Trial_Data |> 
    select(any_of(included)) |>
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
  
  model <- jags.model(textConnection(model_string.Risk),
                      data = initial_values,
                      n.chains=3)
  
  
  ### Sampling from posterior --------------------------------------------------
  
  update(model, Burn_in)
  JAGS_samps  <- coda.samples(model, 
                              variable.names=c("beta0", "beta"), 
                              n.iter= Iterations)
  
  
  ### Calculating posterior means for beta coefficients ------------------------
  
  ### Combining the samples of all three chains
  Samps <- JAGS_samps[[1]] |> rbind(JAGS_samps[[2]], JAGS_samps[[3]])
  
  ### Pulling betas
  beta0 <- Samps[,(p+1)]
  beta <- Samps[,1:p] |> as.data.frame()
  
  colnames(beta) <- colnames(X)
  
  ### Calculating and returning Risk -------------------------------------------
  
  Risk <- data.frame(Risk = mean(beta0) + X %*% as.matrix(colMeans(beta)))
  
  return(Risk)
}


