
########################### Evaluating Heterogeneity ###########################

require(rjags)
require(HDInterval)

Eval_Heterogeneity <- function(Heterogeneity.Data, Burn_in, Iterations, model_string.Heterogeneity) {
  
  ### Checking for Columns -----------------------------------------------------
  
  col_check <- c("PID", "Treatment", "Outcome", "Outcome_BL")
  if(all(col_check %in% colnames(Trial_Data))){} else{
    print("Trial Data should include 'PID', 'Treatment', and 'Outcome'")
    break
  }
  
  
  ### Setting initial values ---------------------------------------------------
  
  ### Design
  X <- Heterogeneity.Data |> 
    select(Treatment, Risk) |>
    mutate(Interaction = Treatment * Risk) |>
    as.matrix()
  
  delta_Y <- Heterogeneity.Data$Outcome - Heterogeneity.Data$Outcome_BL
  
  ### Data dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  ### List of initial values for JAGS samples
  initial_values <- list(Y = delta_Y,
                         n = n,
                         X = X,
                         p = p)
  
  
  ### Initializing model -------------------------------------------------------
  
  model <- jags.model(textConnection(model_string.Heterogeneity), 
                      data = initial_values,
                      n.chains=3)
  
  
  ### Sampling from posterior --------------------------------------------------
  
  update(model, Burn_in)
  JAGS_samps <- coda.samples(model,
                             variable.names=c("beta0", "beta"), 
                             n.iter= Iterations)
  
  
  ### Returning output ---------------------------------------------------------
  
  Samps <- JAGS_samps[[1]] |> rbind(JAGS_samps[[2]], JAGS_samps[[3]])
  
  Post_Means <- data.frame(t(colMeans(Samps)))
  colnames(Post_Means) <- c("beta_Trt", "beta_Risk", "beta_intx", "beta0")
  
  HPDI <- data.frame(hdi(Samps))
  colnames(HPDI) <- c("beta_Trt", "beta_Risk", "beta_intx", "beta0")
  
  Rejection_Het <- (HPDI$beta_intx[1] >= 0) | (HPDI$beta_intx[2] <= 0)
  Heterogeneity <- if_else(Rejection_Het, 
                           "There is evidence of Treatment Effect Heterogeneity",
                           "No evidence of Treatment Effect Heterogeneity" )
  
  R_hat <- gelman.diag(JAGS_samps)
  
  
  
  return(list(Posterior_Means = Post_Means,
              Intervals = HPDI,
              Heterogeneity = Heterogeneity,
              R_hat = R_hat))
  }