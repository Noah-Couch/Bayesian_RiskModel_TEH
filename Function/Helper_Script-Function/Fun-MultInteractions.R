
Opt_MultInt <- function(Trial_Data, included, Burn_in, Iterations,
                        model_string.MultInt) {
  
  ### Checking for Columns -----------------------------------------------------
  
  col_check <- c("PID", "Treatment", "Outcome", "Z3")
  if(all(col_check %in% colnames(Trial_Data))){} else{
    print("Trial Data should include 'PID', 'Treatment', and 'Outcome'")
    break
  }
  
  
  ### Setting initial values ---------------------------------------------------
  
  ### Design
  X <- Trial_Data |> 
    select(any_of(included), Treatment) 
  
  for (i in 1:length(included)) {
    j <- as.numeric(str_remove(included[i], "Z"))
    X <- X |>
      rowwise() |>
      mutate("Interaction_Z{j}":= Treatment * !!sym(included[i]))
  }
  
  X <- X |>
    select(Treatment, contains("Z"), contains("Interaction_Z")) |>
    as.matrix()
  
  ### Outcome
  delta_Y <- Trial_Data$Outcome - Trial_Data$Z3
  
  ### Data dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  ### List of initial values for JAGS samples
  initial_values <- list(Y = delta_Y,
                         n = n,
                         X = X,
                         p = p)
  
  ### Initializing model -------------------------------------------------------
  
  model <- jags.model(textConnection(model_string.MultInt), 
                      data = initial_values,
                      n.chains=3,
                      quiet = TRUE)
  
  
  ### Sampling from posterior --------------------------------------------------
  
  update(model, Burn_in)
  JAGS_samps <- coda.samples(model,
                             variable.names=c("beta0", "beta"), 
                             n.iter= Iterations)
  
  
  ### Returning output ---------------------------------------------------------
  
  Samps <- JAGS_samps[[1]] |> rbind(JAGS_samps[[2]], JAGS_samps[[3]])
  
  Post_Means <- data.frame(t(colMeans(Samps)))
  colnames(Post_Means) <- c("beta_Trt", included, 
                            paste0("Interaction_", included), "beta0")
  
  HPDI <- data.frame(hdi(Samps))
  colnames(HPDI) <- c("beta_Trt", included, 
                      paste0("Interaction_", included), "beta0")
  
  Int_CI <- HPDI |>
    as.data.frame() |>
    select(contains("Interaction"))
  
  Rejection <- if_else(any((Int_CI[1,] >= 0) | (Int_CI[2,] <= 0)),
                       TRUE,
                       FALSE)
    
  return(Rejection)
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  