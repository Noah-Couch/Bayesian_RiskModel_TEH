
### JAGS Scripts for risk modeling approach to evaluate TEH --------------------

### JAGS script for Stochastic Search Variable Selection
model_string.SSVS <- "model{  
  # Likelihood
  for(i in 1:n){
  Y[i]  ~  dnorm(mu[i], tau2)
  mu[i] <- beta0 + inprod(beta[],X[i,])
  }
  #Priors
  for(j in 1:p){
  gamma[j] ~ dnorm(0,0.01)
  delta[j] ~ dbern(prob)
  beta[j]  <- gamma[j] * delta[j]
  }
  beta0 ~ dnorm(0,0.01)
  prob  ~ dunif(0,1)
  tau2 ~ dgamma(0.1,0.1)
  sigma2 <- pow(tau2,-1)}"

### JAGS script for the risk model
model_string.Risk <- "model{
  # Likelihood
  for(i in 1:n){
    Y[i]  ~  dnorm(mu[i], tau2)
    mu[i] <- beta0 + inprod(beta[],X[i,])
  }
  #Priors
   for(j in 1:p){
   beta[j] ~ dnorm(0,0.01)
   }
   beta0 ~ dnorm(0,0.01)
   tau2 ~ dgamma(0.1,0.1)
   sigma2 <- pow(tau2,-1)
 }"

### JAGS Script for the full model
model_string.Heterogeneity <- "model{
  # Likelihood
  for(i in 1:n){
    Y[i] ~ dnorm(mu[i], tau2)
    mu[i] <- beta0 + inprod(beta[],X[i,])
  }
  #Priors
   for(j in 1:p){
      beta[j] ~ dnorm(0,0.01)
   }
   beta0 ~ dnorm(0,0.01)
   tau2 ~ dgamma(0.1,0.1)
   sigma2 <- pow(tau2,-1)
 }"

### Helper function: Stochastic search variable selection

SSVS.function <- function(X, Y) {
  
  ### Data dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  ### List of initial values for JAGS samples
  initial_values <- list(Y = Y,
                         n = n,
                         X = X,
                         p = p)
  
  ### Initializing model
  model <- jags.model(textConnection(model_string.SSVS),
                      data = initial_values,
                      n.chains=3)
  
  ### Running a 10,000 sample burn in
  update(model, 10000)
  
  ### Sampling from the posterior
  JAGS_samps  <- coda.samples(model, 
                              variable.names=c("beta0", "beta","delta"), 
                              n.iter=50000)
  
  Samps <- JAGS_samps[[1]] |>
    rbind(JAGS_samps[[2]], JAGS_samps[[3]])
  
  delta <- Samps[,1:p + p + 1]
  
  ### Defining the column names for which covariates each beta and delta represent
  #colnames(beta.SSVS) <- colnames(X.SSVS)
  colnames(delta.SSVS) <- colnames(X.SSVS) 
  
  ### Pulling the coefficients to include based on SSVS
  included <- delta.SSVS |>
    as.data.frame() |>
    count(across(everything()), sort = TRUE) |>
    slice(1) |>
    select(-n) |>
    pivot_longer(everything(),
                 names_to = "coef",
                 values_to = "value") |>
    filter(value == 1) |>
    pull(coef)
  
  X_included <- X |>
    as.data.frame() |>
    select(any_of(included))
  
  return(X_included)
}

RiskModel.function <- function(X, Y) {
  
  ### Data dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  ### List of initial values for JAGS samples
  initial_values <- list(Y = Y,
                         n = n,
                         X = X,
                         p = p)
  
  ### Initializing model
  model <- jags.model(textConnection(model_string.Risk),
                      data = initial_values,
                      n.chains=3)
  
  ### Running a 10,000 sample burn in
  update(model, 10000)
  
  ### Sampling from the posterior
  JAGS_samps  <- coda.samples(model, 
                              variable.names=c("beta0", "beta"), 
                              n.iter=50000)
  
  Samps <- JAGS_samps[[1]] |>
    rbind(JAGS_samps[[2]], JAGS_samps[[3]])
  
  ### Pulling betas
  beta0 <- Samps[,(p+1)]
  beta <- Samps[,1:p] |> as.data.frame()
  
  colnames(beta) <- colnames(X)
  
  Risk <- data.frame(Risk = X %*% as.matrix(colMeans(beta)))
  
  return(Risk)
}

