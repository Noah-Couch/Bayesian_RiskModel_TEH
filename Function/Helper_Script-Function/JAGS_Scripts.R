
######### JAGS Script for Stochastic Search Variable Selection (SSVS) ##########
################################################################################

model_string.SSVS <- 
  "model{  
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


###################### JAGS Script for Risk Calculations #######################
################################################################################

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


################## JAGS Script for Risk Modeled Heterogeneity ##################
################################################################################

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


############## JAGS Script for Multiple Interactions Heterogeneity #############
################################################################################

model_string.MultInt <- "model{
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