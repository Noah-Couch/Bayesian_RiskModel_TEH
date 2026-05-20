
############################# Generating Outcomes ##############################

calculating_Outcomes <- function(Simulated_Data, Noise.Vars, T_Effect, H_Effect,
                                 Het_Type,
                                 Outcome_SD) {
  
  ### Generating Outcomes ------------------------------------------------------
  
  
  
  
  
  ### Calculating Trial Outcome
  Trial_Data <- Simulated_Data |>
    rowwise() |>
    mutate(mu = if_else(Het_Type == "uniVar",
                        Z3 + 1.00 + 
                          T_Effect * Treatment + 
                          0.75 * Z4 + 
                          H_Effect * Treatment * Z4 + 
                          -0.75 * Z1 +
                          -0.75 * Z2,
                        Het_Type == "linComb", 
                        Z3 + 1.00 +
                          T_Effect * Treatment + 
                          0.75 * Z4 + 
                          H_Effect * Treatment * (0.33 * Z4 + 0.33 * Z2 + 0.33 * Z1) + 
                          -0.75 * Z1 +
                          -0.75 * Z2,),
           Outcome = rnorm(1, mean = mu, sd = Outcome_SD)
           ) |>
    select(-c(mu))
  
  
  ### Generating (weakly) correlated noise variables ---------------------------
  
  if(Noise.Vars >= 1) {
    for(i in 1:Noise.Vars){
      Trial_Data <- Trial_Data |>
        mutate("X{i}":= Z2 * 0.8 + Z3 * 0.8 + rnorm(1, sd = 0.5))
      
    }
    
    
    #Trial_Data <- Trial_Data |>
    #  mutate(X1 = Z2 * 0.8 + Z3 * 0.8 + rnorm(1, sd = 0.5),
    #         X2 = rbinom(1, size = 1, prob = 0.3 + 0.25 * Z1 )
    #         )
  }
  
  return(Trial_Data)
}





