
source("~/Bayesian_RiskModel_TEH/Simulation_Study/Continuous_Outcomes/Helper_Functions/Fun-Label_DF.R")
source("~/Bayesian_RiskModel_TEH/Simulation_Study/Continuous_Outcomes/Helper_Functions/Fun-Plot_Power_DF.R")
source("~/Bayesian_RiskModel_TEH/Simulation_Study/Continuous_Outcomes/Helper_Functions/Fun-Organizing_Plot.R")

Combined_Power_Plot <- function(H_Effects = c(0.00, 0.25,  0.50, 1.00), 
                                T_Effects = c(0.50, 1.00), 
                                Outcome_SD,
                                Noise.Vars,
                                Het_Type, 
                                Z4.type) {
  
  ### Defining Variables -------------------------------------------------------
  
  ### Removing simulations for type I error
  H_Effects <- H_Effects[H_Effects != 0] 
  ### Defining the empty power data frame
  Power <- data.frame(
    N = double(),
    Power.Risk = double(),
    Power.Risk_Noise = double(),
    Power.Intx = double(),
    Power.Intx_adj = double(),
    Power.SDIR = double(),
    T_Effect = double(),
    H_Effect = double()
  )
  ### Defining the base path to pull csv files from
  base_path2 <- paste0("~/Bayesian_RiskModel_TEH/Simulation_Study/Continuous_Outcomes",
                       "/Output/",
                       Het_Type, "/",
                       Z4.type, "/")
  
  
  ### Looping over treatment and heterogeneous effects -------------------------
  
  for (i in 1:length(T_Effects)) {
    
    T_Effect <- T_Effects[i]
    
    for (j in 1:length(H_Effects)) {
      
      H_Effect <- H_Effects[j]
    
      ### Pulling csv files
      file_names <- list.files(paste0(base_path2, "TrtEff_", T_Effect, "-HetEff_", H_Effect, "/"), 
                               pattern = paste0("\\OutcomeSD_",Outcome_SD,".csv$"), 
                               full.names = TRUE)
      list_csvs <- lapply(file_names, read_csv, show_col_types = FALSE)
      names(list_csvs) <- c(100, 150, 200, 30, 50)
      
      ### Calculating and adding rows to the power data frame
      Power_df <- map_dfr(list_csvs, function(df) {
        df |>
          summarise(
            Power.Risk = sum(rejection.Risk == "Reject", na.rm = TRUE) / N_iter,
            Power.Risk_Noise = sum(rejection.Risk_Noise == "Reject", na.rm = TRUE) / N_iter,
            #Power.MultInt = sum(rejection.Opt_MultInt == "Reject", na.rm = TRUE) / N_iter,
            Power.Intx = sum(rejection.interaction == "Reject", na.rm = TRUE) / N_iter,
            Power.Intx_adj = sum(rejection.interaction_adj == "Reject", na.rm = TRUE) / N_iter,
            Power.SDIR = sum(rejection.SDIR == "Reject", na.rm = TRUE) / N_iter
          ) },
        .id = "N") |>
        mutate(N = as.numeric(N),
               T_Effect = T_Effects[i],
               H_Effect = H_Effects[j]) |>
        arrange(N)
      Power <- rbind(Power, Power_df)
      
    } # H Effects
  } # T Effects
  
  ### Cleaning the power data fram for plotting --------------------------------
  Power_Plot.df <- Plot_Power_DF(Power_DF = Power)
  
  ### Creating a data frame for method labels ----------------------------------
  Label_Data <- Label_DF(Power_DF = Power)
  
  Plot <- Power_Plot.df |>
    ggplot(aes(x = N, y = Power, group = Method, linetype = LineType)) %>%
    Organizing_Plot(base_plot = ., Label_Data = Label_Data)
  
  return(Plot)
}

