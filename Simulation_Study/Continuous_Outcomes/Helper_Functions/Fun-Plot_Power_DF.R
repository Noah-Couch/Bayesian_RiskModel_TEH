
Plot_Power_DF <- function(Power_DF) {
  Plot_Power_DF <- Power |>
    mutate(facet = paste0("Treatment Effect: ",T_Effect, ", Heterogeneous Effect: ", H_Effect)) |>
    pivot_longer(
      cols = starts_with("Power"),
      names_to = "Method",
      values_to = "Power") |>
    mutate(Method = str_remove(Method, "Power.")) |>
    mutate(
      Method = case_when(Method == "Risk" ~ "Optimized Model",
                         Method == "Risk_Noise" ~ "Non-optimized Model",
                         #Method == "MultInt" ~ "Optimized Multiple Interactions",
                         Method == "Intx" ~ "Interaction",
                         Method == "Intx_adj" ~ "Interaction (Adjusted)",
                         TRUE ~ Method
      ),
      LineType = case_when(
        Method == "Optimized Model" ~ "solid",
        Method == "Non-optimized Model" ~ "dashed",
        Method == "Interaction" ~ "dotted",
        Method == "Interaction (Adjusted)" ~ "dotdash",
        TRUE ~ "solid"
      )
    )
  
  return(Plot_Power_DF)
}



