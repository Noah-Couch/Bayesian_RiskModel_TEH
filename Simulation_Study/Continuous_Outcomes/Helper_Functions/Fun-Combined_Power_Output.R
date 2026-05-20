
require()




H_Effects <- H_Effects[2:4]
T_Effects


base_path2 <- paste0("~/Bayesian_RiskModel_TEH/Simulation_Study/Continuous_Outcomes",
                     "/Output/",
                    Het_Type, "/",
                    Z4.type, "/")
                    
colnames(Power)
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

for (i in 1:length(T_Effects)) {
  
  T_Effect <- T_Effects[i]

for (j in 1:length(H_Effects)) {
  
  H_Effect <- H_Effects[j]
  
  path <- paste0(base_path2, "TrtEff_", T_Effect, "-HetEff_", H_Effect, "/")
  
  file_names <- list.files(path, pattern = paste0("\\OutcomeSD_",Outcome_SD,".csv$"), full.names = TRUE)
  
  list_csvs <- lapply(file_names, read_csv)
  names(list_csvs) <- c(100, 150, 200, 30, 50)
  
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
}  
    
}
                    
Plot <- Power |>
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
  ) |> 
  mutate(label = if_else(N == max(N), Method, NA_character_)) |>
  ggplot(aes(x = N, y = Power, group = Method, linetype = LineType)) + 
  # Geometric annotations that play the role of grid lines
  geom_vline(
    xintercept = seq(0, 200, by = 25),
    color = "grey91", 
    size = .6
  ) +
  geom_segment(
    data = tibble(y = seq(0, 1, by = 0.25), x1 = 10, x2 = 200),
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    color = "grey91",
    size = .6
  ) +
  geom_segment(
    data = tibble(y = 0, x1 = 10, x2 = 200),
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    color = "grey60",
    size = .8
  ) +
  geom_point(aes(color = Method == "Optimized Model")) +
  geom_line(aes(color = Method == "Optimized Model")) +
  facet_wrap(. ~ facet) +
  geom_text_repel(aes(color = Method == "Optimized Model", label = label, 
                      segment.linetype = LineType),
                  family = "Lato",
                  fontface = "bold",
                  size = 2,
                  direction = "y",
                  xlim = c(220, NA),
                  hjust = 0,
                  segment.size = .7,
                  segment.alpha = .5,
                  #segment.linetype = "dotted",
                  box.padding = .4,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  na.rm = TRUE) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(10, 300), 
    breaks = seq(50, 200, by = 50)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 1), 
    breaks = seq(0, 1, by = 0.25)
  ) +
  scale_linetype_identity(
    aesthetics = c("linetype", "segment.linetype")
  ) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#48494B"))+
  guides(color = "none") +
  xlab("Sample Size, N") +
  ylab("Power")

Plot


