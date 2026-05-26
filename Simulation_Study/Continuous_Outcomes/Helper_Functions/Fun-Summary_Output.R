
require(purrr)
require(ggplot2)
require(ggrepel)
require(ggtext)

theme_set(theme_minimal(base_family = "Lato"))

theme_update(
  # Remove the grid lines that come with ggplot2 plots by default
  panel.grid = element_blank(),
  # Customize margin values (top, right, bottom, left)
  plot.margin = margin(20, 40, 20, 40),
  # Use a light grey color for the background of both the plot and the panel
  plot.background = element_rect(fill = "white", color = "white"),
  panel.background = element_rect(fill = "white", color = "white"),
  # Customize title appearence
  plot.title = element_text(
    color = "grey10", 
    size = 28, 
    face = "bold",
    margin = margin(t = 15)
  ),
  # Customize subtitle appearence
  plot.subtitle = element_markdown(
    color = "grey30", 
    size = 16,
    lineheight = 1.35,
    margin = margin(t = 15, b = 40)
  ),
  # Title and caption are going to be aligned
  plot.title.position = "plot",
  plot.caption.position = "plot",
  plot.caption = element_text(
    color = "grey30", 
    size = 13,
    lineheight = 1.2, 
    hjust = 0,
    margin = margin(t = 40) # Large margin on the top of the caption.
  ),
  axis.title = element_text(size = 18),
  # Remove legend
  legend.position = "none"
)

Summary_Output <- function(base_path, Outcome_SD, N_iter) {
  file_names <- list.files(base_path, pattern = paste0("\\OutcomeSD_",Outcome_SD,".csv$"), full.names = TRUE)
  
  list_csvs <- lapply(file_names, read_csv)
  names(list_csvs) <- c(100, 150, 200, 30, 50)
  
  Power <- map_dfr(list_csvs, function(df) {
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
    mutate(N = as.numeric(N))
  
  Plot <- Power |>
    pivot_longer(
      cols = starts_with("Power"),
      names_to = "Method",
      values_to = "Power"
    ) |>
    mutate(Method = str_remove(Method, "Power.")) |>
    mutate(
      Method = case_when(Method == "Risk" ~ "Optimized Model",
                         Method == "Risk_Noise" ~ "Non-optimized Model",
                         #Method == "MultInt" ~ "Optimized Multiple Interactions",
                         Method == "Intx" ~ "Interaction",
                         Method == "Intx_adj" ~ "Interaction (Adjusted)",
                         TRUE ~ Method
      )
    ) |> 
    mutate(label = if_else(N == max(N), Method, NA_character_)) |>
    ggplot(aes(x = N, y = Power, group = Method)) + 
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
    geom_text_repel(aes(color = Method == "Optimized Model", label = label),,
                    family = "Lato",
                    fontface = "bold",
                    size = 2,
                    direction = "y",
                    xlim = c(220, NA),
                    hjust = 0,
                    segment.size = .7,
                    segment.alpha = .5,
                    segment.linetype = "dotted",
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
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#48494B"))+
    guides(color = "none") +
    xlab("Sample Size, N") +
    ylab("Power")
  
  Var_Selection <- map_dfr(list_csvs, function(df) {
    df |>
      summarise(
        SSVS.X1 = sum(str_detect(SSVS_results, "X1")) / N_iter,
        SSVS.X2 = sum(str_detect(SSVS_results, "X2")) / N_iter,
        SSVS.X3 = sum(str_detect(SSVS_results, "X3")) / N_iter,
        SSVS.X4 = sum(str_detect(SSVS_results, "X4")) / N_iter,
        SSVS.X5 = sum(str_detect(SSVS_results, "X5")) / N_iter,
        SSVS.X6 = sum(str_detect(SSVS_results, "X6")) / N_iter,
        SSVS.X7 = sum(str_detect(SSVS_results, "X7")) / N_iter,
        SSVS.X8 = sum(str_detect(SSVS_results, "X8")) / N_iter,
        SSVS.X9 = sum(str_detect(SSVS_results, "X9")) / N_iter,
        SSVS.X10 = sum(str_detect(SSVS_results, "X10") / N_iter)
      ) },
    .id = "N") |>
    mutate(N = as.numeric(N))
  
  
  
  
  return(list(Plot = Plot,
              Power = Power,
              Accuracy = Accuracy
              Var_Selection = Var_Selection))
}






