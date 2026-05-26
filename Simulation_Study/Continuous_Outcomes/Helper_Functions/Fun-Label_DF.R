
require(purrr)

Label_DF <- function(Power_DF) {
  Label_Data <- Power_DF |>
    mutate(facet = paste0("Treatment Effect: ", T_Effect,
                          ", Heterogeneous Effect: ", H_Effect)) |>
    pivot_longer(
      cols = starts_with("Power"),
      names_to = "Method",
      values_to = "Power"
    ) |>
    mutate(Method = str_remove(Method, "Power.")) |>
    mutate(
      Method = case_when(
        Method == "Risk" ~ "Optimized Model",
        Method == "Risk_Noise" ~ "Non-optimized Model",
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
    filter(H_Effect == max(H_Effect)) |>
    group_by(facet, Method) |>
    filter(N == max(N)) |>
    ungroup() |>
    group_by(facet) |>
    arrange(desc(Power), .by_group = TRUE) |>
    mutate(
      y_raw = Power,
      y_label = purrr::accumulate(
        y_raw,
        ~ min(.y, .x - 0.05)
      ) |> unlist(),
      y_label = pmax(y_label, 0.04),
      x_elbow = 208,
      x_label = 222,
      label = Method
    ) |>
    ungroup()
  
  return(Label_Data)
}



