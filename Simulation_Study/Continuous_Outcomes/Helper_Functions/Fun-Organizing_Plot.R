
Organizing_Plot <- function(base_plot, Label_Data) {
  plot <- base_plot + 
    # Geometric annotations that play the role of grid lines
    geom_vline(
      xintercept = seq(25, 200, by = 25),
      color = "grey91", 
      size = .6
    ) +
    geom_segment(
      data = tibble(y = seq(0, 1, by = 0.25), x1 = 25, x2 = 200),
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      color = "grey91",
      size = .6
    ) +
    geom_segment(
      data = tibble(y = 0, x1 = 25, x2 = 200),
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      color = "grey60",
      size = .8
    ) +
    geom_point(aes(color = Method == "Optimized Model")) +
    geom_line(aes(color = Method == "Optimized Model")) +
    facet_wrap(. ~ facet) +
    geom_segment(
      data = Label_Data,
      aes(
        x = N,
        xend = x_elbow,
        y = Power,
        yend = y_label,
        color = Method == "Optimized Model",
        linetype = LineType
      ),
      inherit.aes = FALSE,
      linewidth = .5,
      alpha = .7
    ) +
    geom_segment(
      data = Label_Data,
      aes(
        x = x_elbow,
        xend = x_label - 1,
        y = y_label,
        yend = y_label,
        color = Method == "Optimized Model",
        linetype = LineType
      ),
      inherit.aes = FALSE,
      linewidth = .5,
      alpha = .7
    ) +
    geom_text(
      data = Label_Data,
      aes(
        x = x_label,
        y = y_label,
        label = label,
        color = Method == "Optimized Model"
      ),
      inherit.aes = FALSE,
      family = "Lato",
      fontface = "bold",
      size = 2.5,
      hjust = 0
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq(50, 200, by = 50)
    ) +
    coord_cartesian(
      xlim = c(25, 200),
      clip = "off"
    ) +
    theme(
      plot.margin = margin(8, 220, 8, 8),
      panel.spacing.x = unit(3.5, "lines"),
      panel.spacing.y = unit(1.2, "lines"),
      aspect.ratio = 0.48,
      axis.title.y = element_text(
        margin = margin(r = 10),
        size = 14
      ),
      axis.title.x = element_text(
        margin = margin(t = 10),
        size = 14
      ),
      strip.text.x = element_text(size = 8, face = "bold", margin = margin(b = 10)),
      strip.text.y = element_text(size = 8, face = "bold")
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
    xlab("Sample Size (N)") +
    ylab("Power")
  
  return(plot)
}







