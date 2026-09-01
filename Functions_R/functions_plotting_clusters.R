# ============================================================
# functions_plotting_clusters.R
# Part B: plotting helpers for the CTRL / CM / LS boxplots
# with pairwise Wilcoxon significance brackets.
# ============================================================

stars_from_p <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}

# Plot builder. Expects `long`, `ann`, `grp_levels`, and `tick` to exist
# in the calling environment (built in the Part B driver chunk of the Rmd).
make_plot <- function(met) {
  d <- filter(long, Metabolite == met)
  a <- filter(ann,  Metabolite == met)

  p <- ggplot(d, aes(class, value)) +
    geom_boxplot(
      aes(fill = class),
      width = 0.6, linewidth = 0.4, alpha = 0.8,
      outlier.shape = 16, outlier.size = 1.2,
      outlier.colour = "black", outlier.alpha = 0.6
    ) +
    scale_fill_manual(
      values = c(CTRL = "grey80", CM = "#D62728", LS = "#00BFC4")
    ) +
    scale_y_continuous(
      breaks = pretty_breaks(n = 6),
      labels = function(x) sprintf("%.1f", x),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    labs(
      x = NULL, y = "log2-metabolite",
      title    = met,
      subtitle = "Pairwise Wilcoxon rank-sum tests, FWER-adjusted"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      aspect.ratio    = 1.2,
      axis.text       = element_text(size = 9,  colour = "black"),
      axis.title      = element_text(size = 10, colour = "black"),
      axis.line       = element_line(linewidth = 0.4,  colour = "black"),
      axis.ticks      = element_line(linewidth = 0.35, colour = "black"),
      plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle   = element_text(size = 10, hjust = 0.5)
    )

  if (nrow(a) > 0) {
    p <- p +
      geom_segment(data = a, aes(x = x1, xend = x2, y = bar_y, yend = bar_y),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "black") +
      geom_segment(data = a, aes(x = x1, xend = x1,
                                 y = bar_y - tick * y_range, yend = bar_y),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "black") +
      geom_segment(data = a, aes(x = x2, xend = x2,
                                 y = bar_y - tick * y_range, yend = bar_y),
                   inherit.aes = FALSE, linewidth = 0.35, colour = "black") +
      geom_text(data = a, aes(x = (x1 + x2) / 2, y = text_y, label = stars),
                inherit.aes = FALSE, size = 3.2, colour = "black")
  }
  p
}
