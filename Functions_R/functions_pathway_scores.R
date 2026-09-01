# ============================================================
# functions_pathway_scores.R
# Part D: per-patient pathway scores (coherent subpathways,
# >=2 validated metabolites moving the same direction) and
# their ranked association with MASH.
# ============================================================

# ---- build a per-patient score for one subpathway ----------------------
# Uses ALL validated metabolites in that subpathway, Z-scored and aligned
# by the sign of their beta so that higher score = more MASH-associated.
# Expects `sig_altast` and `plasma_whole` to exist in the calling
# environment (built in the Part C / Part D driver chunks of the Rmd).
build_score <- function(sub) {
  mets  <- sig_altast$metabolite[sig_altast$Sub_pathway == sub]
  mets  <- intersect(mets, colnames(plasma_whole))
  betas <- sig_altast$beta[match(mets, sig_altast$metabolite)]
  Z <- scale(plasma_whole[, mets, drop = FALSE])   # Z-score each metabolite
  Z <- sweep(Z, 2, sign(betas), `*`)               # align by beta sign -> higher = MASH
  rowMeans(Z, na.rm = TRUE)                         # mean of ALL, per patient
}

# ---- regress each pathway score on MASH (adjusted), rank by significance ----
rank_pathways <- function(score_mat, covars) {
  res <- lapply(colnames(score_mat), function(sub) {
    df <- data.frame(score = score_mat[, sub], covars) |> na.omit()
    fit <- tryCatch(suppressWarnings(
      glm(mash ~ score + Age + sex_bin + BMI + ALT + AST, df, family = binomial())),
      error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    co <- summary(fit)$coefficients["score", ]
    data.frame(
      Sub_pathway = sub,
      n_metabolites = sum(sig_altast$Sub_pathway == sub),
      beta = co["Estimate"],
      OR   = exp(co["Estimate"]),
      CI_low  = exp(co["Estimate"] - 1.96*co["Std. Error"]),
      CI_high = exp(co["Estimate"] + 1.96*co["Std. Error"]),
      pvalue = co["Pr(>|z|)"]
    )
  })
  out <- bind_rows(res)
  # <50 pathways -> FWER (Holm)
  out$p_holm <- p.adjust(out$pvalue, method = "holm")
  out[order(out$pvalue), ]                        # rank by significance
}


# ============================================================
# Forest plot of pathway scores (rank_pathways() output).
#   ranked   : table with Sub_pathway, OR, CI_low, CI_high,
#              n_metabolites, and a ranking column
#   n_show   : how many pathways to plot
#   rank_by  : column used to pick and order them ("p_holm", "pvalue", ...)
#   label_by : optional column naming the rows (default "Sub_pathway")
# ============================================================
forestPathwayFct <- function(ranked,
                             n_show   = 10,
                             rank_by  = "p_holm",
                             label_by = "Sub_pathway",
                             title    = "Sub-pathway scores associated with MASH",
                             subtitle = NULL,
                             point_col = "#185FA5",
                             bar_col   = "#4a4a4a",
                             size_name = "Metabolites\nin score",
                             pad = 1.6) {
  
  d <- ranked %>%
    arrange(.data[[rank_by]]) %>%
    slice_head(n = n_show)
  
  if (nrow(d) == 0) {
    warning("Nothing to plot."); return(invisible(NULL))
  }
  
  if (is.null(subtitle))
    subtitle <- sprintf("Top %d by significance \u00b7 adjusted for age, sex, BMI, ALT, AST",
                        nrow(d))
  
  d <- d %>%
    mutate(
      .row   = factor(.data[[label_by]], levels = rev(.data[[label_by]])),
      .label = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_low, CI_high)
    )
  
  x_max <- max(d$CI_high) * pad      # headroom on the right for the OR labels
  x_min <- min(d$CI_low) * 0.9
  
  ggplot(d, aes(x = OR, y = .row)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                   height = 0.18, colour = bar_col, linewidth = 0.6) +
    geom_point(aes(size = n_metabolites), colour = point_col, shape = 18) +
    geom_text(aes(x = x_max, label = .label), hjust = 1, size = 3.1, colour = "grey25") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8),
                  expand = expansion(mult = c(0.02, 0.02))) +
    scale_size_continuous(range = c(2.5, 6), name = size_name) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(x = "Odds ratio per SD of pathway score (95% CI, log scale)",
         y = NULL, title = title, subtitle = subtitle) +
    theme_minimal(base_size = 11, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(face = "bold", size = 13, margin = ggplot2::margin(b = 2)),
      plot.subtitle = element_text(colour = "grey40", size = 9.5, margin = ggplot2::margin(b = 12)),
      axis.title.x  = element_text(size = 10, margin = ggplot2::margin(t = 10), colour = "grey25"),
      axis.text.y   = element_text(size = 10, colour = "grey15"),
      axis.text.x   = element_text(size = 9,  colour = "grey35"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      legend.position = "bottom",
      legend.title    = element_text(size = 9, colour = "grey30"),
      legend.text     = element_text(size = 9),
      plot.margin     = ggplot2::margin(12, 20, 12, 12)
    )
}