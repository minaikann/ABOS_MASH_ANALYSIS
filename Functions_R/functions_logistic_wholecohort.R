# ============================================================
# functions_logistic_wholecohort.R
# Part C: validation of candidate metabolites in the whole
# biopsy cohort (mash ~ metabolite + Age + sex_bin + BMI + ALT [+ AST]).
# ============================================================

# ---- per-metabolite logistic engine (metabolite standardised -> OR per SD) ----
run_logistic <- function(plasma_mat, covars, adjust, min_n = 30) {
  form <- as.formula(paste("mash ~", paste(c("metabolite", adjust), collapse = " + ")))
  res <- list(); k <- 1
  for (met in colnames(plasma_mat)) {
    df <- data.frame(metabolite = as.numeric(plasma_mat[, met]), covars) |> na.omit()
    if (nrow(df) < min_n || length(unique(df$mash)) < 2 || sd(df$metabolite) == 0) next
    df$metabolite <- as.numeric(scale(df$metabolite))
    fit <- tryCatch(suppressWarnings(glm(form, data = df, family = binomial())),
                    error = function(e) NULL)
    if (is.null(fit)) next
    co <- summary(fit)$coefficients
    if (!"metabolite" %in% rownames(co)) next
    beta <- co["metabolite","Estimate"]; se <- co["metabolite","Std. Error"]
    res[[k]] <- data.frame(
      metabolite = met,
      beta = beta, OR = exp(beta),
      CI_low = exp(beta - 1.96*se), CI_high = exp(beta + 1.96*se),
      pvalue = co["metabolite","Pr(>|z|)"],
      se_max = max(co[,"Std. Error"], na.rm = TRUE),
      n = nrow(df)
    ); k <- k + 1
  }
  out <- bind_rows(res)
  out$FDR <- p.adjust(out$pvalue, method = "BH")
  out[order(out$FDR), ]
}

# ---- annotate significant hits with pathways + direction, 3 dp ----
annotate_tbl <- function(res, chemical_details) {
  res %>%
    filter(FDR <= 0.05, se_max <= 5) %>%
    mutate(direction = ifelse(OR > 1, "\u2191", "\u2193")) %>%
    left_join(
      chemical_details %>%
        select(metabolite    = CHEMICAL_NAME,
               Sub_pathway   = SUB_PATHWAY,
               Super_pathway = SUPER_PATHWAY),
      by = "metabolite"
    ) %>%
    transmute(metabolite, Super_pathway, Sub_pathway, direction,
              beta = round(beta, 3), OR = round(OR, 3),
              `95% CI` = sprintf("%.3f-%.3f", CI_low, CI_high),
              p = signif(pvalue, 3), FDR = signif(FDR, 3), n)
}


# ============================================================
# Volcano plot for run_logistic() output, styled like
# printResultsLimmaFct() (x = log2 OR per SD, y = -log10 raw p,
# cut-offs on FDR).
# ============================================================
volcanoLogitFct = function(resLogit, seuilAdjP = 0.05, seuilLog2OR = 0.26,
                           color = c("blue", "grey", "red"), legendPos = "right",
                           genesToLabel = NULL, seuilSE = 5){
  
  results = as.data.frame(resLogit)
  results = results[order(results$FDR), ]
  row.names(results) = results$metabolite
  
  # near-separation fits: enormous OR with an enormous SE
  nUnstable = length(which(results$se_max > seuilSE))
  if (nUnstable > 0){
    cat("Removed", nUnstable, "metabolite(s) with unstable fits (se_max >", seuilSE, ").\n")
    results = results[results$se_max <= seuilSE, ]
  }
  
  results$log2OR = log2(results$OR)
  
  nSign   = length(which(results$FDR <= seuilAdjP))
  nSignOR = length(which(results$FDR <= seuilAdjP & abs(results$log2OR) >= seuilLog2OR))
  
  cat("There are", nSign, "metabolites associated with MASH considering threshold", seuilAdjP,
      "on FDR, and", nSignOR, "when adding the threshold", seuilLog2OR,
      "on absolute value of log2 Odds Ratio.", "\n")
  
  if (nSign > 0){
    pval_safe = max(results$pvalue[nSign], 1e-20)
    seuilpval = -log10(pval_safe)
  } else {
    seuilpval = -log10(seuilAdjP/nrow(results))
  }
  
  highlight_genes <- c("1-carboxyethyltyrosine","1-carboxyethylleucine",
                       "1-carboxyethylphenylalanine","1-carboxyethylvaline",
                       "1-carboxyethylisoleucine","taurochenodeoxycholate",
                       "taurocholate","cholate","glycocholate")
  
  results$name <- ifelse(row.names(results) %in% highlight_genes, row.names(results), "")
  
  useAutoLabels = TRUE
  
  if (!is.null(genesToLabel)){
    matched = row.names(results) %in% genesToLabel
    if (any(matched)){
      results$name = ifelse(matched, row.names(results), "")
      useAutoLabels = FALSE
    } else {
      cat("None of the genesToLabel were found in the results — falling back to top-5-up/top-5-down labeling.\n")
    }
  }
  
  if (useAutoLabels) {
    cptPos = 0; cptNeg = 0; cpt = 1
    while ((cptPos < 5 | cptNeg < 5) & results[cpt,"FDR"] < seuilAdjP & cpt < nrow(results)){
      if (results[cpt,"log2OR"] < -seuilLog2OR & cptNeg < 5){
        results[cpt, "name"] = row.names(results)[cpt]; cptNeg = cptNeg + 1
      } else if (results[cpt,"log2OR"] > seuilLog2OR & cptPos < 5){
        results[cpt, "name"] = row.names(results)[cpt]; cptPos = cptPos + 1
      }
      cpt = cpt + 1
    }
  }
  
  results$Expression <- ifelse(results$log2OR > seuilLog2OR & results$FDR < seuilAdjP,
                               "Up regulated",
                               ifelse(results$log2OR < -seuilLog2OR & results$FDR < seuilAdjP,
                                      "Down regulated",
                                      "NS"))
  
  ggplot(results, aes(x = log2OR, y = -log10(pvalue), col = Expression)) +
    geom_point() +
    geom_hline(yintercept = seuilpval, color = "grey60", linetype = "dashed") +
    geom_vline(xintercept =  seuilLog2OR, color = "grey60", linetype = "dashed") +
    geom_vline(xintercept = -seuilLog2OR, color = "grey60", linetype = "dashed") +
    scale_color_manual(values = c("Down regulated" = color[1],
                                  "NS"             = color[2],
                                  "Up regulated"   = color[3])) +
    theme_classic() +
    geom_label_repel(data = subset(results, name != ""), aes(label = name),
                     fontface = "italic", size = 3, max.overlaps = 8) +
    xlab("log2 Odds Ratio (per SD)") +
    theme(axis.text   = element_text(size = 12),
          axis.title  = element_text(size = 14),
          legend.position = legendPos,
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 9))
}

# ============================================================
# Forest plot for run_logistic() output.
#   Plots OR (per SD of the metabolite) with 95% CI on a log scale.
#   res      : run_logistic() output, or the subset you want to show
#   n_show   : how many metabolites to plot (ranked by FDR)
#   group_by : optional column to facet on, e.g. "Sub_pathway" or
#              "Super_pathway" - requires an annotated table
#   color    : c(down, up) - OR < 1 and OR > 1
# ============================================================
forestLogitFct <- function(res, title = NULL, subtitle = NULL,
                           n_show   = 20,
                           seuilAdjP = 0.05,
                           seuilSE  = 5,
                           group_by = NULL,
                           color    = c("#2166AC", "#B2182B"),
                           chemical_details = NULL) {
  
  d <- as.data.frame(res)
  
  # unstable fits would blow the x axis apart
  if ("se_max" %in% names(d)) d <- d[!is.na(d$se_max) & d$se_max <= seuilSE, ]
  
  d <- d %>%
    filter(!is.na(OR), OR > 0, FDR <= seuilAdjP) %>%
    arrange(FDR) %>%
    slice_head(n = n_show)
  
  if (nrow(d) == 0) {
    warning("No metabolites pass FDR <= ", seuilAdjP, "; nothing to plot.")
    return(invisible(NULL))
  }
  
  # attach pathways if the caller passed the Metabolon table and wants facets
  if (!is.null(group_by) && !group_by %in% names(d) && !is.null(chemical_details)) {
    d <- d %>%
      left_join(
        chemical_details %>%
          select(metabolite    = CHEMICAL_NAME,
                 Sub_pathway   = SUB_PATHWAY,
                 Super_pathway = SUPER_PATHWAY),
        by = "metabolite"
      )
  }
  
  d <- d %>%
    mutate(
      direction  = ifelse(OR > 1, "Higher in MASH", "Lower in MASH"),
      metabolite = factor(metabolite, levels = rev(unique(metabolite))),
      label      = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_low, CI_high)
    )
  
  x_max <- max(d$CI_high) * 1.9   # headroom on the right for the OR labels
  x_min <- min(d$CI_low) * 0.9
  
  p <- ggplot(d, aes(x = OR, y = metabolite, colour = direction)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey60", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                   height = 0.18, linewidth = 0.6) +
    geom_point(shape = 18, size = 3.2) +
    geom_text(aes(x = x_max, label = label), hjust = 1,
              size = 3.1, colour = "grey25", show.legend = FALSE) +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
    scale_colour_manual(values = c("Lower in MASH"  = color[1],
                                   "Higher in MASH" = color[2]),
                        name = NULL) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    labs(x = "Odds ratio per SD (95% CI, log scale)", y = NULL,
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13,
                                   margin = ggplot2::margin(b = 2)),
      plot.subtitle = element_text(colour = "grey40", size = 9.5,
                                   margin = ggplot2::margin(b = 12)),
      axis.title.x  = element_text(size = 10, colour = "grey25",
                                   margin = ggplot2::margin(t = 10)),
      axis.text.y   = element_text(size = 9.5, colour = "grey15"),
      axis.text.x   = element_text(size = 9,  colour = "grey35"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
      legend.position = "bottom",
      plot.margin     = ggplot2::margin(12, 20, 12, 12)
    )
  
  if (!is.null(group_by) && group_by %in% names(d)) {
    p <- p + facet_grid(rows = vars(.data[[group_by]]),
                        scales = "free_y", space = "free_y",
                        switch = "y") +
      theme(strip.text.y.left = element_text(angle = 0, hjust = 1,
                                             size = 8.5, colour = "grey30"),
            strip.placement = "outside")
  }
  
  p
}