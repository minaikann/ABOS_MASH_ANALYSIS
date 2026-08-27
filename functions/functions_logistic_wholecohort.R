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
