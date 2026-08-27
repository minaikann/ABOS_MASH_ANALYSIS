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
