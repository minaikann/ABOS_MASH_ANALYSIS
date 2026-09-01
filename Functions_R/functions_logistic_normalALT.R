# ============================================================
# functions_logistic_normalALT.R
# Part A: ALT alone vs metabolite alone vs ALT + metabolite,
# restricted to the normal-ALT cohort.
#   - stats helpers (OR/CI, ROC, CV)
#   - per-metabolite worker
#   - driver + selection rule
#   - result table builders
# ============================================================

# ---- stats helpers --------------------------------------------------

# Wald 95% CI on the odds-ratio scale
or_ci <- function(cm, term) {
  est <- cm[term, "Estimate"]
  se  <- cm[term, "Std. Error"]
  c(or = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
}

# ROC on fitted probabilities. direction is fixed, never "auto":
# with direction = "auto" pROC flips the comparison to force AUC > 0.5,
# so a protective marker would be silently reported as predictive.
roc_fit <- function(df, mod) {
  roc(df$mash_lab, predict(mod, type = "response"),
      levels = c("noMASH", "MASH"), direction = "<", quiet = TRUE)
}

# Class-stratified fold assignment: keeps the MASH:noMASH ratio in every
# fold, so no fold ends up with a single outcome class.
stratified_folds <- function(y, K) {
  fid <- integer(length(y))
  for (lv in levels(y)) {
    idx <- which(y == lv)
    fid[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
  }
  fid
}

# Repeated K-fold CV AUC.
#
# Predictions are POOLED across folds and scored once, rather than
# averaging fold-wise AUCs. Each fold holds ~n/K rows, so a fold-wise
# AUC is very noisy and their mean is not the AUC of the pooled
# predictions. `repeats` averages over several random fold splits to
# damp the dependence on any one split.
cv_auc_glm <- function(data, form, K = 10, repeats = 5, seed = 123) {
  set.seed(seed)
  y <- data$mash_lab
  rep_aucs <- rep(NA_real_, repeats)

  for (r in seq_len(repeats)) {
    fold_id <- stratified_folds(y, K)
    p_oof   <- rep(NA_real_, nrow(data))

    for (k in seq_len(K)) {
      tr <- data[fold_id != k, , drop = FALSE]
      te <- data[fold_id == k, , drop = FALSE]

      if (length(unique(tr$mash_lab)) < 2) next
      tr$sex <- droplevels(tr$sex)
      if (nlevels(tr$sex) < 2) next

      # a test row whose sex level never appeared in training cannot
      # be predicted; skip those rows rather than erroring
      keep <- te$sex %in% levels(tr$sex)
      if (!any(keep)) next

      mod <- tryCatch(glm(form, data = tr, family = binomial),
                      error = function(e) NULL)
      if (is.null(mod)) next

      p_oof[which(fold_id == k)[keep]] <-
        predict(mod, newdata = te[keep, , drop = FALSE], type = "response")
    }

    ok <- !is.na(p_oof)
    if (sum(ok) < 10 || length(unique(y[ok])) < 2) next

    rep_aucs[r] <- as.numeric(auc(
      roc(y[ok], p_oof[ok], levels = c("noMASH", "MASH"),
          direction = "<", quiet = TRUE)
    ))
  }

  mean(rep_aucs, na.rm = TRUE)
}

# ---- per-metabolite worker -------------------------------------------

alt_vs_one_met <- function(met_clean, data, sex_var = "sex",
                           K = 10, repeats = 5) {

  df <- data %>%
    mutate(met = .data[[met_clean]],
           sex = factor(.data[[sex_var]])) %>%
    filter(!is.na(met), !is.na(ALT), !is.na(sex), !is.na(mash_lab))

  # ---- guards -------------------------------------------------------
  if (nrow(df) < 20)                             return(NULL)  # too few
  if (length(unique(df$mash_lab)) < 2)           return(NULL)  # one class
  if (sd(df$met) == 0 || !is.finite(sd(df$met))) return(NULL)  # constant

  df$sex <- droplevels(df$sex)
  if (nlevels(df$sex) < 2)                       return(NULL)  # sex const

  df$mash_lab <- factor(df$mash_lab, levels = c("noMASH", "MASH"))
  df$met_z    <- as.numeric(scale(df$met))       # OR per 1 SD

  # ---- models (apparent fit, whole sample) --------------------------
  m_alt  <- glm(mash_lab ~ ALT + sex,         data = df, family = binomial)
  m_met  <- glm(mash_lab ~ met_z + sex,       data = df, family = binomial)
  m_full <- glm(mash_lab ~ ALT + met_z + sex, data = df, family = binomial)

  r_alt  <- roc_fit(df, m_alt)
  r_met  <- roc_fit(df, m_met)
  r_full <- roc_fit(df, m_full)

  # raw, unadjusted marker AUC (no model)
  r_raw <- roc(df$mash_lab, df$met_z,
               levels = c("noMASH", "MASH"), direction = "<", quiet = TRUE)

  auc_alt  <- as.numeric(auc(r_alt))
  auc_met  <- as.numeric(auc(r_met))
  auc_full <- as.numeric(auc(r_full))
  ci_met   <- as.numeric(ci.auc(r_met))          # lower, estimate, upper

  # ---- cross-validated AUCs -----------------------------------------
  # same seed => same fold splits for all three models and all
  # metabolites, so CV AUCs are directly comparable
  cv_alt  <- cv_auc_glm(df, mash_lab ~ ALT + sex,         K, repeats)
  cv_met  <- cv_auc_glm(df, mash_lab ~ met_z + sex,       K, repeats)
  cv_full <- cv_auc_glm(df, mash_lab ~ ALT + met_z + sex, K, repeats)

  # ---- coefficients -------------------------------------------------
  c_alt  <- summary(m_alt)$coefficients
  c_met  <- summary(m_met)$coefficients
  c_full <- summary(m_full)$coefficients

  met_adjsex    <- or_ci(c_met,  "met_z")        # met | sex
  met_adjsexalt <- or_ci(c_full, "met_z")        # met | sex + ALT

  # quasi-separation: huge SE means the MLE is running to infinity
  sep_flag <- max(c_met[ "met_z", "Std. Error"],
                  c_full["met_z", "Std. Error"]) > 5

  # ---- tests --------------------------------------------------------
  # nested -> LRT is the primary test for incremental value over ALT
  lrt_met_over_alt <- lrtest(m_alt, m_full)[2, "Pr(>Chisq)"]

  # non-nested, same data -> DeLong is appropriate here
  delong_met_vs_alt <- roc.test(r_alt, r_met, method = "delong")$p.value
  # nested -> DeLong is anticonservative; reference only
  delong_full_vs_alt <- roc.test(r_alt, r_full, method = "delong")$p.value

  data.frame(
    Metabolite = name_map[[met_clean]],
    n          = nrow(df),
    n_MASH     = sum(df$mash_lab == "MASH"),

    # --- apparent AUCs (in-sample, optimistic) ------------------------
    AUC_ALT_sex     = auc_alt,
    AUC_Met_sex     = auc_met,
    AUC_Met_lo      = ci_met[1],
    AUC_Met_hi      = ci_met[3],
    AUC_Met_raw     = as.numeric(auc(r_raw)),
    AUC_ALT_Met_sex = auc_full,

    # --- cross-validated AUCs (report these) --------------------------
    CV_AUC_ALT_sex     = cv_alt,
    CV_AUC_Met_sex     = cv_met,
    CV_AUC_ALT_Met_sex = cv_full,

    # --- metabolite effect, per 1 SD ----------------------------------
    Met_adjsex_OR = met_adjsex[["or"]],
    Met_adjsex_lo = met_adjsex[["lo"]],
    Met_adjsex_hi = met_adjsex[["hi"]],
    Met_adjsex_p  = c_met["met_z", "Pr(>|z|)"],

    Met_adjsexALT_OR = met_adjsexalt[["or"]],
    Met_adjsexALT_lo = met_adjsexalt[["lo"]],
    Met_adjsexALT_hi = met_adjsexalt[["hi"]],
    Met_adjsexALT_p  = c_full["met_z", "Pr(>|z|)"],

    # --- ALT effect, per 1 U/L ----------------------------------------
    ALT_adjsex_OR    = exp(c_alt[ "ALT", "Estimate"]),
    ALT_adjsex_p     =     c_alt[ "ALT", "Pr(>|z|)"],
    ALT_adjsexMet_OR = exp(c_full["ALT", "Estimate"]),
    ALT_adjsexMet_p  =     c_full["ALT", "Pr(>|z|)"],

    # --- model comparisons --------------------------------------------
    LRT_p_met_over_ALT   = lrt_met_over_alt,
    DeLong_p_met_vs_ALT  = delong_met_vs_alt,
    DeLong_p_full_vs_ALT = delong_full_vs_alt,

    separation_warning = sep_flag,
    row.names = NULL
  )
}

# ---- driver + selection rule ------------------------------------------

## Deltas are created here as real columns (select() cannot compute
## expressions), and FDR is applied across metabolites.
run_alt_vs_met <- function(data, sex_var = "sex", K = 10, repeats = 5) {
  stopifnot(is.data.frame(data))

  res_list <- lapply(clean_metabolites, alt_vs_one_met,
                     data = data, sex_var = sex_var, K = K, repeats = repeats)
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (!length(res_list))
    stop("No metabolite passed the n >= 20 / two-class / two-sex filter.")

  bad <- !vapply(res_list, is.data.frame, logical(1))
  if (any(bad))
    stop("alt_vs_one_met() returned a non-data.frame. ",
         "Check that `name_map` and `clean_metabolites` are not masked.")

  do.call(rbind, res_list) %>%
    mutate(
      # apparent deltas
      Delta_vs_ALT = AUC_Met_sex     - AUC_ALT_sex,
      Delta_added  = AUC_ALT_Met_sex - AUC_ALT_sex,

      # cross-validated deltas
      Delta_CV_vs_ALT = CV_AUC_Met_sex     - CV_AUC_ALT_sex,
      Delta_CV_added  = CV_AUC_ALT_Met_sex - CV_AUC_ALT_sex,

      # optimism: how much the apparent AUC overstates performance
      Optimism_Met = AUC_Met_sex - CV_AUC_Met_sex,

      # multiplicity control across metabolites
      Met_adjsex_FDR        = p.adjust(Met_adjsex_p,        method = "BH"),
      Met_adjsexALT_FDR     = p.adjust(Met_adjsexALT_p,     method = "BH"),
      LRT_FDR_met_over_ALT  = p.adjust(LRT_p_met_over_ALT,  method = "BH"),
      DeLong_FDR_met_vs_ALT = p.adjust(DeLong_p_met_vs_ALT, method = "BH")
    ) %>%
    arrange(desc(CV_AUC_Met_sex))
}

select_metabolites <- function(res, cv_auc_min = 0.70, fdr_max = 0.05) {
  res %>%
    filter(!separation_warning,
           Met_adjsex_FDR < fdr_max,
           CV_AUC_Met_sex > cv_auc_min,
           AUC_Met_lo     > 0.5) %>%   # apparent CI excludes chance
    arrange(desc(CV_AUC_Met_sex))
}

# ---- result table builders --------------------------------------------

# ---- Table 1: metabolite vs ALT, head to head ------------------------
tbl1 <- function(res) {
  datatable(
    res %>%
      arrange(desc(CV_AUC_Met_sex)) %>%
      select(
        Metabolite, n, n_MASH,
        `CV AUC ALT`       = CV_AUC_ALT_sex,
        `CV AUC met`       = CV_AUC_Met_sex,
        `Delta CV AUC`     = Delta_CV_vs_ALT,
        `Apparent AUC met` = AUC_Met_sex,
        `Optimism`         = Optimism_Met,
        `Met OR (per SD)`  = Met_adjsex_OR,
        `Met OR lo`        = Met_adjsex_lo,
        `Met OR hi`        = Met_adjsex_hi,
        `Met FDR`          = Met_adjsex_FDR,
        `DeLong FDR`       = DeLong_FDR_met_vs_ALT
      ),
    caption  = paste("Table 1. Metabolite vs ALT (CV AUC),",
                     "both adjusted for sex - NORMAL ALT"),
    rownames = FALSE,
    options  = list(scrollX = TRUE, pageLength = 25)
  ) %>%
    formatRound(c("CV AUC ALT", "CV AUC met", "Delta CV AUC",
                  "Apparent AUC met", "Optimism",
                  "Met OR (per SD)", "Met OR lo", "Met OR hi"),
                digits = 3) %>%
    formatSignif(c("Met FDR", "DeLong FDR"), digits = 3)
}

# ---- Table 2: ALT + metabolite vs ALT --------------------------------
tbl2 <- function(res) {
  datatable(
    res %>%
      arrange(desc(Delta_CV_added)) %>%
      select(
        Metabolite, n, n_MASH,
        `CV AUC ALT`       = CV_AUC_ALT_sex,
        `CV AUC ALT + met` = CV_AUC_ALT_Met_sex,
        `Delta CV AUC`     = Delta_CV_added,
        `Met OR (per SD, adj. ALT)` = Met_adjsexALT_OR,
        `Met OR lo`          = Met_adjsexALT_lo,
        `Met OR hi`          = Met_adjsexALT_hi,
        `Met FDR (adj. ALT)` = Met_adjsexALT_FDR,
        `LRT FDR`            = LRT_FDR_met_over_ALT
      ),
    caption  = paste("Table 2. ALT + metabolite vs ALT (CV AUC),",
                     "both adjusted for sex - NORMAL ALT"),
    rownames = FALSE,
    options  = list(scrollX = TRUE, pageLength = 25)
  ) %>%
    formatRound(c("CV AUC ALT", "CV AUC ALT + met", "Delta CV AUC",
                  "Met OR (per SD, adj. ALT)", "Met OR lo", "Met OR hi"),
                digits = 3) %>%
    formatSignif(c("Met FDR (adj. ALT)", "LRT FDR"), digits = 3)
}

# ---- Table 3: the shortlist ------------------------------------------
tbl3 <- function(res, cv_auc_min = 0.70, fdr_max = 0.05) {
  datatable(
    select_metabolites(res, cv_auc_min, fdr_max) %>%
      select(
        Metabolite, n, n_MASH,
        `CV AUC met`      = CV_AUC_Met_sex,
        `CV AUC ALT`      = CV_AUC_ALT_sex,
        `Delta CV AUC`    = Delta_CV_vs_ALT,
        `Met OR (per SD)` = Met_adjsex_OR,
        `Met OR lo`       = Met_adjsex_lo,
        `Met OR hi`       = Met_adjsex_hi,
        `Met FDR`         = Met_adjsex_FDR
      ),
    caption = sprintf(
      "Table 3. Selected metabolites (CV AUC > %.2f, Met FDR < %.2f) - NORMAL ALT",
      cv_auc_min, fdr_max),
    rownames = FALSE,
    options  = list(scrollX = TRUE)
  ) %>%
    formatRound(c("CV AUC met", "CV AUC ALT", "Delta CV AUC",
                  "Met OR (per SD)", "Met OR lo", "Met OR hi"),
                digits = 3) %>%
    formatSignif("Met FDR", digits = 3)
}
