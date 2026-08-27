---
title: "Metabolites that identify MASH with normal ALT"
author: "Deborah Mina Ikann"
date: "2026-07-07"


```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo      = TRUE,
  message   = FALSE,
  warning   = FALSE,
  error     = FALSE,
  fig.align = "center"
)
```

# Question

1. **Can plasma metabolites identify MASH when ALT is normal?**
   (MASH vs no-MASH, restricted to normal-ALT patients — **Part A**)

2. **How are these metabolites distributed across the three groups?**
   (Distribution of the metabolites in CTRL, CM and LS — **Part B**)

3. **Are these metabolites also present in the whole cohort?**
   (MASH vs no-MASH, applied to the whole cohort — **Part C**)

4. **Are the biological programs identified in plasma also dysregulated in the liver?**
   (MASH vs no-MASH, applied to the whole cohort — **Part D**)

5. **Which genes are dysregulated with MASH in the liver when ALT is normal?**
   (MASH vs no-MASH, applied to normal-ALT patients — **Part E**)

6. **Which liver metabolites are dysregulated with MASH when ALT is normal?**
   (MASH vs no-MASH, applied to normal-ALT patients — **Part F**)

# Setup

## Libraries

All packages are loaded once here, so no chunk further down needs its own
`library()` call.

```{r libraries}
# Data I/O
library(openxlsx)
library(readxl)

# Wrangling
library(dplyr)
library(tidyr)
library(stringr)

# Differential expression / enrichment
library(limma)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)

# Plotting
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(scales)

# Modelling / evaluation
library(pROC)
library(lmtest)
library(ResourceSelection)

# Tables
library(gtsummary)
library(gt)
library(DT)
library(knitr)

# Flowchart
library(DiagrammeR)
```

## Paths and constants

```{r paths}
data_dir <- "~/Documents/Clustering_ABOS/codes/Data"
fun_dir  <- "~/Documents/Clustering_ABOS/codes/functions"

# Sex coding used throughout (defined once, reused everywhere)
female_codes <- c("F", "Female", "FEMALE", "f", "female", "Femme")
male_codes   <- c("M", "Male",   "MALE",   "m", "male",   "Homme")

# Normal ALT window
alt_lo <- 0
alt_hi <- 35
```

## Data

```{r load-data}
clinical_data    <- read.xlsx(file.path(data_dir, "clinical_data_Mash.xlsx"), rowNames = TRUE)
chemical_details <- read_xlsx(file.path(data_dir, "Debora.xlsx"), sheet = 3)
plasma_data      <- read.xlsx(file.path(data_dir, "plasma_bio_transposed.xlsx"), rowNames = TRUE)
cluster_data     <- read.csv(file.path(data_dir, "cluster_data_mash.csv"), row.names = 1)
transcript_data  <- read.xlsx(file.path(data_dir, "Transcripto_matrix.xlsx"), rowNames = TRUE)
liver_data       <- read.xlsx(file.path(data_dir, "liver_bio_transposed.xlsx"), rowNames = TRUE)
l_chem_details   <- read.xlsx(file.path(data_dir, "Chemical_details_liver.xlsx"))
```

## Helper functions

```{r source-functions}
# Pre-existing helper functions (limma comparison / plotting)
# NOTE: the file name contains a space - kept exactly as on disk.
source(file.path(fun_dir, "functions_ liver_metabo.R"))

# Project function files (each covers one part of the analysis)
source("functions_data_prep.R")            # chemical-name mapping used during data loading
source("functions_logistic_normalALT.R")   # Part A: ALT vs metabolite logistic / CV helpers
source("functions_plotting_clusters.R")    # Part B: cluster boxplot + significance-bracket helpers
source("functions_logistic_wholecohort.R") # Part C: whole-cohort logistic engine + annotation
source("functions_pathway_scores.R")       # Part C: pathway scoring + ranking
source("GSEA_function.R")                  # Parts D/E: gene set enrichment analysis
```

## Name mapping and alignment

```{r map-names}
# Replace COMP_ID row names with human-readable chemical names where available
plasma_data <- map_chemical_names(plasma_data, chemical_details)
liver_data  <- map_metabolite_names(liver_data, l_chem_details)
```

```{r align}
# Patients present in BOTH clinical and plasma tables
common_samples <- intersect(rownames(clinical_data), colnames(plasma_data))
clinical_sub   <- clinical_data[common_samples, , drop = FALSE]
plasma_sub     <- plasma_data[, common_samples, drop = FALSE]   # metabolites x samples

cat("Common (analysed) patients:", length(common_samples), "\n")
```

## Working dataset (built once)

```{r work-table}
# Normal ALT threshold (upper limit only)
alt_hi <- 35

work <- data.frame(
  patient = common_samples,
  mash    = clinical_sub$Mash,
  sex     = clinical_sub$sexe,
  ALT     = clinical_sub$TGPUL,
  Age     = clinical_sub$AgeJourIntervention,
  BMI     = clinical_sub$BMI,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(mash), !is.na(ALT)) %>%
  mutate(
    sex_label = case_when(
      as.character(sex) %in% female_codes ~ "Female",
      as.character(sex) %in% male_codes   ~ "Male",
      TRUE ~ NA_character_
    ),
    alt_normal = ALT <= alt_hi,
    alt_high   = ALT >  alt_hi,
    mash_lab   = factor(ifelse(mash == 1, "MASH", "noMASH"),
                        levels = c("noMASH", "MASH")),
    alt_level  = ifelse(alt_normal, "normal", "high")
  ) %>%
  filter(!is.na(sex_label))

work_normal_alt <- work %>% filter(alt_level == "normal")

cat("Patients with normal ALT (<= ", alt_hi, "): ",
    nrow(work_normal_alt), "\n", sep = "")
print(table(work_normal_alt$mash_lab, useNA = "ifany"))
```

# Part A — Normal ALT

## A1. Flowchart

```{r flow-counts}
n_total  <- nrow(work)
n_high   <- sum(work$alt_level == "high",   na.rm = TRUE)
n_normal <- sum(work$alt_level == "normal", na.rm = TRUE)

norm_df    <- work %>% filter(alt_level == "normal")
n_mash     <- sum(norm_df$mash_lab == "MASH")
n_nomash   <- sum(norm_df$mash_lab == "noMASH")
n_mash_f   <- sum(norm_df$mash_lab == "MASH"   & norm_df$sex_label == "Female")
n_mash_m   <- sum(norm_df$mash_lab == "MASH"   & norm_df$sex_label == "Male")
n_nomash_f <- sum(norm_df$mash_lab == "noMASH" & norm_df$sex_label == "Female")
n_nomash_m <- sum(norm_df$mash_lab == "noMASH" & norm_df$sex_label == "Male")

cat("Total:", n_total,
    "| normal ALT:", n_normal,
    "| high ALT (excluded):", n_high, "\n")
```

```{r flowchart}
grViz(sprintf("
digraph cohort {
  graph [rankdir = TB, splines = ortho, nodesep = 0.6]
  node  [shape = box, style = filled, fillcolor = '#f7f7f7',
         fontname = Helvetica, fontsize = 11, width = 2.6]
  total     [label = 'Patients with plasma + clinical data\\n(n = %d)']
  excl_high [label = 'Excluded: high ALT (>35)\\n(n = %d)',
             fillcolor = '#ffffff', style = 'filled,dashed', width = 2.0]
  normal    [label = 'Normal ALT (<= 35)\\n(n = %d)', fillcolor = '#ffffff']
  mash      [label = 'MASH\\n(n = %d)\\nFemale %d / Male %d', fillcolor = '#ffffff']
  nomash    [label = 'no-MASH\\n(n = %d)\\nFemale %d / Male %d', fillcolor = '#ffffff']
  # invisible node keeps the exclusion arrow horizontal
  node [shape = point, width = 0.01, height = 0.01, label = '']
  gap
  edge [arrowhead = normal, color = '#555555']
  total -> gap [arrowhead = none]
  gap   -> normal
  gap   -> excl_high [style = dashed]
  normal -> mash
  normal -> nomash
  { rank = same; gap; excl_high }
}
", n_total, n_high, n_normal,
   n_mash,   n_mash_f,   n_mash_m,
   n_nomash, n_nomash_f, n_nomash_m))
```

## A2. Summary (normal ALT)

Group comparison test: Wilcoxon for continuous variables, Fisher for
categorical ones, printed via `kable()` rather than an interactive widget.

```{r cohort-table-normal}
clinical_table <- work %>%
  filter(alt_level == "normal") %>%
  select(mash_lab, sex_label, Age, BMI, ALT) %>%
  tbl_summary(
    by   = mash_lab,
    type = all_continuous() ~ "continuous2",
    statistic = list(
      all_continuous()  ~ c("{median} ({p25} - {p75})", "{min} - {max}"),
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      sex_label ~ "sex",
      Age       ~ "Age",
      BMI       ~ "BMI",
      ALT       ~ "ALT"
    ),
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test",
      all_categorical() ~ "fisher.test"
    )
  ) %>%
  modify_header(label = "**Characteristic**")

clinical_table %>%
  as_kable(caption = "Normal-ALT cohort: MASH vs no-MASH")
```

## A3. Differential abundance — limma (MASH vs no-MASH)

```{r limma-normal}
# Metadata + metabolite matrix (metabolites x samples).
#
# NOTE ON THE COHORT USED HERE:
# the filter below removes only HIGH-ALT patients, so low-ALT patients (<10)
# are still included. That is why the group sizes printed here are slightly
# larger than the normal-ALT counts from Part A1.
# To restrict strictly to normal ALT, replace the filter line with:
#   m <- m[m$alt_level == "normal", ]
m <- work
rownames(m) <- m$patient
m <- m[!(m$alt_high %in% TRUE), ]     # NA alt_high treated as "not high"
m$grp <- factor(m$mash_lab, levels = c("noMASH", "MASH"))

cat("Group sizes entering limma:\n")
print(table(m$grp, useNA = "ifany"))

x <- as.matrix(plasma_sub[, rownames(m), drop = FALSE])

design <- model.matrix(~ -1 + grp, data = m)
fit    <- lmFit(x, design)
comp   <- comparisonsLimmaFct(fit, "grpMASH - grpnoMASH", design, nrow(x))

res_alt_norm <- printResultsLimmaFct(comp, topPrintHist = FALSE, topPrintVolc = FALSE)

sig_alt_norm <- rownames(
  res_alt_norm$results[
    !is.na(res_alt_norm$results$adj.P.Val) &
      res_alt_norm$results$adj.P.Val < 0.05, , drop = FALSE]
)
cat("MASH vs noMASH, NORMAL ALT — significant (FDR < 0.05):",
    length(sig_alt_norm), "\n")

datatable(
  res_alt_norm$results[order(res_alt_norm$results$logFC, decreasing = TRUE), ],
  caption = "MASH vs no-MASH within normal ALT",
  options = list(scrollX = TRUE)
)

res_alt_norm$volcano
```

### Metabolon annotation

```{r annotate-normal}
sig <- res_alt_norm$results
idx <- match(rownames(sig), chemical_details$CHEMICAL_NAME)
cat("Unannotated:", sum(is.na(idx)), "/", nrow(sig), "\n")

sig$Super_pathway <- chemical_details$SUPER_PATHWAY[idx]
sig$Sub_pathway   <- chemical_details$SUB_PATHWAY[idx]
sig$metabolite    <- rownames(sig)

sig_show <- sig[, c("Super_pathway", "Sub_pathway", "logFC", "FC", "adj.P.Val")]

datatable(sig_show,
          caption = "MASH vs no-MASH, normal ALT — Metabolon annotation",
          options = list(scrollX = TRUE))
```

### Summary by sub-pathway

```{r subpathway-summary-normal}
# Sorted by |logFC| first, so `Strongest` really lists the strongest hits.
sub_sum <- sig %>%
  arrange(desc(abs(logFC))) %>%
  group_by(Super_pathway, Sub_pathway) %>%
  summarise(
    No_Sig    = n(),
    No_Up     = sum(logFC > 0),
    No_Down   = sum(logFC < 0),
    Direction = case_when(
      all(logFC > 0) ~ "All \u2191",
      all(logFC < 0) ~ "All \u2193",
      TRUE ~ paste0(sum(logFC > 0), "\u2191 / ", sum(logFC < 0), "\u2193")
    ),
    Strongest = paste(head(metabolite, 3), collapse = "; "),
    .groups   = "drop"
  ) %>%
  arrange(desc(No_Sig))

datatable(sub_sum, options = list(scrollX = TRUE),
          caption = "Significant metabolites by Metabolon sub-pathway")
```

## A4. Logistic regression: ALT alone vs metabolite alone vs ALT + metabolite

### Data prep (RINT)

```{r logistic-data-prep}
## RINT is applied to every metabolite BEFORE any outcome is examined,
## so it cannot leak outcome information.

work_aligned   <- work
plasma_aligned <- as.data.frame(t(plasma_sub))          # samples x metabolites
plasma_aligned <- plasma_aligned[work_aligned$patient, , drop = FALSE]

# Rank-based inverse normal transform
rint <- function(x) {
  n <- sum(!is.na(x))
  if (n < 3) return(x)
  qnorm((rank(x, na.last = "keep") - 0.5) / n)
}
plasma_aligned <- as.data.frame(lapply(plasma_aligned, rint))
rownames(plasma_aligned) <- work_aligned$patient        # lapply() drops row names

# Metabolite names: clean for modelling, keep a map back to the originals
orig_metabolites  <- rownames(plasma_sub)
clean_metabolites <- make.names(orig_metabolites)
colnames(plasma_aligned) <- clean_metabolites
name_map <- setNames(orig_metabolites, clean_metabolites)

# Modelling frames: built AFTER the transform
df_all <- work_aligned %>%
  select(patient, mash_lab, ALT, alt_level, sex) %>%
  cbind(plasma_aligned[work_aligned$patient, , drop = FALSE]) %>%
  filter(!is.na(mash_lab), !is.na(ALT))

df_norm <- df_all %>% filter(alt_level == "normal")

# Verify the transform reached the modelling frame
stopifnot(max(abs(df_norm[[clean_metabolites[1]]]), na.rm = TRUE) < 4)

cat("All patients:", nrow(df_all),
    "| Normal-ALT patients:", nrow(df_norm), "\n")
cat("Max |z| in df_norm:",
    round(max(abs(as.matrix(df_norm[, clean_metabolites])), na.rm = TRUE), 2), "\n")
```

### Run

Currently disabled (`eval=FALSE`). Set `eval=TRUE` on this chunk to run it.

```{r logistic-run, eval=FALSE}
res_norm <- run_alt_vs_met(df_norm)

if (any(res_norm$separation_warning)) {
  warning(sum(res_norm$separation_warning),
          " metabolite(s) show quasi-separation (SE > 5); ORs unstable.")
}

tbl1(res_norm)
tbl2(res_norm)
tbl3(res_norm, cv_auc_min = 0.70, fdr_max = 0.05)
```

# Part B — Selected metabolites across clusters (CTRL vs CM vs LS)

Selected metabolites compared across the CTRL / CM / LS clusters: pairwise
Wilcoxon rank-sum tests (Holm-adjusted within each metabolite), shown as
per-metabolite boxplots with significance brackets.

Plotting helpers (`stars_from_p`, `make_plot`) live in
`functions_plotting_clusters.R`, sourced above.

## B1. Metabolites of interest

```{r cluster-metabolite-list}
met_list <- c(
  "1-carboxyethylvaline",
  "glycochenodeoxycholate glucuronide (1)",
  "palmitoyl-arachidonoyl-glycerol (16:0/20:4) [2]*",
  "palmitoyl-arachidonoyl-glycerol (16:0/20:4) [1]*",
  "homocitrulline",
  "cholic acid glucuronide",
  "palmitoyl-oleoyl-glycerol (16:0/18:1) [2]*",
  "palmitoyl-oleoyl-glycerol (16:0/18:1) [1]*",
  "deoxycholic acid glucuronide",
  "1-carboxyethylisoleucine",
  "1,5-anhydroglucitol (1,5-AG)"
)

missing_mets <- setdiff(met_list, rownames(plasma_data))
if (length(missing_mets) > 0) {
  stop("Not in rownames(plasma_data): ", paste(missing_mets, collapse = ", "))
}
```

## B2. Long-format table

```{r cluster-long}
plasma_t        <- as.data.frame(t(plasma_data[met_list, , drop = FALSE]))
plasma_t$sample <- rownames(plasma_t)

cluster_data$sample <- rownames(cluster_data)

# Column of cluster_data holding CTRL / CM / LS
class_col <- "class"
stopifnot(class_col %in% colnames(cluster_data))

long <- plasma_t %>%
  left_join(
    cluster_data %>% select(sample, class = all_of(class_col)),
    by = "sample"
  ) %>%
  filter(!is.na(class)) %>%
  pivot_longer(
    cols      = all_of(met_list),
    names_to  = "Metabolite",
    values_to = "value"
  ) %>%
  mutate(class = factor(class, levels = c("CTRL", "CM", "LS")))

cat("long built:", nrow(long), "rows |",
    dplyr::n_distinct(long$sample), "samples |",
    dplyr::n_distinct(long$Metabolite), "metabolites\n")
print(table(unique(long[, c("sample", "class")])$class))   # group sizes

grp_levels <- levels(long$class)
```

## B3. Pairwise tests and bracket positions

```{r cluster-tests}
pw <- long %>%
  group_by(Metabolite) %>%
  group_modify(~ {
    res <- pairwise.wilcox.test(.x$value, .x$class,
                                p.adjust.method = "holm", exact = FALSE)
    as.data.frame(as.table(res$p.value)) %>%
      setNames(c("g2", "g1", "p_adj")) %>%
      filter(!is.na(p_adj))
  }) %>%
  ungroup() %>%
  mutate(
    g1    = as.character(g1),
    g2    = as.character(g2),
    stars = stars_from_p(p_adj),
    x1    = match(g1, grp_levels),
    x2    = match(g2, grp_levels)
  ) %>%
  mutate(                                   # ensure x1 < x2
    lo = pmin(x1, x2), hi = pmax(x1, x2),
    x1 = lo, x2 = hi
  ) %>%
  select(-lo, -hi)

# Drop non-significant brackets (comment out to keep them)
pw <- filter(pw, stars != "ns")

# Bracket heights, anchored to the plotted data
rng <- long %>%
  group_by(Metabolite) %>%
  summarise(
    ymin    = min(value, na.rm = TRUE),
    ymax    = max(value, na.rm = TRUE),
    y_range = ymax - ymin,
    .groups = "drop"
  )

ann <- pw %>%
  left_join(rng, by = "Metabolite") %>%
  group_by(Metabolite) %>%
  arrange(x2 - x1, x1, .by_group = TRUE) %>%   # shortest span sits lowest
  mutate(
    bar_y  = ymax + y_range * (0.06 + 0.09 * (row_number() - 1)),
    text_y = bar_y + y_range * 0.015
  ) %>%
  ungroup()

tick <- 0.015   # bracket end-tick length, as a fraction of y_range
```

## B4. One plot per metabolite

```{r cluster-boxplots, fig.width=5.5, fig.height=5.5}
mets  <- sort(unique(long$Metabolite))
plots <- lapply(mets, make_plot)

# pdf("metabolites_per_page.pdf", width = 5.5, height = 5.5)
for (p in plots) print(p)
# dev.off()
```

# Part C — Validation in the whole biopsy cohort

## C1. Data prep

Helper functions for this section (`run_logistic`, `annotate_tbl`) live in
`functions_logistic_wholecohort.R`, sourced above.

```{r whole-cohort-prep}
### mash ~ metabolite + Age + sex_bin + BMI + ALT (+ AST)

plasma_whole <- t(plasma_sub)
pts          <- intersect(rownames(plasma_whole), rownames(clinical_sub))
plasma_whole <- plasma_whole[pts, , drop = FALSE]

covars <- data.frame(
  Age     = clinical_sub[pts, "AgeJourIntervention"],
  sex_bin = ifelse(as.character(clinical_sub[pts, "sexe"]) %in% female_codes, 0, 1),
  BMI     = clinical_sub[pts, "BMI"],
  ALT     = clinical_sub[pts, "TGPUL"],
  AST     = clinical_sub[pts, "AST"],
  mash    = clinical_sub[pts, "Mash"]
)

cat("Patients entering the whole-cohort models:", nrow(plasma_whole), "\n")
```

## C2. Model 1 — adjusted for ALT

```{r logistic-alt}
res_alt <- run_logistic(plasma_whole, covars,
                        adjust = c("Age", "sex_bin", "BMI", "ALT"))
sig_alt <- annotate_tbl(res_alt, chemical_details)

cat("Significant, ALT-adjusted:", nrow(sig_alt), "\n")
cat("Total number of sub-pathways, ALT-adjusted:",
    length(unique(sig_alt$Sub_pathway)), "\n")

datatable(sig_alt, caption = "MASH ~ metabolite + Age + sex + BMI + ALT",
          options = list(scrollX = TRUE))

# Sub-pathways represented by more than one / exactly one metabolite
multi_alt      <- names(which(table(sig_alt$Sub_pathway) > 1))
single_alt     <- names(which(table(sig_alt$Sub_pathway) == 1))

sig_alt_multi <- sig_alt %>%
  filter(Sub_pathway %in% multi_alt) %>%
  arrange(Sub_pathway) %>%
  select(Sub_pathway, metabolite)

sig_alt_single <- sig_alt %>%
  filter(Sub_pathway %in% single_alt) %>%
  arrange(Sub_pathway) %>%
  select(Sub_pathway, metabolite)
```

## C3. Model 2 — adjusted for ALT + AST

```{r logistic-alt-ast}
res_altast <- run_logistic(plasma_whole, covars,
                           adjust = c("Age", "sex_bin", "BMI", "ALT", "AST"))
sig_altast <- annotate_tbl(res_altast, chemical_details)

cat("Significant, ALT+AST-adjusted:", nrow(sig_altast), "\n")
cat("Total number of sub-pathways, ALT+AST-adjusted:",
    length(unique(sig_altast$Sub_pathway)), "\n")

datatable(sig_altast, caption = "MASH ~ metabolite + Age + sex + BMI + ALT + AST",
          options = list(scrollX = TRUE))
```

## C4. Metabolites robust to both models

```{r logistic-robust}
robust   <- intersect(sig_alt$metabolite, sig_altast$metabolite)
mash_ash <- sig_altast %>% filter(metabolite %in% robust)

cat("Significant in BOTH models:", length(robust), "\n")

datatable(mash_ash,
          caption = "MASH-associated independent of ALT and AST",
          options = list(scrollX = TRUE))
```

## C5. Singleton sub-pathways (ALT + AST model)

```{r logistic-singletons}
single_altast <- names(which(table(sig_altast$Sub_pathway) == 1))

sig_altast_m1 <- sig_altast %>%
  filter(Sub_pathway %in% single_altast) %>%
  arrange(Sub_pathway) %>%
  transmute(
    Super_pathway,
    Sub_pathway,
    metabolite,
    direction,
    OR       = round(OR, 3),
    `95% CI` = `95% CI`,
    p        = signif(p, 3),
    FDR      = signif(FDR, 3)
  )

cat("Sub-pathways with 1 significant metabolite each:", nrow(sig_altast_m1), "\n")

datatable(
  sig_altast_m1,
  rownames = FALSE,
  options  = list(scrollX = TRUE),
  caption  = paste("Singleton sub-pathways: MASH-associated metabolites represented",
                   "by a single metabolite (adjusted for age, sex, BMI, ALT, AST)")
)
```

## C6. Model diagnostics (optional)

Calibration (Hosmer-Lemeshow) and discrimination (AUC) per metabolite.
Disabled by default because it refits several hundred models; set
`eval=TRUE` on these chunks to run them.

```{r fit-diagnostics, eval=FALSE}
# Refit each significant metabolite in a given model and return HL p + AUC
fit_diagnostics <- function(sig_tbl, adjust, plasma_mat, covars) {
  form <- as.formula(paste("mash ~", paste(c("metabolite", adjust), collapse = " + ")))
  out <- lapply(sig_tbl$metabolite, function(met) {
    df <- data.frame(metabolite = as.numeric(plasma_mat[, met]), covars) |> na.omit()
    df$metabolite <- as.numeric(scale(df$metabolite))
    fit  <- suppressWarnings(glm(form, data = df, family = binomial()))
    ynum <- as.numeric(as.character(df$mash))
    prob <- fitted(fit)
    data.frame(
      metabolite = met,
      AUC  = round(as.numeric(auc(roc(ynum, prob, direction = "<", quiet = TRUE))), 3),
      HL_p = round(tryCatch(hoslem.test(ynum, prob, g = 10)$p.value,
                            error = function(e) NA), 3)
    )
  })
  bind_rows(out)
}

diag_alt    <- fit_diagnostics(sig_alt,    c("Age", "sex_bin", "BMI", "ALT"),
                               plasma_whole, covars)
diag_altast <- fit_diagnostics(sig_altast, c("Age", "sex_bin", "BMI", "ALT", "AST"),
                               plasma_whole, covars)

cat("== Model 1 (ALT) fit diagnostics ==\n")
cat("Mean AUC:", round(mean(diag_alt$AUC), 3),
    "| calibrated models (HL p > 0.05):", sum(diag_alt$HL_p > 0.05, na.rm = TRUE),
    "of", nrow(diag_alt), "\n")
datatable(diag_alt, caption = "Model 1 (ALT): AUC + Hosmer-Lemeshow per metabolite")

cat("\n== Model 2 (ALT+AST) fit diagnostics ==\n")
cat("Mean AUC:", round(mean(diag_altast$AUC), 3),
    "| calibrated models (HL p > 0.05):", sum(diag_altast$HL_p > 0.05, na.rm = TRUE),
    "of", nrow(diag_altast), "\n")
datatable(diag_altast, caption = "Model 2 (ALT+AST): AUC + Hosmer-Lemeshow per metabolite")
```

```{r model-comparison, eval=FALSE}
# Did adding AST change anything?
alt_only    <- setdiff(sig_alt$metabolite,    sig_altast$metabolite)  # dropped when AST added
altast_only <- setdiff(sig_altast$metabolite, sig_alt$metabolite)     # emerged when AST added

cat("Robust to both        :", length(robust), "\n")
cat("ALT only (lost w/ AST):", length(alt_only), "\n")
cat("ALT+AST only (gained) :", length(altast_only), "\n\n")

comparison <- full_join(
  res_alt    %>% select(metabolite, OR_alt    = OR, FDR_alt    = FDR),
  res_altast %>% select(metabolite, OR_altast = OR, FDR_altast = FDR),
  by = "metabolite"
) %>%
  filter(FDR_alt <= 0.05 | FDR_altast <= 0.05) %>%
  mutate(status = case_when(
    FDR_alt <= 0.05 & FDR_altast <= 0.05 ~ "both",
    FDR_alt <= 0.05                      ~ "ALT only",
    TRUE                                 ~ "ALT+AST only"
  )) %>%
  transmute(metabolite, status,
            OR_alt    = round(OR_alt, 3),    FDR_alt    = signif(FDR_alt, 3),
            OR_altast = round(OR_altast, 3), FDR_altast = signif(FDR_altast, 3)) %>%
  arrange(status, FDR_altast)

datatable(comparison,
          caption = "Model 1 (ALT) vs Model 2 (ALT+AST): per-metabolite comparison")
```

```{r poor-calibration, eval=FALSE}
# Metabolites whose model was NOT well calibrated (HL p < 0.05)
poorly_cal_alt    <- diag_alt$metabolite[!is.na(diag_alt$HL_p)       & diag_alt$HL_p       < 0.05]
poorly_cal_altast <- diag_altast$metabolite[!is.na(diag_altast$HL_p) & diag_altast$HL_p    < 0.05]

cat("== Poorly calibrated (HL p < 0.05) ==\n")
cat("Model 1 (ALT):", length(poorly_cal_alt), "metabolites\n")
print(poorly_cal_alt)
cat("\nModel 2 (ALT+AST):", length(poorly_cal_altast), "metabolites\n")
print(poorly_cal_altast)

datatable(
  diag_altast %>% filter(metabolite %in% poorly_cal_alt) %>% arrange(HL_p),
  caption = "Model 1 metabolites failing calibration (HL p < 0.05)"
)

datatable(
  diag_altast %>% filter(metabolite %in% poorly_cal_altast) %>% arrange(HL_p),
  caption = "Model 2 metabolites failing calibration (HL p < 0.05)"
)
```

```{r calibration-plot, eval=FALSE}
# Visual calibration check for one poorly calibrated metabolite
met <- poorly_cal_alt[15]
met

df  <- data.frame(metabolite = scale(as.numeric(plasma_whole[, met])), covars) |> na.omit()
fit <- glm(mash ~ metabolite + Age + sex_bin + BMI + ALT, df, family = binomial())
df$pred <- fitted(fit)
df$obs  <- as.numeric(as.character(df$mash))

df %>%
  mutate(bin = ntile(pred, 10)) %>%
  group_by(bin) %>%
  summarise(pred = mean(pred), obs = mean(obs)) %>%
  {
    plot(.$pred, .$obs, xlim = c(0, 1), ylim = c(0, 1))
    abline(0, 1, lty = 2)
  }
```

## C7. Pathway scores — coherent sub-pathways

Coherent = at least two validated metabolites moving in the same direction.

For each such sub-pathway:

1. take its validated metabolites;
2. Z-score each across patients;
3. flip metabolites with a negative beta, so higher = MASH-associated;
4. score = unweighted mean of the aligned Z-scores (per patient);
5. regress the score on MASH, adjusted for age / sex / BMI / ALT / AST.

Pathways are then ranked by the association of their score with MASH.
Helpers (`build_score`, `rank_pathways`) live in `functions_pathway_scores.R`.

```{r coherent-pathways}
coherent <- sig_altast %>%
  group_by(Sub_pathway) %>%
  summarise(
    n    = n(),
    n_up = sum(direction == "\u2191"),
    n_dn = sum(direction == "\u2193"),
    .groups = "drop"
  ) %>%
  filter(n_up >= 2 | n_dn >= 2)     # at least 2 in the same direction

cat("Sub-pathways to score:", nrow(coherent), "\n")
print(coherent$Sub_pathway)
```

```{r pathway-scores}
# Score using all metabolites in the sub-pathway (aligned by beta)
score_mat <- sapply(coherent$Sub_pathway, build_score)
rownames(score_mat) <- rownames(plasma_whole)

# Regress each score on MASH (adjusted), collect stats
pathway_ranked <- rank_pathways(score_mat, covars) %>%
  arrange(pvalue) %>%
  mutate(rank = row_number())

datatable(
  pathway_ranked %>% transmute(
    rank,
    Sub_pathway, n_metabolites,
    OR       = round(OR, 3),
    `95% CI` = sprintf("%.3f-%.3f", CI_low, CI_high),
    p        = signif(pvalue, 3),
    `p (Holm)` = signif(p_holm, 3)
  ),
  rownames = FALSE,
  options  = list(scrollX = TRUE),
  caption  = paste("Coherent sub-pathway scores ranked by association with MASH",
                   "(adjusted for age, sex, BMI, ALT, AST)")
)
```

```{r pathway-forest, fig.width=7.5, fig.height=5}
top10 <- pathway_ranked %>%
  arrange(p_holm) %>%
  slice_head(n = 10) %>%
  mutate(
    Sub_pathway = factor(Sub_pathway, levels = rev(Sub_pathway)),  # rank order top -> bottom
    label       = sprintf("%.2f (%.2f\u2013%.2f)", OR, CI_low, CI_high)
  )

# Headroom on the right for the OR labels (log scale)
x_max <- max(top10$CI_high) * 1.6

p_forest <- ggplot(top10, aes(x = OR, y = Sub_pathway)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high),
                 height = 0.18, colour = "#4a4a4a", linewidth = 0.6) +
  geom_point(aes(size = n_metabolites), colour = "#185FA5", shape = 18) +
  geom_text(aes(x = x_max, label = label), hjust = 1, size = 3.1, colour = "grey25") +
  scale_x_log10(
    breaks = c(0.5, 1, 2, 4, 8),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_size_continuous(range = c(2.5, 6), name = "Metabolites\nin score") +
  coord_cartesian(xlim = c(min(top10$CI_low) * 0.9, x_max), clip = "off") +
  labs(
    x        = "Odds ratio per SD of pathway score (95% CI, log scale)",
    y        = NULL,
    title    = "Sub-pathway scores associated with MASH",
    subtitle = "Top 10 by significance \u00b7 adjusted for age, sex, BMI, ALT, AST"
  ) +
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

p_forest
```

# Part D — Liver transcriptomics (whole cohort)

## D1. Align transcriptomic data to the MASH cohort

```{r tx-align}
# All patients with non-missing MASH status in clinical_data
# (no plasma / ALT restriction here - this is the full transcriptomic cohort)
common_samples_tx <- intersect(rownames(clinical_data), colnames(transcript_data))

clinical_tx <- clinical_data[common_samples_tx, , drop = FALSE]
tx_sub      <- transcript_data[, common_samples_tx, drop = FALSE]  # genes x samples

cat("Common (transcriptomic) patients:", length(common_samples_tx), "\n")
```

```{r tx-work-table}
work_tx <- data.frame(
  patient = common_samples_tx,
  mash    = clinical_tx$Mash,
  sex     = clinical_tx$sexe,
  ALT     = clinical_tx$TGPUL,
  AST     = clinical_tx$AST,
  Age     = clinical_tx$AgeJourIntervention,
  BMI     = clinical_tx$BMI,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(mash), !is.na(sex), !is.na(ALT), !is.na(AST),
         !is.na(Age), !is.na(BMI)) %>%
  mutate(
    sex_label = ifelse(as.character(sex) %in% female_codes, "Female", "Male"),
    mash_lab  = factor(ifelse(mash == 1, "MASH", "noMASH"),
                       levels = c("noMASH", "MASH"))
  )

cat("Patients entering limma (complete covariates):", nrow(work_tx), "\n")
print(table(work_tx$mash_lab, useNA = "ifany"))
```

## D2. Differential expression — limma

MASH vs no-MASH, adjusted for Age, sex, BMI, ALT and AST.

```{r tx-limma}
mt <- work_tx
rownames(mt) <- mt$patient
mt$grp     <- factor(mt$mash_lab, levels = c("noMASH", "MASH"))
mt$sex_bin <- ifelse(mt$sex_label == "Female", 0, 1)

cat("Group sizes entering limma:\n")
print(table(mt$grp, useNA = "ifany"))

xt <- as.matrix(tx_sub[, rownames(mt), drop = FALSE])

design_tx <- model.matrix(~ -1 + grp + Age + sex_bin + BMI + ALT + AST, data = mt)
fit_tx    <- lmFit(xt, design_tx)
comp_tx   <- comparisonsLimmaFct(fit_tx, "grpMASH - grpnoMASH", design_tx, nrow(xt))
res_tx    <- printResultsLimmaFct(comp_tx, topPrintHist = FALSE, topPrintVolc = FALSE)

sig_tx <- rownames(
  res_tx$results[
    !is.na(res_tx$results$adj.P.Val) &
      res_tx$results$adj.P.Val < 0.05, , drop = FALSE]
)
cat("MASH vs noMASH, whole cohort, adjusted for Age/sex/BMI/ALT/AST —",
    "significant (FDR < 0.05):", length(sig_tx), "\n")

res_tx$volcano
```

## D3. Gene set resources (built once)

The SYMBOL-to-ENTREZ map and the MSigDB GO:BP table are built once here and
reused in both Part D and Part E.

```{r gsea-resources}
# SYMBOL -> ENTREZID
symbEntrezid <- bitr(rownames(transcript_data),
                     fromType = "SYMBOL",
                     toType   = "ENTREZID",
                     OrgDb    = "org.Hs.eg.db")
symbEntrezid <- symbEntrezid[order(symbEntrezid$SYMBOL), ]
symbEntrezid <- symbEntrezid[!duplicated(symbEntrezid$SYMBOL), ]

nMatch <- nrow(symbEntrezid)
pMatch <- round(nMatch / nrow(transcript_data) * 100, 2)
cat("Genes matched to ENTREZ identifiers:", nMatch, "of",
    nrow(transcript_data), "(", pMatch, "%)\n")

# MSigDB GO:BP - works with both the old (category/subcategory) and the
# new (collection/subcollection) msigdbr argument names.
gobp_human <- tryCatch(
  msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  error = function(e)
    msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
)
gobp_human <- gobp_human[, c("gs_name", "entrez_gene")]
gobp_human$gs_name <- str_to_lower(substr(gobp_human$gs_name, 6, nchar(gobp_human$gs_name)))
```

## D4. Pathway enrichment analysis (whole cohort)

```{r tx-enrichment}
gobp_tx_mash <- enrichMsigdbFct(results      = comp_tx$mstat,
                                dataBase     = gobp_human,
                                symbEntrezid = symbEntrezid)

printResultsEnrchHistFct(enrchResults = gobp_tx_mash$resultHist,
                         xlab = "MASH vs no-MASH")

rownames(gobp_tx_mash$resultEnrich) <- NULL
sig_pathway <- gobp_tx_mash$resultEnrich[gobp_tx_mash$resultEnrich$p.adjust <= 0.05, ]

datatable(
  sig_pathway[, c("ID", "setSize", "NES", "rank", "p.adjust",
                  "n_leading", "leading_edge")],
  options = list(scrollX = TRUE),
  caption = "GO-BP enrichment, whole cohort (FDR <= 0.05)"
) %>%
  formatRound(columns = "NES", digits = 6) %>%
  formatSignif(columns = "p.adjust", digits = 6)
```

# Part E — Liver transcriptomics, normal ALT

## E1. Restrict to normal-ALT patients

```{r tx-normal-align}
common_samples_tx_nrm <- intersect(rownames(clinical_data), colnames(transcript_data))

clinical_tx_nrm <- clinical_data[common_samples_tx_nrm, , drop = FALSE]
tx_sub_nrm      <- transcript_data[, common_samples_tx_nrm, drop = FALSE]

cat("Common transcriptomic patients:", length(common_samples_tx_nrm), "\n")
```

```{r tx-normal-work-table}
work_tx_normal_alt <- data.frame(
  patient = common_samples_tx_nrm,
  mash    = clinical_tx_nrm$Mash,
  sex     = clinical_tx_nrm$sexe,
  ALT     = clinical_tx_nrm$TGPUL,
  AST     = clinical_tx_nrm$AST,
  Age     = clinical_tx_nrm$AgeJourIntervention,
  BMI     = clinical_tx_nrm$BMI,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(mash), !is.na(ALT)) %>%
  mutate(
    sex_label = case_when(
      as.character(sex) %in% female_codes ~ "Female",
      as.character(sex) %in% male_codes   ~ "Male",
      TRUE ~ NA_character_
    ),
    alt_normal = ALT >= alt_lo & ALT <= alt_hi,
    mash_lab   = factor(ifelse(mash == 1, "MASH", "noMASH"),
                        levels = c("noMASH", "MASH")),
    alt_level  = ifelse(alt_normal, "normal", "abnormal")
  ) %>%
  filter(!is.na(sex_label), alt_level == "normal")

cat("Patients entering the normal-ALT transcriptomic analysis:",
    nrow(work_tx_normal_alt), "\n")
print(table(work_tx_normal_alt$mash_lab, useNA = "ifany"))
```

## E2. Differential expression — limma

MASH vs no-MASH, adjusted for Age, sex and BMI. ALT and AST are **not** in the
model here: the cohort is already restricted to normal ALT.

```{r tx-normal-limma}
mt_normal_alt <- work_tx_normal_alt
rownames(mt_normal_alt) <- mt_normal_alt$patient

mt_normal_alt$grp     <- factor(mt_normal_alt$mash_lab, levels = c("noMASH", "MASH"))
mt_normal_alt$sex_bin <- ifelse(mt_normal_alt$sex_label == "Female", 0, 1)

cat("Group sizes entering limma — normal ALT:\n")
print(table(mt_normal_alt$grp, useNA = "ifany"))

# Expression matrix: genes x normal-ALT samples
xt_normal_alt <- as.matrix(tx_sub_nrm[, rownames(mt_normal_alt), drop = FALSE])

design_tx_normal_alt <- model.matrix(~ -1 + grp + Age + sex_bin + BMI,
                                     data = mt_normal_alt)
fit_tx_normal_alt    <- lmFit(xt_normal_alt, design_tx_normal_alt)
comp_tx_normal_alt   <- comparisonsLimmaFct(fit_tx_normal_alt, "grpMASH - grpnoMASH",
                                            design_tx_normal_alt, nrow(xt_normal_alt))
res_tx_normal_alt    <- printResultsLimmaFct(comp_tx_normal_alt,
                                             topPrintHist = FALSE, topPrintVolc = FALSE)

sig_tx_normal_alt <- rownames(
  res_tx_normal_alt$results[
    !is.na(res_tx_normal_alt$results$adj.P.Val) &
      res_tx_normal_alt$results$adj.P.Val < 0.05, , drop = FALSE]
)

cat("MASH vs noMASH, normal ALT, adjusted for Age/sex/BMI — ",
    "significant genes, FDR < 0.05: ", length(sig_tx_normal_alt), "\n", sep = "")

res_tx_normal_alt$volcano
```

## E3. Pathway enrichment analysis (normal ALT)

```{r tx-normal-enrichment}
gobp_tx_normal_alt <- enrichMsigdbFct(results      = comp_tx_normal_alt$mstat,
                                      dataBase     = gobp_human,
                                      symbEntrezid = symbEntrezid)

printResultsEnrchHistFct(enrchResults = gobp_tx_normal_alt$resultHist,
                         xlab = "MASH vs no-MASH — normal ALT")

rownames(gobp_tx_normal_alt$resultEnrich) <- NULL
sig_pathway_normal_alt <- gobp_tx_normal_alt$resultEnrich[
  gobp_tx_normal_alt$resultEnrich$p.adjust <= 0.05, ]

datatable(
  sig_pathway_normal_alt[, c("ID", "setSize", "NES", "rank", "p.adjust",
                             "n_leading", "leading_edge")],
  options = list(scrollX = TRUE),
  caption = "GO-BP enrichment, normal ALT (FDR <= 0.05)"
) %>%
  formatRound(columns = c("NES", "p.adjust"), digits = 6)

```

## E4. Comparison with the whole cohort

```{r tx-enrichment-comparison}
whole <- sig_pathway %>%
  transmute(ID,
            NES_whole = round(NES, 3),
            FDR_whole = signif(p.adjust, 3))

normal <- sig_pathway_normal_alt %>%
  transmute(ID,
            NES_normal = round(NES, 3),
            FDR_normal = signif(p.adjust, 3))

comp <- full_join(whole, normal, by = "ID") %>%
  mutate(
    in_whole  = !is.na(NES_whole),
    in_normal = !is.na(NES_normal),
    category = case_when(
       in_whole &  in_normal ~ "Shared biology",
      !in_whole &  in_normal ~ "Newly enriched in normal ALT",
       in_whole & !in_normal ~ "Only in whole cohort",
      TRUE ~ NA_character_
    ),
    interpretation = case_when(
      category == "Shared biology" & abs(NES_normal) > abs(NES_whole) ~
        "Shared biology (more prominent in normal ALT)",
      category == "Shared biology" ~ "Shared biology",
      TRUE ~ category
    )
  ) %>%
  arrange(desc(in_normal), FDR_normal, FDR_whole)

datatable(
  comp %>% select(
    Pathway            = ID,
    `NES (whole)`      = NES_whole,
    `FDR (whole)`      = FDR_whole,
    `NES (normal ALT)` = NES_normal,
    `FDR (normal ALT)` = FDR_normal,
    Interpretation     = interpretation
  ),
  rownames = FALSE,
  caption  = "GO-BP enrichment: whole cohort vs normal ALT — shared and distinct pathways",
  options  = list(pageLength = 15, scrollX = TRUE)
)
```

# Part F — Liver metabolomics, normal ALT

## F1. Restrict to normal-ALT patients with complete covariates

```{r liver-align}
common_liver <- intersect(work$patient, colnames(liver_data))

m_liver <- work %>%
  filter(patient %in% common_liver, alt_level == "normal",
         !is.na(Age), !is.na(BMI), !is.na(sex_label), !is.na(mash_lab)) %>%
  mutate(
    grp     = factor(mash_lab, levels = c("noMASH", "MASH")),
    sex_bin = ifelse(sex_label == "Female", 0, 1)
  )
rownames(m_liver) <- m_liver$patient

cat("Liver + clinical patients:", length(common_liver),
    "| normal ALT with complete covariates:", nrow(m_liver), "\n")
print(table(m_liver$grp))
```

## F2. Differential abundance — limma (adjusted for Age, sex, BMI)

```{r liver-limma}
x_liver <- as.matrix(liver_data[, m_liver$patient, drop = FALSE])

design_liver <- model.matrix(~ -1 + grp + Age + sex_bin + BMI, data = m_liver)
fit_liver    <- lmFit(x_liver, design_liver)
comp_liver   <- comparisonsLimmaFct(fit_liver, "grpMASH - grpnoMASH",
                                    design_liver, nrow(x_liver))
res_liver    <- printResultsLimmaFct(comp_liver,
                                     topPrintHist = FALSE, topPrintVolc = FALSE)

res_liver$volcano

datatable(
  res_liver$results[order(res_liver$results$logFC, decreasing = TRUE), ],
  caption = "Liver: MASH vs no-MASH within normal ALT",
  options = list(scrollX = TRUE)
)
```

## F3. Metabolon annotation

```{r liver-annotate}
sig_liver <- res_liver$results
idx <- match(rownames(sig_liver), l_chem_details$BIOCHEMICAL)
cat("Unannotated:", sum(is.na(idx)), "/", nrow(sig_liver), "\n")

sig_liver$Super_pathway <- l_chem_details$SUPER.PATHWAY[idx]
sig_liver$Sub_pathway   <- l_chem_details$SUB.PATHWAY[idx]
sig_liver$metabolite    <- rownames(sig_liver)
```

## F4. Summary by sub-pathway

```{r liver-subpathway-summary}
sub_sum_liver <- sig_liver %>%
  arrange(desc(abs(logFC))) %>%
  group_by(Super_pathway, Sub_pathway) %>%
  summarise(
    No_sig    = n(),
    No_up     = sum(logFC > 0),
    No_down   = sum(logFC < 0),
    Direction = case_when(
      all(logFC > 0) ~ "all up",
      all(logFC < 0) ~ "all down",
      TRUE ~ paste0(sum(logFC > 0), " up / ", sum(logFC < 0), " down")
    ),
    median_logFC = round(median(logFC), 3),
    median_FC    = round(2^median(logFC), 2),
    best_FDR     = signif(min(adj.P.Val), 3),
    strongest    = paste(head(metabolite, 3), collapse = "; "),
    .groups      = "drop"
  ) %>%
  arrange(desc(No_sig), best_FDR)

datatable(sub_sum_liver, options = list(scrollX = TRUE),
          caption = "Liver: significant metabolites by Metabolon sub-pathway (normal ALT)")
```

# Session info

```{r session-info}
#sessionInfo()
```
