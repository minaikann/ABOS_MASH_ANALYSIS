---
title: "Metabolites that identify MASH with normal ALT"
author: "Deborah Mina Ikann"
date: "2026-07-07"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo      = TRUE,
  message   = FALSE,
  warning   = FALSE,
  error     = FALSE,
  fig.align = "center"
)
```

# Questions

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

Structure of the document:

| Section | Content |
|---|---|
| Setup | libraries, constants, data, shared helpers |
| Cohort overview | counts + flowcharts for the whole cohort and the normal-ALT subset |
| Part A | plasma metabolomics, normal ALT |
| Part B | selected metabolites across CTRL / CM / LS |
| Part C | plasma metabolomics, whole cohort (validation) |
| Part D | liver transcriptomics, whole cohort |
| Part E | liver transcriptomics, normal ALT |
| Part F | liver metabolomics, normal ALT |

Within every part the order is the same: **align → describe (flowchart / table)
→ model → annotate → summarise**.

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

Every constant used anywhere in the document is defined **here and only here**

```{r paths}
data_dir <- "~/Documents/Clustering_ABOS/codes/Data"
fun_dir <- "~/Documents/Clustering_ABOS/codes/MASH/Functions_R"

# Sex coding used throughout
female_codes <- c("F", "Female", "FEMALE", "f", "female", "Femme")
male_codes   <- c("M", "Male",   "MALE",   "m", "male",   "Homme")

# Normal ALT window (flat cut-off, applied identically in Parts A, E and F)
alt_lo <- 0
alt_hi <- 35

# Covariates required for a patient to be considered "complete"
covar_cols <- c("sexe", "AgeJourIntervention", "BMI", "TGPUL", "AST")
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

## External helper functions

```{r source-functions}
# Pre-existing helper functions (limma comparison / plotting)
# NOTE: the file name contains a space - kept exactly as on disk.
source(file.path(fun_dir,"functions_ liver_metabo.R"))       # Diffrential Analysis function 
source(file.path(fun_dir,"functions_data_prep.R"))            # chemical-name mapping used during data loading
source(file.path(fun_dir,"functions_logistic_normalALT.R"))   # Part A: ALT vs metabolite logistic / CV helpers
source(file.path(fun_dir,"functions_plotting_clusters.R"))    # Part B: cluster boxplot + significance-bracket helpers
source(file.path(fun_dir,"functions_logistic_wholecohort.R")) # Part C: whole-cohort logistic engine + annotation
source(file.path(fun_dir,"functions_pathway_scores.R"))       # Part C: pathway scoring + ranking
source(file.path(fun_dir,"GSEA_function.R"))                  # Parts D/E: gene set enrichment analysis
source(file.path(fun_dir,"functions_flowchart.R"))            # 
```




## Name mapping and alignment

```{r map-names}
# Replace COMP_ID row names with human-readable chemical names where available
plasma_data <- map_chemical_names(plasma_data, chemical_details)
liver_data  <- map_metabolite_names(liver_data, l_chem_details)
```




# Cohort overview

Accounting for the full ABOS cohort and for the normal-ALT subset, across the
three omics layers. Counts first, flowcharts immediately after.

## Counts

```{r overview-counts}
n_total <- nrow(clinical_data)

# Variables a participant must have to enter the analysis
keep_cols <- c(covar_cols, "Mash")

complete_covars <- complete_rows(clinical_data, keep_cols)
complete_ids    <- rownames(complete_covars)
n_complete      <- length(complete_ids)
n_excl          <- n_total - n_complete

# Normal- vs high-ALT among complete-covariate patients
normal_ids <- complete_ids[complete_covars$TGPUL >= alt_lo &
                           complete_covars$TGPUL <= alt_hi]
n_normal   <- length(normal_ids)
n_high     <- n_complete - n_normal

# Omics availability, whole cohort and normal-ALT subset
omics <- list(plasma = colnames(plasma_data),
              tx     = colnames(transcript_data),
              liver  = colnames(liver_data))

n_omics      <- sapply(omics, function(ids) length(intersect(complete_ids, ids)))
n_omics_norm <- sapply(omics, function(ids) length(intersect(normal_ids,   ids)))

# Per-variable missingness, with display names for the exclusion panel
miss_cols <- c("ALT" = "TGPUL", "AST" = "AST", "MASH status" = "Mash")
n_miss <- sapply(miss_cols, function(v) sum(is.na(clinical_data[[v]])))
```

## Flowchart — whole cohort

```{r overview-flow-whole}
flow_chart(
  title = "ABOS cohort",
  start = sprintf("%s participants with clinical data", n_fmt(n_total)),
  steps = list(
    list(keep = sprintf("%s participants with complete Age, Sex, BMI, ALT and AST",
                        n_fmt(n_complete)),
         excl = c("Exclusion:",
                  sprintf("- Missing %s (n=%s)", names(n_miss), n_fmt(n_miss)),
                  sprintf("  (categories overlap; %s excluded in total)", n_fmt(n_excl))))
  ),
  leaves = c(
    sprintf("%s participants with plasma metabolomic data available",   n_fmt(n_omics["plasma"])),
    sprintf("%s participants with liver transcriptomic data available", n_fmt(n_omics["tx"])),
    sprintf("%s participants with liver metabolomic data available",    n_fmt(n_omics["liver"]))
  )
)
```



## Flowchart — normal-ALT subset

```{r overview-flow-normal}
flow_chart(
  title = "ABOS cohort, normal ALT",
  start = sprintf("%s participants with complete Age, Sex, BMI, ALT, AST and MASH status ",
                  n_fmt(n_complete)),
  steps = list(
    list(
      keep = sprintf("%s participants with normal ALT (%d-%d U/L)",
                     n_fmt(n_normal), alt_lo, alt_hi),
      excl = c("Exclusion:",
               sprintf("- ALT outside %d-%d U/L \n (n=%s)", alt_lo, alt_hi, n_fmt(n_high)))
    )
  ),
  leaves = c(
    sprintf("%s participants with plasma metabolomic data available",   n_fmt(n_omics_norm["plasma"])),
    sprintf("%s participants with liver transcriptomic data available", n_fmt(n_omics_norm["tx"])),
    sprintf("%s participants with liver metabolomic data available",    n_fmt(n_omics_norm["liver"]))
  )
)
```

# Part A — Plasma metabolomics, normal ALT



## Working dataset (reused by Parts A, B, C)

```{r align}
# Patients present in BOTH clinical and plasma tables
common_samples <- intersect(rownames(clinical_data), colnames(plasma_data))
clinical_sub   <- clinical_data[common_samples, , drop = FALSE]
plasma_sub     <- plasma_data[, common_samples, drop = FALSE]   # metabolites x samples

cat("Common (analysed) patients:", length(common_samples), "\n")
```

```{r work-table}
work <- data.frame(
  patient = common_samples,
  mash    = clinical_sub$Mash,
  sex     = clinical_sub$sexe,
  ALT     = clinical_sub$TGPUL,
  AST     = clinical_sub$AST,
  Age     = clinical_sub$AgeJourIntervention,
  BMI     = clinical_sub$BMI,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(mash), !is.na(ALT), !is.na(AST), !is.na(Age), !is.na(BMI)) %>%
  mutate(
    sex_label = case_when(
      as.character(sex) %in% female_codes ~ "Female",
      as.character(sex) %in% male_codes   ~ "Male",
      TRUE ~ NA_character_
    ),
    mash_lab  = factor(ifelse(mash == 1, "MASH", "noMASH"),
                       levels = c("noMASH", "MASH")),
    alt_level = ifelse(ALT >= alt_lo & ALT <= alt_hi, "normal", "high")
  ) %>%
  filter(!is.na(sex_label))

stopifnot(!anyNA(work))    # no missing values anywhere in the table

# Normal-ALT subset, derived once and reused everywhere below
work_normal_alt <- work %>% filter(alt_level == "normal")

cat("Plasma \u2229 clinical:", length(common_samples),
    "| retained in `work`:", nrow(work),
    "| normal ALT (", alt_lo, "-", alt_hi, " U/L):", nrow(work_normal_alt), "\n")
print(table(work_normal_alt$mash_lab, useNA = "ifany"))
```

## A1. Cohort accounting

```{r A-counts}
plasma_ids      <- colnames(plasma_data)
plasma_clin_ids <- intersect(plasma_ids, rownames(clinical_data))
nA_plasma_total <- length(plasma_clin_ids)

clinical_plasma <- clinical_data[plasma_clin_ids, , drop = FALSE]

# Complete covariates AND known MASH status, among patients with plasma data
A_complete <- complete_rows(clinical_plasma, c(covar_cols, "Mash")) %>%
  mutate(
    sex_label = ifelse(as.character(sexe) %in% female_codes, "Female", "Male"),
    mash_lab  = ifelse(Mash == 1, "MASH", "noMASH")
  )

nA_complete    <- nrow(A_complete)
nA_excl_covars <- nA_plasma_total - nA_complete

nA_cols <- c("MASH status" = "Mash", "AST" = "AST")
nA_miss <- sapply(nA_cols, function(v) sum(is.na(clinical_plasma[[v]])))

# Normal-ALT subset (from `work`, built once in Setup)
nA_high   <- sum(work$alt_level == "high")
nA_normal <- nrow(work_normal_alt)

cat("Plasma patients:", nA_plasma_total,
    "| complete covariates:", nA_complete, "(excluded", nA_excl_covars, ")\n",
    "MASH:", sum(A_complete$mash_lab == "MASH"),
    "| no-MASH:", sum(A_complete$mash_lab == "noMASH"), "\n",
    "Normal ALT:", nA_normal, "| high ALT (excluded):", nA_high, "\n")
```

### Flowchart — all plasma patients

```{r A-flow-all}
flow_chart(
  title = "Plasma metabolomics: Whole cohort",
  start = sprintf("Participants with plasma metabolomic data  \n (n = %s)", n_fmt(nA_plasma_total)),
  steps = list(
    list(keep = sprintf("Participants with complete  plasma metabolomi data \n (n = %s)",
                        n_fmt(nA_complete)),
         excl = c("Exclusion:",
                  sprintf("- Missing %s (n=%s)", names(nA_miss), n_fmt(nA_miss)),
                  sprintf("  (categories overlap; %s excluded in total)",
                          n_fmt(nA_excl_covars))))
  ),
  leaves = c(
    sprintf("Participants with MASH \n (n = %s)",
            n_fmt(sum(A_complete$mash_lab == "MASH"))),
    sprintf("Participants without MASH \n (n = %s)",
            n_fmt(sum(A_complete$mash_lab == "noMASH")))
  )
)
```

### Flowchart — normal-ALT plasma cohort

```{r A-flow-normal}
flow_chart(
  title = "Plasma metabolomics: Normal ALT",
  start = sprintf("Participants with plasma metabolomic data \n (n = %s)", n_fmt(nA_plasma_total)),
  steps = list(
    list(keep = sprintf("Participants with complete data \n (n = %s)",
                        n_fmt(nA_complete)),
         excl = c("Exclusion:",
                  sprintf("- Missing %s (n=%s)", names(nA_miss), n_fmt(nA_miss)),
                  sprintf("  (categories overlap; %s excluded in total)",
                          n_fmt(nA_excl_covars)))),
    list(keep = sprintf("Participants with normal ALT (%d-%d U/L) \n (n = %s)",
                        alt_lo, alt_hi, n_fmt(nA_normal)),
         excl = c("Exclusion:",
                  sprintf("- ALT outside %d-%d U/L (n=%s)",
                          alt_lo, alt_hi, n_fmt(nA_high))))
  ),
  leaves = c(
    sprintf("Participants with MASH \n (n = %s)",
            n_fmt(sum(work_normal_alt$mash_lab == "MASH"))),
    sprintf("Participants without MASH \n (n = %s)",
            n_fmt(sum(work_normal_alt$mash_lab == "noMASH")))
  )
)
```

## A2. Summary table (normal ALT)

Group comparison test: Wilcoxon for continuous variables, Fisher for
categorical ones.

```{r A-cohort-table}
clinical_table <- work_normal_alt %>%
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

The model is fitted on `work_normal_alt`, i.e. exactly the cohort counted in
A1. (The previous version filtered on `!alt_high`, which let a slightly
different set of patients through and produced group sizes that did not match
the flowchart.)

```{r A-limma}
m <- work_normal_alt
rownames(m) <- m$patient
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

## A4. Metabolon annotation

```{r A-annotate}
sig <- res_alt_norm$results
idx <- match(rownames(sig), chemical_details$CHEMICAL_NAME)
cat("Unannotated:", sum(is.na(idx)), "/", nrow(sig), "\n")

sig$Super_pathway <- chemical_details$SUPER_PATHWAY[idx]
sig$Sub_pathway   <- chemical_details$SUB_PATHWAY[idx]
sig$metabolite    <- rownames(sig)

datatable(sig[, c("Super_pathway", "Sub_pathway", "logFC", "FC", "adj.P.Val")],
          caption = "MASH vs no-MASH, normal ALT — Metabolon annotation",
          options = list(scrollX = TRUE))
```

## A5. Summary by sub-pathway

```{r A-subpathway-summary}
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

## A6. Logistic regression: ALT alone vs metabolite alone vs ALT + metabolite

### Data prep (RINT)

```{r A-logistic-data-prep}
## RINT is applied to every metabolite BEFORE any outcome is examined,
## so it cannot leak outcome information.

plasma_aligned <- as.data.frame(t(plasma_sub))          # samples x metabolites
plasma_aligned <- plasma_aligned[work$patient, , drop = FALSE]

# Rank-based inverse normal transform
rint <- function(x) {
  n <- sum(!is.na(x))
  if (n < 3) return(x)
  qnorm((rank(x, na.last = "keep") - 0.5) / n)
}
plasma_aligned <- as.data.frame(lapply(plasma_aligned, rint))
rownames(plasma_aligned) <- work$patient                # lapply() drops row names

# Metabolite names: clean for modelling, keep a map back to the originals
orig_metabolites  <- rownames(plasma_sub)
clean_metabolites <- make.names(orig_metabolites)
colnames(plasma_aligned) <- clean_metabolites
name_map <- setNames(orig_metabolites, clean_metabolites)

# Modelling frames: built AFTER the transform
df_all <- work %>%
  select(patient, mash_lab, ALT, alt_level, sex) %>%
  cbind(plasma_aligned[work$patient, , drop = FALSE]) %>%
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

```{r A-logistic-run, eval=FALSE}
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

```{r B-metabolite-list}
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

```{r B-long}
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

```{r B-tests}
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

```{r B-boxplots, fig.width=5.5, fig.height=5.5}
mets  <- sort(unique(long$Metabolite))
plots <- lapply(mets, make_plot)

# pdf("metabolites_per_page.pdf", width = 5.5, height = 5.5)
for (p in plots) print(p)
# dev.off()
```

# Part C — Validation in the whole biopsy cohort

Helper functions for this section (`run_logistic`, `annotate_tbl`) live in
`functions_logistic_wholecohort.R`; the pathway-score helpers (`build_score`,
`rank_pathways`) live in `functions_pathway_scores.R`.

## C1. Data prep

```{r C-prep}
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

```{r C-logistic-alt}
res_alt <- run_logistic(plasma_whole, covars,
                        adjust = c("Age", "sex_bin", "BMI", "ALT"))
sig_alt <- annotate_tbl(res_alt, chemical_details)

cat("Significant, ALT-adjusted:", nrow(sig_alt), "\n")
cat("Total number of sub-pathways, ALT-adjusted:",
    length(unique(sig_alt$Sub_pathway)), "\n")

datatable(sig_alt, caption = "MASH ~ metabolite + Age + sex + BMI + ALT",
          options = list(scrollX = TRUE))

# Sub-pathways represented by more than one / exactly one metabolite
multi_alt  <- names(which(table(sig_alt$Sub_pathway) > 1))
single_alt <- names(which(table(sig_alt$Sub_pathway) == 1))

sig_alt_multi <- sig_alt %>%
  filter(Sub_pathway %in% multi_alt) %>%
  arrange(Sub_pathway) %>%
  select(Sub_pathway, metabolite)

sig_alt_single <- sig_alt %>%
  filter(Sub_pathway %in% single_alt) %>%
  arrange(Sub_pathway) %>%
  select(Sub_pathway, metabolite)

volcanoLogitFct(res_alt)
forestLogitFct(res_alt, n_show = 10,
               title    = "Metabolites associated with MASH",
               subtitle = "Top 10 by FDR \u00b7 adjusted for Age, Sex, BMI, ALT")
```

## C3. Model 2 — adjusted for ALT + AST

```{r C-logistic-alt-ast}
res_altast <- run_logistic(plasma_whole, covars,
                           adjust = c("Age", "sex_bin", "BMI", "ALT", "AST"))
sig_altast <- annotate_tbl(res_altast, chemical_details)

cat("Significant, ALT+AST-adjusted:", nrow(sig_altast), "\n")
cat("Total number of sub-pathways, ALT+AST-adjusted:",
    length(unique(sig_altast$Sub_pathway)), "\n")

datatable(sig_altast, caption = "MASH ~ metabolite + Age + sex + BMI + ALT + AST",
          options = list(scrollX = TRUE))

volcanoLogitFct(res_altast)

forestLogitFct(res_altast, n_show = 10,
               title    = "Metabolites associated with MASH",
               subtitle = "Top 10 by FDR \u00b7 adjusted for Age, Sex, BMI, ALT, AST")
```

## C4. Metabolites robust to both models

```{r C-logistic-robust}
robust   <- intersect(sig_alt$metabolite, sig_altast$metabolite)
mash_ash <- sig_altast %>% filter(metabolite %in% robust)

cat("Significant in BOTH models:", length(robust), "\n")

datatable(mash_ash,
          caption = "MASH-associated independent of ALT and AST",
          options = list(scrollX = TRUE))
```
## C5. Normal-ALT hits: replication in the whole cohort

```{r C-normalalt-replication}
# Part A hits (limma, normal ALT). `sig` was annotated in A4.
a_hits <- sig %>%
  filter(!is.na(adj.P.Val), adj.P.Val < 0.05) %>%
  transmute(
    metabolite,
    Super_pathway, Sub_pathway,
    logFC_normalALT = round(logFC, 3),
    FDR_normalALT   = signif(adj.P.Val, 3),
    dir_normalALT   = ifelse(logFC > 0, "\u2191", "\u2193")
  )

# Direction + FDR from each whole-cohort model, for ALL metabolites
# (using the unfiltered run_logistic() output, so "absent" can be told apart
#  from "not tested")
m1 <- res_alt %>%
  transmute(metabolite,
            OR_alt  = round(OR, 3),
            FDR_alt = signif(FDR, 3),
            dir_alt = ifelse(OR > 1, "\u2191", "\u2193"))

m2 <- res_altast %>%
  transmute(metabolite,
            OR_altast  = round(OR, 3),
            FDR_altast = signif(FDR, 3),
            dir_altast = ifelse(OR > 1, "\u2191", "\u2193"))

replication <- a_hits %>%
  left_join(m1, by = "metabolite") %>%
  left_join(m2, by = "metabolite") %>%
  mutate(
    in_alt = case_when(
      is.na(FDR_alt)                            ~ "not tested",
      metabolite %in% sig_alt$metabolite        ~ "present",
      TRUE                                      ~ "absent"
    ),
    in_altast = case_when(
      is.na(FDR_altast)                         ~ "not tested",
      metabolite %in% sig_altast$metabolite     ~ "present",
      TRUE                                      ~ "absent"
    ),
    # direction agreement is only meaningful where the metabolite replicated
    concord_alt = case_when(
      in_alt != "present"            ~ NA_character_,
      dir_alt == dir_normalALT       ~ "same",
      TRUE                           ~ "opposite"
    ),
    concord_altast = case_when(
      in_altast != "present"         ~ NA_character_,
      dir_altast == dir_normalALT    ~ "same",
      TRUE                           ~ "opposite"
    ),
    summary = case_when(
      in_alt == "present" & in_altast == "present" &
        concord_alt == "same" & concord_altast == "same" ~ "Replicated in both, same direction",
      in_alt == "present" & in_altast == "present"       ~ "Replicated in both, direction differs",
      in_altast == "present"                             ~ "Replicated with ALT+AST only",
      in_alt == "present"                                ~ "Replicated with ALT only (lost when AST added)",
      in_alt == "not tested" | in_altast == "not tested" ~ "Not testable in whole cohort",
      TRUE                                               ~ "Not replicated"
    )
  ) %>%
  arrange(summary, FDR_normalALT)

cat("Normal-ALT hits:", nrow(replication), "\n")
print(table(replication$summary))
cat("\nDirection agreement where replicated (ALT+AST model):\n")
print(table(replication$concord_altast, useNA = "no"))

datatable(
  replication %>% transmute(
    Metabolite            = metabolite,
    `Super pathway`       = Super_pathway,
    `Sub pathway`         = Sub_pathway,
    `Normal ALT: logFC`   = logFC_normalALT,
    `Normal ALT: FDR`     = FDR_normalALT,
    `Normal ALT`          = dir_normalALT,
    `Model 1 (ALT)`       = in_alt,
    `M1 OR`               = OR_alt,
    `M1 FDR`              = FDR_alt,
    `M1 direction`        = dir_alt,
    `M1 vs normal ALT`    = concord_alt,
    `Model 2 (ALT+AST)`   = in_altast,
    `M2 OR`               = OR_altast,
    `M2 FDR`              = FDR_altast,
    `M2 direction`        = dir_altast,
    `M2 vs normal ALT`    = concord_altast,
    Outcome               = summary
  ),
  rownames = FALSE,
  options  = list(pageLength = 20, scrollX = TRUE),
  caption  = paste("Normal-ALT metabolites (limma, FDR < 0.05) and their status in the",
                   "whole-cohort logistic models adjusted for ALT, and for ALT + AST")
)
```

## C5. Singleton sub-pathways (ALT + AST model)

```{r C-logistic-singletons}
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

## C6. Pathway scores — coherent sub-pathways

Coherent = at least two validated metabolites moving in the same direction.

For each such sub-pathway:

1. take its validated metabolites;
2. Z-score each across patients;
3. flip metabolites with a negative beta, so higher = MASH-associated;
4. score = unweighted mean of the aligned Z-scores (per patient);
5. regress the score on MASH, adjusted for age / sex / BMI / ALT / AST.

```{r C-coherent-pathways}
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

```{r C-pathway-scores}
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
    OR         = round(OR, 3),
    `95% CI`   = sprintf("%.3f-%.3f", CI_low, CI_high),
    p          = signif(pvalue, 3),
    `p (Holm)` = signif(p_holm, 3)
  ),
  rownames = FALSE,
  options  = list(scrollX = TRUE),

  caption  = paste("Coherent sub-pathway scores ranked by association with MASH",
                   "(adjusted for age, sex, BMI, ALT, AST)")
)


forestPathwayFct(pathway_ranked)
```




# Part D — Liver transcriptomics, whole cohort

## D1. Align transcriptomic data

`common_samples_tx`, `clinical_tx` and `tx_sub` are built once here and reused
by Part E (the original code rebuilt an identical copy under a different name).

```{r D-align}
common_samples_tx <- intersect(rownames(clinical_data), colnames(transcript_data))

clinical_tx <- clinical_data[common_samples_tx, , drop = FALSE]
tx_sub      <- transcript_data[, common_samples_tx, drop = FALSE]  # genes x samples

cat("Common (transcriptomic) patients:", length(common_samples_tx), "\n")
```

## D2. Working table (complete covariates)

```{r D-work-table}
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

## D3. Flowchart

```{r D-flow}
nD_total    <- length(common_samples_tx)
nD_retained <- nrow(work_tx)
nD_excl     <- nD_total - nD_retained

nD_cols <- c("Mash", "sexe", "AgeJourIntervention", "BMI", "TGPUL", "AST")
nD_miss <- sapply(nD_cols, function(v) sum(is.na(clinical_tx[[v]])))
```

```{r}
flow_chart(
  title = "Liver transcriptomics: whole cohort",
  start = sprintf("Participants with liver transcriptomics data \n (n = %s)", n_fmt(nD_total)),
  steps = list(
    list(keep = sprintf("Participants with data \n (n = %s)", n_fmt(nD_retained)),
         excl = c("Exclusion:",
                  sprintf("- Missing MASH status (n=%s)", n_fmt(nD_miss["Mash"])),
                  sprintf("- Missing ALT (n=%s)",         n_fmt(nD_miss["TGPUL"])),
                  sprintf("- Missing AST (n=%s)",         n_fmt(nD_miss["AST"])),
                  sprintf("  (categories overlap; %s excluded in total)", n_fmt(nD_excl))))
  ),
  leaves = c(sprintf("Participants with MASH \n (n = %s)",    n_fmt(sum(work_tx$mash_lab == "MASH"))),
             sprintf("Participants with no MASH \n (n = %s)", n_fmt(sum(work_tx$mash_lab == "noMASH"))))
)
```


## D4. Differential expression — limma

MASH vs no-MASH, adjusted for Age, sex, BMI, ALT and AST.

```{r D-limma}
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

## D5. Gene set resources (built once, reused in Part E)

```{r D-gsea-resources}
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

## D6. Pathway enrichment (whole cohort)

```{r D-enrichment}
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

Reuses `clinical_tx` / `tx_sub` from D1 and the gene-set resources from D5.

## E1. Working table (normal ALT)

```{r E-work-table}
work_tx_normal_alt <- data.frame(
  patient = common_samples_tx,
  mash    = clinical_tx$Mash,
  sex     = clinical_tx$sexe,
  ALT     = clinical_tx$TGPUL,
  AST     = clinical_tx$AST,
  Age     = clinical_tx$AgeJourIntervention,
  BMI     = clinical_tx$BMI,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(mash), !is.na(ALT), !is.na(AST), !is.na(Age), !is.na(BMI)) %>%
  mutate(
    sex_label = case_when(
      as.character(sex) %in% female_codes ~ "Female",
      as.character(sex) %in% male_codes   ~ "Male",
      TRUE ~ NA_character_
    ),
    mash_lab  = factor(ifelse(mash == 1, "MASH", "noMASH"),
                       levels = c("noMASH", "MASH")),
    alt_level = ifelse(ALT >= alt_lo & ALT <= alt_hi, "normal", "high")
  ) %>%
  filter(!is.na(sex_label))

stopifnot(!anyNA(work_tx_normal_alt))

# Step counts for the flowchart
nE_total <- length(common_samples_tx)
nE_step1 <- nrow(work_tx_normal_alt)          # complete covariates

work_tx_normal_alt <- work_tx_normal_alt %>% filter(alt_level == "normal")
nE_step2 <- nrow(work_tx_normal_alt)

# Step 1: missingness among the transcriptomic patients
nE_cols <- c("MASH status" = "Mash", "sex" = "sexe", "age" = "AgeJourIntervention",
             "BMI" = "BMI", "ALT" = "TGPUL", "AST" = "AST")
nE_miss <- sapply(nE_cols, function(v) sum(is.na(clinical_tx[[v]])))
nE_miss <- nE_miss[nE_miss > 0]

# Step 2: high ALT among the complete-covariate patients
nE_high <- nE_step1 - nE_step2

cat("Transcriptomic patients:", nE_total,
    "| complete covariates:", nE_step1,
    "| normal ALT:", nE_step2, "\n")
print(table(work_tx_normal_alt$mash_lab, useNA = "ifany"))
```

## E2. Flowchart

```{r E-flow}
flow_chart(
  title = "Liver transcriptomics: normal ALT",
  start = sprintf("Participants with liver transcriptomic data \n (n = %s)", n_fmt(nE_total)),
  steps = list(
    list(keep = sprintf("Participants with MASH & ALT data \n (n = %s)", n_fmt(nE_step1)),
         excl = c("Exclusion:",
                  sprintf("- Missing %s (n=%s)", names(nE_miss), n_fmt(nE_miss)),
                  sprintf("  (categories overlap; %s excluded in total)",
                          n_fmt(nE_total - nE_step1)))),
    list(keep = sprintf("Participants with Normal-ALT \n (n = %s)", n_fmt(nE_step2)),
         excl = c("Exclusion:",
                  sprintf("- ALT outside %d-%d U/L (n=%s)", alt_lo, alt_hi, n_fmt(nE_high))))
         ),
         
  leaves = c(sprintf("Participants with MASH \n (n = %s)",    n_fmt(sum(work_tx_normal_alt$mash_lab == "MASH"))),
             sprintf("Participants with no MASH \n (n = %s)", n_fmt(sum(work_tx_normal_alt$mash_lab == "noMASH"))))
)
```

## E3. Differential expression — limma

MASH vs no-MASH, adjusted for Age, sex and BMI. ALT and AST are **not** in the
model here: the cohort is already restricted to normal ALT.

```{r E-limma}
mt_normal_alt <- work_tx_normal_alt
rownames(mt_normal_alt) <- mt_normal_alt$patient

mt_normal_alt$grp     <- factor(mt_normal_alt$mash_lab, levels = c("noMASH", "MASH"))
mt_normal_alt$sex_bin <- ifelse(mt_normal_alt$sex_label == "Female", 0, 1)

cat("Group sizes entering limma — normal ALT:\n")
print(table(mt_normal_alt$grp, useNA = "ifany"))

# Expression matrix: genes x normal-ALT samples
xt_normal_alt <- as.matrix(tx_sub[, rownames(mt_normal_alt), drop = FALSE])

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

## E4. Pathway enrichment (normal ALT)

```{r E-enrichment}
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

## E5. Comparison with the whole cohort

```{r E-enrichment-comparison}
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

```{r}
## F1. Liver metabolomics working table
liver_ids      <- intersect(colnames(liver_data), rownames(clinical_data))
clinical_liver <- clinical_data[liver_ids, , drop = FALSE]

work_liver <- data.frame(
  patient = liver_ids,
  mash    = clinical_liver$Mash,
  sex     = clinical_liver$sexe,
  ALT     = clinical_liver$TGPUL,
  AST     = clinical_liver$AST,
  Age     = clinical_liver$AgeJourIntervention,
  BMI     = clinical_liver$BMI,
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(mash), !is.na(ALT), !is.na(AST), !is.na(Age), !is.na(BMI)) %>%
  mutate(
    sex_label = case_when(
      as.character(sex) %in% female_codes ~ "Female",
      as.character(sex) %in% male_codes   ~ "Male",
      TRUE ~ NA_character_
    ),
    mash_lab  = factor(ifelse(mash == 1, "MASH", "noMASH"),
                       levels = c("noMASH", "MASH")),
    grp       = mash_lab,
    sex_bin   = ifelse(sex_label == "Female", 0, 1),
    alt_level = ifelse(ALT >= alt_lo & ALT <= alt_hi, "normal", "high")
  ) %>%
  filter(!is.na(sex_label))

stopifnot(!anyNA(work_liver))

# Normal-ALT subset: the cohort Part F actually models
m_liver <- work_liver %>% filter(alt_level == "normal")
rownames(m_liver) <- m_liver$patient

cat("Liver + clinical patients:", length(liver_ids),
    "| complete covariates:", nrow(work_liver),
    "| normal ALT:", nrow(m_liver), "\n")
print(table(m_liver$grp))
```


```{r F-align}
nF_total    <- length(liver_ids)
nF_retained <- nrow(m_liver)

nF_cols <- c("MASH status" = "Mash", "sex" = "sexe", "age" = "AgeJourIntervention",
             "BMI" = "BMI", "ALT" = "TGPUL")
nF_miss <- sapply(nF_cols, function(v) sum(is.na(clinical_liver[[v]])))
nF_miss <- nF_miss[nF_miss > 0]

nF_high <- sum(work_liver$alt_level == "high")

stopifnot(nF_high + sum(nF_miss) >= nF_total - nF_retained)
```

## F2. Flowchart

```{r F-flow}

flow_chart(
  title = "Liver metabolomics: normal ALT",
  start = sprintf("Participants with liver metabolomic data\n (n = %s)", n_fmt(nF_total)),
  steps = list(
    list(keep = sprintf("Participants with normal ALT and complete covariates \n (n = %s)",
                        n_fmt(nF_retained)),
         excl = c("Exclusion:",
                  sprintf("- ALT outside %d-%d U/L \n (n=%s)", alt_lo, alt_hi, n_fmt(nF_high)),
                  sprintf("- Missing %s \n (n=%s)", names(nF_miss), n_fmt(nF_miss))))
  ),
  leaves = c(
    sprintf("Participants with MASH \n (n = %s)",    n_fmt(sum(m_liver$grp == "MASH"))),
    sprintf("Participants without MASH \n (n = %s)", n_fmt(sum(m_liver$grp == "noMASH")))
  )
)
```






## F3. Differential abundance — limma (adjusted for Age, sex, BMI)

```{r F-limma}
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

## F4. Metabolon annotation

```{r F-annotate}
sig_liver <- res_liver$results
idx <- match(rownames(sig_liver), l_chem_details$BIOCHEMICAL)
cat("Unannotated:", sum(is.na(idx)), "/", nrow(sig_liver), "\n")

sig_liver$Super_pathway <- l_chem_details$SUPER.PATHWAY[idx]
sig_liver$Sub_pathway   <- l_chem_details$SUB.PATHWAY[idx]
sig_liver$metabolite    <- rownames(sig_liver)
```

## F5. Summary by sub-pathway

```{r F-subpathway-summary}
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
    strongest    = paste(head(metabolite, 3), collapse = "; "),
    .groups      = "drop"
  ) %>%
  arrange(desc(No_sig))

datatable(sub_sum_liver, options = list(scrollX = TRUE),
          caption = "Liver: significant metabolites by Metabolon sub-pathway (normal ALT)")
```

# Session info

```{r session-info}
# sessionInfo()
```
