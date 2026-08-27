## =============================================================================
## functions_normalALT_addendum.R
##
## Helper functions for the normal-ALT follow-up analyses:
##   - GSEA on the full ranked limma gene list (transcriptomics)
##   - mapping GSEA gene sets onto named biological programs
##   - Metabolon super/sub-pathway summaries of limma results
##   - side-by-side comparison of whole cohort vs normal ALT
##
## Sourced from MASH_metabolites_analysis.Rmd
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(msigdbr)
})


## -----------------------------------------------------------------------------
## 1. MSigDB retrieval (works with msigdbr < 10 and >= 10)
## -----------------------------------------------------------------------------
## msigdbr >= 10 renamed category/subcategory to collection/subcollection.
## This wrapper tries the new API first and falls back to the old one.

get_msigdb <- function(collection,
                       subcollection = NULL,
                       species = "Homo sapiens",
                       id_col  = "gene_symbol") {

  args_new <- list(species = species, collection = collection)
  if (!is.null(subcollection)) args_new$subcollection <- subcollection

  args_old <- list(species = species, category = collection)
  if (!is.null(subcollection)) args_old$subcategory <- subcollection

  df <- tryCatch(
    do.call(msigdbr::msigdbr, args_new),
    error = function(e) do.call(msigdbr::msigdbr, args_old)
  )

  df <- as.data.frame(df)
  if (!id_col %in% colnames(df)) {
    stop("Column '", id_col, "' not returned by msigdbr(). Available: ",
         paste(colnames(df), collapse = ", "))
  }

  out <- df[, c("gs_name", id_col)]
  colnames(out) <- c("gs_name", "gene")
  dplyr::distinct(out)
}


## -----------------------------------------------------------------------------
## 2. Ranked gene list from a limma result
## -----------------------------------------------------------------------------
## GSEA needs EVERY gene, ranked - not just the FDR-significant ones.
## Default statistic is the moderated t (as requested).
## Duplicated gene symbols are collapsed to the most extreme statistic.

rank_from_limma <- function(res, stat_col = "t") {

  tab <- if (is.data.frame(res)) res else res$results
  if (is.null(tab)) stop("Could not find a results table in the object supplied.")

  if (!stat_col %in% colnames(tab)) {
    alt <- intersect(c("t", "mstat", "moderated.t", "B", "logFC"), colnames(tab))
    if (!length(alt)) {
      stop("No usable ranking statistic. Columns present: ",
           paste(colnames(tab), collapse = ", "))
    }
    warning("'", stat_col, "' not found; ranking on '", alt[1], "' instead.")
    stat_col <- alt[1]
  }

  r <- stats::setNames(as.numeric(tab[[stat_col]]), rownames(tab))
  r <- r[!is.na(r) & !is.na(names(r)) & names(r) != ""]

  # collapse duplicated symbols -> keep most extreme statistic
  if (any(duplicated(names(r)))) {
    r <- tapply(r, names(r), function(z) z[which.max(abs(z))])
    r <- stats::setNames(as.numeric(r), names(r))
  }

  sort(r, decreasing = TRUE)
}


## -----------------------------------------------------------------------------
## 3. GSEA wrapper
## -----------------------------------------------------------------------------
## pvalueCutoff = 1 so that NON-significant programs are still returned - we
## need them for the comparison table (a program that is absent in one cohort
## is as informative as one that is present).

run_gsea <- function(ranks, t2g,
                     minGSSize = 10,
                     maxGSSize = 500,
                     seed = 42) {

  stopifnot(!is.unsorted(rev(ranks)))   # must be decreasing

  set.seed(seed)
  res <- clusterProfiler::GSEA(
    geneList      = ranks,
    TERM2GENE     = t2g,
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    eps           = 0,
    seed          = TRUE,
    verbose       = FALSE
  )

  out <- as.data.frame(res)
  rownames(out) <- NULL
  out
}


## -----------------------------------------------------------------------------
## 4. Biological programs requested in the email
## -----------------------------------------------------------------------------
## Regexes are matched against lower-cased MSigDB set names
## (e.g. "gobp_leukocyte_migration", "reactome_bile_acid_and_bile_salt_metabolism").
## Edit / extend freely - this is the only place the mapping lives.

program_patterns <- list(
  "Immune-cell migration" =
    "leukocyte_migration|leukocyte_chemotaxis|myeloid_leukocyte_migration|granulocyte_chemotaxis|cell_chemotaxis|leukocyte_cell_cell_adhesion",

  "Lymphocyte / T-cell migration" =
    "lymphocyte_migration|lymphocyte_chemotaxis|t_cell_migration|t_cell_chemotaxis",

  "IL-2 signalling" =
    "interleukin_2",

  "IL-6 signalling" =
    "interleukin_6",

  "TNF signalling" =
    "tumor_necrosis_factor|tnfa_signaling",

  "Sterol / cholesterol biosynthesis" =
    "sterol_biosynthetic_process|cholesterol_biosynthetic_process|secondary_alcohol_biosynthetic_process|cholesterol_biosynthesis",

  "Sterol / cholesterol metabolism" =
    "sterol_metabolic_process|cholesterol_metabolic_process|cholesterol_homeostasis|sterol_homeostasis",

  "Acetyl-CoA metabolism" =
    "acetyl_coa",

  "Isoprenoid metabolism" =
    "isoprenoid|terpenoid|farnesyl|prenyl",

  "Prostanoid metabolism" =
    "prostanoid|prostaglandin|icosanoid|eicosanoid|arachidonic_acid",

  "Extracellular matrix organization" =
    "extracellular_matrix_organization|extracellular_structure_organization|collagen_fibril|external_encapsulating_structure_organization",

  "Oxidative stress / detoxification" =
    "response_to_oxidative_stress|reactive_oxygen_species_metabolic|cellular_detoxification|detoxification|xenobiotic_metabolic_process|glutathione_metabolic",

  "Bile acid synthesis / transport" =
    "bile_acid"
)


## -----------------------------------------------------------------------------
## 5. Collapse GSEA output onto the named programs
## -----------------------------------------------------------------------------
## For each program: how many matching sets were tested, how many reached FDR,
## and the details of the best-scoring set (NES, FDR, leading edge).

.leading_genes <- function(core, n = 15) {
  vapply(strsplit(as.character(core), "/", fixed = TRUE),
         function(g) paste(utils::head(g, n), collapse = ", "),
         character(1))
}

.empty_program_row <- function(prog) {
  tibble::tibble(
    program       = prog,
    n_sets_tested = 0L,
    n_sets_FDR    = 0L,
    gene_set      = NA_character_,
    setSize       = NA_integer_,
    NES           = NA_real_,
    pvalue        = NA_real_,
    FDR           = NA_real_,
    n_leading     = NA_integer_,
    leading_edge  = NA_character_
  )
}

summarise_programs <- function(gsea_df,
                               patterns   = program_patterns,
                               fdr_cut    = 0.05,
                               n_sets     = 1,
                               n_leading  = 15) {

  if (!nrow(gsea_df)) {
    return(dplyr::bind_rows(lapply(names(patterns), .empty_program_row)))
  }

  g <- gsea_df %>% dplyr::mutate(.key = stringr::str_to_lower(ID))

  rows <- lapply(names(patterns), function(prog) {

    hits <- g %>% dplyr::filter(stringr::str_detect(.key, patterns[[prog]]))
    if (!nrow(hits)) return(.empty_program_row(prog))

    n_tested <- nrow(hits)
    n_fdr    <- sum(hits$p.adjust < fdr_cut, na.rm = TRUE)

    hits %>%
      dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) %>%
      dplyr::slice_head(n = n_sets) %>%
      dplyr::transmute(
        program       = prog,
        n_sets_tested = n_tested,
        n_sets_FDR    = n_fdr,
        gene_set      = ID,
        setSize       = as.integer(setSize),
        NES           = round(NES, 3),
        pvalue        = signif(pvalue, 3),
        FDR           = signif(p.adjust, 3),
        n_leading     = lengths(strsplit(core_enrichment, "/", fixed = TRUE)),
        leading_edge  = .leading_genes(core_enrichment, n_leading)
      )
  })

  dplyr::bind_rows(rows)
}


## -----------------------------------------------------------------------------
## 6. Arrow notation + whole-cohort vs normal-ALT comparison
## -----------------------------------------------------------------------------

arrowize <- function(NES, FDR, fdr_cut = 0.05, strong = 1.75) {
  ifelse(
    is.na(NES) | is.na(FDR), "?",
    ifelse(FDR >= fdr_cut, "\u00b1",                     # tested, not significant
      ifelse(NES > 0,
             ifelse(abs(NES) >= strong, "\u2191\u2191", "\u2191"),
             ifelse(abs(NES) >= strong, "\u2193\u2193", "\u2193"))))
}

.interpret <- function(nes_w, fdr_w, nes_n, fdr_n, fdr_cut, delta = 0.30) {

  sig_w <- !is.na(fdr_w) & fdr_w < fdr_cut
  sig_n <- !is.na(fdr_n) & fdr_n < fdr_cut

  dplyr::case_when(
    is.na(nes_w) & is.na(nes_n)              ~ "Not represented in either cohort",
    is.na(nes_n)                             ~ "Not represented in normal ALT",
    is.na(nes_w)                             ~ "Not represented in whole cohort",
    sig_w & sig_n & sign(nes_w) != sign(nes_n) ~ "Discordant direction",
    sig_w & sig_n & abs(nes_n) - abs(nes_w) >  delta ~ "Preserved - more prominent in normal ALT",
    sig_w & sig_n & abs(nes_w) - abs(nes_n) >  delta ~ "Preserved - attenuated in normal ALT",
    sig_w & sig_n                            ~ "Shared biology",
    sig_w & !sig_n                           ~ "Reduced / absent in normal ALT",
    !sig_w & sig_n                           ~ "Newly enriched in normal ALT",
    TRUE                                     ~ "Not enriched in either cohort"
  )
}

compare_programs <- function(prog_whole, prog_normal,
                             fdr_cut = 0.05,
                             delta   = 0.30) {

  dplyr::full_join(
    prog_whole  %>% dplyr::select(program, NES, FDR, gene_set, leading_edge),
    prog_normal %>% dplyr::select(program, NES, FDR, gene_set, leading_edge),
    by = "program", suffix = c("_whole", "_normal")
  ) %>%
    dplyr::mutate(
      `Whole cohort`   = sprintf("%s  (NES %s / FDR %s)",
                                 arrowize(NES_whole, FDR_whole, fdr_cut),
                                 ifelse(is.na(NES_whole), "-", format(round(NES_whole, 2), nsmall = 2)),
                                 ifelse(is.na(FDR_whole), "-", format(signif(FDR_whole, 2)))),
      `Normal ALT`     = sprintf("%s  (NES %s / FDR %s)",
                                 arrowize(NES_normal, FDR_normal, fdr_cut),
                                 ifelse(is.na(NES_normal), "-", format(round(NES_normal, 2), nsmall = 2)),
                                 ifelse(is.na(FDR_normal), "-", format(signif(FDR_normal, 2)))),
      Interpretation   = .interpret(NES_whole, FDR_whole, NES_normal, FDR_normal, fdr_cut, delta),
      dNES             = round(abs(NES_normal) - abs(NES_whole), 2)
    ) %>%
    dplyr::select(`Biological program` = program,
                  `Whole cohort`, `Normal ALT`, dNES, Interpretation,
                  gene_set_whole, gene_set_normal)
}


## -----------------------------------------------------------------------------
## 7. Metabolon annotation of a limma result table
## -----------------------------------------------------------------------------
## Detects the Metabolon column names rather than hard-coding them, because the
## sheet layout varies between exports.

metabolon_cols <- function(chem) {
  nm  <- colnames(chem)
  low <- stringr::str_to_lower(nm)
  pick <- function(pat) {
    i <- which(stringr::str_detect(low, pat))
    if (!length(i)) NA_character_ else nm[i[1]]
  }
  list(
    name  = pick("chemical.?name|^biochemical|^metabolite$|^short.?name"),
    super = pick("super.?pathway"),
    sub   = pick("sub.?pathway"),
    comp  = pick("^comp.?id")
  )
}

annotate_metabolites <- function(res_tab, chemical_details, data_is_log2 = TRUE) {

  tab <- if (is.data.frame(res_tab)) res_tab else res_tab$results
  cn  <- metabolon_cols(chemical_details)

  if (is.na(cn$name) || is.na(cn$sub)) {
    stop("Could not detect chemical-name / sub-pathway columns in chemical_details. ",
         "Columns are: ", paste(colnames(chemical_details), collapse = ", "))
  }

  map <- as.data.frame(chemical_details) %>%
    dplyr::transmute(
      metabolite    = as.character(.data[[cn$name]]),
      Super_pathway = if (is.na(cn$super)) NA_character_ else as.character(.data[[cn$super]]),
      Sub_pathway   = as.character(.data[[cn$sub]])
    ) %>%
    dplyr::distinct(metabolite, .keep_all = TRUE)

  tab %>%
    tibble::rownames_to_column("metabolite") %>%
    dplyr::left_join(map, by = "metabolite") %>%
    dplyr::mutate(
      Super_pathway = dplyr::coalesce(Super_pathway, "Unannotated"),
      Sub_pathway   = dplyr::coalesce(Sub_pathway,   "Unannotated"),
      direction     = ifelse(logFC > 0, "\u2191", "\u2193"),
      FC            = if (data_is_log2) 2^logFC else NA_real_
    )
}


## -----------------------------------------------------------------------------
## 8. Summarise significant metabolites by Metabolon pathway
## -----------------------------------------------------------------------------
## No enrichment test - just the descriptive summary requested:
## n significant, direction, strongest metabolites, median fold change.

summarise_metabolon <- function(annot,
                                level   = c("Sub_pathway", "Super_pathway"),
                                fdr_cut = 0.05,
                                n_top   = 3,
                                data_is_log2 = TRUE) {

  level <- match.arg(level)

  tested <- annot %>%
    dplyr::group_by(.pw = .data[[level]]) %>%
    dplyr::summarise(n_tested = dplyr::n(), .groups = "drop")

  sig <- annot %>% dplyr::filter(!is.na(adj.P.Val), adj.P.Val < fdr_cut)

  if (!nrow(sig)) {
    warning("No metabolite reached FDR < ", fdr_cut)
    return(tibble::tibble())
  }

  sig %>%
    dplyr::group_by(.pw = .data[[level]]) %>%
    dplyr::arrange(dplyr::desc(abs(logFC)), .by_group = TRUE) %>%
    dplyr::summarise(
      Super_pathway = dplyr::first(Super_pathway),
      n_sig         = dplyr::n(),
      n_up          = sum(logFC > 0),
      n_down        = sum(logFC < 0),
      direction     = dplyr::case_when(
                        all(logFC > 0) ~ "\u2191 all up",
                        all(logFC < 0) ~ "\u2193 all down",
                        sum(logFC > 0) > sum(logFC < 0) ~ "mostly \u2191",
                        sum(logFC < 0) > sum(logFC > 0) ~ "mostly \u2193",
                        TRUE ~ "mixed"),
      median_logFC  = round(stats::median(logFC), 3),
      median_FC     = if (data_is_log2) round(2^stats::median(logFC), 2) else NA_real_,
      max_abs_logFC = round(logFC[which.max(abs(logFC))], 3),
      best_FDR      = signif(min(adj.P.Val, na.rm = TRUE), 3),
      strongest     = paste(sprintf("%s (%s%.2f)",
                                    utils::head(metabolite, n_top),
                                    utils::head(direction, n_top),
                                    utils::head(logFC, n_top)),
                            collapse = "; "),
      .groups = "drop"
    ) %>%
    dplyr::left_join(tested, by = ".pw") %>%
    dplyr::mutate(pct_sig = round(100 * n_sig / n_tested, 1)) %>%
    dplyr::arrange(dplyr::desc(n_sig), best_FDR) %>%
    dplyr::rename(!!level := .pw) %>%
    dplyr::select(dplyr::any_of("Super_pathway"), dplyr::all_of(level),
                  n_sig, n_tested, pct_sig, n_up, n_down, direction,
                  median_logFC, median_FC, max_abs_logFC, best_FDR, strongest)
}


## -----------------------------------------------------------------------------
## 9. Metabolomics comparison: whole cohort vs normal ALT, by pathway
## -----------------------------------------------------------------------------

compare_metabolon <- function(annot_whole, annot_normal,
                              level   = "Sub_pathway",
                              fdr_cut = 0.05) {

  side <- function(annot, tag) {
    annot %>%
      dplyr::group_by(.pw = .data[[level]]) %>%
      dplyr::summarise(
        n_tested = dplyr::n(),
        n_sig    = sum(!is.na(adj.P.Val) & adj.P.Val < fdr_cut),
        med_logFC_sig = ifelse(n_sig > 0,
                               round(stats::median(logFC[!is.na(adj.P.Val) &
                                                           adj.P.Val < fdr_cut]), 3),
                               NA_real_),
        best_FDR = ifelse(all(is.na(adj.P.Val)), NA_real_,
                          signif(min(adj.P.Val, na.rm = TRUE), 3)),
        .groups = "drop"
      ) %>%
      dplyr::rename_with(~ paste0(.x, "_", tag), -.pw)
  }

  dplyr::full_join(side(annot_whole, "whole"),
                   side(annot_normal, "normal"), by = ".pw") %>%
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("n_sig_"), ~ tidyr::replace_na(.x, 0L)),
      `Whole cohort` = dplyr::case_when(
        n_sig_whole == 0 ~ "\u00b1",
        med_logFC_sig_whole > 0 & n_sig_whole >= 3 ~ "\u2191\u2191",
        med_logFC_sig_whole > 0 ~ "\u2191",
        n_sig_whole >= 3 ~ "\u2193\u2193",
        TRUE ~ "\u2193"),
      `Normal ALT` = dplyr::case_when(
        n_sig_normal == 0 ~ "\u00b1",
        med_logFC_sig_normal > 0 & n_sig_normal >= 3 ~ "\u2191\u2191",
        med_logFC_sig_normal > 0 ~ "\u2191",
        n_sig_normal >= 3 ~ "\u2193\u2193",
        TRUE ~ "\u2193"),
      Interpretation = dplyr::case_when(
        n_sig_whole > 0 & n_sig_normal > 0 &
          abs(med_logFC_sig_normal) > abs(med_logFC_sig_whole) ~ "Shared - stronger in normal ALT",
        n_sig_whole > 0 & n_sig_normal > 0 ~ "Shared biology",
        n_sig_whole > 0 & n_sig_normal == 0 ~ "Lost in normal ALT",
        n_sig_whole == 0 & n_sig_normal > 0 ~ "Emerges in normal ALT",
        TRUE ~ "Not significant in either")
    ) %>%
    dplyr::arrange(dplyr::desc(n_sig_normal + n_sig_whole)) %>%
    dplyr::rename(`Metabolic pathway` = .pw)
}
