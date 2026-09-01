comparisonsLimmaFct = function(fitLM, contr, design, nbGenes){
  
  contr <<- contr
  contrastMatrix = makeContrasts(contr, levels = design)
  
  fit2 = contrasts.fit(fitLM, contrastMatrix)
  
  fit2 = eBayes(fit2)
  
  mstat = topTable(fit2, coef = 1, number = nbGenes)
  mstat <- mstat[order(mstat$adj.P.Val),]
  
  return(list(mstat = mstat, contr = contr))
}

printResultsLimmaFct = function(resultComp, seuilAdjP = 0.05, seuilLogFC = 0.26,
                                topPrintHist = TRUE, topPrintVolc = TRUE,
                                color = c("blue", "grey", "red"), legendPos = "right", 
                                topMetabo = FALSE, genesToLabel = NULL){
  
  results = resultComp[["mstat"]]
  results = results[order(results$adj.P.Val),]
  
  nSign = length(which(results$adj.P.Val<=seuilAdjP))
  nSignFC = length(which(results$adj.P.Val<=seuilAdjP & abs(results$logFC) >= seuilLogFC))
  
  cat("There are", nSign, "differentially expressed features considering threshold", seuilAdjP, "on adjusted P-values, and", nSignFC, "when adding the threshold", seuilLogFC, "on absolute value of log2 Fold Change." ,"\n")
  
  if (nSign > 0){
    pval_safe <- max(results$P.Value[nSign], 1e-20)
    seuilpval = -log10(pval_safe)
  } else {
    seuilpval = -log10(seuilAdjP/nrow(results))
  }
  
  highlight_genes <- c("CHI3L1", "AKR1B10", "CCL20", "SPP1", "LPL",
                       "PTGDS", "CXCL10", "IL32", "VCAN", "LUM", "FAM151A", "AASS", "CYP2C19", "GRTP1",
                       "1-carboxyethyltyrosine","1-carboxyethylleucine","1-carboxyethylphenylalanine",
                       "1-carboxyethylvaline","1-carboxyethylisoleucine", "taurochenodeoxycholate", "taurocholate",
                       "cholate", "glycocholate")
  
  results$name <- ifelse(row.names(results) %in% highlight_genes, row.names(results), "")
  
  useAutoLabels = TRUE
  
  if (!is.null(genesToLabel)){
    # ---- Try labeling the user-supplied set of gene names ----
    matched = row.names(results) %in% genesToLabel
    if (any(matched)){
      results$name = ifelse(matched, row.names(results), "")
      useAutoLabels = FALSE
    } else {
      cat("None of the genesToLabel were found in the results — falling back to top-5-up/top-5-down labeling.\n")
    }
  }
  
  if (useAutoLabels) {
    # ---- Original behavior: auto top-5-up / top-5-down by adj.P.Val order ----
    cptPos = 0
    cptNeg = 0
    cpt = 1
    while ((cptPos < 5 | cptNeg < 5) & results[cpt,"adj.P.Val"] < seuilAdjP & cpt < nrow(results)){
      if (results[cpt,"logFC"] < -seuilLogFC & cptNeg < 5){
        results[cpt, "name"] = row.names(results)[cpt]
        cptNeg = cptNeg + 1
      } else if (results[cpt,"logFC"] > seuilLogFC & cptPos < 5){
        results[cpt, "name"] = row.names(results)[cpt]
        cptPos = cptPos + 1
      }
      cpt = cpt + 1
    }
  }
  
  results$Expression <- ifelse(results$logFC > seuilLogFC & results$adj.P.Val < seuilAdjP,
                               "Up regulated",
                               ifelse(results$logFC < -seuilLogFC & results$adj.P.Val < seuilAdjP,
                                      "Down regulated",
                                      "NS"))
  
  if (topPrintHist) hist(results$P.Value, main = 'Histogram of raw p-values', xlab = 'p-values')
  
  p = ggplot(results, aes(x = logFC, y = -log10(P.Value), col = Expression)) +
    geom_point() +
    geom_hline(yintercept = seuilpval, color = "grey60", linetype = "dashed") +
    geom_vline(xintercept = seuilLogFC, color = "grey60", linetype = "dashed") +
    geom_vline(xintercept = -seuilLogFC, color = "grey60", linetype = "dashed") +
    scale_color_manual(values = c("Down regulated" = color[1], "NS" = color[2], "Up regulated" = color[3])) +
    theme_classic() +
    geom_label_repel(data = subset(results, name != ""), aes(label = name),
                     fontface = "italic", size = 3, max.overlaps = 8) +
    xlab("log2FoldChange") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.position = legendPos,
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))

  
  if (topPrintVolc) print(p)
  
  # NEW: Add Fold Change column
  results$FC <- 2^results$logFC
  results$logFC = round(results$logFC, 6)
  results$FC = round(results$FC, 2)  # Round FC to 2 decimals
  results$adj.P.Val = round(results$adj.P.Val, 6)
  
  # NEW: keep all rows before filtering
  results_all <- results
  
  # UPDATED: Include FC in output
  if (topMetabo){
    results_sig = results[results$adj.P.Val < seuilAdjP & abs(results$logFC) >= seuilLogFC,
                          c("CHEMICAL_NAME", "logFC", "FC", "adj.P.Val")]
  } else {
    results_sig = results[results$adj.P.Val < seuilAdjP & abs(results$logFC) >= seuilLogFC,
                          c("logFC", "FC", "adj.P.Val")]
  }
  
  return(list(results = results_sig, results_all = results_all, volcano = p))
}

# UPDATED BAR PLOT with FC labels
create_top_plot <- function(results_obj, comparison_name) {
  top <- head(results_obj$results[order(results_obj$results$logFC, decreasing = TRUE), ], 10)
  top$label_pos <- ifelse(top$logFC > 0, 1.3, -0.3)
  
  # Color for upregulated (logFC > 0): blueviolet only for LS vs CTRL
  pos_color <- ifelse(comparison_name == "LS vs CTRL", "darkcyan", 
                      ifelse(comparison_name == "CM vs CTRL","blue" ,"blueviolet"))
  
  ggplot(top, aes(x = reorder(rownames(top), logFC), y = logFC, fill = logFC > 0)) +
    geom_col(width = 0.9) +
    geom_text(aes(label = round(logFC, 2), hjust = label_pos),
              position = position_dodge(width = 0.9),
              color = "white", size = 3.5) +
    coord_flip() +
    scale_fill_manual(values = c("FALSE" = "forestgreen", "TRUE" = pos_color),
                      guide = "none") +
    labs(x = "Metabolite", y = "log2Fold Change") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 14))
}


map_metabolite_names <- function(mstat_object, sample_info){ 
  # Map row names to chemical names safely
  sample_info$COMP.ID <- as.character(sample_info$COMP.ID)
  indices <- match(row.names(mstat_object), sample_info$COMP.ID)
  chemical_names <- sample_info$BIOCHEMICAL[indices]
  
  # Replace row names (fallback to ID if missing)
  row.names(mstat_object) <- ifelse(is.na(chemical_names) | chemical_names == "", 
                                    row.names(mstat_object), chemical_names)
  return(mstat_object)
}

map_chemical_names <- function(mstat, chemical_details) {
  indices <- match(row.names(mstat), chemical_details$COMP_ID)
  chemical_names <- chemical_details$CHEMICAL_NAME[indices]
  
  row.names(mstat) <- ifelse(
    is.na(chemical_names) | chemical_names == "", 
    row.names(mstat), 
    chemical_names
  )
  return(mstat)
}