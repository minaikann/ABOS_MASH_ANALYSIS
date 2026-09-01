# Enrichment analysis (GSEA on the ranked |t| list)
enrichMsigdbFct = function(results, dataBase, symbEntrezid,
                           seuilAdjP = 0.05, seuilLogFC = 0.26, seuilEnrch = 0.05){
  
  ## Prepare data — recup entrezid
  dataEnrcht = merge(results, symbEntrezid, by.x = 0, by.y = "SYMBOL", all.x = TRUE)
  names(dataEnrcht)[1] = "SYMBOL"
  
  ## Compter combien on en perd d'importants
  isSign = dataEnrcht$adj.P.Val < seuilAdjP & abs(dataEnrcht$logFC) > seuilLogFC
  nLose  = sum(isSign & is.na(dataEnrcht$ENTREZID))
  nSign  = sum(isSign)
  
  ## Vecteur ranke
  ok = !is.na(dataEnrcht$ENTREZID)
  genesSorted = setNames(abs(dataEnrcht[ok, "t"]), dataEnrcht[ok, "ENTREZID"])
  genesSorted = sort(genesSorted, decreasing = TRUE)
  
  ## Enrichment
  set.seed(24)
  gseaEnrich = GSEA(genesSorted, TERM2GENE = dataBase, minGSSize = 30, maxGSSize = 500,
                    pvalueCutoff = seuilEnrch, scoreType = "pos", verbose = FALSE)
  resultEnrich = gseaEnrich@result
  
  ## Leading-edge : ENTREZID -> SYMBOL
  entrez2symb = setNames(as.character(symbEntrezid$SYMBOL), as.character(symbEntrezid$ENTREZID))
  leIds = strsplit(as.character(resultEnrich$core_enrichment), "/", fixed = TRUE)
  resultEnrich$n_leading = lengths(leIds)
  resultEnrich$leading_edge = vapply(leIds, function(x){   # vapply : OK si 0 gene set
    s = entrez2symb[x]; s[is.na(s)] = x[is.na(s)]          # garder les IDs non mappes
    paste(s, collapse = ", ")
  }, character(1))
  
  ## Table histogramme — [1:15,] donnerait des lignes NA si < 15 gene sets
  nHist = min(15, nrow(resultEnrich))
  resultHist = resultEnrich[seq_len(nHist), c("Description", "p.adjust", "pvalue")]
  if (nHist == 0) warning("Aucun gene set ne passe pvalueCutoff = ", seuilEnrch)
  resultHist$col = ifelse(resultHist$p.adjust < 0.05, "0",
                          ifelse(resultHist$pvalue < 0.05, "1", "2"))
  resultHist$cat = substr(resultHist$Description, 1, 35)   # alleger les etiquettes
  
  list(resultHist = resultHist, resultEnrich = resultEnrich)
}


# Plot the enrichment results (histogram)
printResultsEnrchHistFct = function(enrchResults, xlab,
                                    color = c("0" = "darkblue", "1" = "cyan", "2" = "grey70")){
  
  ggplot(enrchResults, aes(x = reorder(cat, -pvalue), y = -log10(pvalue), fill = col)) +
    geom_bar(stat = "identity") + coord_flip() +
    theme(legend.position = "none",
          axis.text  = element_text(size = 14),
          axis.title = element_text(size = 18)) +
    xlab(xlab) +
    scale_fill_manual(values = color)
}

