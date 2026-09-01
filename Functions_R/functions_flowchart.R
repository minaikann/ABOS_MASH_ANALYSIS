## ---------------------------------------------------------------------------
## Flowchart style: one source of truth for every chart in the document
## ---------------------------------------------------------------------------
flow_style <- list(
  font     = "Helvetica",
  fontsize = 10,
  title_fs = 13,
  excl_fs  = 8.5,
  top_col  = "#1a1a1a",
  keep_col = "#2E7D32",
  excl_col = "#9e9e9e",
  edge_col = "#1a1a1a",
  box_w    = 2.0,    # was 2.5 - narrower boxes, text wraps to more lines
  excl_w   = 2.6,    # was 3.5
  char_w   = 0.55
)
## Thousands separator for counts inside labels: n_fmt(1545) -> "1,545"
n_fmt <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")


complete_rows <- function(df, cols) {
  df[stats::complete.cases(df[, cols, drop = FALSE]), , drop = FALSE]
}


wrap_lab <- function(x, width_in, fs, indent = "", char_w = flow_style$char_w) {
  n_chars <- max(8, floor(width_in * 72 / (fs * char_w)))
  lines   <- unlist(strsplit(paste(x, collapse = "\n"), "\n", fixed = TRUE))
  unlist(lapply(lines, function(l) {
    if (!nzchar(trimws(l))) return(" ")
    strwrap(l, width = n_chars, prefix = indent, initial = "")
  }))
}


flow_chart <- function(start, steps = list(), leaves = list(),
                       title = NULL, style = flow_style) {
  
  s <- style
  
  ## centred label: "\n" = line break
  lab <- function(x) {
    lines <- unlist(strsplit(paste(x, collapse = "\n"), "\n", fixed = TRUE))
    paste(lines, collapse = "\\n")
  }
  
  ## flush-left HTML label, one line per vector element (exclusion panels)
  lab_left <- function(x) {
    x <- unlist(strsplit(paste(x, collapse = "\n"), "\n", fixed = TRUE))
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;",  x, fixed = TRUE)
    x <- gsub(">", "&gt;",  x, fixed = TRUE)
    paste0("<", paste0(x, '<br align="left"/>', collapse = ""), ">")
  }
  
  
  box <- function(id, text, col, width = s$box_w, fs = s$fontsize) {
    sprintf("  %s [label = \"%s\", color = \"%s\", width = %s, fontsize = %s]",
            id, lab(wrap_lab(text, width, fs, char_w = s$char_w)), col, width, fs)
  }
  
  dot <- c(
    "digraph flow {",
    sprintf("  graph [rankdir = TB, splines = ortho, nodesep = 0.35, ranksep = 0.45%s]",
            if (is.null(title)) "" else
              sprintf(paste0(", label = \"%s\", labelloc = t, labeljust = c,",
                             " fontsize = %s, fontname = \"%s-Bold\""),
                      lab(title), s$title_fs, s$font)),
    sprintf(paste0("  node [shape = box, style = filled, fillcolor = white,",
                   " penwidth = 0.9, fontname = \"%s\", fontsize = %s]"),
            s$font, s$fontsize),
    sprintf(paste0("  edge [arrowhead = normal, arrowsize = 0.7,",
                   " color = \"%s\", penwidth = 0.8]"), s$edge_col),
    box("n0", start, s$top_col)
  )
  
  prev      <- "n0"
  rank_same <- character(0)
  
  for (i in seq_along(steps)) {
    keep <- sprintf("keep%d", i)
    excl <- sprintf("excl%d", i)
    gap  <- sprintf("gap%d",  i)
    
    dot <- c(dot,
             box(keep, steps[[i]]$keep, s$keep_col),
             sprintf("  %s [label = %s, color = \"%s\", width = %s, fontsize = %s]",
                     excl,
                     lab_left(wrap_lab(steps[[i]]$excl, s$excl_w, s$excl_fs,
                                       indent = "  ", char_w = s$char_w)),
                     s$excl_col, s$excl_w, s$excl_fs),
             sprintf("  %s [shape = point, width = 0.01, height = 0.01, label = \"\"]", gap),
             sprintf("  %s -> %s [arrowhead = none]", prev, gap),
             sprintf("  %s -> %s", gap, keep),
             sprintf("  %s -> %s [color = \"%s\"]", gap, excl, s$excl_col)
    )
    rank_same <- c(rank_same, sprintf("  { rank = same; %s; %s }", gap, excl))
    prev <- keep
  }
  
  if (length(leaves)) {
    if (is.character(leaves)) leaves <- lapply(leaves, function(l) list(label = l))
    ids <- sprintf("leaf%d", seq_along(leaves))
    
    ## a hub point turns the split into a square bracket rather than a fan
    if (length(leaves) > 1) {
      dot <- c(dot,
               "  hub [shape = point, width = 0.01, height = 0.01, label = \"\"]",
               sprintf("  %s -> hub [arrowhead = none]", prev))
      parent <- "hub"
    } else {
      parent <- prev
    }
    
    for (j in seq_along(leaves)) {
      dot <- c(dot, box(ids[j], leaves[[j]]$label, s$keep_col),
               sprintf("  %s -> %s", parent, ids[j]))
      
      kids <- leaves[[j]]$children          # boxes hanging under this leaf
      for (k in seq_along(kids)) {
        kid <- sprintf("leaf%d_%d", j, k)
        dot <- c(dot, box(kid, kids[k], s$keep_col),
                 sprintf("  %s -> %s", ids[j], kid))
      }
    }
    rank_same <- c(rank_same,
                   sprintf("  { rank = same; %s }", paste(ids, collapse = "; ")))
  }
  
  grViz(paste(c(dot, rank_same, "}"), collapse = "\n"))
}

## MASH / no-MASH leaf label with the sex split
grp_leaf <- function(df, grp) {
  d <- df[df$mash_lab == grp, , drop = FALSE]
  sprintf("%s (n = %s)\nFemale %s / Male %s",
          grp, n_fmt(nrow(d)),
          n_fmt(sum(d$sex_label == "Female")), n_fmt(sum(d$sex_label == "Male")))
}