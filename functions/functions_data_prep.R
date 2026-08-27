# ============================================================
# functions_data_prep.R
# Helpers used while loading / aligning the raw data tables.
# ============================================================

# Replace COMP_ID row names with human-readable chemical names where available
map_chemical_names <- function(mstat, chemical_details) {
  idx <- match(row.names(mstat), chemical_details$COMP_ID)
  nm  <- chemical_details$CHEMICAL_NAME[idx]
  row.names(mstat) <- ifelse(is.na(nm) | nm == "", row.names(mstat), nm)
  mstat
}
