library(r4projects)
setwd(get_project_wd())
rm(list = ls())
# ---- 1. Read MetaCyc data ----
# Adjust the paths to your MetaCyc files
enzrxn_lines <- readLines("~/compound_database_shenlab/metacyc/data/data/enzrxns.dat")
rxn_lines <- readLines("~/compound_database_shenlab/metacyc/data/data/reactions.dat")

# ---- 2. Helper function to split records ----
split_records <- function(lines) {
  idx <- which(grepl("^//", lines))
  starts <- c(1, idx[-length(idx)] + 1)
  records <- mapply(function(s, e) lines[s:e], starts, idx, SIMPLIFY = FALSE)
  records
}

# ---- 3. Parse reactions.dat ----
rxn_df <- do.call(rbind, lapply(rxn_records, function(r) {
  id   <- paste(sub("UNIQUE-ID - ", "", r[grepl("^UNIQUE-ID", r)]), collapse = ";")
  left <- paste(sub("LEFT - ", "", r[grepl("^LEFT", r)]), collapse = ";")
  right<- paste(sub("RIGHT - ", "", r[grepl("^RIGHT", r)]), collapse = ";")
  
  if (id   == "") id   <- NA
  if (left == "") left <- NA
  if (right== "") right<- NA
  
  data.frame(
    Reaction_ID = id,
    Left  = left,
    Right = right,
    stringsAsFactors = FALSE
  )
}))

# ---- 4. Parse enzrxns.dat ----
enz_df <- do.call(rbind, lapply(enz_records, function(r) {
  id   <- paste(sub("REACTION - ", "", r[grepl("^REACTION", r)]), collapse = ";")
  enz  <- paste(sub("ENZYME - ", "", r[grepl("^ENZYME", r)]), collapse = ";")
  ec   <- paste(sub("EC-NUMBER - ", "", r[grepl("^EC-NUMBER", r)]), collapse = ";")
  name <- paste(sub("COMMON-NAME - ", "", r[grepl("^COMMON-NAME", r)]), collapse = ";")
  
  if (id   == "") id   <- NA
  if (enz  == "") enz  <- NA
  if (ec   == "") ec   <- NA
  if (name == "") name <- NA
  
  data.frame(
    Reaction_ID = id,
    Enzyme_ID   = enz,
    EC_Number   = ec,
    Enzyme_Name = name,
    stringsAsFactors = FALSE
  )
}))

dim(enz_df)
head(enz_df)
# ---- 5. Merge both datasets ----
metacyc_ezmrxn_merged_df <- merge(rxn_df, enz_df, by = "Reaction_ID", all.x = TRUE)

head(metacyc_ezmrxn_merged_df)
save(metacyc_ezmrxn_merged_df, file = "~/compound_database_shenlab/metacyc/metacyc_ezmrxn_merged_df.rds")
