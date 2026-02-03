load("~/Xu_project/25-11-30/2_data/amide_plasma/amide_pos_plasma/Result/Annotation/object_amide_pos3")
load("~/Xu_project/25-11-30/2_data/amide_plasma/amide_neg_plasma/Result/Annotation/object_amide_neg3")
load("~/Xu_project/25-11-30/2_data/c18_plasma/c18_plasma_pos/Result/Annotation/object_rplc_pos3")
load("~/Xu_project/25-11-30/2_data/c18_plasma/c18_plasma_neg/Result/Annotation/object_rplc_neg3")

object_amide_pos4 <- 
  object_amide_pos3 %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2)
object_amide_neg4 <- 
  object_amide_neg3 %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2)
object_rplc_pos4 <- 
  object_rplc_pos3 %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2)
object_rplc_neg4 <- 
  object_rplc_neg3 %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2)

hilic_pos <- object_amide_pos4
hilic_neg <- object_amide_neg4
rplc_pos <- object_rplc_pos4
rplc_neg <- object_rplc_neg4

head(colnames(hilic_pos), 19)
# #[1] "amide pos blank 41 A" "amide pos blank 42 A" "amide pos blank 43 A" "amide pos blank 44 A"
# [5] "Amide pos FAB2 10 A"  "Amide pos FAB2 11 A"  "Amide pos FAB2 12 A"  "Amide pos FAB2 13 A" 
# [9] "Amide pos FAB2 6 A"   "Amide pos FAB2 7 A"   "Amide pos FAB2 8 A"   "Amide pos FAFA 5 A"  
# [13] "Amide pos FAFA 6 A"   "amide QC 1"           "amide QC 2"           "Amide QC 3"          
# [17] "Amide QC 4"           "Amide QC 5"           "Amide QC 6"
head(colnames(hilic_neg), 19)
# [1] "amide neg blank 41" "amide neg blank 42" "amide neg blank 43" "amide neg blank 44"
# [5] "Amide neg FAB2 10"  "Amide neg FAB2 11"  "Amide neg FAB2 12"  "Amide neg FAB2 13" 
# [9] "Amide neg FAB2 6"   "Amide neg FAB2 7"   "Amide neg FAB2 8"   "Amide neg FAFA 5"  
# [13] "Amide neg FAFA 6"   "amide QC 1"         "amide QC 2"         "Amide QC 3"        
# [17] "Amide QC 4"         "Amide QC 5"         "Amide QC 6"
head(colnames(rplc_pos), 19)
#  [1] "C18 pos blank 41" "C18 pos blank 42" "C18 pos blank 43" "C18 pos blank 44" "C18 pos FAB2 10" 
# [6] "C18 pos FAB2 11"  "C18 pos FAB2 12"  "C18 pos FAB2 13"  "C18 pos FAB2 6"   "C18 pos FAB2 7"  
# [11] "C18 pos FAB2 8"   "C18 pos FAFA 5"   "C18 pos FAFA 6"   "C18 pos QC 1"     "C18 pos QC 2"    
# [16] "C18 pos QC 3"     "C18 pos QC 4"     "C18 pos QC 5"     "C18 pos QC 6" 
head(colnames(rplc_neg), 19)
# [1] "C18 neg blank 38" "C18 neg blank 42" "C18 neg blank 43" "C18 neg blank 44" "C18 neg FAB2 10" 
# [6] "C18 neg FAB2 11"  "C18 neg FAB2 12"  "C18 neg FAB2 13"  "C18 neg FAB2 6"   "C18 neg FAB2 7"  
# [11] "C18 neg FAB2 8"   "C18 neg FAFA 5"   "C18 neg FAFA 6"   "C18 neg QC 1"     "C18 neg QC 2"    
# [16] "C18 neg QC 3"     "C18 neg QC 4"     "C18 neg QC 5"     "C18 neg QC 6" 

head(colnames(hilic_pos@expression_data), 19)
head(hilic_pos@sample_info$sample_id, 19)


#make sure all sample_ids are the same across
clean_sample_id <- function(x) {
  
  # fix known typo
  x <- gsub("^C18 neg blank 38$", "C18 neg blank 41", x, ignore.case = TRUE)
  
  # remove platform labels at start
  x <- gsub("(?i)^(amide|C18)\\s+", "", x, perl = TRUE)
  
  # remove polarity at start OR in middle
  x <- gsub("(?i)^(pos|neg)\\s+|\\s+(pos|neg)\\s+", " ", x, perl = TRUE)
  
  # remove trailing " A"
  x <- gsub("\\s+A$", "", x)
  
  # clean extra spaces
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  
  x
}


#to make it safe, do this

#apply fxn to all 4 datasets
#HILIC_POS
old_ids <- colnames(hilic_pos@expression_data)
new_ids <- clean_sample_id(old_ids)

stopifnot(
  length(new_ids) == length(old_ids),
  !any(is.na(new_ids)),
  !any(new_ids == ""),
  length(unique(new_ids)) == length(new_ids)
)

colnames(hilic_pos@expression_data) <- new_ids
hilic_pos@sample_info$sample_id     <- new_ids

#HILIC_NEG
old_ids <- colnames(hilic_neg@expression_data)
new_ids <- clean_sample_id(old_ids)

stopifnot(
  length(new_ids) == length(old_ids),
  !any(is.na(new_ids)),
  !any(new_ids == ""),
  length(unique(new_ids)) == length(new_ids)
)

colnames(hilic_neg@expression_data) <- new_ids
hilic_neg@sample_info$sample_id     <- new_ids

#C18 POS
old_ids <- colnames(rplc_pos@expression_data)
new_ids <- clean_sample_id(old_ids)

stopifnot(
  length(new_ids) == length(old_ids),
  !any(is.na(new_ids)),
  !any(new_ids == ""),
  length(unique(new_ids)) == length(new_ids)
)

colnames(rplc_pos@expression_data) <- new_ids
rplc_pos@sample_info$sample_id     <- new_ids

#rplc_neg
old_ids <- colnames(rplc_neg@expression_data)
new_ids <- clean_sample_id(old_ids)

stopifnot(
  length(new_ids) == length(old_ids),
  !any(is.na(new_ids)),
  !any(new_ids == ""),
  length(unique(new_ids)) == length(new_ids)
)

colnames(rplc_neg@expression_data) <- new_ids
rplc_neg@sample_info$sample_id     <- new_ids

#CHECK BEFORE MERGING
stopifnot(
  identical(colnames(hilic_pos@expression_data),
            hilic_pos@sample_info$sample_id),
  
  identical(colnames(hilic_neg@expression_data),
            hilic_neg@sample_info$sample_id),
  
  identical(colnames(rplc_pos@expression_data),
            rplc_pos@sample_info$sample_id),
  
  identical(colnames(rplc_neg@expression_data),
            rplc_neg@sample_info$sample_id),
  
  identical(colnames(hilic_pos@expression_data),
            colnames(hilic_neg@expression_data)),
  
  identical(colnames(hilic_pos@expression_data),
            colnames(rplc_pos@expression_data)),
  
  identical(colnames(hilic_pos@expression_data),
            colnames(rplc_neg@expression_data))
)

colnames(hilic_pos@expression_data)


#merge datasets
#remember to keep polarity data before merging so that this information is not lost
hilic_pos@variable_info$polarity <- "positive"
hilic_neg@variable_info$polarity <- "negative"
rplc_pos@variable_info$polarity <- "positive"
rplc_neg@variable_info$polarity <- "negative"

library(massdataset)
hilic_pos@sample_info <- 
  dplyr::mutate(hilic_pos@sample_info, dplyr::across(everything(), as.character))

hilic_neg@sample_info <- 
  dplyr::mutate(hilic_neg@sample_info, dplyr::across(everything(), as.character))
keep_vars <- c("variable_id", "mz", "rt")

hilic_pos@variable_info <- hilic_pos@variable_info[, keep_vars, drop = FALSE]
hilic_neg@variable_info <- hilic_neg@variable_info[, keep_vars, drop = FALSE]

length(hilic_pos@sample_info) #5
ncol(hilic_pos@expression_data) #19
nrow(hilic_pos@sample_info) #19
length(hilic_neg@sample_info) #5
ncol(hilic_neg@expression_data) #19
nrow(hilic_neg@sample_info) #19

ncol(hilic_pos@sample_info) #5
nrow(hilic_pos@sample_info_note) #4

ncol(hilic_neg@sample_info) #5
nrow(hilic_neg@sample_info_note) #4
#mismatch found, therefore need to fix this before matching
hilic_pos@sample_info_note <- data.frame(
  name    = colnames(hilic_pos@sample_info),
  meaning = colnames(hilic_pos@sample_info),
  stringsAsFactors = FALSE
)

hilic_neg@sample_info_note <- data.frame(
  name    = colnames(hilic_neg@sample_info),
  meaning = colnames(hilic_neg@sample_info),
  stringsAsFactors = FALSE
)
#check if this passes 
stopifnot(
  ncol(hilic_pos@sample_info) == nrow(hilic_pos@sample_info_note),
  ncol(hilic_neg@sample_info) == nrow(hilic_neg@sample_info_note)
)

#another error after trying to merge, check varaible_info column no. shld be the same as variable_info_note's row number
hilic_pos@variable_info_note <- data.frame(
  name    = colnames(hilic_pos@variable_info),
  meaning = colnames(hilic_pos@variable_info),
  stringsAsFactors = FALSE
)

hilic_neg@variable_info_note <- data.frame(
  name    = colnames(hilic_neg@variable_info),
  meaning = colnames(hilic_neg@variable_info),
  stringsAsFactors = FALSE
)
#check again if sanity check passes
stopifnot(
  ncol(hilic_pos@variable_info) == nrow(hilic_pos@variable_info_note),
  ncol(hilic_neg@variable_info) == nrow(hilic_neg@variable_info_note)
)


hilic_merged <- merge_mass_dataset(
  x = hilic_pos,
  y = hilic_neg,
  sample_direction = "inner",       # keep only samples that exist in both
  variable_direction = "full",      # keep all features
  sample_by = "sample_id",
  variable_by = c("variable_id", "mz", "rt")
)

save(hilic_merged, file = "~/Xu_project/25-11-30/2_data/plasma_merged/hilic_merged")

#rplc_merged
rplc_pos@sample_info <- 
  dplyr::mutate(rplc_pos@sample_info, dplyr::across(everything(), as.character))

rplc_neg@sample_info <- 
  dplyr::mutate(rplc_neg@sample_info, dplyr::across(everything(), as.character))
keep_vars <- c("variable_id", "mz", "rt")

rplc_pos@variable_info <- rplc_pos@variable_info[, keep_vars, drop = FALSE]
rplc_neg@variable_info <- rplc_neg@variable_info[, keep_vars, drop = FALSE]

length(rplc_pos@sample_info) #5
ncol(rplc_pos@expression_data) #19
nrow(rplc_pos@sample_info) #19
length(rplc_neg@sample_info) #5
ncol(rplc_neg@expression_data) #19
nrow(rplc_neg@sample_info) #19

ncol(rplc_pos@sample_info) #5
nrow(rplc_pos@sample_info_note) #4

ncol(rplc_neg@sample_info) #5
nrow(rplc_neg@sample_info_note) #4
#mismatch found, therefore need to fix this before matching
rplc_pos@sample_info_note <- data.frame(
  name    = colnames(rplc_pos@sample_info),
  meaning = colnames(rplc_pos@sample_info),
  stringsAsFactors = FALSE
)

rplc_neg@sample_info_note <- data.frame(
  name    = colnames(rplc_neg@sample_info),
  meaning = colnames(rplc_neg@sample_info),
  stringsAsFactors = FALSE
)
#check if this passes 
stopifnot(
  ncol(rplc_pos@sample_info) == nrow(rplc_pos@sample_info_note),
  ncol(rplc_neg@sample_info) == nrow(rplc_neg@sample_info_note)
)

#another error after trying to merge, check varaible_info column no. shld be the same as variable_info_note's row number
rplc_pos@variable_info_note <- data.frame(
  name    = colnames(rplc_pos@variable_info),
  meaning = colnames(rplc_pos@variable_info),
  stringsAsFactors = FALSE
)

rplc_neg@variable_info_note <- data.frame(
  name    = colnames(rplc_neg@variable_info),
  meaning = colnames(rplc_neg@variable_info),
  stringsAsFactors = FALSE
)
#check again if sanity check passes
stopifnot(
  ncol(rplc_pos@variable_info) == nrow(rplc_pos@variable_info_note),
  ncol(rplc_neg@variable_info) == nrow(rplc_neg@variable_info_note)
)


rplc_merged <- merge_mass_dataset(
  x = rplc_pos,
  y = rplc_neg,
  sample_direction = "inner",       # keep only samples that exist in both
  variable_direction = "full",      # keep all features
  sample_by = "sample_id",
  variable_by = c("variable_id", "mz", "rt")
)

save(rplc_merged, file = "~/Xu_project/25-11-30/2_data/plasma_merged/rplc_merged")

#Merge hilic_merged and rplc_merged
load("~/Xu_project/25-11-30/2_data/plasma_merged/rplc_merged")
load("~/Xu_project/25-11-30/2_data/plasma_merged/hilic_merged")
#unique variable_id
hilic_merged@variable_info$variable_id <-
  paste0("HILIC_", hilic_merged@variable_info$variable_id)

rplc_merged@variable_info$variable_id <-
  paste0("RPLC_", rplc_merged@variable_info$variable_id)

rownames(hilic_merged@expression_data) <-
  hilic_merged@variable_info$variable_id

rownames(rplc_merged@expression_data) <-
  rplc_merged@variable_info$variable_id

stopifnot(
  identical(
    rownames(hilic_merged@expression_data),
    hilic_merged@variable_info$variable_id
  ),
  identical(
    rownames(rplc_merged@expression_data),
    rplc_merged@variable_info$variable_id
  )
)
library(massdataset)

all_merged <- merge_mass_dataset(
  x = rplc_merged,
  y = hilic_merged,
  sample_direction   = "inner",   # same samples
  variable_direction = "full",    # keep all features
  sample_by   = "sample_id",
  variable_by = c("variable_id", "mz", "rt")
)
#final validation
stopifnot(
  ncol(all_merged@expression_data) ==
    nrow(all_merged@sample_info),
  
  length(unique(all_merged@variable_info$variable_id)) ==
    nrow(all_merged@variable_info)
)

table(grepl("^HILIC_", all_merged@variable_info$variable_id))
table(grepl("^RPLC_",  all_merged@variable_info$variable_id))
save(
  all_merged,
  file = "~/Xu_project/25-11-30/2_data/plasma_merged/all_merged"
)
