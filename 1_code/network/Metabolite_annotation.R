library(metid)
load("compound_database_shenlab/MS1/merged_ms1_library.rda")
class("merged_ms1")
slotNames(merged_ms1)
merged_ms1@database.info
names(merged_ms1@spectra.data) # only mz
names(merged_ms1@spectra.info) # metadata and experimental info
names(merged_ms1@database.info) #Version, source, link, creater, email, RT
merged_ms1@spectra.info

library(tidymass)
load("MS2Library_ChiFang/In_Vitro/RPLC/POS/Result/data_cleaning/object_pos")
load("data_cleaning/NEG/object_pos2")
object_pos <-
  mutate_ms1(
    object = object_pos,
    column = "rp",
    polarity = "positive",
    ms1.match.mz.tol = 15,
    ms1.match.rt.tol = 30,
    path = "MS2Library_ChiFang/In_Vitro/HILIC/POS"
  )

extract_sample_info(object_pos2) %>%
  select(sample_id, batch, class, file_name) %>%
  head()
    
#Annotation
object_pos2
extract_ms1_database(object_pos2)
load("compound_database_shenlab/MS1/hmdb_ms1.rda")
hmdb_ms1
object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = hmdb_ms1,
    ## based_on = c("ms1", "rt"),
    polarity = "positive",
    column = "hilic",
    #adduct.table = NULL
  )

load("compound_database_shenlab/MS1/bloodexposome_ms1.rda")
bloodexposome_ms1
#bloodexposome_ms1 does not have RT and lab IDs, therefore create dummy columns to prevent crashing
# Add Lab.ID 
if (!"Lab.ID" %in% colnames(bloodexposome_ms1@spectra.info)) {
  bloodexposome_ms1@spectra.info$Lab.ID <- paste0("bloodex_", seq_len(nrow(bloodexposome_ms1@spectra.info)))
}

# Add RT
if (!"RT" %in% colnames(bloodexposome_ms1@spectra.info)) {
  bloodexposome_ms1@spectra.info$RT <- NA
}

object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = bloodexposome_ms1,
    polarity = "positive",
    column = "rp",
    
  )

load("compound_database_shenlab/MS1/t3db_ms1.rda")
t3db_ms1
object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = t3db_ms1,
    ## based_on = c("ms1", "rt"),
    polarity = "positive",
    column = "rp",
    #adduct.table = NULL
  )

load("compound_database_shenlab/MS1/kegg_ms1.rda")
kegg_ms1
object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = kegg_ms1,
    ## based_on = c("ms1", "rt"),
    polarity = "positive",
    column = "rp",
    #adduct.table = NULL
  )

save(object_pos2, file = "object_hilic_pos2.rda")
head(object_pos2@annotation_table)


#Run NEG for RPLC
load("MS2Library_ChiFang/In_Vitro/RPLC/NEG/Result/data_cleaning/object_pos")
load("data_cleaning/NEG/object_pos2")

#Annotation
object_pos2
extract_ms1_database(object_pos2)
load("compound_database_shenlab/MS1/hmdb_ms1.rda")
hmdb_ms1
object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = hmdb_ms1,
    ## based_on = c("ms1", "rt"),
    polarity = "posative",
    column = "rp",
    #adduct.table = NULL
  )

load("compound_database_shenlab/MS1/bloodexposome_ms1.rda")
bloodexposome_ms1
#bloodexposome_ms1 does not have RT and lab IDs, therefore create dummy columns to prevent crashing
# Add Lab.ID 
if (!"Lab.ID" %in% colnames(bloodexposome_ms1@spectra.info)) {
  bloodexposome_ms1@spectra.info$Lab.ID <- paste0("bloodex_", seq_len(nrow(bloodexposome_ms1@spectra.info)))
}

# Add RT
if (!"RT" %in% colnames(bloodexposome_ms1@spectra.info)) {
  bloodexposome_ms1@spectra.info$RT <- NA
}

object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = bloodexposome_ms1,
    polarity = "posative",
    column = "rp",
  )

load("compound_database_shenlab/MS1/t3db_ms1.rda")
t3db_ms1
object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = t3db_ms1,
    ## based_on = c("ms1", "rt"),
    polarity = "posative",
    column = "rp",
    #adduct.table = NULL
  )

load("compound_database_shenlab/MS1/kegg_ms1.rda")
kegg_ms1
object_pos2 <-
  metid::annotate_metabolites_mass_dataset(
    object = object_pos2,
    database = kegg_ms1,
    ## based_on = c("ms1", "rt"),
    polarity = "posative",
    column = "rp",
    #adduct.table = NULL
  )


save(object_neg2, file = "object_rplc_neg2.rda")

#Merge RPLC NEG and POS "object", and HILIC together, since I ran annotation without merging them
#Use rbind
#Items to merge: "object_pos2" and "object_neg2"

#library(massdataset)
# annotated_object_merged <- rbind_mass_dataset(object_pos2, object_neg2)
#Failed because diff columns, move on to version-agnostic manual approach
library(massdataset)
library(dplyr)
var_pos <- extract_variable_info(object_pos2)
var_neg <- extract_variable_info(object_neg2)
var_merged <- bind_rows(var_pos, var_neg) %>% distinct
sample_pos <- extract_sample_info(object_pos2)
sample_neg <- extract_sample_info(object_neg2)
sample_pos_character <- sample_pos %>% mutate(across(everything(), as.character))
sample_neg_character <- sample_neg %>% mutate(across(everything(), as.character))
sample_merged <- bind_rows(sample_pos_character, sample_neg_character) %>% distinct


#Merge expression data
expr_pos <- object_pos2@expression_data
expr_neg <- object_neg2@expression_data
all_vars <- union(rownames(expr_pos), rownames(expr_neg))
all_samples <- union(colnames(expr_pos), colnames(expr_neg))
expr_merged <- matrix(NA, nrow = length(all_vars), ncol = length(all_samples),
                      dimnames = list(all_vars, all_samples))
expr_merged[rownames(expr_pos), colnames(expr_pos)] <- expr_pos
#expr_merged[rownames(expr_neg), colnames(expr_neg)] <- expr_neg #error
expr_pos_mat <- as.matrix(extract_expression_data(object_pos2))
expr_neg_mat <- as.matrix(extract_expression_data(object_neg2))
all_vars <- union(rownames(expr_pos_mat), rownames(expr_neg_mat))
all_samples <- union(colnames(expr_pos_mat), colnames(expr_neg_mat))
expr_merged <- matrix(NA, nrow = length(all_vars), ncol = length(all_samples),
                      dimnames = list(all_vars, all_samples))
#Create merged object
# POS
expr_merged[match(rownames(expr_pos_mat), all_vars),
            match(colnames(expr_pos_mat), all_samples)] <- expr_pos_mat

# NEG
expr_merged[match(rownames(expr_neg_mat), all_vars),
            match(colnames(expr_neg_mat), all_samples)] <- expr_neg_mat

expr_merged_df <- as.data.frame(expr_merged)

object_merged <- object_pos2
object_merged@expression_data <- expr_merged_df
object_merged@sample_info <- sample_merged
object_merged@variable_info <- bind_rows(
  extract_variable_info(object_pos2),
  extract_variable_info(object_neg2)
) %>% distinct()

#Check that everything is in place
dim(object_merged@expression_data)
dim(object_merged@sample_info)
dim(object_merged@variable_info)
head(colnames(object_merged@expression_data))
head(rownames(object_merged@expression_data))
head(object_merged@variable_info$variable_id)

annotated_object_merged <- object_merged
save(annotated_object_merged, file = "annotated_object_merged.rda")

#Annotation results
head(extract_annotation_table(object = object_pos2))
variable_info_pos <-
  extract_variable_info(object = object_pos2)
head(variable_info_pos)
table(variable_info_pos$Level)
table(variable_info_pos$Database)