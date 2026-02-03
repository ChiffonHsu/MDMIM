#Annotate using ms1_database
load("2_data/ms1_database/ms1_database.rda")
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/data_cleaning/object_pos2")
object_pos2 <-
  annotate_metabolites(
    object = object_pos2,
    ms1.match.ppm = 15,
    polarity = "positive",
    column = "rp",
    database = ms1_database,
    based_on = c("ms1"),
    adduct.table = NULL
  )
variable_info_pos <-
  extract_variable_info(object = object_pos2)
head(variable_info_pos)
table(variable_info_pos$Level)
table(variable_info_pos$Database)
sum(!is.na(variable_info_pos$Compound.name)) #24821
sum(is.na(variable_info_pos$Compound.name)) #323
object_pos3 <- object_pos2
save(object_pos3, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/ms1_annotation/object_pos3")
########################################################################################################################
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/ms1_database.rda")
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/data_cleaning/object_neg2")
object_neg2 <-
  annotate_metabolites(
    object = object_neg2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = ms1_database,
    based_on = c("ms1"),
    adduct.table = NULL
  )
variable_info_neg <-
  extract_variable_info(object = object_neg2)
head(variable_info_neg)
table(variable_info_neg$Level)
table(variable_info_neg$Database)
sum(!is.na(variable_info_neg$Compound.name)) #1100
sum(is.na(variable_info_neg$Compound.name)) #214
object_neg3 <- object_neg2
save(object_neg3, file = "2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/ms1_annotation/object_neg3")
