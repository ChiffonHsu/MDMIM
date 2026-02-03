#Annotate using ms2_database
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/ms1_annotation/object_pos3")
object_pos3_ms2 <-
  mutate_ms2(
    object = object_pos3,
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/rawdata/"
  )

load("~/compound_database_shenlab/MS2/gnps_ms2.rda")
load("~/compound_database_shenlab/MS2/hmdb_ms2.rda")
load("~/compound_database_shenlab/MS2/massbank_ms2.rda")
load("~/compound_database_shenlab/MS2/metlin_ms2.rda")
load("~/compound_database_shenlab/MS2/mona_ms2.rda")
load("~/compound_database_shenlab/MS2/mpsnyder_rplc_ms2.rda")
load("~/compound_database_shenlab/MS2/nist_ms2.rda")

load("~/compound_database_shenlab/MS2/gnps_ms2.rda")
object_pos3_ms2 <-
  annotate_metabolites(
    object = object_pos3_ms2,
    database = gnps_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

hmdb_ms2
object_pos3_ms2 <-
  annotate_metabolites(
    object = object_pos3_ms2,
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

massbank_ms2
object_pos3_ms2 <-
  annotate_metabolites(
    object = object_pos3_ms2,
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

metlin_ms2
object_pos3_ms2 <-
  annotate_metabolites(
    object = object_pos3_ms2,
    database = metlin_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

mona_ms2
object_pos3_ms2 <-
  annotate_metabolites(
    object = object_pos3_ms2,
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

mpsnyder_rplc_ms2
object_pos3_ms2 <-
  annotate_metabolites(
    object = object_pos3_ms2,
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

nist_ms2
object_pos3_ms2 <-
  annotate_metabolites(
    object = object_pos3_ms2,
    database = nist_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

variable_info_pos <-
  extract_variable_info(object = object_pos3_ms2)
head(variable_info_pos)
table(variable_info_pos$Level)
# 1     2     3 
# 39   443 24339 
table(variable_info_pos$Database)
# GNPS_2022-04-29              HMDB_2022-04-11          MassBank_2022-04-27              METLIN_20220425 
# 54                           30                           23                           61 
# Michael_Snyder_RPLC_20220424              MoNA_2022-04-27                NIST_20220425               Shen-Lab_0.0.1 
# 39                           27                          248                        24339 
sum(!is.na(variable_info_pos$Compound.name)) #24821
sum(is.na(variable_info_pos$Compound.name)) #323
save(object_pos3_ms2, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/object_pos3_ms2")
########################################################################################################################

load("2_data/preliminary_determination_bacterial_metabolomics/MS2/NEG/object_neg4")
object_neg4_ms2 <-
  annotate_metabolites(
    object = object_neg4,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = gnps_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )
object_neg4_ms2 <-
  annotate_metabolites(
    object = object_neg4_ms2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

object_neg4_ms2 <-
  annotate_metabolites(
    object = object_neg4_ms2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )
object_neg4_ms2 <-
  annotate_metabolites(
    object = object_neg4_ms2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = metlin_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )
object_neg4_ms2 <-
  annotate_metabolites(
    object = object_neg4_ms2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )
object_neg4_ms2 <-
  annotate_metabolites(
    object = object_neg4_ms2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )
object_neg4_ms2 <-
  annotate_metabolites(
    object = object_neg4_ms2,
    ms1.match.ppm = 15,
    polarity = "negative",
    column = "rp",
    database = nist_ms2,
    based_on = c("ms1", "ms2"),
    adduct.table = NULL
  )

variable_info_neg <-
  extract_variable_info(object = object_neg4_ms2)
head(variable_info_neg)
table(variable_info_neg$Level)
table(variable_info_neg$Database)
sum(!is.na(variable_info_neg$Compound.name)) #1101
sum(is.na(variable_info_neg$Compound.name)) #213
save(object_neg4_ms2, file = "2_data/preliminary_determination_bacterial_metabolomics/MS2/NEG/object_neg4_ms2")
