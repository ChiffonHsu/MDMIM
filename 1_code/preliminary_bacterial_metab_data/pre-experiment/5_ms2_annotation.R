#Annotate using ms2_database
load("2_data/preliminary_determination_bacterial_metabolomics/MS2/POS/object_pos4")
#add dummy adduct column
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/data_cleaning//object_pos2")
object_pos4_ms2 <-
  mutate_ms2(
    object = object_pos4,
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/preliminary_determination_bacterial_metabolomics/MS2/POS"
  )

load("~/compound_database_shenlab/MS2/gnps_ms2.rda")
load("~/compound_database_shenlab/MS2/hmdb_ms2.rda")
load("~/compound_database_shenlab/MS2/massbank_ms2.rda")
load("~/compound_database_shenlab/MS2/metlin_ms2.rda")
load("~/compound_database_shenlab/MS2/mona_ms2.rda")
load("~/compound_database_shenlab/MS2/mpsnyder_rplc_ms2.rda")
load("~/compound_database_shenlab/MS2/nist_ms2.rda")

load("~/compound_database_shenlab/MS2/gnps_ms2.rda")
object_pos4_ms2 <-
  annotate_metabolites(
    object = object_pos4_ms2,
    database = gnps_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

hmdb_ms2
object_pos4_ms2 <-
  annotate_metabolites(
    object = object_pos4_ms2,
    database = hmdb_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

massbank_ms2
object_pos4_ms2 <-
  annotate_metabolites(
    object = object_pos4_ms2,
    database = massbank_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

metlin_ms2
object_pos4_ms2 <-
  annotate_metabolites(
    object = object_pos4_ms2,
    database = metlin_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

mona_ms2
object_pos4_ms2 <-
  annotate_metabolites(
    object = object_pos4_ms2,
    database = mona_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

mpsnyder_rplc_ms2
object_pos4_ms2 <-
  annotate_metabolites(
    object = object_pos4_ms2,
    database = mpsnyder_rplc_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

nist_ms2
object_pos4_ms2 <-
  annotate_metabolites(
    object = object_pos4_ms2,
    database = nist_ms2,
    based_on = c("ms1", "ms2"),
    polarity = "positive",
    column = "rp",
    adduct.table = NULL
  )

variable_info_pos <-
  extract_variable_info(object = object_pos4_ms2)
head(variable_info_pos)
table(variable_info_pos$Level)
table(variable_info_pos$Database)
sum(!is.na(variable_info_pos$Compound.name)) #3452
sum(is.na(variable_info_pos$Compound.name)) #19
save(object_pos4_ms2, file = "2_data/preliminary_determination_bacterial_metabolomics/MS2/POS/object_pos4_ms2")
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
