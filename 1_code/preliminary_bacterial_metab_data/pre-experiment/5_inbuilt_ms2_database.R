#built ms2 database using MS2 DATA, Then get ms2 matching plots
 #add ms2 data into massdataset
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/ms1_annotation/object_pos3")
object_pos3 <-
  mutate_ms2(
    object = object_pos3,
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/preliminary_determination_bacterial_metabolomics/MS2/POS"
  )
object_pos4 <- object_pos3
save(object_pos4, file = "2_data/preliminary_determination_bacterial_metabolomics/MS2/POS/object_pos4")
###########################################################################################################load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/ms1_annotation/object_pos3")
object_neg3 <-
  mutate_ms2(
    object = object_neg3,
    column = "rp",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "2_data/preliminary_determination_bacterial_metabolomics/MS2/NEG"
  )
object_neg4 <- object_neg3
save(object_neg4, file = "2_data/preliminary_determination_bacterial_metabolomics/MS2/NEG/object_neg4")
