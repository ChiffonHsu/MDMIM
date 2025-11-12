#annotation
load("/public5/qifang/MS2Library_ChiFang/ms1_database.rda")
object_ecoli_sig <- readRDS("/public5/qifang/MS2Library_ChiFang/ecoli/object_ecoli_sig.rds")
library(metid)
library(tidymass)
library(readr)
load("/public5/qifang/MS2Library_ChiFang/ecoli/object_ecoli_sig.rds")

object_ecoli_sig_rplcneg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_ecoli_sig_rplcneg,
    database = ms1_database,
    polarity = "negative",
    column = "rp",
    #adduct.table = NULL
  )
object_ecoli_sig_annotated<- object_ecoli_sig
save(object_ecoli_sig, file = "object_ecoli_sig_annotated.rda")
saveRDS(object_ecoli_sig, file = "object_ecoli_sig_rplcneg_annot.rds")

object_ecoli_sig_rplcneg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_ecoli_sig_rplcneg,
    database = ms1_database,
    polarity = "negative",
    column = "rp",
    #adduct.table = NULL
  )
object_ecoli_sig_annotated<- object_ecoli_sig
save(object_ecoli_sig, file = "object_ecoli_sig_annotated.rda")
saveRDS(object_ecoli_sig, file = "object_ecoli_sig_annotated.rds")