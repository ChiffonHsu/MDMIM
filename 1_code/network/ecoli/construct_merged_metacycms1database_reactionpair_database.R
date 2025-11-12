metacyc_ms1_nodes_all <- read.csv("~/compound_database_shenlab/MS1/metacyc_ms1_nodes_all.csv")
library(dplyr)
#Add in RT column with all NA,
#Add in Lab.ID. column as Shen_lab01, all the way till the end
#Change Compoiund_name to Compound.name
metacyc_ms1_nodes_all <- metacyc_ms1_nodes_all %>%
    mutate(Lab.ID = paste0("Shen_Lab", sprintf("%02d", row_number()) ) ) 

head(metacyc_ms1_nodes_all)
colnames(metacyc_ms1_nodes_all)
sum(is.na(metacyc_ms1_nodes_all$INCHIKEY_ID)) #20452
sum(is.na(metacyc_ms1_nodes_all$INCHIKEY_ID_all)) #20452

metacyc_ms1_nodes_all <- metacyc_ms1_nodes_all %>%
  dplyr::select(
    Lab.ID, Compound.name, mz, RT, CAS.ID, everything()
  )
colnames(metacyc_ms1_nodes_all)
# Export the spectra info as CSV
write.csv(
  metacyc_ms1_nodes_all, 
  file = "~/compound_database_shenlab/MS1/metabolite.info.csv",
  row.names = FALSE
)


library(tidymass)
library(metid)
new.path <- file.path("~/compound_database_shenlab/MS1")
metacyc_ms1_nodes_all <- construct_database(
  path = new.path,
  version = "0.0.1",
  metabolite.info.name = "metabolite.info.csv",
  source = "Shen-Lab",
  link = "https://www.shen-lab.org/",
  creater = "Chi-Fang Hsu",
  email = "shenxt1990@163.com",
  rt = TRUE,
  mz.tol = 15,
  rt.tol = 30,
  threads = 3
)
class(metacyc_ms1_nodes_all) #databaseClass
slotNames(metacyc_ms1_nodes_all)
save(metacyc_ms1_nodes_all, file = "~/compound_database_shenlab/MS1/metacyc_ms1_nodes_all.rda")
load("~/compound_database_shenlab/MS1/metacyc_ms1_nodes_all.rda")

library(dplyr)

metacyc_edges_filtered <- metacyc_edges_all %>%
  dplyr::select(from_inchi, to_inchi)
saveRDS(metacyc_edges_filtered, file = "~/compound_database_shenlab/metacyc/metacyc_edges_filtered.rds")
