#load massdataset ms1_database
load("~/compound_database_shenlab/MS1/ms1_database.rda")
#check compounds
ms1_database@spectra.info
#CHECK FOR missing values in ms1_databse  (532488 obs. of 68 variables)
#Lab.ID          Compound.name                     mz                     RT 
#0                      1                   8378                 532488 
#CAS.ID                KEGG.ID                Formula               Synonyms 
#493835                 518221                   2673                  16313 
#Formula_all    Monoisotopic_weight   Compound_description        Database_source 
#2673                 253976                 292422                      0 
#HMDB.ID            HMDB_ID_all            KEGG_ID_all             CAS_ID_all 
#313910                 313910                 518221                 494053 
#INCHI_ID           INCHI_ID_all            INCHIKEY_ID        INCHIKEY_ID_all 
#19941                  19941                  22200                  22200 
#SMILES_ID          SMILES_ID_all           KEGG_DRUG_ID          CHEMSPIDER_ID 
#25185                  25185                 520408                 499594 
#DRUGBANK_ID               FOODB_ID    PUBCHEM_COMPOUND_ID   PUBCHEM_SUBSTANCE_ID 
#522705                 436942                 409951                 530791 
#CHEBI_ID           CHEBI_ID_all              CHEMBL_ID             PDB_CCD_ID 
#486303                 486303                 519995                 529443 
#X3DMET_ID             NIKKAJI_ID            KNAPSACK_ID           LIPIDMAPS_ID 
#527192                 518055                 527924                 526199 
#LIPIDBANK_ID              BIOCYC_ID                BIGG_ID     BIGG_IDENTIFIER_ID 
#532488                 528995                 529857                 525468 
#WIKIPEDIA_ID              METLIN_ID                T3DB_ID            REACTOME_ID 
#523443                 530907                 528861                 529757 
#MODELSEED_ID              MIMEDB_ID               LOTUS_ID             from_human 
#517250                 504979                 256716                      0 
#from_which_part          from_bacteria    from_which_bacteria       bacteria_ncbi_id 
#0                      0                     13                      0 
#bacteria_phylum         bacteria_class         bacteria_order        bacteria_family 
#0                      0                      0                      0 
#bacteria_genus       bacteria_species             from_plant       from_which_plant 
#0                      0                      0                     29 
#from_animal      from_which_animal       from_environment from_which_environment 
#0                     40                      0                      0 
#from_drug        from_which_drug              from_food        from_which_food 
#0                      0                      0                      0 

rowsum(ms1_database@spectra.info) #532488
nrow(metacyc_nodes_all) #26087

   #merge with metacyc reaction pairs using inchi_all as the primary compound
#library(dplyr)
#merge_inchi_from <- metacyc_nodes_all %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("from_", .)),
#            by = c("from_inchi" = "from_INCHI_ID_all"))
#merge_inchi_to <- metacyc_nodes_all %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("to_", .)),
#            by = c("to_inchi" = "to_INCHI_ID_all"))
     #join metacyc and ms1_databse by inchi first and then extract remaining rows for smiles later
#merged_edges_from_inchi <- metacyc_nodes_all %>%
#  left_join(ms1_database@spectra.info %>% 
#              dplyr::select(INCHI_ID_all, everything()) %>%
#              rename_with(~paste0("from_", .), -INCHI_ID_all),
#            by = c("from_inchi" = "INCHI_ID_all")) 

#sum(is.na(merged_edges_from_inchi$from_INCHI_ID_all))  #16800494
#sum(!is.na(merged_edges_from_inchi$from_INCHI_ID_all)) #3716124
    
    #merge with metacyc reaction pairs using inchi as the primary compound instead
#merge_inchi_unique_from <- metacyc_nodes_all %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("from_", .)),
#            by = c("from_inchi" = "from_INCHI_ID"))
#nrow(merge_inchi_unique_from) # 20516752
      #merge with metacyc reaction pairs using inchi and smiles togetheras the primary compound instead
#merge_inchismiles_unique_from <- metacyc_nodes_all %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("from_", .)),
#            by = c("from_inchi" = "from_INCHI_ID",
#                   "from_smiles" = "from_SMILES_ID"))
#nrow(merge_inchismiles_unique_from) # 7162953
#merge_inchismiles_unique_to <- metacyc_nodes_all %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("to_", .)),
#            by = c("to_inchi" = "to_INCHI_ID",
#                   "to_smiles" = "to_SMILES_ID"))
#nrow(merge_inchismiles_unique_to) # 7703561
       #merge using the unique pairing in metacyc
          #merge "from" into metacyc first, then add "to" into the edited metacyc with the "from" info
metacyc_ms1_nodes <- metacyc_nodes_all %>%
  left_join(ms1_database@spectra.info, 
            by = c("inchi" = "INCHI_ID"))

            nrow(metacyc_ms1_nodes) # 46446538
            sum(is.na(metacyc_ms1_nodes$inchi)) # 46422648
            metacyc_ms1_nodes_filtered <- metacyc_ms1_nodes %>%
              filter(!is.na(inchi))
            nrow(metacyc_ms1_nodes_filtered) #23890
            colnames(metacyc_ms1_nodes_filtered)
            merged_reaction_pairs <- metacyc_ms1_nodes_filtered
            
            saveRDS(merged_reaction_pairs, file = "~/compound_database_shenlab/MS1/merged_reaction_pairs.rds")
            write.csv(merged_reaction_pairs,
                      file = "~/compound_database_shenlab/MS1/merged_reaction_pairs.csv",
                      row.names = FALSE)
            colnames(merged_reaction_pairs)
            
            
             #           [1] "compound_id"            "compound_name"          "smiles"                 "inchi"                  "chebi"                  "Lab.ID"                
#            [7] "Compound.name"          "mz"                     "RT"                     "CAS.ID"                 "KEGG.ID"                "Formula"               
#            [13] "Synonyms"               "Formula_all"            "Monoisotopic_weight"    "Compound_description"   "Database_source"        "HMDB.ID"               
#            [19] "HMDB_ID_all"            "KEGG_ID_all"            "CAS_ID_all"             "INCHI_ID_all"           "INCHIKEY_ID"            "INCHIKEY_ID_all"       
#            [25] "SMILES_ID"              "SMILES_ID_all"          "KEGG_DRUG_ID"           "CHEMSPIDER_ID"          "DRUGBANK_ID"            "FOODB_ID"              
#            [31] "PUBCHEM_COMPOUND_ID"    "PUBCHEM_SUBSTANCE_ID"   "CHEBI_ID"               "CHEBI_ID_all"           "CHEMBL_ID"              "PDB_CCD_ID"            
#            [37] "X3DMET_ID"              "NIKKAJI_ID"             "KNAPSACK_ID"            "LIPIDMAPS_ID"           "LIPIDBANK_ID"           "BIOCYC_ID"             
#            [43] "BIGG_ID"                "BIGG_IDENTIFIER_ID"     "WIKIPEDIA_ID"           "METLIN_ID"              "T3DB_ID"                "REACTOME_ID"           
#            [49] "MODELSEED_ID"           "MIMEDB_ID"              "LOTUS_ID"               "from_human"             "from_which_part"        "from_bacteria"         
#            [55] "from_which_bacteria"    "bacteria_ncbi_id"       "bacteria_phylum"        "bacteria_class"         "bacteria_order"         "bacteria_family"       
#            [61] "bacteria_genus"         "bacteria_species"       "from_plant"             "from_which_plant"       "from_animal"            "from_which_animal"     
#            [67] "from_environment"       "from_which_environment" "from_drug"              "from_which_drug"        "from_food"              "from_which_food" 
           
            
             #check how heavily duplicated merged_from is
sum(is.na(metacyc_nodes_all$from_inchi)) #1026
length(unique(metacyc_ms1_nodes_filtered$inchi)) #8718
length(ms1_database@spectra.info$INCHI_ID) #532488
#too many duplicates, remove all na rows that didnt match any ms1 info
merged_from_filtered <- merged_from %>%
  filter(!is.na(from_inchi))
nrow(merged_from_filtered) #56951
merged_both <- merged_from_filtered %>%
  left_join(ms1_database@spectra.info %>% rename_with(~paste0("to_", .)),
            by = c("to_inchi" = "to_INCHI_ID",
                   "to_smiles" = "to_SMILES_ID"))
nrow(merged_both) #5038227
sum(is.na(merged_both$from)) #0
sum(is.na(merged_both$to)) #0
merged_both_filtered <- merged_both %>%
  filter(!is.na(to_inchi))
nrow(merged_both_filtered) #56078

merged_reaction_pairs <- merged_both_filtered

saveRDS(merged_reaction_pairs, file = "~/compound_database_shenlab/MS1/merged_reaction_pairs.rds")
write.csv(merged_reaction_pairs,
          file = "~/compound_database_shenlab/MS2/merged_reaction_pairs.csv",
          row.names = FALSE)
colnames(merged_reaction_pairs)

 #merged_edges_inchi <- merged_edges_from_inchi %>%
#  left_join(ms1_database@spectra.info %>%
#              dplyr::select(INCHI_ID_all, everything()) %>%
#              rename_with(~paste0("to_", .), -INCHI_ID_all),
#            by = c("to_inchi" = "INCHI_ID_all"))
#
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("to_", .)),
#            by = c("to_inchi" = "INCHI_ID_all"))
#for rows with missing inchi, merge by smiles
#nrow(merge_inchi_to %>% filter(is.na(to_INCHIKEY_ID_all))) #19182172
#nrow(merge_inchi_from %>% filter(is.na(from_INCHIKEY_ID_all))) #16800494
#from_missing <- merge_inchi_from %>% filter(is.na(from_INCHIKEY_ID_all))

#merge the remaining using smiles
#from_merged_smiles <- from_missing %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("from_", .)),
#            by = c( "from_smiles" = "from_SMILES_ID_all"),
#            suffix = c("", "_from_SMILES"))
    #check for duplicates of smiles id
#n_distinct(from_missing$from_smiles) #9289
#nrow(from_missing) #16800494
#n_distinct(ms1_database@spectra.info$SMILES_ID_all) #506538
#nrow(ms1_database@spectra.info) #532488
#             #check how many times each smiles appear in ms1_database
#smiles_dup_stats <- ms1_database@spectra.info %>% 
#  group_by(SMILES_ID_all) %>%
#  summarise(n = n()) %>%
#  arrange(desc(n))
#total_unique_smiles <- n_distinct(ms1_database@spectra.info$SMILES_ID_all)
#duplicated_smiles <- smiles_dup_stats %>% filter(n > 1)
#cat("Total unique SMILES:", total_unique_smiles, "\n") #3506538
#cat("Duplicated SMILES:", nrow(duplicated_smiles), "\n") #755 (~0.15% of all SMILES)
#head(duplicated_smiles, 20)
   #remove all NAs
#from_missing_nonNA <- from_missing %>% filter(!is.na(from_smiles))
   #only merge unique SMILES
#unique_smiles <- smiles_dup_stats %>% filter(n == 1) %>% pull(SMILES_ID_all)
#from_missing_unique <- from_missing_nonNA %>%
#  filter(from_smiles %in% unique_smiles)
#from_merged_smiles <- from_missing_unique %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("from_", .)),
#            by = c( "from_smiles" = "from_SMILES_ID_all"))
#for (col in names(from_merged_smiles)[endsWith(names(from_merged_smiles), "_from_SMILES")]) {
#  orig_col <- sub("_from_SMILES$", "", col)
#  from_merged_smiles[[orig_col]] <- coalesce(from_merged_smiles[[orig_col]], from_merged_smiles[[col]])
#}

#from_merged_smiles <- from_merged_smiles[, !endsWith(names(from_merged_smiles), "_from_SMILES")]
#merged_from_final <- merge_inchi_from %>%
#  filter(!is.na(from_INCHIKEY_ID_all)) %>%
#  bind_rows(from_merged_smiles)
#nrow(from_merged_smiles) #416536


 #repeat for to_INCHIKEY_ID_all
#to_missing <- merge_inchi_to %>% filter(is.na(to_INCHIKEY_ID_all))

#merge the remaining using smiles
#to_merged_smiles <- to_missing %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("to_", .)),
#            by = c( "to_smiles" = "to_SMILES_ID_all"),
#            suffix = c("", "_to_SMILES")) #failed cos too many rows
#check for duplicates of smiles id
#n_distinct(to_missing$to_smiles) #11386
#nrow(to_missing) #19182172
#n_distinct(ms1_database@spectra.info$SMILES_ID_all) #506538
#nrow(ms1_database@spectra.info) #532488
#check how many times each smiles appear in ms1_database
#smiles_dup_stats <- ms1_database@spectra.info %>% 
#  group_by(SMILES_ID_all) %>%
#  summarise(n = n()) %>%
#  arrange(desc(n))
#total_unique_smiles <- n_distinct(ms1_database@spectra.info$SMILES_ID_all)
#duplicated_smiles <- smiles_dup_stats %>% filter(n > 1)
#cat("Total unique SMILES:", total_unique_smiles, "\n") #506538
#cat("Duplicated SMILES:", nrow(duplicated_smiles), "\n") #755 (~0.15% of all SMILES)
#head(duplicated_smiles, 20)
#remove all NAs
#to_missing_nonNA <- to_missing %>% filter(!is.na(to_smiles))
#only merge unique SMILES
#unique_smiles <- smiles_dup_stats %>% filter(n == 1) %>% pull(SMILES_ID_all)
#to_missing_unique <- to_missing_nonNA %>%
#  filter(to_smiles %in% unique_smiles)
#to_merged_smiles <- to_missing_unique %>%
#  left_join(ms1_database@spectra.info %>% rename_with(~paste0("to_", .)),
#            by = c( "to_smiles" = "to_SMILES_ID_all"))

#for (col in names(to_merged_smiles)[endsWith(names(to_merged_smiles), "_to_SMILES")]) {
#  orig_col <- sub("_to_SMILES$", "", col)
#  to_merged_smiles[[orig_col]] <- coalesce(to_merged_smiles[[orig_col]], to_merged_smiles[[col]])
#}

#to_merged_smiles <- to_merged_smiles[, !endsWith(names(to_merged_smiles), "_to_SMILES")]
#merged_to_final <- merge_inchi_to %>%
#  filter(!is.na(to_INCHIKEY_ID_all)) %>%
#  bind_rows(to_merged_smiles)
#nrow(to_merged_smiles) #448996
 #FINALLY, merged to and from sets according to reaction pairs in metacyc_nodes_all
#ms1_unique_to <- ms1_database@spectra.info %>%
#  filter(!is.na(INCHI_ID_all)) %>%
#  distinct(INCHI_ID_all, .keep_all = TRUE) %>%
#  rename_with(~paste0("to_", .))
#merged_final_inchi <- merged_edges %>%
#  left_join(ms1_unique_to, by = c("to_inchi" = "to_INCHI_ID_all"))
 #identify missing rows
# Make sure columns align
#merged_final_edges <- merged_final_inchi

#nrow(merged_final_edges) #23427779
#colnames(merged_final_edges)
#saveRDS(merged_final_edges, file = "~/compound_database_shenlab/metacyc/merged_final_edges.rds")
