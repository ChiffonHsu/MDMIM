library(dplyr)
library(stringr)
library(purrr)

#1. Parse metacyc 
###1.  ---- Parse compounds.dat ----
metacyc_cpd <- "~/compound_database_shenlab/metacyc/data/data/compounds.dat"

lines <- readLines(metacyc_cpd, encoding="latin1")
lines <- lines[!grepl("^#", lines)]

blocks <- split(lines, cumsum(lines == "//"))

parse_compound <- function(block){
  block <- block[block != "//"]
  
  out <- list(
    compound_id   = NA,
    compound_name = NA,
    smiles        = NA,
    inchi         = NA,
    chebi         = NA,
    kegg          = NA
  )
  
  for(line in block){
    
    if(startsWith(line, "UNIQUE-ID"))
      out$compound_id <- sub("UNIQUE-ID - ", "", line)
    
    if(startsWith(line, "COMMON-NAME"))
      out$compound_name <- sub("COMMON-NAME - ", "", line)
    
    if(startsWith(line, "SMILES"))
      out$smiles <- sub("SMILES - ", "", line)
    
    if(startsWith(line, "NON-STANDARD-INCHI"))
      out$inchi <- sub("NON-STANDARD-INCHI - ", "", line)
    
    # ChEBI mapping inside DBLINKS
    if(grepl("DBLINKS", line) && grepl("CHEBI", line))
      out$chebi <- sub('.*"([^"]+)".*', "\\1", line)
    
    # KEGG (LIGAND-CPD)
    if(grepl("DBLINKS", line) && grepl("LIGAND-CPD", line))
      out$kegg <- sub('.*"([^"]+)".*', "\\1", line)
  }
  
  as.data.frame(out, stringsAsFactors=FALSE)
}


metacyc_nodes_all <- bind_rows(lapply(blocks, parse_compound))
saveRDS(metacyc_nodes_all, "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_nodes_all.rds")
colnames(metacyc_nodes_all)
###2.  ---- Parse reactions.dat ----
reaction_file <- "~/compound_database_shenlab/metacyc/data/data/reactions.dat"

rx_lines <- readLines(reaction_file, encoding="latin1")
rx_lines <- rx_lines[!grepl("^#", rx_lines)]
rx_blocks <- split(rx_lines, cumsum(rx_lines == "//"))

parse_reaction_block <- function(block){
  block <- block[block != "//"]
  
  rxn_id <- NA
  LEFT <- c()
  RIGHT <- c()
  REV <- FALSE
  
  for(line in block){
    
    if(startsWith(line, "UNIQUE-ID"))
      rxn_id <- sub("UNIQUE-ID - ", "", line)
    
    if(startsWith(line, "LEFT"))
      LEFT <- c(LEFT, sub("LEFT - ", "", line))
    
    if(startsWith(line, "RIGHT"))
      RIGHT <- c(RIGHT, sub("RIGHT - ", "", line))
    
    if(startsWith(line, "REVERSIBLE?"))
      REV <- grepl("T", line)
  }
  
  if(is.na(rxn_id) || length(LEFT)==0 || length(RIGHT)==0)
    return(NULL)
  
  edges <- expand.grid(
    reaction_id = rxn_id,
    substrate   = LEFT,
    product     = RIGHT,
    stringsAsFactors = FALSE
  )
  
  if(REV){
    edges <- rbind(
      edges,
      data.frame(
        reaction_id = rxn_id,
        substrate   = RIGHT,
        product     = LEFT,
        stringsAsFactors = FALSE
      )
    )
  }
  
  edges
}

# parse_reaction_block <- function(block){
#   block <- block[block != "//"]
#   
#   LEFT <- c()
#   RIGHT <- c()
#   REV <- FALSE
#   
#   for(line in block){
#     if(startsWith(line, "LEFT"))
#       LEFT <- c(LEFT, sub("LEFT - ", "", line))
#     
#     if(startsWith(line, "RIGHT"))
#       RIGHT <- c(RIGHT, sub("RIGHT - ", "", line))
#     
#     if(startsWith(line, "REVERSIBLE?"))
#       REV <- grepl("T", line)
#   }
#   
#   if(length(LEFT)==0 || length(RIGHT)==0) return(NULL)
#   
#   edges <- expand.grid(
#     substrate = LEFT,
#     product  = RIGHT,
#     stringsAsFactors = FALSE
#   )
#   
#   if(REV){
#     edges <- rbind(edges,
#                    data.frame(substrate=RIGHT, product=LEFT))
#   }
#   
#   edges
# }

reaction_edges_df <- bind_rows(lapply(rx_blocks, parse_reaction_block))
###3.  ---- Map compound names-KEGG.ID ----
edges_kegg <- reaction_edges_df %>%
  left_join(metacyc_nodes_all %>% 
              dplyr::select(compound_id, kegg),
            by = c("substrate" = "compound_id")) %>%
  dplyr::rename(substrate_kegg = kegg) %>%
  left_join(metacyc_nodes_all %>% 
              dplyr::select(compound_id, kegg),
            by = c("product" = "compound_id")) %>%
  dplyr::rename(product_kegg = kegg)
saveRDS(edges_kegg,
        "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_edges_kegg.rds")
###4.  ---- Map compound names-smiles ----
edges_smiles <- reaction_edges_df %>%
  left_join(metacyc_nodes_all %>% 
              dplyr::select(compound_id, smiles),
            by = c("substrate" = "compound_id")) %>%
  dplyr::rename(substrate_smiles = smiles) %>%
  left_join(metacyc_nodes_all %>% 
              dplyr::select(compound_id, smiles),
            by = c("product" = "compound_id")) %>%
  dplyr::rename(product_smiles = smiles)
saveRDS(edges_smiles,
        "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_edges_smiles.rds")
###5.  ---- Map compound names-chebi ----
edges_chebi <- reaction_edges_df %>%
  left_join(metacyc_nodes_all %>% 
              dplyr::select(compound_id, chebi),
            by = c("substrate" = "compound_id")) %>%
  dplyr::rename(substrate_chebi = chebi) %>%
  left_join(metacyc_nodes_all %>% 
              dplyr::select(compound_id, chebi),
            by = c("product" = "compound_id")) %>%
  dplyr::rename(product_chebi = chebi)
saveRDS(edges_chebi,
        "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_edges_chebi.rds")
###5.  ---- Map compound names-inchi ----
#extract only valid inchi
valid_idx <- which(!is.na(metacyc_nodes_all$inchi) &
                     grepl("^InChI=", metacyc_nodes_all$inchi))

valid_inchi <- metacyc_nodes_all$inchi[valid_idx]

length(valid_inchi) #10570
writeLines(valid_inchi, "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_inchi_valid.txt")

    #change to correct inchi format
# Path to the converted InChIKeys
inchikey_path <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_inchikey_valid.txt"

# Read keys
inchikeys <- readLines(inchikey_path)

# Ensure equal length
length(valid_inchi)
length(inchikeys)

#filter only riws with inchi
metacyc_valid <- metacyc_nodes_all %>%
  filter(inchi %in% valid_inchi)
nrow(metacyc_valid) #10570
# Add column to your MetaCyc nodes
metacyc_valid$inchikey <- inchikeys
#check
colnames(metacyc_valid)
head(metacyc_valid[, c("inchi", "inchikey")])
# Save final MetaCyc nodes table
saveRDS(metacyc_valid, 
        "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_valid.rds")
edges_inchikey <- reaction_edges_df %>%
  left_join(metacyc_valid %>% 
              dplyr::select(compound_id, inchikey),
            by = c("substrate" = "compound_id")) %>%
  dplyr::rename(substrate_inchikey = inchikey) %>%
  left_join(metacyc_valid %>% 
              dplyr::select(compound_id, inchikey),
            by = c("product" = "compound_id")) %>%
  dplyr::rename(product_inchikey = inchikey)
saveRDS(edges_inchikey,
        "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_edges_inchikey.rds")
