#parse enzyme reactions
# 1. Path to enzrxns.dat
enzrxns_file <- "~/compound_database_shenlab/metacyc/data/data/enzrxns.dat"

  # Read file
    lines <- readLines(enzrxns_file, encoding="latin1")
    lines <- lines[!grepl("^#", lines)]
    blocks <- split(lines, cumsum(lines == "//"))

# 2. Parse function
parse_enzrxn <- function(block){
  block <- block[block != "//"]
  
  out <- list(
    reaction_id = NA,
    enzyme_id   = NA,
    ec_number   = NA
  )
  
  for(line in block){
    if(startsWith(line, "REACTION"))
      out$reaction_id <- sub("REACTION - ", "", line)
    
    if(startsWith(line, "ENZYME"))
      out$ec_number <- sub("ENZYME - ", "", line)
    
    if(startsWith(line, "UNIQUE-ID"))
      out$enzyme_id <- sub("UNIQUE-ID - ", "", line)
  }
  
  as.data.frame(out, stringsAsFactors=FALSE)
}

enzrxns_all <- bind_rows(lapply(blocks, parse_enzrxn))
head(enzrxns_all)
  # 1 REACTION-DIRECTION - PHYSIOL-LEFT-TO-RIGHT ENZRXN-24513      MONOMER-19692
  # 2                                   RXN-1824 ENZRXN-18829      MONOMER-16303
  # 3                                  RXN-17029 ENZRXN-24283      MONOMER-19574
  # 4                                  RXN-10079 ENZRXN-16139      MONOMER-14434
  # 5                                   RXN-8166 ENZRXN-12985 G1G01-3828-MONOMER
  # 6 REACTION-DIRECTION - PHYSIOL-RIGHT-TO-LEFT ENZRXN-25525          CPLX-9193

# 3. parse proteins.dat to get ec numbers for monomer-
enzymes_file <- "~/compound_database_shenlab/metacyc/data/data/reactions.dat"
rxn_lines <- readLines(enzymes_file, encoding="latin1", warn=FALSE)
rxn_lines <- rxn_lines[!grepl("^#", rxn_lines)]
rxn_blocks <- split(rxn_lines, cumsum(rxn_lines == "//"))
parse_rxn_ec <- function(block){
  block <- block[block != "//"]
  rxn_id <- NA_character_
  ecs <- character(0)
  
  for(line in block){
    if(startsWith(line, "UNIQUE-ID - "))
      rxn_id <- sub("^UNIQUE-ID - ", "", line)
    
    if(startsWith(line, "EC-NUMBER - "))
      ecs <- c(ecs, sub("^EC-NUMBER - ", "", line))
  }
  
  if(is.na(rxn_id)) return(NULL)
  
  tibble(
    reaction_id = rxn_id,
    EC_number = if(length(ecs)==0) NA_character_ else unique(ecs)
  )
}

rxn_ec_map <- bind_rows(lapply(rxn_blocks, parse_rxn_ec)) %>%
  unnest(EC_number) %>%
  filter(!is.na(EC_number) & nzchar(EC_number)) %>%
  distinct()

head(rxn_ec_map)

# 4. Join reactin_id and ec_number into edges_inchikey
library(dplyr)
edges_annot <- edges_inchikey %>%
  left_join(enzrxns_all, by="reaction_id")

    #check
      colnames(edges_annot)
      # [1] "reaction_id"        "substrate"          "product"            "substrate_inchikey" "product_inchikey"  
      # [6] "ec_number"
      table(is.na(edges_annot$substrate_inchikey))
      # FALSE  TRUE 
      # 80182 84562 
      table(is.na(edges_annot$product_inchikey))
      # FALSE  TRUE 
      # 90685 74059 
      table(is.na(edges_annot$product_inchikey) & is.na(edges_annot$substrate_inchikey))
      # FALSE   TRUE 
      # 118831  45913 
# 5. Remove duplicate rows
edges_annot <- edges_annot %>%
      distinct()
nrow(edges_annot) #164534
save(edges_annot, file= "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/edges_annot")
head(edges_annot)
    # reaction_id                substrate   product substrate_inchikey product_inchikey     ec_number
    # 1   RXN-16092                CPD-17341 CPD-17342               <NA>             <NA> MONOMER-19055
    # 2   RXN-16092                 Donor-H2 CPD-17342               <NA>             <NA> MONOMER-19055
    # 3   RXN-16092                CPD-17341  Acceptor               <NA>             <NA> MONOMER-19055
    # 4   RXN-16092                 Donor-H2  Acceptor               <NA>             <NA> MONOMER-19055
    # 5   RXN-23105 3-hydroxy-pentanoyl-SpnA  SpnA-PKS               <NA>             <NA>          <NA>
    #   6   RXN-23105                 SpnB-PKS  SpnA-PKS               <NA>             <NA>          <NA>

#6. error in ec_number format, parse reactions.dat to get the correct ec_number
library(dplyr)
library(tidyr)

enzymes_file <- "~/compound_database_shenlab/metacyc/data/data/reactions.dat"
rxn_lines <- readLines(enzymes_file, encoding="latin1", warn=FALSE)
rxn_lines <- rxn_lines[!grepl("^#", rxn_lines)]
rxn_blocks <- split(rxn_lines, cumsum(rxn_lines == "//"))

parse_rxn_ec <- function(block){
  block <- block[block != "//"]
  rxn_id <- NA_character_
  ecs <- character(0)
  
  for(line in block){
    if(startsWith(line, "UNIQUE-ID - "))
      rxn_id <- sub("^UNIQUE-ID - ", "", line)
    
    if(startsWith(line, "EC-NUMBER - "))
      ecs <- c(ecs, sub("^EC-NUMBER - ", "", line))
  }
  
  if(is.na(rxn_id)) return(NULL)
  
  tibble(
    reaction_id = rxn_id,
    EC_number = if(length(ecs)==0) NA_character_ else unique(ecs)
  )
}

rxn_ec_map <- bind_rows(lapply(rxn_blocks, parse_rxn_ec)) %>%
  unnest(EC_number) %>%
  filter(!is.na(EC_number) & nzchar(EC_number)) %>%
  distinct()

head(rxn_ec_map)
# A tibble: 6 × 2
    # reaction_id                     EC_number    
    # <chr>                           <chr>        
    #   1 RXN-16092                       EC-1.1.99    
    # 2 RXN-20782                       EC-2.1.1.M65 
    # 3 ACETYL-COA-CARBOXYLTRANSFER-RXN EC-6.4.1.2   
    # 4 RXN-22227                       |EC-7.5.99.a|
    #   5 RXN-23229                       EC-1.1.1.81  
    # 6 RXN-13301                       EC-1.1.1.330 

#7. Correct the ec_number in edges_annot
  #check if the column types for reaction_id are the same
rxn_ec_map$reaction_id <- as.character(rxn_ec_map$reaction_id)
edges_annot$reaction_id <- as.character(edges_annot)
  #join the two data frames using reaction_id as the common column
edges_annot_updated <- edges_annot %>%
  left_join(rxn_ec_map, 
            by = "reaction_id") 
  #check and drop unnecessary labels in each box
    head(edges_annot_updated)
    # reaction_id                substrate   product substrate_inchikey product_inchikey     ec_number EC_number
    # 1   RXN-16092                CPD-17341 CPD-17342               <NA>             <NA> MONOMER-19055 EC-1.1.99
    # 2   RXN-16092                 Donor-H2 CPD-17342               <NA>             <NA> MONOMER-19055 EC-1.1.99
    # 3   RXN-16092                CPD-17341  Acceptor               <NA>             <NA> MONOMER-19055 EC-1.1.99
    # 4   RXN-16092                 Donor-H2  Acceptor               <NA>             <NA> MONOMER-19055 EC-1.1.99
    # 5   RXN-23105 3-hydroxy-pentanoyl-SpnA  SpnA-PKS               <NA>             <NA>          <NA>      <NA>
    #   6   RXN-23105                 SpnB-PKS  SpnA-PKS               <NA>             <NA>          <NA>      <NA>

    #drop the unnecessary EC- in EC_number
    edges_annot_updated <- edges_annot_updated %>%
      mutate(ec_id = ec_number) %>%
      mutate(EC_number = sub("^EC-", "", EC_number)) %>%
      select(-ec_number)

head(edges_annot_updated) #167632 rows
    # reaction_id                substrate   product substrate_inchikey product_inchikey EC_number         ec_id
    # 1   RXN-16092                CPD-17341 CPD-17342               <NA>             <NA>    1.1.99 MONOMER-19055
    # 2   RXN-16092                 Donor-H2 CPD-17342               <NA>             <NA>    1.1.99 MONOMER-19055
    # 3   RXN-16092                CPD-17341  Acceptor               <NA>             <NA>    1.1.99 MONOMER-19055
    # 4   RXN-16092                 Donor-H2  Acceptor               <NA>             <NA>    1.1.99 MONOMER-19055
    # 5   RXN-23105 3-hydroxy-pentanoyl-SpnA  SpnA-PKS               <NA>             <NA>      <NA>          <NA>
    #   6   RXN-23105                 SpnB-PKS  SpnA-PKS               <NA>             <NA>      <NA>          <NA>

save(edges_annot_updated, file= "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/edges_annot_updated")
