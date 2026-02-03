# Use KEGGREST for biochemical reaction annotation for reaction-level annotation
library(KEGGREST)
library(dplyr)
library(tidyr)

#1. Prepare EC numbers
gb_only_ec_tidy <- gb_only_ec %>%
  tidyr::separate_rows(EC_number, sep = ",") %>%   # split multiple ECs
  mutate(EC_number = trimws(EC_number)) %>%       # remove white space(blanks)
  filter(!EC_number %in% c("", "NA", "2.3.-.-"))
#2. Map EC number to kegg
ec_to_kegg <- gb_only_ec_tidy %>%
  rowwise() %>%
  mutate(
    kegg_reactions = list(tryCatch(keggLink("rn", paste0("ec:", EC_number)), 
                                   error = function(e) NA))
  ) %>%
  ungroup()
  #Check
  head(ec_to_kegg)
#3. Map reactions to KEGG pathways
ec_to_kegg <- ec_to_kegg %>%
    rowwise() %>%
    mutate(
      kegg_pathways = list(
        if (length(kegg_reactions) > 0 && !all(is.na(kegg_reactions))) {
          unlist(lapply(kegg_reactions, function(r) {
            tryCatch(keggLink("pathway", r), error = function(e) NA)
          }))
        } else {
          NA_character_
        }
      )
    ) %>%
    ungroup()
#4. Clean it such that one row per EC per pathway/reaction
# ec_kegg_flat <- ec_to_kegg %>%
#     tidyr::unnest(cols = c(kegg_reactions, kegg_pathways)) %>%
#     distinct()
ec_reactions_flat <- ec_to_kegg %>%
  select(protein_id, EC_number, product, kegg_reactions) %>%
  tidyr::unnest(cols = c(kegg_reactions)) %>%
  distinct()
ec_pathways_flat <- ec_reactions_flat %>%
  rowwise() %>%
  mutate(kegg_pathways = list(
    if (!is.na(kegg_reactions) && kegg_reactions != "") {
      tryCatch(keggLink("pathway", kegg_reactions), error = function(e) NA)
    } else {
      NA_character_
    }
  )) %>%
  ungroup() %>%
  tidyr::unnest(cols = c(kegg_pathways)) %>%
  distinct()
  #check
  head(ec_reactions_flat)
  # A tibble: 6 × 4
  # protein_id     EC_number product                                                                 kegg_reactions
  # <chr>          <chr>     <chr>                                                                   <chr>         
  #   1 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synthetase GshB rn:R10993     
  # 2 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synthetase GshB rn:R00894     
  # 3 WP_002287901.1 6.3.2.3   bifunctional glutamate--cysteine ligaseGshA/glutathione synthetase GshB rn:R10994     
  # 4 WP_002287901.1 6.3.2.3   bifunctional glutamate--cysteine ligaseGshA/glutathione synthetase GshB rn:R00497     
  # 5 WP_002297133.1 2.7.1.201 PTS system trehalose-specific EIIBC component                           rn:R02780     
  # 6 WP_002297133.1 2.7.1.201 PTS system trehalose-specific EIIBC component                           rn:R06229   
  colnames(ec_reactions_flat)
  #[1] "protein_id"     "EC_number"      "product"        "kegg_reactions"
  
  head(ec_pathways_flat)
  # A tibble: 6 × 5
  # protein_id     EC_number product                                                        kegg_reactions kegg_pathways
  # <chr>          <chr>     <chr>                                                          <chr>          <chr>        
  #   1 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synth… rn:R10993      path:map00270
  # 2 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synth… rn:R10993      path:rn00270 
  # 3 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synth… rn:R10993      path:map01100
  # 4 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synth… rn:R10993      path:rn01100 
  # 5 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synth… rn:R00894      path:map00480
  # 6 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA/glutathione synth… rn:R00894      path:rn00480 
  colnames(ec_pathways_flat) 
  #[1] "protein_id"     "EC_number"      "product"        "kegg_reactions" "kegg_pathways" 
  save(ec_pathways_flat, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/ec_pathways_flat")
  save(ec_reactions_flat, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/ec_reactions_flat")
  