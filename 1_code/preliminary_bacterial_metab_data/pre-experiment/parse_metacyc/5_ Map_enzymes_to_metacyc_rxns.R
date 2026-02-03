#1. Load datasets
#collapsed_results: results of differential features and metabolites
#enzyme_catalogue: enzyme data derived from bacterial genome
#edges_annot_updated: edges of reaction pairs derived from metacyc
load("2_data/preliminary_determination_bacterial_metabolomics/20251010/MS2/differential_analysis_each_bacterium/collapsed_results")
load("2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/enzyme_catalogue")
load("2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/edges_annot_updated")

#2. Load differential analysis data
#collapsed_results does not have the kegg id or inchikey in the annotation for further identifying of pathways and reaction_ids
#add in this data using the file "annotation_lookup"
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/MS2/differential_analysis_each_bacterium/annotation_lookup")
colnames(collapsed_results$`Enterococcus faecium`) #"variable_id"   "class"         "fc_final"      "pval_final"    "Compound.name"
colnames(annotation_lookup) #[1] "variable_id" "KEGG.ID"     "HMDB.ID"     "CAS.ID"      "INCHIKEY_ID"
library(dplyr)
Enterococcus_faecium_withid <- collapsed_results$`Enterococcus faecium` %>%
  left_join(
    annotation_lookup %>%
      select(variable_id, INCHIKEY_ID),
    by = "variable_id"
  )
colnames(Enterococcus_faecium_withid) 
#[1] "variable_id"   "class"         "fc_final"      "pval_final"    "Compound.name" "INCHIKEY_ID"  


Enterococcus_faecium_withid_keggid <- collapsed_results$`Enterococcus faecium` %>%
  left_join(
    annotation_lookup %>%
      select(variable_id, KEGG.ID),
    by = "variable_id"
  )
colnames(Enterococcus_faecium_withid_keggid) 
#[1] "variable_id"   "class"         "fc_final"      "pval_final"    "Compound.name" "KEGG.ID"

      #this section is re-organised into 6_reaction_edges_inchikey
      # #3. According to variable_id, put in ec_number and kegg
      # load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/edges_annot_updated")
      # colnames(edges_annot_updated) 
      # #[1] "reaction_id"  "substrate"  "product" "substrate_inchikey" "product_inchikey"   "EC_number". "ec_id     
      # reaction_pairs_hits <- reaction_pairs_ef %>%
      #   filter(!is.na(substrate) | !is.na(product))
      # 
      # nrow(reaction_pairs_hits)
      # head(reaction_pairs_hits, 20)
      # reaction_pairs_both <- reaction_pairs_ef %>%
      #   filter(!is.na(substrate) & !is.na(product))
      # 
      # nrow(reaction_pairs_both)
      # head(reaction_pairs_both, 20)
      # 
      #     #check the table as an overview
      #     overview_pairs <- reaction_pairs_ef %>%
      #       summarise(
      #         total_pairs = n(),
      #         substrate_matched = sum(!is.na(substrate)),
      #         product_matched   = sum(!is.na(product)),
      #         both_matched      = sum(!is.na(substrate) & !is.na(product)),
      #         neither_matched   = sum(is.na(substrate) & is.na(product))
      #       )
      #     overview_pairs
      #     
      #   #merge substrate and products  
      #   overview_rxn <- reaction_pairs_ef %>%
      #       mutate(
      #         sub_hit = !is.na(substrate),
      #         prod_hit = !is.na(product),
      #         both_hit = sub_hit & prod_hit
      #       ) %>%
      #       group_by(reaction_id) %>%
      #       summarise(
      #         any_substrate_hit = any(sub_hit),
      #         any_product_hit   = any(prod_hit),
      #         any_both_hit      = any(both_hit),
      #         .groups = "drop"
      #       ) %>%
      #       summarise(
      #         total_reactions = n(),
      #         reactions_with_any_substrate = sum(any_substrate_hit),
      #         reactions_with_any_product   = sum(any_product_hit),
      #         reactions_with_both_in_same_pair = sum(any_both_hit)
      #       )
      #     
      #   overview_rxn
      #     #check
      #     reaction_pairs_ef %>%
      #       arrange(desc(!is.na(substrate) | !is.na(product))) %>%
      #       head(20)
      #       # reaction_id          substrate_inchikey        substrate            product_inchikey     product  EC_number
      #       # 1          NITRATE-REDUCTASE-NADPH-RXN IOVCWXUNBOPUCH-UHFFFAOYSA-M             <NA> NHNBFGGVMKEFGY-UHFFFAOYSA-N     Nitrate    1.7.1.3
      #       # 2          NITRATE-REDUCTASE-NADPH-RXN XLYOFNOQVPJJNP-UHFFFAOYSA-N             <NA> NHNBFGGVMKEFGY-UHFFFAOYSA-N     Nitrate    1.7.1.3
      #       # 3          NITRATE-REDUCTASE-NADPH-RXN                        <NA>             <NA> NHNBFGGVMKEFGY-UHFFFAOYSA-N     Nitrate    1.7.1.3
      #       # 4                            RXN-19913                        <NA>             <NA> WPYMKLBDIGXBTP-UHFFFAOYSA-M    benzoate   1.2.1.30
      #       # 5                            RXN-19913 XPPKVPWEQAFLFU-UHFFFAOYSA-K             <NA> WPYMKLBDIGXBTP-UHFFFAOYSA-M    benzoate   1.2.1.30
      #       # 6                            RXN-19913 HUMNYLRZRPPJDN-UHFFFAOYSA-N             <NA> WPYMKLBDIGXBTP-UHFFFAOYSA-M    benzoate   1.2.1.30
      #       # 7  BENZALDEHYDE-DEHYDROGENASE-NAD+-RXN HUMNYLRZRPPJDN-UHFFFAOYSA-N             <NA> WPYMKLBDIGXBTP-UHFFFAOYSA-M    benzoate   1.2.1.28
      #       # 8  BENZALDEHYDE-DEHYDROGENASE-NAD+-RXN XLYOFNOQVPJJNP-UHFFFAOYSA-N             <NA> WPYMKLBDIGXBTP-UHFFFAOYSA-M    benzoate   1.2.1.28
      #       # 9  BENZALDEHYDE-DEHYDROGENASE-NAD+-RXN                        <NA>             <NA> WPYMKLBDIGXBTP-UHFFFAOYSA-M    benzoate   1.2.1.28
      #       # 10                            RXN-4444 XFRVVPUIAFSTFO-UHFFFAOYSA-N       Tridecanol BGEHHAVMRVXCGR-UHFFFAOYSA-N        <NA>   1.1.3.20
      #       # 11                            RXN-4444 XFRVVPUIAFSTFO-UHFFFAOYSA-N       Tridecanol MHAJPDPJQMAIIY-UHFFFAOYSA-N        <NA>   1.1.3.20
      #       # 12                           RXN-24428 COJBGNAUUSNXHX-UHFFFAOYSA-M             <NA> ZDXPYRJPNDTMRX-VKHMYHEASA-N L-Glutamine      2.6.1
      #       # 13                           RXN-24428                        <NA>             <NA> ZDXPYRJPNDTMRX-VKHMYHEASA-N L-Glutamine      2.6.1
      #       # 14                            R247-RXN NHNBFGGVMKEFGY-UHFFFAOYSA-N          Nitrate                        <NA>        <NA>      1.7.2
      #       # 15                            R247-RXN NHNBFGGVMKEFGY-UHFFFAOYSA-N          Nitrate IOVCWXUNBOPUCH-UHFFFAOYSA-M        <NA>      1.7.2
      #       # 16                            R247-RXN NHNBFGGVMKEFGY-UHFFFAOYSA-N          Nitrate XLYOFNOQVPJJNP-UHFFFAOYSA-N        <NA>      1.7.2
      #       # 17                           RXN-10118 DXOSJQLIRGXWCF-UHFFFAOYSA-N 3-Fluorocatechol PTYIDKBTGREHJW-HSFFGMMNSA-L        <NA> 1.13.11.M6
      #       # 18                           RXN-10118 DXOSJQLIRGXWCF-UHFFFAOYSA-N 3-Fluorocatechol GPRLSGONYQIRFK-UHFFFAOYSA-N        <NA> 1.13.11.M6
      #       # 19                           RXN-23849                        <NA>             <NA> ZDXPYRJPNDTMRX-VKHMYHEASA-N L-Glutamine      1.5.1
      #       # 20                           RXN-23849 XLYOFNOQVPJJNP-UHFFFAOYSA-N             <NA> ZDXPYRJPNDTMRX-VKHMYHEASA-N L-Glutamine      1.5.1
      #     reaction_pairs_ef %>%
      #       arrange(desc(!is.na(substrate) & !is.na(product))) %>%
      #       head(20)
      #     # Enterococcus_faecium_ec_filtered <- reaction_pairs_ef %>%
      #     #   filter(!is.na(substrate) & !is.na(product))
      #     #The above command gave redundant rows
      #     #The next step is to remove redundant rows 
      #     Enterococcus_faecium_ec_filtered <- reaction_pairs_ef %>%
      #       filter(substrate_inchikey != product_inchikey)
      #     nrow(Enterococcus_faecium_ec_filtered) #39780
      # #Match genome-derived ec_numbers with filtered metabolite-derived ec_number
      #     #check structure of enzyme_catalogue
      #     colnames(enzyme_catalogue)
      #     #protein_id"     "EC_number"      "gene_id"        "product"        "sources"        "n_sources"      "confidence"   "kegg_reactions" "kegg_supported"
      # 
      #     #check structure of Enterococcus_faecium_ec_filtered
      #     colnames(Enterococcus_faecium_ec_filtered)
      #     #"reaction_id"        "substrate_inchikey" "substrate"          "product_inchikey"   "product"            "EC_number"    
      #     
      #   #1. clean the EC_number columns for both tables
      #      met_ec <- Enterococcus_faecium_ec_filtered %>%
      #       mutate(
      #         EC_number = str_trim(as.character(EC_number))
      #       )
      #     
      #     gen_ec <- enzyme_catalogue %>%
      #       mutate(
      #         EC_number = str_trim(as.character(EC_number))
      #       )
      #   #2. met_gen_merged_long <- met_ec_nonNA %>%
      #   met_gen_merged_long <- met_ec %>%
      #       left_join(
      #         gen_ec %>%
      #           transmute(
      #             EC_number,
      #             protein_id,
      #             gene_id,
      #             gene_confidence = confidence,
      #             enzyme_product = product   # rename here
      #           ),
      #         by = "EC_number"
      #       )
      #  
      #     nrow(met_ec_merged_long) #40478
      #   #3. Filter out NA rows for subtrate and product columns
      #     met_gen_merged_long_filtered <- met_gen_merged_long %>%
      #       filter(!is.na(substrate) & !is.na(product))
      #     #the columns substrate and product are empty again, add them in again using edges_annot_updated
      #     inchikey_name_map <- bind_rows(
      #       edges_annot_updated %>%
      #         transmute(inchikey = substrate_inchikey, metacyc_name = substrate),
      #       edges_annot_updated %>%
      #         transmute(inchikey = product_inchikey,   metacyc_name = product)
      #     ) %>%
      #       filter(!is.na(inchikey), !is.na(metacyc_name)) %>%
      #       distinct(inchikey, .keep_all = TRUE)
      #     met_gen_merged_filled <- met_gen_merged_long %>%
      #       left_join(
      #         inchikey_name_map %>% rename(substrate_inchikey = inchikey,
      #                                      substrate_metacyc = metacyc_name),
      #         by = "substrate_inchikey"
      #       ) %>%
      #       left_join(
      #         inchikey_name_map %>% rename(product_inchikey = inchikey,
      #                                      product_metacyc = metacyc_name),
      #         by = "product_inchikey"
      #       ) %>%
      #       mutate(
      #         substrate = coalesce(substrate, substrate_metacyc),
      #         product   = coalesce(product,   product_metacyc)
      #       ) %>%
      #       select(-substrate_metacyc, -product_metacyc)
      #     met_gen_merged_filled %>%
      #       summarise(
      #         n = n(),
      #         substrate_filled = sum(!is.na(substrate)),
      #         product_filled   = sum(!is.na(product))
      #       )
      #     
      #     head(met_gen_merged_filled)
      #     #       n substrate_filled product_filled
      #     # 1 40478            40478          40478
      #     
      # #Remove redundant metabolites (called currency metabolites)
      #     currency_metabolites <- c(
      #       "PROTON", "H+", "WATER", "H2O",
      #       "CARBON-DIOXIDE", "CO2",
      #       "OXYGEN", "O2",
      #       "NAD", "NADH", "NADP", "NADPH",
      #       "ATP", "ADP", "AMP",
      #       "PHOSPHATE", "ORTHOPHOSPHATE",
      #       "HCO3", "BICARBONATE"
      #     )
      #     met_gen_filtered <- met_gen_merged_filled %>%
      #       filter(
      #         !substrate %in% currency_metabolites,
      #         !product   %in% currency_metabolites
      #       )
      # # Add in annotation_confidence
      #     #1. Define the function for annotation confidence
      #     library(stringr)
      #     add_annotation_confidence <- function(df) {
      #       df %>%
      #         mutate(
      #           gene_conf_score = case_when(
      #             is.na(gene_confidence) ~ 0L,
      #             str_to_lower(str_trim(gene_confidence)) == "low" ~ 1L,
      #             str_to_lower(str_trim(gene_confidence)) == "medium" ~2L,
      #             str_to_lower(str_trim(gene_confidence)) == "medium-high" ~ 3L,
      #             str_to_lower(str_trim(gene_confidence)) == "high" ~ 4L,
      #             TRUE ~0L
      #           ),
      #           subprod_score = case_when(
      #             !is.na(substrate) & !is.na(product) ~ 2L,
      #             xor(!is.na(substrate), !is.na(product)) ~1L,
      #             TRUE ~ 0L
      #           ),
      #           ec_score = if_else(!is.na(EC_number) & str_trim(as.character(EC_number)) != "", 1L, 0L),
      #           protein_score = if_else(!is.na(protein_id) & str_trim(as.character(protein_id)) != "", 1L, 0L),
      #           gene_id_score = if_else(!is.na(gene_id) & str_trim(as.character(gene_id)) != "", 1L, 0L),
      #           annotation_confidence = gene_conf_score + subprod_score + ec_score + protein_score + gene_id_score
      #         ) %>%
      #         arrange(desc(annotation_confidence))
      #     }
      #     #2. Apply function to the object
      #     met_gen_scored <- add_annotation_confidence(met_gen_filtered)
      #     #3. Remove the columns gene_conf_score, subprod_score, ec_score, protein_score, gene_id_score
      #     met_gen_scored <- met_gen_scored %>%
      #       select(
      #         -gene_conf_score,
      #         -subprod_score,
      #         -ec_score,
      #         -protein_score,
      #         -gene_id_score
      #       )
      #     nrow(met_gen_scored) #16867
      #     #4. Remove redundant metabolites again
      #     bad_patterns <- c(
      #       "Primary-|Secondary-|Monocarboxylates|Alcohols",
      #       "Orthophosphoric",
      #       "Glycerophosphodiesters",
      #       "Uracil-tRNAs|tRNA",
      #       "Protein|Peptidoglycan|Ubiquitin|carrier",
      #       "THF-GLU|Formyl-THF",
      #       "^CPD-(110|6442)"  # known non-metabolite CPDs
      #     )
      #     
      #     met_gen_clean <- met_gen_scored %>%
      #       filter(
      #         !str_detect(substrate, paste(bad_patterns, collapse = "|")),
      #         !str_detect(product,   paste(bad_patterns, collapse = "|"))
      #       )
      #     nrow(met_gen_clean) #14478
      #     nrow(distinct(met_gen_clean)) #14478
      #     
      #     ef_inchikeys <- Enterococcus_faecium_withid %>%
      #       transmute(INCHIKEY_ID = toupper(trimws(INCHIKEY_ID))) %>%
      #       filter(!is.na(INCHIKEY_ID)) %>%
      #       distinct()
      #     
      #     met_gen_clean_tagged <- met_gen_clean %>%
      #       mutate(
      #         substrate_hit = toupper(trimws(substrate_inchikey)) %in% ef_inchikeys$INCHIKEY_ID,
      #         product_hit   = toupper(trimws(product_inchikey))   %in% ef_inchikeys$INCHIKEY_ID,
      #         hit_type = case_when(
      #           substrate_hit & product_hit ~ "both",
      #           substrate_hit & !product_hit ~ "substrate_only",
      #           !substrate_hit & product_hit ~ "product_only",
      #           TRUE ~ "neither"
      #         )
      #       )
      #     
      #     met_gen_clean_tagged %>% count(hit_type)
      #     
      #     
      #     met_gen_clean_tagged %>% count(hit_type)
      #     
      # #Add in pvalue and fc values
      #     #1. Select necessary information from Enterococcus_faecium_ec
      #       #1A. Start with joining substrate info first
      #     ef_met_stats <- Enterococcus_faecium_ec %>%
      #       transmute(
      #         inchikey = toupper(trimws(INCHIKEY_ID)),
      #         substrate_fc_final = fc_final,
      #         substrate_pval_final = pval_final
      #       ) %>%
      #       filter(!is.na(inchikey)) %>%
      #       group_by(inchikey) %>%
      #       arrange(substrate_pval_final, desc(abs(substrate_fc_final))) %>%
      #       slice(1) %>%   # keep the most significant if duplicates exist
      #       ungroup()
      #         #Join stats value to reaction pair info
      #     met_gen_with_met <- met_gen_clean %>%
      #       mutate(
      #         substrate_inchikey = toupper(trimws(substrate_inchikey)),
      #         product_inchikey   = toupper(trimws(product_inchikey))
      #       ) %>%
      #       left_join(
      #         ef_met_stats,
      #         by = c("substrate_inchikey" = "inchikey")
      #       )
      #     #1B. Join for product
      #     met_gen_with_met <- met_gen_with_met %>%
      #       left_join(
      #         ef_met_stats %>%
      #           rename(
      #             product_inchikey = inchikey,
      #             product_fc_final = substrate_fc_final,
      #             product_pval_final = substrate_pval_final,
      #           ),
      #         by = "product_inchikey"
      #       )
      #     #2. Check successful pairs
      #     sum(!is.na(met_gen_with_met$substrate_fc_final) & !is.na(met_gen_with_met$substrate_pval_final)) #166
      #     sum(!is.na(met_gen_with_met$product_fc_final) & !is.na(met_gen_with_met$product_pval_final)) #97
      #     sum(!is.na(met_gen_with_met$substrate_fc_final) & !is.na(met_gen_with_met$substrate_pval_final) & 
      #           !is.na(met_gen_with_met$product_fc_final) & !is.na(met_gen_with_met$product_pval_final)) #1
      #     met_gen_with_met %>%
      #       filter(
      #         !is.na(substrate_fc_final),
      #         !is.na(substrate_pval_final),
      #         !is.na(product_fc_final),
      #         !is.na(product_pval_final)
      #       )
    