############################################################
# Enterococcus faecium metabolomics ↔ MetaCyc ↔ genome EC
# Versioned pipeline:
#  - Redundant metabolite removal (currency/generic/macro) is versioned
#  - Annotation confidence scoring is versioned
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
})

############################################################
# 0) Load datasets
############################################################
load("2_data/preliminary_determination_bacterial_metabolomics/20251010/MS2/differential_analysis_each_bacterium/collapsed_results")
load("2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/enzyme_catalogue")
load("2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/edges_annot_updated")
load("2_data/preliminary_determination_bacterial_metabolomics/20251010/MS2/differential_analysis_each_bacterium/annotation_lookup")


#collapsed_results: results of differential features and metabolites
#enzyme_catalogue: enzyme data derived from bacterial genome
#edges_annot_updated: edges of reaction pairs derived from metacyc

############################################################
# 1) Build Enterococcus faecium metabolite table with InChIKeys
############################################################
Enterococcus_faecium_withid <- collapsed_results$`Enterococcus faecium` %>%
  left_join(
    annotation_lookup %>% select(variable_id, INCHIKEY_ID),
    by = "variable_id"
  ) %>%
  mutate(
    INCHIKEY_ID = toupper(str_trim(INCHIKEY_ID))
  )

# Map InChIKey -> metabolite stats (fc/p) (deduplicate by best p, then biggest |fc|)
ef_met_stats <- Enterococcus_faecium_withid %>%
  transmute(
    inchikey = INCHIKEY_ID,
    fc_final = fc_final,
    pval_final = pval_final,
    Compound.name = Compound.name
  ) %>%
  filter(!is.na(inchikey), inchikey != "") %>%
  group_by(inchikey) %>%
  arrange(pval_final, desc(abs(fc_final))) %>%
  slice(1) %>%
  ungroup()

ef_keys <- ef_met_stats %>% distinct(inchikey)

############################################################
# 2) Build MetaCyc edge table (reaction pairs) and keep those
#    where either substrate or product InChIKey hits EF metabolites
############################################################
edges_clean <- edges_annot_updated %>%
  transmute(
    reaction_id,
    substrate_inchikey = toupper(str_trim(substrate_inchikey)),
    product_inchikey   = toupper(str_trim(product_inchikey)),
    substrate = substrate,
    product   = product,
    EC_number = str_trim(as.character(EC_number))
  )

# Keep only edges where at least one end is in EF metabolites
reaction_pairs_ef <- edges_clean %>%
  filter(
    (!is.na(substrate_inchikey) & substrate_inchikey %in% ef_keys$inchikey) |
      (!is.na(product_inchikey) & product_inchikey %in% ef_keys$inchikey)
  ) %>%
  # drop self-loops (same compound on both ends)
  filter(is.na(substrate_inchikey) | is.na(product_inchikey) | substrate_inchikey != product_inchikey) %>%
  distinct(reaction_id, substrate_inchikey, product_inchikey, .keep_all = TRUE)

# Fill names from MetaCyc mapping by InChIKey (robust even if substrate/product columns are missing earlier)
inchikey_name_map <- bind_rows(
  edges_clean %>% transmute(inchikey = substrate_inchikey, metacyc_name = substrate),
  edges_clean %>% transmute(inchikey = product_inchikey,   metacyc_name = product)
) %>%
  filter(!is.na(inchikey), inchikey != "", !is.na(metacyc_name), metacyc_name != "") %>%
  distinct(inchikey, .keep_all = TRUE)

reaction_pairs_ef <- reaction_pairs_ef %>%
  left_join(
    inchikey_name_map %>% rename(substrate_inchikey = inchikey, substrate_metacyc = metacyc_name),
    by = "substrate_inchikey"
  ) %>%
  left_join(
    inchikey_name_map %>% rename(product_inchikey = inchikey, product_metacyc = metacyc_name),
    by = "product_inchikey"
  ) %>%
  mutate(
    substrate = coalesce(substrate, substrate_metacyc),
    product   = coalesce(product,   product_metacyc)
  ) %>%
  select(-substrate_metacyc, -product_metacyc)

############################################################
# 3) Join genome enzyme catalogue (by EC_number)
############################################################
gen_ec <- enzyme_catalogue %>%
  mutate(EC_number = str_trim(as.character(EC_number))) %>%
  transmute(
    EC_number,
    protein_id,
    gene_id,
    gene_confidence = confidence,
    enzyme_product  = product
  )

met_gen_merged <- reaction_pairs_ef %>%
  left_join(gen_ec, by = "EC_number")

############################################################
# 4) Attach EF metabolomics stats to substrate/product (by InChIKey)
############################################################
met_gen_with_met <- met_gen_merged %>%
  left_join(
    ef_met_stats %>%
      transmute(
        substrate_inchikey = inchikey,
        substrate_fc_final = fc_final,
        substrate_pval_final = pval_final,
        substrate_compound_name = Compound.name
      ),
    by = "substrate_inchikey"
  ) %>%
  left_join(
    ef_met_stats %>%
      transmute(
        product_inchikey = inchikey,
        product_fc_final = fc_final,
        product_pval_final = pval_final,
        product_compound_name = Compound.name
      ),
    by = "product_inchikey"
  )

############################################################
# 5) Versioned "redundant metabolite" handling (v1)
#    - Prefer FLAGS first (for audit), then filter to "clean_v1"
############################################################
currency_metabolites_v1 <- c(
  "PROTON", "H+", "WATER", "H2O",
  "CARBON-DIOXIDE", "CO2",
  "OXYGEN", "O2",
  "NAD", "NADH", "NADP", "NADPH",
  "ATP", "ADP", "AMP",
  "PHOSPHATE", "ORTHOPHOSPHATE",
  "HCO3", "BICARBONATE", "AMMONIUM"
)

generic_class_patterns_v1 <- c(
  "Primary-|Secondary-|Monocarboxylates|Alcohols",
  "Orthophosphoric",
  "Glycerophosphodiesters",
  "Quinones|Flavanones|Monocarboxylates"
)

macromolecule_patterns_v1 <- c(
  "Uracil-tRNAs|tRNA",
  "Protein|Peptidoglycan|Ubiquitin|carrier"
)

met_gen_flagged_v1 <- met_gen_with_met %>%
  mutate(
    substrate_u = toupper(str_trim(substrate)),
    product_u   = toupper(str_trim(product)),
    
    is_currency_v1 =
      substrate_u %in% currency_metabolites_v1 |
      product_u   %in% currency_metabolites_v1,
    
    is_generic_class_v1 =
      str_detect(substrate, paste(generic_class_patterns_v1, collapse = "|")) |
      str_detect(product,   paste(generic_class_patterns_v1, collapse = "|")),
    
    is_macromolecule_v1 =
      str_detect(substrate, paste(macromolecule_patterns_v1, collapse = "|")) |
      str_detect(product,   paste(macromolecule_patterns_v1, collapse = "|")),
    
    is_metacyc_placeholder_v1 =
      str_detect(substrate, "^CPD-") | str_detect(product, "^CPD-")
  ) %>%
  select(-substrate_u, -product_u)

# "clean" table for v1 (remove currency + generic + macro)
met_gen_clean_v1 <- met_gen_flagged_v1 %>%
  filter(
    !is_currency_v1,
    !is_generic_class_v1,
    !is_macromolecule_v1
  )

############################################################
# 6) Versioned annotation confidence scoring (v1)
#    Note: your stated components sum to max 9; we keep that.
#          If you want max 8, drop gene_id_score OR protein_score.
############################################################
score_annotation_v1 <- function(df) {
  df %>%
    mutate(
      gene_conf_score_v1 = case_when(
        is.na(gene_confidence) ~ 0L,
        str_to_lower(str_trim(gene_confidence)) == "low" ~ 1L,
        str_to_lower(str_trim(gene_confidence)) == "medium" ~ 2L,
        str_to_lower(str_trim(gene_confidence)) %in% c("medium-high", "medium high", "medium_high") ~ 3L,
        str_to_lower(str_trim(gene_confidence)) == "high" ~ 4L,
        TRUE ~ 0L
      ),
      
      subprod_score_v1 = case_when(
        !is.na(substrate) & !is.na(product) ~ 2L,
        xor(!is.na(substrate), !is.na(product)) ~ 1L,
        TRUE ~ 0L
      ),
      
      ec_score_v1 = if_else(!is.na(EC_number) & str_trim(as.character(EC_number)) != "", 1L, 0L),
      protein_score_v1 = if_else(!is.na(protein_id) & str_trim(as.character(protein_id)) != "", 1L, 0L),
      gene_id_score_v1 = if_else(!is.na(gene_id) & str_trim(as.character(gene_id)) != "", 1L, 0L),
      
      annotation_confidence_v1 =
        gene_conf_score_v1 + subprod_score_v1 + ec_score_v1 + protein_score_v1 + gene_id_score_v1
    ) %>%
    arrange(desc(annotation_confidence_v1))
}

met_gen_scored_v1 <- score_annotation_v1(met_gen_clean_v1)

############################################################
# 7) Optional: de-duplicate rows (choose your rule)
#    A) strict row-level distinct:
############################################################
met_gen_scored_v1_distinct <- met_gen_scored_v1 %>%
  distinct(reaction_id, substrate_inchikey, product_inchikey, EC_number, protein_id, gene_id, .keep_all = TRUE)

############################################################
# 8) Diagnostics / sanity checks
############################################################
# How many edges had substrate/product that are EF metabolites?
met_gen_scored_v1 %>%
  summarise(
    n_rows = n(),
    substrate_is_EF = sum(!is.na(substrate_inchikey) & substrate_inchikey %in% ef_keys$inchikey),
    product_is_EF   = sum(!is.na(product_inchikey) & product_inchikey %in% ef_keys$inchikey),
    both_are_EF     = sum(!is.na(substrate_inchikey) & substrate_inchikey %in% ef_keys$inchikey &
                            !is.na(product_inchikey) & product_inchikey %in% ef_keys$inchikey),
    max_conf_v1 = max(annotation_confidence_v1, na.rm = TRUE)
  )

# Count hit types (for interpretation)
met_gen_scored_v1 %>%
  mutate(
    substrate_hit = !is.na(substrate_inchikey) & substrate_inchikey %in% ef_keys$inchikey,
    product_hit   = !is.na(product_inchikey) & product_inchikey %in% ef_keys$inchikey,
    hit_type = case_when(
      substrate_hit & product_hit ~ "both_EF",
      substrate_hit & !product_hit ~ "substrate_only_EF",
      !substrate_hit & product_hit ~ "product_only_EF",
      TRUE ~ "neither" # should be ~0 because we filtered earlier
    )
  ) %>%
  count(hit_type)
    # hit_type   n
    # 1           both_EF   1
    # 2   product_only_EF 146
    # 3 substrate_only_EF 229
Enterococcus_faecium_gen_scored_v1 <- met_gen_scored_v1
save(Enterococcus_faecium_gen_scored_v1, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/reaction_pairs/Enterococcus_faecium_gen_scored_v1")
############################################################
# OUTPUT OBJECTS YOU CARE ABOUT
#
# - reaction_pairs_ef        : MetaCyc edges anchored to EF metabolites (pre-genome join)
# - met_gen_with_met         : edges + genome join + EF stats (pre-filter)
# - met_gen_clean_v1         : versioned redundant-metabolite filtered table
# - met_gen_scored_v1  & Enterococcus_faecium_gen_scored_v1  : v1 confidence scored (KEEP components for later edits)
# - met_gen_scored_v1_distinct : optional deduplicated view
############################################################


