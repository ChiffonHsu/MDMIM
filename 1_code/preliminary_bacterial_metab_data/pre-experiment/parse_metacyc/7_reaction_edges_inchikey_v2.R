#This script includes the addition of reaction confidence level. For example, having a negative expression of metabolite for substrate,
#having a positive expressinon of metabolite for product, adds on to the confidence level.
#2nd step is to clear redundant metabolites

load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/reaction_pairs/Enterococcus_faecium_gen_scored_v1")

#1. add reaction confidence level
score_reaction_v1 <- function(df) {
  df %>%
    mutate(
      # 1) Collapse substrate/product into reaction-level metrics
      reaction_pval = pmin(
        substrate_pval_final,
        product_pval_final,
        na.rm = TRUE
      ),
      
      reaction_fc = pmax(
        abs(substrate_fc_final),
        abs(product_fc_final),
        na.rm = TRUE
      ),
      
      # Handle rows where both were NA
      reaction_pval = if_else(is.infinite(reaction_pval), NA_real_, reaction_pval),
      reaction_fc   = if_else(is.infinite(reaction_fc), NA_real_, reaction_fc),
      
      # 2) Reaction score (v1)
      reaction_score_v1 = case_when(
        # No metabolite evidence
        is.na(reaction_pval) & is.na(reaction_fc) ~ 0L,
        
        # Not significant
        reaction_pval > 0.05 ~ 1L,
        
        # Significant but small effect
        reaction_pval <= 0.05 & reaction_fc < 1 ~ 2L,
        
        # Moderate changes
        reaction_pval <= 0.05 & reaction_fc >= 1 & reaction_fc < 2 ~ 3L,
        
        reaction_pval <= 0.05 & reaction_fc >= 2 & reaction_fc < 5 ~ 4L,
        
        # Large changes
        reaction_pval <= 0.05 & reaction_fc >= 5 & reaction_fc < 10 ~ 5L,
        
        # Extreme changes
        reaction_pval <= 0.05 & reaction_fc >= 10 ~ 6L,
        
        TRUE ~ 0L
      )
    )
}

met_gen_scored_v2 <- met_gen_scored_v1 %>%
  score_reaction_v1()
nrow(met_gen_scored_v2) #376
Enterococcus_faecium_gen_scored_v2 <- met_gen_scored_v2
save(Enterococcus_faecium_gen_scored_v2, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/reaction_pairs/Enterococcus_faecium_gen_scored_v2")

met_gen_scored_v2_df <- met_gen_scored_v2 %>%
  filter(met_gen_scored_v2$substrate_pval_final < 0.05 | met_gen_scored_v2$product_pval_final < 0.05)
nrow(met_gen_scored_v2_df) #121
Enterococcus_faecium_gen_scored_v2_df <- met_gen_scored_v2_df
save(Enterococcus_faecium_gen_scored_v2_df, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/reaction_pairs/Enterococcus_faecium_gen_scored_v2_df")

#2. Remove redundant metabolites
normalize_name <- function(x) {
  x %>%
    str_trim() %>%
    str_to_upper() %>%
    str_replace_all("\\s+", " ")
}
  #2A. Expand currency metabolite list
currency_metabolites_v2 <- normalize_name(c(
  # core
  "PROTON","H+","WATER","H2O","CARBON-DIOXIDE","CO2","OXYGEN","O2",
  "NAD","NADH","NADP","NADPH",
  "ATP","ADP","AMP","GTP","GDP","GMP",
  "PHOSPHATE","ORTHOPHOSPHATE","PI","PPI","PPi",
  "HCO3","BICARBONATE",
  "AMMONIA","AMMONIUM",
  "SULFATE","SO4",
  "HYDROGEN-PEROXIDE",
  "CO-A","COA"
))
  #2B. Remove non metabolite patterns
non_metabolite_patterns_v2 <- c(
  # placeholders / MetaCyc generic objects
  "^CPD0?-",            # CPD-, CPD0- placeholders
  "^ORSEL",             # ORSEL... (often pathway/PKS objects)
  "PKS",                # polyketide synthase objects (enzyme-like nodes)
  
  # generic electron / acceptor/donor pools
  "^E-$|^E\\-$",
  "DONOR", "ACCEPTOR",
  "OXIDIZED", "REDUCED",
  "FERREDOXIN", "CYTOCHROME", "NAPC",
  "HEMOPROTEIN", "REDUCTASES",
  
  # macromolecules / conjugates you usually don't want in small-molecule network
  "PROTEIN", "PEPTIDOGLYCAN", "UBIQUITIN", "CARRIER",
  "TRNA|tRNA"
)
generic_class_patterns_v2 <- c(
  # broad classes / bins
  "LONG[- ]CHAIN", "FATTY ACIDS", "ALKANES", "ALDEHYDES", "ALCOHOLS",
  "DIACYLGLYCEROL", "TRIACYLGLYCEROL|\\bTGS\\b",
  "QUINONES|UBIQUINONES|MENAQUINONES",
  "MONOCARBOXYLATES", "FLAVANONES|FLAVONES",
  
  # inorganic / polymeric / mixtures
  "GLUCURONANS", "GLYCEROPHOSPHODIESTERS",
  "ORTHOPHOSPHORIC"
)

    #2C. Remove
met_gen_flagged_v2 <- met_gen_with_met %>%
  mutate(
    substrate_n = normalize_name(substrate),
    product_n   = normalize_name(product),
    
    is_currency_v2 =
      substrate_n %in% currency_metabolites_v2 |
      product_n   %in% currency_metabolites_v2,
    
    is_non_metabolite_v2 =
      str_detect(substrate_n, paste(non_metabolite_patterns_v2, collapse="|")) |
      str_detect(product_n,   paste(non_metabolite_patterns_v2, collapse="|")),
    
    is_generic_class_v2 =
      str_detect(substrate_n, paste(generic_class_patterns_v2, collapse="|")) |
      str_detect(product_n,   paste(generic_class_patterns_v2, collapse="|")),
    
    # helpful “why removed” label
    removal_reason_v2 = case_when(
      is_currency_v2 ~ "currency",
      is_non_metabolite_v2 ~ "non_metabolite_entity",
      is_generic_class_v2 ~ "generic_class",
      TRUE ~ "keep"
    )
  ) %>%
  select(-substrate_n, -product_n)

# Inspect what you would remove
met_gen_flagged_v2 %>%
  count(removal_reason_v2, sort = TRUE)

met_gen_flagged_v2 %>%
  filter(removal_reason_v2 != "keep") %>%
  count(substrate, sort = TRUE) %>%
  head(50)

met_gen_flagged_v2 %>%
  filter(removal_reason_v2 != "keep") %>%
  count(product, sort = TRUE) %>%
  head(50)

# Apply filter (hard-remove all three categories)
met_gen_clean_v2 <- met_gen_flagged_v2 %>%
  filter(removal_reason_v2 == "keep")
    #2D. Remove not so important ones or penalize them
library(tidyr)

hub_table <- met_gen_with_met %>%
  select(substrate, product) %>%
  pivot_longer(cols = c(substrate, product), values_to = "metabolite") %>%
  filter(!is.na(metabolite), metabolite != "") %>%
  mutate(metabolite = normalize_name(metabolite)) %>%
  count(metabolite, name = "degree_proxy") %>%
  arrange(desc(degree_proxy))

# Choose a threshold (example: remove metabolites appearing in > 200 edges)
hub_cutoff <- 200

hub_metabolites <- hub_table %>%
  filter(degree_proxy > hub_cutoff) %>%
  pull(metabolite)

met_gen_clean_v2_hub <- met_gen_clean_v2 %>%
  mutate(
    substrate_n = normalize_name(substrate),
    product_n   = normalize_name(product)
  ) %>%
  filter(
    !(substrate_n %in% hub_metabolites),
    !(product_n %in% hub_metabolites)
  ) %>%
  select(-substrate_n, -product_n)

Enterococcus_faecium_gen_scored_v2_cleaned <- met_gen_clean_v2_hub
save(met_gen_clean_v2_hub, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/reaction_pairs/Enterococcus_faecium_gen_scored_v2_cleaned")
