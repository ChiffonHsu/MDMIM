#1. Extract and Match Enterococcus faecium with MetaCyc
   #a. Extract data for Enterococcus faecium from collapsed_results
enterococcus_data <- collapsed_results$`Enterococcus faecium`

   #b. Get the InChIKeys from the bacterium's data
enterococcus_inchi_keys <- enterococcus_data$INCHIKEY_ID  # Adjust the column name if necessary

   #c. Match with MetaCyc edges for valid pairings (substrate and product)
matched_edges_enterococcus <- edges_inchikey %>%
  filter(substrate_inchikey %in% enterococcus_inchi_keys | product_inchikey %in% enterococcus_inchi_keys)

   #d. Check the number of matched pairs
nrow(matched_edges_enterococcus)  #94227
      # Filter rows where both substrate and product are not NA
matched_edges_enterococcus_valid <- matched_edges_enterococcus %>%
  filter(!is.na(substrate_inchikey) & !is.na(product_inchikey))

     # Check the number of valid matched pairs
nrow(matched_edges_enterococcus_valid)  #438

     # Check the first few rows to ensure the filtering worked
head(matched_edges_enterococcus_valid)
   #e. Check number of unique rows, no duplicates of rows
     # Get unique rows based on all columns
unique_matched_edges_enterococcus_valid <- matched_edges_enterococcus_valid %>%
  distinct()

     # Check the number of unique rows
nrow(unique_matched_edges_enterococcus_valid)  # This gives the number of unique (non-duplicate) rows

     # Check the first few rows of the unique data
head(unique_matched_edges_enterococcus_valid) #238

#2. Add in directions
reg_lookup <- collapsed_results$`Enterococcus faecium` %>%
  dplyr::select(INCHIKEY_ID, fc_final, pval_final) %>%
  dplyr::mutate(
    actual_reg = case_when(
      fc_final > 0 ~ "up",
      fc_final < 0 ~ "down",
      TRUE ~ "none"
    )
  ) %>%
  dplyr::rename(inchikey = INCHIKEY_ID)

annot_edges <- annot_edges %>%
  mutate(
    substrate_match_type = case_when(
      substrate_actual == "down" ~ "matched",
      substrate_actual == "up"   ~ "opposite",
      TRUE ~ "none"
    ),
    producct_match_type = case_when(
      product_actual == "down" ~ "opposite",
      product_actual == "up"   ~ "matched",
      TRUE ~ "none"
    )
  )

annot_edges <- unique_matched_edges_enterococcus_valid %>%
  # Add SUBSTRATE regulation
  left_join(reg_lookup,
            by = c("substrate_inchikey" = "inchikey"),
            suffix = c("", "_sub")) %>%
  # Add PRODUCT regulation
  left_join(reg_lookup,
            by = c("product_inchikey" = "inchikey"),
            suffix = c("_sub", "_prod"))
annot_edges <- annot_edges %>%
  dplyr::rename(
    substrate_actual = actual_reg_sub,
    product_actual   = actual_reg_prod
  )
annot_edges <- annot_edges %>%
  dplyr::mutate(
    substrate_match_type = case_when(
      substrate_actual == "down" ~ "matched",
      substrate_actual == "up"  ~ "opposite",
      TRUE ~ "none"),
    product_match_type = case_when(  
      product_actual == "down" ~ "opposite",
      product_actual == "up" ~ "matched",
      TRUE ~ "none")
    )
sum(annot_edges$substrate_match_type == "matched" &
      annot_edges$product_match_type == "matched") #NA
sum(annot_edges$substrate_match_type == "matched") #81
sum(annot_edges$product_match_type == "matched") #32
sum(annot_edges$substrate_match_type == "opposite" &
      annot_edges$product_match_type == "opposite") #0
sum(annot_edges$substrate_match_type == "opposite" |
      annot_edges$product_match_type == "opposite") #129

#3. Eliminate common compounds for stricter analysis
currency_names <- c(
  "WATER",
  "OXYGEN",
  "CARBON DIOXIDE",
  "AMMONIUM",
  "AMMONIA",
  "HYDROGEN ION",
  "PROTON",
  "PHOSPHATE",
  "DIPHOSPHATE",
  "ATP",
  "ADP",
  "AMP",
  "NAD",
  "NADH",
  "NADP",
  "NADPH",
  "FAD",
  "FADH2",
  "COENZYME A"
  # add more as needed
)

edges_filtered <- annot_edges %>%
  dplyr::filter(
    !toupper(substrate) %in% currency_names,
    !toupper(product)   %in% currency_names
  )
edges_filtered <- edges_filtered %>%
  filter(substrate_inchikey != product_inchikey)

cols_needed <- c(
  "fc_final_sub",
  "pval_final_sub",
  "substrate_actual",
  "fc_final_prod",
  "pval_final_prod",
  "product_actual"
)
count_valid <- edges_filtered %>%
  filter(if_all(all_of(cols_needed), ~ !is.na(.x))) %>%
  nrow()

count_valid #1
rows_valid <- edges_filtered %>%
  filter(if_all(all_of(cols_needed), ~ !is.na(.x)))

rows_valid
