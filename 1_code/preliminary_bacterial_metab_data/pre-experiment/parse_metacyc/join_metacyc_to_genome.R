#Join metacyc reaction data with the ec numbers to genome enzyme table
reaction_protein_map <- edges_annot %>%
  filter(!is.na(ec_number)) %>%
  left_join(
    ec_protein_table,
    by = c("ec_number" = "EC_number")
  )
head(reaction_protein_map)
colnames(reaction_protein_map)
# [1] "reaction_id"        "substrate"          "product.x"          "substrate_inchikey" "product_inchikey"  
# [6] "ec_number"          "protein_id"         "ec_number.y"        "product.y"          "genome_accession"  
# [11] "organism"