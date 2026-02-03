#Merge the annotations together
gb_ec <- gb_only_ec %>%
  transmute(
    gene_id = NA_character_,       
    protein_id,
    EC_number = ec_number,
    product,
    source = "GenBank"
  )

prokka_ec2 <- prokka_ec %>%
  filter(!is.na(EC_number)) %>%
  transmute(
    gene_id = locus_tag,
    protein_id = locus_tag,         
    EC_number,
    product,
    source = "Prokka"
  )
#2. Tidy ec numbers
enzyme_raw <- bind_rows(gb_ec, prokka_ec2) %>%
  separate_rows(EC_number, sep = ",") %>%
  mutate(EC_number = trimws(EC_number)) %>%
  filter(
    !is.na(EC_number),
    EC_number != "",
    EC_number != "NA"
  )
#3. Assign confidence scores
  #1) Count evidence per protein–EC pair
enzyme_conf <- enzyme_raw %>%
  group_by(protein_id, EC_number) %>%
  summarise(
    gene_id = first(gene_id),
    product = first(product),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sources = n_distinct(source),
    .groups = "drop"
  )

  #2) Set the confidence rules
enzyme_catalogue <- enzyme_conf %>%
  mutate(
    confidence = case_when(
      n_sources >= 2                          ~ "high",
      grepl("\\.-\\.-", EC_number)            ~ "low",
      n_sources == 1                          ~ "medium",
      TRUE                                    ~ "low"
    )
  )

  #3) KEGG validation
library(KEGGREST)

enzyme_catalogue <- enzyme_catalogue %>%
  rowwise() %>%
  mutate(
    kegg_reactions = list(
      tryCatch(
        keggLink("rn", paste0("ec:", EC_number)),
        error = function(e) NA
      )
    ),
    kegg_supported = !all(is.na(kegg_reactions))
  ) %>%
  ungroup()

  #4) Upgrade confidence
enzyme_catalogue <- enzyme_catalogue %>%
  mutate(
    confidence = ifelse(
      confidence == "medium" & kegg_supported,
      "medium-high",
      confidence
    )
  )

  #5) Check
enzyme_catalogue %>%
  select(
    gene_id,
    protein_id,
    EC_number,
    product,
    sources,
    confidence
  ) %>%
  arrange(confidence, EC_number)
save(enzyme_catalogue, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/enzyme_catalogue")
