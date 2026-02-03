#1. Ensure comparison_results have inchikey_id for annotation
annotation_lookup <- object_annot_best_inchi %>%
  dplyr::select(
    variable_id,
    KEGG.ID,
    HMDB.ID,
    CAS.ID,
    INCHIKEY_ID
  )

comparison_results <- map(
  comparison_results,
  ~ .x %>%
    dplyr::select(-any_of(c("KEGG.ID", "HMDB.ID", "CAS.ID", "INCHIKEY_ID"))) %>%
    left_join(annotation_lookup, by = "variable_id")
)

#2. Pathway enrichment
 #A. Define objects
df <- comparison_results$`Enterococcus faecium`
 #B. Extract background (ALL) metabolites
lcms_to_metacyc <- df %>%
  filter(!is.na(KEGG.ID)) %>%
  left_join(
    metacyc_valid %>% 
      dplyr::select(kegg, compound_id),
    by = c("KEGG.ID" = "kegg")
  )
 #CHECK
mean(!is.na(lcms_to_metacyc$compound_id)) #0.2857143

 #C. Extract significant metabolites 
df_sig <- df %>% filter(pval_final < 0.05)

sig_kegg <- df_sig$KEGG.ID
bg_kegg  <- df$KEGG.ID[!is.na(df$KEGG.ID)]

#3. Build metacyc compound
lcms_to_compound <- df %>%
  filter(!is.na(KEGG.ID)) %>%
  left_join(metacyc_nodes, by = c("KEGG.ID" = "kegg"))

df_pathway_map <- lcms_to_compound %>%
  inner_join(mc_compound_pathway, by = c("KEGG.ID" = "kegg"))
#check
head(df_pathway_map)
nrow(df_pathway_map) #1128

#4. Conduct 

# all LC-MS annotated → pathway mapped
bg <- df_pathway_map

# significant LC-MS metabolites
df_sig <- df %>% filter(pval_final < 0.05)

sig <- df_pathway_map %>%
  filter(variable_id %in% df_sig$variable_id)
#check
nrow(bg)
nrow(sig)
library(dplyr)
library(purrr)

pwy_enrich <- bg %>%
  group_by(pathway_id) %>%
  summarise(
    bg_count  = n(),
    sig_count = sum(variable_id %in% df_sig$variable_id)
  ) %>%
  mutate(
    pval = map2_dbl(
      sig_count, bg_count,
      ~ fisher.test(
        matrix(c(.x,
                 nrow(sig) - .x,
                 .y,
                 nrow(bg)  - .y),
               nrow = 2)
      )$p.value
    ),
    fdr = p.adjust(pval, method = "fdr")
  ) %>%
  arrange(fdr)


#inspect
pwy_enrich %>% filter(fdr < 0.1)





library(dplyr)
library(purrr)
library(tidyr)

# 1. Prepare mapping table
pwy_map <- df_pathway_map %>%
  dplyr::select(variable_id, INCHIKEY_ID, pathway_id)

# 2. Background table (all metabolites)
bg <- pwy_map %>%
  filter(INCHIKEY_ID %in% bg_kegg)

# 3. Significant metabolites
sig <- pwy_map %>%
  filter(INCHIKEY_ID %in% df_sig$INCHIKEY_ID)

# 4. Enrichment
pwy_enrich <- bg %>%
  group_by(pathway_id) %>%
  summarise(
    bg_count = n(),
    sig_count = sum(INCHIKEY_ID %in% df_sig$INCHIKEY_ID)
  ) %>%
  mutate(
    pval = map2_dbl(sig_count, bg_count, ~ fisher.test(
      matrix(c(.x, length(sig_kegg) - .x,
               .y, length(bg_kegg) - .y),
             nrow = 2)
    )$p.value),
    fdr = p.adjust(pval, method = "fdr")
  ) %>%
  arrange(fdr)

pwy_enrich #0

#5. Try to remove currency metabolites before enrichment
currency <- c("C00001","C00002","C00003","C00004","C00008","C00009",
              "C00010","C00011","C00020","C00080")

df_pathway_map_clean <- df_pathway_map %>%
  filter(!KEGG.ID %in% currency)
bg <- df_pathway_map_clean
sig <- df_pathway_map_clean %>%
  filter(variable_id %in% df_sig$variable_id)
library(dplyr)
library(purrr)

pwy_enrich <- bg %>%
  group_by(pathway_id) %>%
  summarise(
    bg_count  = n(),
    sig_count = sum(variable_id %in% df_sig$variable_id),
    .groups = "drop"
  ) %>%
  mutate(
    pval = map2_dbl(
      sig_count, bg_count,
      ~ fisher.test(
        matrix(c(.x,
                 nrow(sig) - .x,
                 .y,
                 nrow(bg)  - .y),
               nrow = 2)
      )$p.value
    ),
    fdr = p.adjust(pval, method = "fdr")
  ) %>%
  arrange(fdr)
pwy_enrich %>% filter(fdr < 0.1)

#Run using tidymass
#build massdataset
bacterium_name <- "Enterococcus faecium"
df_bac <- comparison_results[[bacterium_name]]
classes <- unique(df_bac$class)
classes
sample_info_all <- object_annot_best@sample_info
samples_enterococcus <- sample_info_all %>%
  dplyr::filter(grepl("10h", sample_id)) %>%   # selects sample IDs containing "10"
  dplyr::pull(sample_id)
samples_enterococcus
   #get expression data
samples_keep <- sample_info_enterococcus$sample_id
expr_enterococcus_clean <-
  object_annot_best@expression_data %>%
  tibble::rownames_to_column(var = "variable_id") %>%
  dplyr::select(variable_id, all_of(samples_keep))
  #get sample_info
sample_info_medium <- object_annot_best@sample_info %>%
  dplyr::filter(sample_id %in% samples_enterococcus) 
sample_info_enterococcus <- object_annot_best@sample_info %>%
  dplyr::filter(class %in% classes) 
sample_info_enterococcus <- dplyr::bind_rows(
  sample_info_medium,
  sample_info_enterococcus)
  #get variable_info
variable_info <- comparison_results$`Enterococcus faecium`
  #make massdataset for kegg pathway enrichment with tidymass
expr_enterococcus_mat <-
  expr_enterococcus_clean %>%
  tibble::column_to_rownames("variable_id")   # makes variable_id the rownames

variable_info_enterococcus <- 
  object_annot_best@variable_info %>%
  dplyr::filter(variable_id %in% rownames(expr_enterococcus_mat))
nrow(variable_info_enterococcus)
all(rownames(expr_enterococcus_mat) == variable_info_enterococcus$variable_id)

variable_info_enterococcus <- comparison_results$`Enterococcus faecium` %>%
  left_join(
    object_annot_best@variable_info %>% 
      dplyr::select(variable_id, everything()),
    by = "variable_id"
  ) %>%
  dplyr::select(variable_id, everything())

#check
nrow(variable_info_enterococcus)       # should be 4555
nrow(expr_enterococcus_mat)            # 4555
all(rownames(expr_enterococcus_mat) == variable_info_enterococcus$variable_id)  # TRUE


library(massdataset)
md_enterococcus <- create_mass_dataset(
  expression_data = expr_enterococcus_mat,
  sample_info     = sample_info_enterococcus,
  variable_info   = variable_info_enterococcus
)


#Run enrichment analysis using tidymass
library(tidyverse)
library(tidymass)
dir.create(path = "3_data_analysis/preliminary_determination_bacterial_metabolomics/pathway_enrichment", showWarnings = FALSE)

diff_metabolites <-
  md_enterococcus %>% 
  activate_mass_dataset(what = "variable_info") %>% 
  filter(pval_final < 0.05) %>% 
  extract_variable_info()
head(diff_metabolites)
 #load kegg pathway
data("kegg_hsa_pathway", package = "metpath")
kegg_hsa_pathway
get_pathway_class(kegg_hsa_pathway)
 #remove disease pathways
kegg_hsa_pathway_clean <- kegg_hsa_pathway %>%
  dplyr::filter(!stringr::str_starts(pathway_id, "hsa05"))
head(kegg_hsa_pathway_clean)
kegg_id <-
  diff_metabolites$KEGG.ID 
kegg_id <-
  kegg_id[!is.na(kegg_id)]
kegg_id
result <- 
  enrich_kegg(query_id = kegg_id, 
              query_type = "compound", 
              id_type = "KEGG",
              pathway_database = kegg_hsa_pathway_clean, 
              p_cutoff = 0.05, 
              p_adjust_method = "BH", 
              threads = 3)
#> 206 pathways.
#> 17  pathways p-values < 0.05
# Pyrimidine metabolism
# Amino sugar and nucleotide sugar metabolism
# ABC transporters
# Regulation of lipolysis in adipocytes
# Mineral absorption ... (only top 5 shows)
  #plot
enrich_bar_plot(object = result,
                x_axis = "p_value",
                cutoff = 0.05)
enrich_scatter_plot(object = result, y_axis = "p_value", y_axis_cutoff = 0.05)
enrich_network(
  object = result,
  point_size = "p_value",
  p_cutoff = 0.05,
  only_significant_pathway = TRUE
)
