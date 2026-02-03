# 1. Download .fna and .gbff and refseq files
# Packages
library(biomaRt)
library(rentrez)
exists("getGenome") #TRUE

bacterial_genome_dir <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/bacterial_genome"
#1. Extract bacteria names from the object hybrid_sig_summary
# 1. Extract the names from your tibble
bacteria_to_download <- hybrid_sig_summary$bacterium

# 2. Start the Loop
if (!dir.exists(bacterial_genome_dir)) dir.create(bacterial_genome_dir, recursive = TRUE)

bacteria_to_download <- hybrid_sig_summary$bacterium

for (species in bacteria_to_download) {
  
  # Create species-specific folder inside your defined directory
  clean_name <- gsub(" ", "_", species)
  species_folder <- file.path(bacterial_genome_dir, clean_name)
  
  if (!dir.exists(species_folder)) dir.create(species_folder)
  
  message("\n>>> Target: ", species)
  
  tryCatch({
    # 1. SEARCH: Find the RefSeq 'Reference' assembly (The Green Tick)
    # This query filters for the latest version and the 'reference' status
    query <- paste0(species, "[Organism] AND (latest[filter] AND \"reference genome\"[filter] AND all[filter] NOT atypical[filter])")
    search_res <- entrez_search(db = "assembly", term = query)
    
    # Fallback: If no strict 'Reference' exists, search for 'Representative'
    if (length(search_res$ids) == 0) {
      message("   No strict 'Reference' found, searching for 'Representative'...")
      query <- paste0(species, "[Organism] AND (latest[filter] AND \"representative genome\"[filter] AND all[filter] NOT atypical[filter])")
      search_res <- entrez_search(db = "assembly", term = query)
    }
    
    if (length(search_res$ids) == 0) stop("No RefSeq assembly found for this name.")
    
    # 2. EXTRACT: Get the FTP path for the RefSeq version
    sum_res <- entrez_summary(db = "assembly", id = search_res$ids[1])
    ftp_path <- sum_res$ftppath_refseq
    
    if (ftp_path == "") stop("FTP path is empty for this assembly.")
    
    # 3. CONSTRUCT: Build the exact URLs for .fna and .gbff
    base_name <- basename(ftp_path)
    fna_url  <- paste0(ftp_path, "/", base_name, "_genomic.fna.gz")
    gbff_url <- paste0(ftp_path, "/", base_name, "_genomic.gbff.gz")
    
    # 4. DOWNLOAD: Save directly into your folder
    message("   Downloading FNA...")
    download.file(fna_url, destfile = file.path(species_folder, paste0(base_name, ".fna.gz")), mode = "wb", quiet = TRUE)
    
    message("   Downloading GBFF...")
    download.file(gbff_url, destfile = file.path(species_folder, paste0(base_name, ".gbff.gz")), mode = "wb", quiet = TRUE)
    
    message("   Successfully saved to: ", species_folder)
    
  }, error = function(e) {
    message("!! Failed for ", species, ": ", e$message)
  })
}




# 
# 
# 
# 
# #3. For .gbff files, unzip to get
# library(gggenomes)
# 
# gbff_file <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/GCF_009734005.1_ASM973400v2_genomic.gbff"
# 
# gb_data <- read_gbk(gbff_file)
# # Features read
# # A tibble: 10 × 3
# # source type             n
# # <chr>  <chr>        <int>
# # 1 NA     CDS           2808
# # 2 NA     gene          2900
# # 3 NA     misc_binding     9
# # 4 NA     misc_feature     4
# # 5 NA     ncRNA            3
# # 6 NA     rRNA            18
# # 7 NA     region           2
# # 8 NA     regulatory      12
# # 9 NA     tRNA            70
# # 10 NA    tmRNA            1
# 
# cds <- gb_data[gb_data$type == "CDS", ]
# 
# # Extract the columns you care about
# gbff_ec <- data.frame(
#   locus_tag  = cds$locus_tag,
#   protein_id = cds$protein_id,
#   product    = cds$product,
#   EC_number  = cds$ec_number,   # lowercase 'ec_number'
#   stringsAsFactors = FALSE
# )
# 
# # Remove rows where EC_number is NA or empty
# gbff_ec <- gbff_ec[!is.na(gbff_ec$EC_number) & nzchar(gbff_ec$EC_number), ]
# 
# # Inspect
# head(gbff_ec)
# # locus_tag     protein_id                                              product   EC_number
# # 2  E6A31_RS00010 WP_002288361.1               DNA polymerase III subunit beta   2.7.7.7
# # 5  E6A31_RS00025 WP_002288364.1 DNA topoisomerase (ATP-hydrolyzing) subunit B   5.6.2.2
# # 6  E6A31_RS00030 WP_002321575.1                          DNA gyrase subunit A   5.6.2.2
# # 12 E6A31_RS00060 WP_002288373.1                      replicative DNA helicase  3.6.4.12
# # 13 E6A31_RS00065 WP_002288375.1                     adenylosuccinate synthase   6.3.4.4
# # 17 E6A31_RS00085 WP_002288384.1                              pyruvate oxidase   1.2.3.3
# 
# 
# #For gbff, get ec number and genome accessin number
# ec_protein_table <- cds %>%
#   as.data.frame() %>%
#   dplyr::select(protein_id, ec_number, product) %>%
#   dplyr::filter(!is.na(ec_number)) %>%
#   tidyr::separate_rows(ec_number, sep = ",") %>%  # split multiple ECs per protein
#   dplyr::mutate(
#     EC_number = trimws(ec_number),
#     genome_accession = "GCF_009734005.1",
#     organism = "Enterococcus faecium"
#   ) %>%
#   dplyr::distinct()
# #check output
# head(ec_protein_table)
# # protein_id     ec_number product                                      EC_number genome_accession organism
# # <chr>          <chr>     <chr>                                        <chr>     <chr>            <chr>   
# #   1 WP_002288361.1 2.7.7.7   DNA polymerase III subunit beta              2.7.7.7   GCF_009734005.1  Enteroc…
# # 2 WP_002288364.1 5.6.2.2   DNA topoisomerase (ATP-hydrolyzing) subunit… 5.6.2.2   GCF_009734005.1  Enteroc…
# # 3 WP_002321575.1 5.6.2.2   DNA gyrase subunit A                         5.6.2.2   GCF_009734005.1  Enteroc…
# # 4 WP_002288373.1 3.6.4.12  replicative DNA helicase                     3.6.4.12  GCF_009734005.1  Enteroc…
# # 5 WP_002288375.1 6.3.4.4   adenylosuccinate synthase                    6.3.4.4   GCF_009734005.1  Enteroc…
# # 6 WP_002288384.1 1.2.3.3   pyruvate oxidase                             1.2.3.3   GCF_009734005.1  Enteroc…
# 
# # 4. Compare .gbff to .fna results
# #1) Read .fna result which is tsv file
# library(dplyr)
# library(tidyr)
# prokka_tsv <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/prokka_Enterococcus_faecium/Enterococcus_faecium.tsv"
# 
# # read Prokka output
# prokka_ec <- read.delim(prokka_tsv, stringsAsFactors = FALSE)
# 
# # Check column names
# colnames(prokka_ec)
# # [1] "locus_tag" "ftype"     "length_bp" "gene"      "EC_number" "COG"       "product"  
# 
# #2) tidy prokka ec numbers
# prokka_ec_tidy <- prokka_ec %>%
#   dplyr::select(protein_id = locus_tag, EC_number, product) %>%
#   dplyr::filter(!is.na(EC_number) & nzchar(EC_number)) %>%
#   tidyr::separate_rows(EC_number, sep = ",") %>%
#   dplyr::mutate(EC_number = trimws(EC_number)) %>%
#   dplyr::distinct()
# 
# head(prokka_ec_tidy)
# # protein_id     EC_number product                             
# # <chr>          <chr>     <chr>                               
# # 1 PAANNGHF_00005 5.6.2.2   DNA gyrase subunit B                
# # 2 PAANNGHF_00010 3.1.4.59  Cyclic-di-AMP phosphodiesterase GdpP
# # 3 PAANNGHF_00012 5.6.2.3   Replicative DNA helicase DnaC       
# # 4 PAANNGHF_00013 6.3.4.4   Adenylosuccinate synthetase         
# # 5 PAANNGHF_00017 1.2.3.3   Pyruvate oxidase                    
# # 6 PAANNGHF_00029 5.4.2.6   Beta-phosphoglucomutase    
# 
# #3) Find common ec numbers in both gbff and tsv (fna)
# common_ec <- inner_join(
#   ec_protein_table,      # from GenBank
#   prokka_ec_tidy,        # from Prokka
#   by = "EC_number"
# )
# 
# head(common_ec)
# 
# #4) Find EC numbers in GenBank but missing in Prokka
# gb_only_ec <- anti_join(
#   ec_protein_table,
#   prokka_ec_tidy,
#   by = "EC_number"
# )
# head(gb_only_ec)
# # protein_id     ec_number product                                      EC_number genome_accession organism
# # <chr>          <chr>     <chr>                                        <chr>     <chr>            <chr>   
# # 1 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA… 6.3.2.2   GCF_009734005.1  Enteroc…
# # 2 WP_002287901.1 6.3.2.3   bifunctional glutamate--cysteine ligaseGshA… 6.3.2.3   GCF_009734005.1  Enteroc…
# # 3 WP_002297133.1 2.7.1.201 PTS system trehalose-specific EIIBC compone… 2.7.1.201 GCF_009734005.1  Enteroc…
# # 4 WP_002294234.1 5.4.99.5  chorismate mutase                            5.4.99.5  GCF_009734005.1  Enteroc…
# # 5 WP_002296548.1 2.3.-.-   GNAT family N-acetyltransferase              2.3.-.-   GCF_009734005.1  Enteroc…
# # 6 WP_002296567.1 1.2.1.10  bifunctional acetaldehyde-CoA/alcoholdehydr… 1.2.1.10  GCF_009734005.1  Enteroc…
# 
# #5) Find EC numbers in Prokka but missing in GenBank
# prokka_only_ec <- anti_join(
#   prokka_ec_tidy,
#   ec_protein_table,
#   by = "EC_number"
# )
# 
# head(prokka_only_ec)
# # A tibble: 6 × 3
# # protein_id     EC_number product                                        
# # <chr>          <chr>     <chr>                                          
# #   1 PAANNGHF_00010 3.1.4.59  Cyclic-di-AMP phosphodiesterase GdpP           
# # 2 PAANNGHF_00012 5.6.2.3   Replicative DNA helicase DnaC                  
# # 3 PAANNGHF_00043 3.2.1.26  Sucrose-6-phosphate hydrolase                  
# # 4 PAANNGHF_00045 2.7.1.4   Putative fructokinase                          
# # 5 PAANNGHF_00048 2.1.1.172 Ribosomal RNA small subunit methyltransferase C
# # 6 PAANNGHF_00066 3.6.5.3   Elongation factor Tu   
# #6) Summarize overlap
# data.frame(
#   total_gb = nrow(ec_protein_table),
#   total_prokka = nrow(prokka_ec_tidy),
#   common = nrow(common_ec),
#   gb_only = nrow(gb_only_ec),
#   prokka_only = nrow(prokka_only_ec)
# )
# #      total_gb total_prokka common  gb_only  prokka_only
# # 1     747          861      2258     115         272
# save(gb_only_ec, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/gb_only_ec")
# save(prokka_only_ec, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/prokka_only_ec")
