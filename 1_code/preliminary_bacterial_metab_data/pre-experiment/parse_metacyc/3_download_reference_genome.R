# 1. Download .fna and .gbff and refseq files
    # URLs
fna_url  <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/734/005/GCF_009734005.1_ASM973400v2/GCF_009734005.1_ASM973400v2_genomic.fna.gz"
gbff_url <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/734/005/GCF_009734005.1_ASM973400v2/GCF_009734005.1_ASM973400v2_genomic.gbff.gz"

    # Set destination directory
dest_dir <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence"

    # Must ensure that the directory exists
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

    # Set destination paths
fna_dest  <- file.path(dest_dir, "GCF_009734005.1_ASM973400v2_genomic.fna.gz")
gbff_dest <- file.path(dest_dir, "GCF_009734005.1_ASM973400v2_genomic.gbff.gz")

    # Download
download.file(fna_url,  destfile = fna_dest,  mode = "wb")
download.file(gbff_url, destfile = gbff_dest, mode = "wb")

    # Unzip in R (Use R.utils)
R.utils::gunzip(fna_dest, remove = FALSE)  
R.utils::gunzip(gbff_dest, remove = FALSE)

#2. Use Prokka for gene calling for .fna files
#Run this entire step 2 at local r studio with the following command:
library(Biostrings)

fna_file <- "/Volumes/ClareHsu/converted_rawdata/1_genome_sequence/GCF_009734005.1_ASM973400v2_genomic.fna"

# Read the sequences
dna <- readDNAStringSet(fna_file)

# Write to a new FASTA
temp_fasta <- "ASM973400v2_modified.fna"
writeXStringSet(dna, temp_fasta)

# Run Prokka
library(Biostrings)
outdir <- "prokka_Enterococcus_faecium"
prefix <- "Enterococcus_faecium"

prokka_path <- "/Users/chifang/Applications/miniconda/envs/prokka/bin/prokka"

cmd <- sprintf(
  "source /Users/chifang/Applications/miniconda/bin/activate prokka && prokka --outdir %s --prefix %s --genus Enterococcus --species faecium --strain ASM973400v2 %s",
  outdir, prefix, temp_fasta
)

system(cmd)
#Output of step 2: should have the following files all uploaded to the files
# [1] "Enterococcus_faecium.err"   "Enterococcus_faecium.faa"   "Enterococcus_faecium.ffn"  
# [4] "Enterococcus_faecium.fna"   "Enterococcus_faecium.fsa"   "Enterococcus_faecium.gbf-r"
# [7] "Enterococcus_faecium.gbk"   "Enterococcus_faecium.gff"   "Enterococcus_faecium.log"  
# [10] "Enterococcus_faecium.sqn"   "Enterococcus_faecium.tbl"   "Enterococcus_faecium.tsv"  
# [13] "Enterococcus_faecium.txt"

#3. For .gbff files, unzip to get
library(gggenomes)

gbff_file <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/GCF_009734005.1_ASM973400v2_genomic.gbff"

gb_data <- read_gbk(gbff_file)
# Features read
# A tibble: 10 × 3
# source type             n
# <chr>  <chr>        <int>
# 1 NA     CDS           2808
# 2 NA     gene          2900
# 3 NA     misc_binding     9
# 4 NA     misc_feature     4
# 5 NA     ncRNA            3
# 6 NA     rRNA            18
# 7 NA     region           2
# 8 NA     regulatory      12
# 9 NA     tRNA            70
# 10 NA    tmRNA            1

    cds <- gb_data[gb_data$type == "CDS", ]

    # Extract the columns you care about
    gbff_ec <- data.frame(
      locus_tag  = cds$locus_tag,
      protein_id = cds$protein_id,
      product    = cds$product,
      EC_number  = cds$ec_number,   # lowercase 'ec_number'
      stringsAsFactors = FALSE
    )

    # Remove rows where EC_number is NA or empty
    gbff_ec <- gbff_ec[!is.na(gbff_ec$EC_number) & nzchar(gbff_ec$EC_number), ]

    # Inspect
    head(gbff_ec)
# locus_tag     protein_id                                              product   EC_number
# 2  E6A31_RS00010 WP_002288361.1               DNA polymerase III subunit beta   2.7.7.7
# 5  E6A31_RS00025 WP_002288364.1 DNA topoisomerase (ATP-hydrolyzing) subunit B   5.6.2.2
# 6  E6A31_RS00030 WP_002321575.1                          DNA gyrase subunit A   5.6.2.2
# 12 E6A31_RS00060 WP_002288373.1                      replicative DNA helicase  3.6.4.12
# 13 E6A31_RS00065 WP_002288375.1                     adenylosuccinate synthase   6.3.4.4
# 17 E6A31_RS00085 WP_002288384.1                              pyruvate oxidase   1.2.3.3


     #For gbff, get ec number and genome accessin number
    ec_protein_table <- cds %>%
      as.data.frame() %>%
      dplyr::select(protein_id, ec_number, product) %>%
      dplyr::filter(!is.na(ec_number)) %>%
      tidyr::separate_rows(ec_number, sep = ",") %>%  # split multiple ECs per protein
      dplyr::mutate(
        EC_number = trimws(ec_number),
        genome_accession = "GCF_009734005.1",
        organism = "Enterococcus faecium"
      ) %>%
      dplyr::distinct()
    #check output
    head(ec_protein_table)
        # protein_id     ec_number product                                      EC_number genome_accession organism
        # <chr>          <chr>     <chr>                                        <chr>     <chr>            <chr>   
        #   1 WP_002288361.1 2.7.7.7   DNA polymerase III subunit beta              2.7.7.7   GCF_009734005.1  Enteroc…
        # 2 WP_002288364.1 5.6.2.2   DNA topoisomerase (ATP-hydrolyzing) subunit… 5.6.2.2   GCF_009734005.1  Enteroc…
        # 3 WP_002321575.1 5.6.2.2   DNA gyrase subunit A                         5.6.2.2   GCF_009734005.1  Enteroc…
        # 4 WP_002288373.1 3.6.4.12  replicative DNA helicase                     3.6.4.12  GCF_009734005.1  Enteroc…
        # 5 WP_002288375.1 6.3.4.4   adenylosuccinate synthase                    6.3.4.4   GCF_009734005.1  Enteroc…
        # 6 WP_002288384.1 1.2.3.3   pyruvate oxidase                             1.2.3.3   GCF_009734005.1  Enteroc…

# 4. Compare .gbff to .fna results
    #1) Read .fna result which is tsv file
    library(dplyr)
    library(tidyr)
    prokka_tsv <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/prokka_Enterococcus_faecium/Enterococcus_faecium.tsv"
    
    # read Prokka output
    prokka_ec <- read.delim(prokka_tsv, stringsAsFactors = FALSE)
    
    # Check column names
    colnames(prokka_ec)
    # [1] "locus_tag" "ftype"     "length_bp" "gene"      "EC_number" "COG"       "product"  
    
    #2) tidy prokka ec numbers
    prokka_ec_tidy <- prokka_ec %>%
      dplyr::select(protein_id = locus_tag, EC_number, product) %>%
      dplyr::filter(!is.na(EC_number) & nzchar(EC_number)) %>%
      tidyr::separate_rows(EC_number, sep = ",") %>%
      dplyr::mutate(EC_number = trimws(EC_number)) %>%
      dplyr::distinct()
    
    head(prokka_ec_tidy)
    # protein_id     EC_number product                             
    # <chr>          <chr>     <chr>                               
    # 1 PAANNGHF_00005 5.6.2.2   DNA gyrase subunit B                
    # 2 PAANNGHF_00010 3.1.4.59  Cyclic-di-AMP phosphodiesterase GdpP
    # 3 PAANNGHF_00012 5.6.2.3   Replicative DNA helicase DnaC       
    # 4 PAANNGHF_00013 6.3.4.4   Adenylosuccinate synthetase         
    # 5 PAANNGHF_00017 1.2.3.3   Pyruvate oxidase                    
    # 6 PAANNGHF_00029 5.4.2.6   Beta-phosphoglucomutase    
    
    #3) Find common ec numbers in both gbff and tsv (fna)
    common_ec <- inner_join(
      ec_protein_table,      # from GenBank
      prokka_ec_tidy,        # from Prokka
      by = "EC_number"
    )
    
    head(common_ec)
    
    #4) Find EC numbers in GenBank but missing in Prokka
    gb_only_ec <- anti_join(
      ec_protein_table,
      prokka_ec_tidy,
      by = "EC_number"
    )
    head(gb_only_ec)
    # protein_id     ec_number product                                      EC_number genome_accession organism
    # <chr>          <chr>     <chr>                                        <chr>     <chr>            <chr>   
    # 1 WP_002287901.1 6.3.2.2   bifunctional glutamate--cysteine ligaseGshA… 6.3.2.2   GCF_009734005.1  Enteroc…
    # 2 WP_002287901.1 6.3.2.3   bifunctional glutamate--cysteine ligaseGshA… 6.3.2.3   GCF_009734005.1  Enteroc…
    # 3 WP_002297133.1 2.7.1.201 PTS system trehalose-specific EIIBC compone… 2.7.1.201 GCF_009734005.1  Enteroc…
    # 4 WP_002294234.1 5.4.99.5  chorismate mutase                            5.4.99.5  GCF_009734005.1  Enteroc…
    # 5 WP_002296548.1 2.3.-.-   GNAT family N-acetyltransferase              2.3.-.-   GCF_009734005.1  Enteroc…
    # 6 WP_002296567.1 1.2.1.10  bifunctional acetaldehyde-CoA/alcoholdehydr… 1.2.1.10  GCF_009734005.1  Enteroc…
    
    #5) Find EC numbers in Prokka but missing in GenBank
    prokka_only_ec <- anti_join(
      prokka_ec_tidy,
      ec_protein_table,
      by = "EC_number"
    )
    
    head(prokka_only_ec)
    # A tibble: 6 × 3
    # protein_id     EC_number product                                        
    # <chr>          <chr>     <chr>                                          
    #   1 PAANNGHF_00010 3.1.4.59  Cyclic-di-AMP phosphodiesterase GdpP           
    # 2 PAANNGHF_00012 5.6.2.3   Replicative DNA helicase DnaC                  
    # 3 PAANNGHF_00043 3.2.1.26  Sucrose-6-phosphate hydrolase                  
    # 4 PAANNGHF_00045 2.7.1.4   Putative fructokinase                          
    # 5 PAANNGHF_00048 2.1.1.172 Ribosomal RNA small subunit methyltransferase C
    # 6 PAANNGHF_00066 3.6.5.3   Elongation factor Tu   
    #6) Summarize overlap
    data.frame(
      total_gb = nrow(ec_protein_table),
      total_prokka = nrow(prokka_ec_tidy),
      common = nrow(common_ec),
      gb_only = nrow(gb_only_ec),
      prokka_only = nrow(prokka_only_ec)
    )
    #      total_gb total_prokka common  gb_only  prokka_only
    # 1     747          861      2258     115         272
  save(gb_only_ec, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/gb_only_ec")
  save(prokka_only_ec, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/1_genome_sequence/prokka_only_ec")
  