#Try for one bacterial genome first before doing for all bacterial genomes
# Recommended: Use 'biomartr' or 'genbankr' to read the gbff
# This extracts the 'EC_number' tag from the features
library(biomartr)
library(seqinr)
library(genbankr)
# Example for one bacterium
file_path <- "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/bacterial_genome/Citrobacter_youngae/GCF_030294585.1_ASM3029458v1.gbff.gz"
  #unpack .gz compression
con <- gzfile(file_path, "r")
gb_lines <- readLines(con)
close(con)

  #1. Use grep to find the protein functions (products) and locus tags
    # We look for lines starting with whitespace + /product or /locus_tag
product_lines <- grep("/product=", gb_lines, value = TRUE)
locus_lines   <- grep("/locus_tag=", gb_lines, value = TRUE)

  #2. Clean up the text (remove the tags and quotes)
clean_products <- gsub('.*/product="([^"]+)".*', "\\1", product_lines)
clean_locus    <- gsub('.*/locus_tag="([^"]+)".*', "\\1", locus_lines)

  #3. Create a data frame (ensuring lengths match)
n <- min(length(clean_products), length(clean_locus))
cy_metabolic_profile <- data.frame(
  locus_tag = clean_locus[1:n],
  function_description = clean_products[1:n],
  stringsAsFactors = FALSE
)

  #4. View your results
head(cy_metabolic_profile, 10)



# Filter for CDS (Coding Sequences) that have an EC number
ec_list <- annotation[annotation$type == "CDS" & !is.na(Citrobacter_youngae_annot$EC_number), ]


#Now try to parse .fna files
# Use the microseq package to find ORFs in your downloaded .fna files
library(microseq)

# Example for one genome
genome <- readFasta("Citrobacter_youngae/GCF_xxx_genomic.fna")
orf_table <- findOrfs(genome, trans.tab = 11) # table 11 is for bacteria
# This gives you Start, End, and the Sequence of every possible gene.