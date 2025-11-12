library(dplyr)
library(stringr)
library(readr)

### ---- Parse compounds.dat ----
metacyc_cpd <- "~/compound_database_shenlab/metacyc/data/data/compounds.dat"
compound_lines <- readLines(metacyc_cpd, encoding = "latin1")
# Show the first 100 non-comment lines
 compound_lines_no_comments <- compound_lines[!str_starts(compound_lines, "#")]
 head(compound_lines_no_comments, 100)
  #Each block is separated by "//" and has multiple fields (key-value pairs).
  #split the file into blocks and then extract fields from each block.
 # remove comments first
 compound_lines_no_comments <- compound_lines[!grepl("^#", compound_lines)]
 
 # split blocks by "//"
 compound_blocks <- split(compound_lines_no_comments, 
                          cumsum(compound_lines_no_comments == "//"))
 
 # function to parse a single block
 parse_block <- function(block) {
   block <- block[block != "//"]  # drop the //
   fields <- strsplit(block, " - ")
   df <- data.frame(
     key = sapply(fields, `[`, 1),
     value = sapply(fields, `[`, 2),
     stringsAsFactors = FALSE
   )
   # return as named list
   out <- as.list(setNames(df$value, df$key))
   return(out)
 }
 
 # parse all blocks
 compound_list <- lapply(compound_blocks, parse_block)
 
 # convert to dataframe (one compound per row)
 parse_compound_block <- function(block) {
   block <- block[block != "//"]
   x <- list(
     compound_id = NA,
     compound_name = NA,
     smiles = NA,
     inchi = NA,
     chebi = NA
   )
   
   for (line in block) {
     if (str_starts(line, "UNIQUE-ID")) x$compound_id <- str_split_fixed(line, " - ", 2)[,2]
     if (str_starts(line, "COMMON-NAME")) x$compound_name <- str_split_fixed(line, " - ", 2)[,2]
     if (str_starts(line, "SMILES")) x$smiles <- str_split_fixed(line, " - ", 2)[,2]
     if (str_starts(line, "NON-STANDARD-INCHI")) x$inchi <- str_split_fixed(line, " - ", 2)[,2]
     if (str_starts(line, "DBLINKS") & str_detect(line, "CHEBI")) {
       x$chebi <- sub('.*"(.*)".*', '\\1', line)
     }
   }
   x
 }
 
 compound_list <- lapply(compound_blocks, parse_compound_block)
 compound_df <- bind_rows(lapply(compound_list, as.data.frame))
   #check for missing values in each columns
 colSums(is.na(compound_df[, c("compound_id", "compound_name", "smiles", "inchi", "chebi")]))
 
 
### ---- Parse reactions.dat ----
  reaction_file <- "~/compound_database_shenlab/metacyc/data/data/reactions.dat"
  reaction_lines <- readLines(reaction_file, encoding = "latin1")
  # Remove comment lines
  reaction_lines_no_comments <- reaction_lines[!str_starts(reaction_lines, "#")]
  
  # Split reactions by "//"
  reaction_blocks <- split(reaction_lines_no_comments, cumsum(reaction_lines_no_comments == "//"))
  
  # Function to parse a single reaction block
  parse_reaction_block <- function(block, compound_map) {
    block <- block[block != "//"]  # drop the separator
    entry <- list(LEFT = character(), RIGHT = character(), REVERSIBLE = FALSE)
    
    for (line in block) {
      if (str_starts(line, "LEFT")) {
        entry$LEFT <- c(entry$LEFT, str_split_fixed(line, " - ", 2)[,2])
      } else if (str_starts(line, "RIGHT")) {
        entry$RIGHT <- c(entry$RIGHT, str_split_fixed(line, " - ", 2)[,2])
      } else if (str_starts(line, "REVERSIBLE?")) {
        entry$REVERSIBLE <- str_detect(line, "T")
      }
    }
    
    edges <- list()
    for (l in entry$LEFT) {
      for (r in entry$RIGHT) {
        if (l %in% compound_map$compound_id & r %in% compound_map$compound_id) {
          lname <- compound_map$compound_name[compound_map$compound_id == l]
          rname <- compound_map$compound_name[compound_map$compound_id == r]
          edges <- append(edges, list(data.frame(substrate = lname, product = rname)))
          if (entry$REVERSIBLE) {
            edges <- append(edges, list(data.frame(substrate = rname, product = lname)))
          }
        }
      }
    }
    
    if (length(edges) > 0) bind_rows(edges) else NULL
  }
  
  # Parse all reactions
  reaction_edges <- lapply(reaction_blocks, parse_reaction_block, compound_map = compound_df)
  
   # Combine into a single dataframe
  reaction_edges_df <- bind_rows(reaction_edges)
  
  # Quick check
  head(reaction_edges_df)
  
  #filter out missing ids
  reaction_edges_df <- reaction_edges_df %>% 
    filter(!is.na(substrate) & !is.na(product))
  
  #filter out exact duplicates 
  reaction_edges_df <- reaction_edges_df %>% distinct()
  
 ############################################################make into node and edge table
  library(dplyr)
  
  # ---- Prepare edges ----
  # edges = reaction pairs
  edges_df <- reaction_edges_df %>%
    filter(!is.na(substrate) & !is.na(product)) %>%
    distinct() %>%
    dplyr::rename(from = substrate,
                  to   = product)
  
  # Quick check
  head(edges_df)
 
  #merge compound_df to edges_df
  colnames(compound_df)
  library(dplyr)
  nodes_df <- data.frame(name = unique(c(edges_df$from, edges_df$to)),
                         stringsAsFactors = FALSE)
  nodes_df <- nodes_df %>%
    left_join(compound_df, by = c("name" = "compound_name"))
  
  # join metadata for "from" compounds
 edges_with_from <- edges_df %>%
   left_join(nodes_df, by = c("from" = "name")) 
 colnames(edges_with_from)
 edges_with_from <- edges_df %>%
   left_join(nodes_df, by = c("from" = "name")) %>%
 dplyr::rename(
          from_compound_id = compound_id,
          from_smiles      = smiles,
          from_inchi      = inchi,
          from_chebi       = chebi,
          from_cas        =)
colnames(edges_with_from)
  
  # join metadata for "to" compounds
edges_with_meta <- edges_with_from %>%
  left_join(nodes_df, by = c("to" = "name")) 
colnames(edges_with_meta)
edges_with_meta <- edges_with_from %>%
    left_join(nodes_df, by = c("to" = "name")) %>%
    dplyr::rename(
           to_compound_id = compound_id,
           to_smiles      = smiles,
           to_inchi       = inchi,
           to_chebi       = chebi)
  
head(edges_with_meta)
colnames(edges_with_meta)
metacyc_edges_all <- edges_with_meta
saveRDS(metacyc_edges_all, file = "~/compound_database_shenlab/metacyc/metacyc_edges_all.rds")

#remove redundant reaction pairs => Zero missing values for two and fro, 
metacyc_edges_filtered_unique <- unique(metacyc_edges_filtered)
saveRDS(metacyc_edges_filtered_unique, file = "~/compound_database_shenlab/metacyc/metacyc_edges_filtered_unique.rds")
count(metacyc_edges_filtered_unique) #24052

     # ---- Prepare nodes ----
  # nodes = unique metabolites
compound_df <- as.data.frame(compound_df)
nodes_df <- compound_df %>%
  dplyr::select(compound_id,
         compound_name,
         smiles,
         inchi,
         chebi) %>%
  distinct()
 colnames(nodes_df)   
 metacyc_nodes_all <- nodes_df
 saveRDS(metacyc_nodes_all, file = "~/compound_database_shenlab/metacyc/metacyc_nodes_all.rds") 

 