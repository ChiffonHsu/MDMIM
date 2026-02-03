#1. Buid MetaCyc compound → kegg map
metacyc_valid <- readRDS("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_valid.rds")
metacyc_nodes <- metacyc_valid %>%
  dplyr::select(compound_id, kegg) %>%
  filter(!is.na(kegg)) %>%
  distinct()

  #check
  #nrow(metacyc_nodes)     # 10570
  #length(unique(metacyc_nodes$kegg)) #8717
#2. Parse metacyc reactions
parse_reactions <- function(file) {
  lines <- readLines(file)
  
  rxn_id <- NULL
  left <- c()
  right <- c()
  
  out <- list()   # store raw parsed reactions
  
  for (line in lines) {
    
    # New reaction starts
    if (startsWith(line, "UNIQUE-ID -")) {
      # save previous reaction
      if (!is.null(rxn_id)) {
        out[[rxn_id]] <- list(substrates = left, products = right)
      }
      
      # reset for new reaction
      rxn_id <- sub("UNIQUE-ID - ", "", line)
      left <- c()
      right <- c()
    }
    
    # Capture LEFT (substrates)
    if (startsWith(line, "LEFT -")) {
      left <- c(left, sub("LEFT - ", "", line))
    }
    
    # Capture RIGHT (products)
    if (startsWith(line, "RIGHT -")) {
      right <- c(right, sub("RIGHT - ", "", line))
    }
  }
  
  # save final reaction
  if (!is.null(rxn_id)) {
    out[[rxn_id]] <- list(substrates = left, products = right)
  }
  
  # Expand each reaction into substrate–product pairs
  df <- purrr::map_dfr(
    names(out),
    function(rid) {
      subs <- out[[rid]]$substrates
      prods <- out[[rid]]$products
      
      # If no substrates OR no products → still output so we can inspect later
      if (length(subs) == 0) subs <- NA_character_
      if (length(prods) == 0) prods <- NA_character_
      
      # Cartesian expansion (one row per substrate-product pair)
      tidyr::expand_grid(
        rxn_id = rid,
        substrate_id = subs,
        product_id = prods
      )
    }
  )
  
  return(df)
}

mc_reactions <- parse_reactions(
  "~/compound_database_shenlab/metacyc/data/data/reactions.dat"
)

#3. Add kegg_id into reactions
mc_reactions_kegg <- mc_reactions %>%
  left_join(metacyc_nodes, by = c("substrate_id" = "compound_id")) %>%
  dplyr::rename(substrate_kegg = kegg) %>%
  left_join(metacyc_nodes, by = c("product_id" = "compound_id")) %>%
  dplyr::rename(product_kegg = kegg)
#check
mean(!is.na(mc_reactions_kegg$substrate_kegg)) # 0.4052576
mean(!is.na(mc_reactions_kegg$product_kegg)) # 0.4503094

#4. Parse metacyc pathways
parse_reaction_pathways <- function(file) {
  lines <- readLines(file)
  
  rxn <- NULL
  pathways <- c()
  out <- list()
  
  for (line in lines) {
    
    # New reaction frame
    if (startsWith(line, "UNIQUE-ID -")) {
      # save previous reaction's pathway list
      if (!is.null(rxn)) {
        out[[rxn]] <- pathways
      }
      rxn <- sub("UNIQUE-ID - ", "", line)
      pathways <- c()
    }
    
    # pathway links
    if (startsWith(line, "IN-PATHWAY -")) {
      pathways <- c(pathways, sub("IN-PATHWAY - ", "", line))
    }
  }
  
  # save last reaction
  if (!is.null(rxn)) {
    out[[rxn]] <- pathways
  }
  
  tibble(
    rxn_id = rep(names(out), times = sapply(out, length)),
    pathway_id = unlist(out)
  )
}

mc_rxn_pathway <- parse_reaction_pathways(
  "~/compound_database_shenlab/metacyc/data/data/reactions.dat"
)

#5. Combine everything
mc_compound_pathway <- mc_reactions_kegg %>%
  left_join(mc_rxn_pathway, by = "rxn_id") %>%
  tidyr::pivot_longer(
    cols = c(substrate_kegg, product_kegg),
    names_to = "role",
    values_to = "kegg"
  ) %>%
  filter(!is.na(kegg), !is.na(pathway_id)) %>%
  dplyr::select(kegg, pathway_id) %>%
  distinct()
#check
head(mc_compound_pathway)
#kegg                    pathway_id
# <chr>                       <chr>     
#   1 CURLTUGMZLYLDI-UHFFFAOYSA-N PWY-7100  
# 2 XLYOFNOQVPJJNP-UHFFFAOYSA-N PWY-7100  
nrow(mc_compound_pathway) #25656
length(unique(mc_compound_pathway$kegg)) #4246
length(unique(mc_compound_pathway$pathway_id)) #4376
