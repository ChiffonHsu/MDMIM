load("2_data/preliminary_determination_bacterial_metabolomics/MS2/object_annot_best")
load("2_data/preliminary_determination_bacterial_metabolomics/MS2/differential_analysis_each_bacterium/collapsed_results")
#1. Create simple and cleaned tibble for each individual bacterium
network_input <- lapply(collapsed_results, function(df){
  df %>%
    mutate(
      regulation = case_when(
        pval_final < 0.05 & fc_final > 1 ~ "up",
        pval_final < 0.05 & fc_final < -1 ~ "down",
        TRUE ~ "none"
      )
    ) %>%
    filter(regulation != "none") %>%
    dplyr::select(class, Compound.name, regulation, fc_final, pval_final)
})
names(network_input) <- bacteria_classes

#2. build network edges from downregulated to upregulated
  #A. ADD in inchi_key and KEGG, CAS.ID
#add column "inchi" back to the annotation tabler of the object
load("~/compound_database_shenlab/MS1/ms1_database.rda")
colnames(ms1_database@spectra.info)#  "Lab.ID"  "INCHI_ID"  "INCHI_ID_all"  "INCHIKEY_ID" "INCHIKEY_ID_all"  
annotation_table <- object_annot_best@annotation_table
ms1_database_df <- ms1_database@spectra.info
inchi_info <- ms1_database_df[, c("INCHIKEY_ID", "Lab.ID")]
object_annot_best_inchi <- merge(
  annotation_table, 
  inchi_info, 
  by = "Lab.ID", 
  all.x = TRUE
)
annotation_lookup <- object_annot_best_inchi %>%
  dplyr::select(
    variable_id,
    KEGG.ID,
    HMDB.ID,
    CAS.ID,
    INCHIKEY_ID
  )

network_input <- lapply(collapsed_results, function(df){
  df %>%
    mutate(
      regulation = case_when(
        pval_final < 0.05 & fc_final > 1 ~ "up",
        pval_final < 0.05 & fc_final < -1 ~ "down",
        TRUE ~ "none"
      )
    ) %>%
    filter(regulation != "none") %>%
    dplyr::select(class, Compound.name, regulation, fc_final, pval_final, KEGG.ID, HMDB.ID, CAS.ID, INCHIKEY_ID)
})
names(network_input) <- bacteria_classes

  #B. Prepare edges
network_nodes <- lapply(network_input, function(df){
    df %>%
      mutate(
        direction = case_when(
          regulation == "down" ~ "from",
          regulation == "up"   ~ "to",
          TRUE ~ NA_character_
        )
      )
  })
  
  #C. Check for pairing with metacyc, annotate using kegg.id
# per bacterium
metacyc_edges_all <- readRDS("~/compound_database_shenlab/metacyc/metacyc_edges_all.rds")
metacyc_nodes_all <- readRDS("~/compound_database_shenlab/metacyc/metacyc_nodes_all.rds")

metacyc_pairs_kegg <- lapply(names(network_nodes), function(bact){
  
  df <- collapsed_results[[bact]]
  
  kegg_ids <- df$KEGG.ID %>% na.omit() %>% unique()
  
  # reaction pairs where either substrate OR product matches
  matched <- metacyc_edges_all %>%
    filter(
      from_kegg %in% kegg_ids |
        to_kegg   %in% kegg_ids
    )
  
  matched
})

names(metacyc_pairs_all) <- names(metacyc_pairs_kegg)



library(tidymass)
library(dplyr)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)

# -------------------------------
# 1. Prepare nodes
# -------------------------------

# Extract info
df <- as.data.frame(comparison_results$`Enterococcus faecium`)

# Make sure info are correct
df <- df %>%
  mutate(
    variable_id = as.character(variable_id),
    Compound.name = ifelse(is.na(Compound.name) | Compound.name == "",
                           variable_id, Compound.name)
  )

# build nodes
nodes <- df %>%
  group_by(variable_id) %>%
  summarise(
    label     = paste(unique(Compound.name), collapse = ";"),
    mean_fc   = mean(fc_final, na.rm = TRUE),
    direction = case_when(
      mean_fc >  4.644 ~ "up",    # tweak thresholds to taste
      mean_fc < -4.644 ~ "down",
      TRUE        ~ "neutral"
    ),
    fc_abs = abs(mean_fc),
    .groups = "drop"
  ) %>%
  filter(direction != "neutral")   # keep only up/down

# Add central node
central_name <- "Enterococcus faecium"
central_node <- tibble(
  variable_id = central_name,
  label       = central_name,
  direction   = "central",
  fc_abs      = NA_real_
)
# Combine nodes
nodes_all <- bind_rows(central_node, nodes) %>%
  mutate(
    size_scaled = if_else(direction == "central",
                          20,
                          scales::rescale(fc_abs, to = c(5, 15))),
    # igraph expects a "name" column for vertex IDs
    name = variable_id,
    # create a 'size' column to use directly in igraph/ggraph
    size = size_scaled
  )
nodes_all <- nodes_all %>%
  mutate(
    fontface = if_else(direction == "central", "bold", "plain")
  )

# --------------------------------------------------------
# Pick top 10 up and top 10 down regulated nodes for labels
# --------------------------------------------------------
top_up <- nodes_all %>%
  filter(direction == "up") %>%
  arrange(desc(fc_abs)) %>%
  slice_head(n = 10) %>%
  pull(variable_id)

top_down <- nodes_all %>%
  filter(direction == "down") %>%
  arrange(desc(fc_abs)) %>%
  slice_head(n = 10) %>%
  pull(variable_id)

# Label only these nodes + the central node
nodes_all <- nodes_all %>%
  mutate(
    label_to_plot = case_when(
      variable_id == central_name      ~ label,     # always keep bacterium label
      variable_id %in% top_up          ~ label,     # top 10 up
      variable_id %in% top_down        ~ label,     # top 10 down
      TRUE                              ~ ""         # hide all others
    )
  )

# Scale node sizes
size_range <- c(5, 15)
nodes <- nodes %>%
  mutate(size_scaled = ifelse(direction != "central",
                              scales::rescale(fc_abs, to = size_range),
                              20))

# -------------------------------
# 2. Create edges
# -------------------------------
edges <- nodes_all %>%
  filter(direction != "central") %>%
  transmute(from = variable_id, to = central_name)

# -------------------------------
# 3. Build igraph object
# -------------------------------
g <- graph_from_data_frame(edges, vertices = nodes_all, directed = TRUE)
layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)

V(g)$label_to_plot <- nodes_all$label_to_plot
V(g)$fontface <- nodes_all$fontface


central_idx <- which(V(g)$direction == "central")
up_idx      <- which(V(g)$direction == "up")
down_idx    <- which(V(g)$direction == "down")

compute_angles <- function(idx_vec, start_angle, end_angle) {
  n <- length(idx_vec)
  if (n == 1) return((start_angle + end_angle) / 2)
  seq(start_angle, end_angle, length.out = n + 1)[-1]
}

# Up-regulated nodes: top semicircle
if (length(up_idx) > 0) {
  angles_up <- compute_angles(up_idx, 0, pi)
  radii_up  <- 1.5 * (1 - scales::rescale(V(g)$size[up_idx], to = c(0, 1))) + 0.5
  layout_matrix[up_idx, ] <- cbind(radii_up * cos(angles_up),
                                   radii_up * sin(angles_up))
}

# Down-regulated nodes: bottom semicircle
if (length(down_idx) > 0) {
  angles_down <- compute_angles(down_idx, pi, 2 * pi)
  radii_down  <- 1.5 * (1 - scales::rescale(V(g)$size[down_idx], to = c(0, 1))) + 0.5
  layout_matrix[down_idx, ] <- cbind(radii_down * cos(angles_down),
                                     radii_down * sin(angles_down))
}

# Central node at origin
layout_matrix[central_idx, ] <- c(0, 0)

g_layout <- create_layout(
  g,
  layout = "manual",
  x = layout_matrix[, 1],
  y = layout_matrix[, 2]
)

#plot
ggraph(g_layout) +
  geom_edge_link(arrow = arrow(length = unit(1, "mm")), alpha = 0.6) +
  geom_node_point(aes(size = size, color = direction)) +
  geom_node_text(
    aes(label = label_to_plot, fontface = fontface),
    repel = TRUE,
    size = 4
  )+
  scale_color_manual(values = c(up = "red",
                                down = "blue",
                                central = "yellow")) +
  theme_void()


png("D:/NTU/Research/Projects/Microbiome_driven_metabolite_interaction_mapping/MDMIM/3_data_analysis/rplc_pos/ecoli/annotatedecoli_starnetwork_fc2.png",
    width = 6000, height = 6000, res = 300)

ggraph(g_layout) +
  geom_edge_link(arrow = arrow(length = unit(4, "mm")), alpha = 0.6) +
  geom_node_point(aes(size = size, color = direction)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 4) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "central" = "yellow")) +
  scale_size_continuous(range = c(5, 15), name = "|FC|") +
  theme_void()

dev.off()