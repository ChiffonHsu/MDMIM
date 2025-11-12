.libPaths("D:/APPS/R sTUDIO/r-2025/Rlibs")
library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

#annotation
object_ecoli_sig <- readRDS("D:/NTU/Research/Projects/Microbiome_driven_metabolite_interaction_mapping/MDMIM/2_data/differential_analysis/rplc_pos/ecoli/object_ecoli_sig.rds")
library(metid)
library(tidymass                                                                                                                                                                                                                                                                                                                                                                                                                         )
library(readr)

  load("~/MS2Library_ChiFang/ms1_database.rda")

object_ecoli_sig_rplcneg_annot <-
  metid::annotate_metabolites_mass_dataset(
    object = object_ecoli_sig_rplcneg,
    database = ms1_database,
    polarity = "negative",
    column = "rp",
    #adduct.table = NULL
  )
saveRDS(object_ecoli_sig_rplcneg, file = "~/MS2Library_ChiFang/In_Vitro/R_objects/object_ecoli_sig_rplcneg_annot.rds")


object_ecoli_hpos_annot <-
  metid::annotate_metabolites_mass_dataset(
    object = object_ecoli_hpos_annot,
    database = ms1_database,
    polarity = "positive",
    column = "hilic",
    #adduct.table = NULL
  )
saveRDS(object_ecoli_sig_rplcneg, file = "~/MS2Library_ChiFang/In_Vitro/R_objects/object_ecoli_sig_hpos_annot.rds")

#merge annotated objects together
load("~/MS2Library_ChiFang/In_Vitro/RPLC/POS/annotation/object_ecoli_sig_annotated.rda")
head(colnames(object_ecoli_hpos_annot))
head(colnames(object_ecoli_sig_rplcneg_annot))
object_ecoli_sig_rplc_annot_merged <- 
  merge_mass_dataset(x = object_ecoli_hpos_annot, 
                     y = object_ecoli_sig_rplcneg_annot, 
                     sample_direction = "inner",
                     variable_direction = "full", 
                     sample_by = "sample_id", 
                     variable_by = c("variable_id", "mz", "rt"))

head(colnames(object_ecoli_sig_annotated))
head(colnames(object_ecoli_sig_rplc_annot_merged))
object_ecoli_sig_annot_allnerged <- 
  merge_mass_dataset(x = object_ecoli_sig_annotated, 
                     y = object_ecoli_sig_rplc_annot_merged, 
                     sample_direction = "inner",
                     variable_direction = "full", 
                     sample_by = "sample_id", 
                     variable_by = c("variable_id", "mz", "rt"))



  #construct volcano map for annotated metabolites
library(tidymass)
library(ggplot2)

p <- volcano_plot(
  object = object_ecoli_sig_rplcneg_annot,
  p_value_cutoff = 0.05,
  add_text = TRUE,
  text_from = "Compound.name", 
  point_size_scale = "p_value"
) +
  scale_size_continuous(range = c(0.5, 10))

   # Save as png
ggsave(
  filename = "D:/NTU/Research/Projects/Microbiome_driven_metabolite_interaction_mapping/MDMIM/3_data_analysis/rplc_pos/ecoli/annotatedecoli_network_volc.png",
  plot = p,
  device = "png",
  width = 20,  # in inches
  height = 20
)


#Radial arrangement
library(tidymass)
library(dplyr)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)

# -------------------------------
# 1. Prepare nodes
# -------------------------------

# Extract variable info
df <- as.data.frame(object_ecoli_sig_annotated@variable_info)
df$variable_id <- as.character(df$variable_id)

# Extract annotation table
anno <- as.data.frame(extract_annotation_table(object_ecoli_sig_annotated))
anno$variable_id <- as.character(anno$variable_id)

# Merge compound names
df <- df %>%
  left_join(anno[, c("variable_id", "Compound.name")], by = "variable_id") %>%
  mutate(Compound.name = ifelse(is.na(Compound.name) | Compound.name == "", variable_id, Compound.name))

# Collapse nodes by variable_id, filter strong regulation only
nodes_collapsed <- df %>%
  group_by(variable_id) %>%
  summarize(
    label = paste(unique(Compound.name), collapse = ";"),
    mean_fc = mean(fc, na.rm = TRUE),
    direction = case_when(
      mean_fc > 2 ~ "up",
      mean_fc < 0.5 ~ "down",
      TRUE ~ "neutral"
    ),
    fc_abs = abs(mean_fc),
    .groups = "drop"
  ) %>%
  filter(direction != "neutral")   # keep only strong regulation

# Add central node
central_node <- tibble(
  variable_id = "Escherichia coli MG-1655 ATCC 7000926",
  label = "Escherichia coli MG-1655 ATCC 7000926",
  direction = "central",
  fc_abs = NA
)

# Combine nodes
nodes <- bind_rows(central_node, nodes_collapsed)

# Scale node sizes
size_range <- c(5, 15)
nodes <- nodes %>%
  mutate(size_scaled = ifelse(direction != "central",
                              scales::rescale(fc_abs, to = size_range),
                              20))

# -------------------------------
# 2. Create edges
# -------------------------------
edges <- nodes %>%
  filter(direction != "central") %>%
  transmute(from = variable_id, to = central_node$variable_id)

# -------------------------------
# 3. Build igraph object
# -------------------------------
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# -------------------------------
# 4. Radial layout
# -------------------------------
vcount_g <- vcount(g)
layout_matrix <- matrix(0, nrow = vcount_g, ncol = 2)
central_idx <- which(V(g)$direction == "central")
up_idx <- which(V(g)$direction == "up")
down_idx <- which(V(g)$direction == "down")

# Assign node sizes and attributes
V(g)$size <- nodes$size_scaled
V(g)$label <- nodes$label
V(g)$direction <- nodes$direction

# Function to compute angles evenly along arc
compute_angles <- function(idx_vec, start_angle, end_angle){
  n <- length(idx_vec)
  if(n == 1) return((start_angle + end_angle)/2)
  seq(start_angle, end_angle, length.out = n + 1)[-1]
}

# Up nodes: top half
if(length(up_idx) > 0){
  angles_up <- compute_angles(up_idx, 0, pi)
  radii_up <- 1.5 * (1 - (V(g)$size[up_idx] - min(V(g)$size[up_idx])) /
                       (max(V(g)$size[up_idx]) - min(V(g)$size[up_idx]) + 1e-6)) + 0.5
  layout_matrix[up_idx, 1] <- radii_up * cos(angles_up)
  layout_matrix[up_idx, 2] <- radii_up * sin(angles_up)
}

# Down nodes: bottom half
if(length(down_idx) > 0){
  angles_down <- compute_angles(down_idx, pi, 2*pi)
  radii_down <- 1.5 * (1 - (V(g)$size[down_idx] - min(V(g)$size[down_idx])) /
                         (max(V(g)$size[down_idx]) - min(V(g)$size[down_idx]) + 1e-6)) + 0.5
  layout_matrix[down_idx, 1] <- radii_down * cos(angles_down)
  layout_matrix[down_idx, 2] <- radii_down * sin(angles_down)
}

# Central node at origin
layout_matrix[central_idx, ] <- c(0, 0)

# Create layout object for ggraph
g_layout <- create_layout(g, layout = "manual",
                          x = layout_matrix[,1], y = layout_matrix[,2])

# -------------------------------
# 5. Plot and save
# -------------------------------
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

############for fc>5, <1/5
# ===========================================
# Radial star network for E. coli metabolites
# Filter: fc > 5 or fc < 1/5
# Symmetry: fc and 1/fc have same distance & node size
# ===========================================

library(tidymass)
library(dplyr)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)

# -------------------------------
# 1. Prepare nodes
# -------------------------------

# Extract variable info
df <- as.data.frame(object_ecoli_sig_annotated@variable_info)
df$variable_id <- as.character(df$variable_id)

# Extract annotation table
anno <- as.data.frame(extract_annotation_table(object_ecoli_sig_annotated))
anno$variable_id <- as.character(anno$variable_id)

# Merge compound names
df <- df %>%
  left_join(anno[, c("variable_id", "Compound.name")], by = "variable_id") %>%
  mutate(Compound.name = ifelse(is.na(Compound.name) | Compound.name == "", variable_id, Compound.name))

# Collapse nodes by variable_id, compute log2 FC
nodes_collapsed <- df %>%
  group_by(variable_id) %>%
  summarize(
    label = paste(unique(Compound.name), collapse = ";"),
    mean_fc = mean(fc, na.rm = TRUE),
    log2_fc = mean(log2(fc), na.rm = TRUE),   # log2 transform
    direction = case_when(
      mean_fc > 5 ~ "up",
      mean_fc < 0.2 ~ "down",
      TRUE ~ "neutral"
    ),
    fc_abs = abs(log2_fc),   # symmetric magnitude
    .groups = "drop"
  ) %>%
  filter(direction != "neutral")

# Add central node
central_node <- tibble(
  variable_id = "Escherichia coli MG-1655 ATCC 7000926",
  label = "Escherichia coli MG-1655 ATCC 7000926",
  direction = "central",
  fc_abs = NA
)

nodes <- bind_rows(central_node, nodes_collapsed)

# Scale node sizes (based on |log2FC|)
size_range <- c(5, 15)
nodes <- nodes %>%
  mutate(size_scaled = ifelse(direction != "central",
                              scales::rescale(fc_abs, to = size_range),
                              20))

# -------------------------------
# 2. Create edges
# -------------------------------
edges <- nodes %>%
  filter(direction != "central") %>%
  transmute(from = variable_id, to = central_node$variable_id)

# -------------------------------
# 3. Build igraph object
# -------------------------------
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# Assign attributes to vertices
V(g)$label <- nodes$label
V(g)$direction <- nodes$direction
V(g)$fc_abs <- nodes$fc_abs
V(g)$size <- nodes$size_scaled

# -------------------------------
# 4. Radial layout
# -------------------------------
vcount_g <- vcount(g)
layout_matrix <- matrix(0, nrow = vcount_g, ncol = 2)

central_idx <- which(V(g)$direction == "central")
up_idx <- which(V(g)$direction == "up")
down_idx <- which(V(g)$direction == "down")

# Function to compute evenly spaced angles
compute_angles <- function(idx_vec, start_angle, end_angle){
  n <- length(idx_vec)
  if(n == 1) return((start_angle + end_angle)/2)
  seq(start_angle, end_angle, length.out = n + 1)[-1]
}

# Up nodes: top half
if(length(up_idx) > 0){
  angles_up <- compute_angles(up_idx, 0, pi)
  radii_up <- 1.5 * (1 - (V(g)$fc_abs[up_idx] - min(V(g)$fc_abs[up_idx])) /
                       (max(V(g)$fc_abs[up_idx]) - min(V(g)$fc_abs[up_idx]) + 1e-6)) + 0.5
  layout_matrix[up_idx, 1] <- radii_up * cos(angles_up)
  layout_matrix[up_idx, 2] <- radii_up * sin(angles_up)
}

# Down nodes: bottom half
if(length(down_idx) > 0){
  angles_down <- compute_angles(down_idx, pi, 2*pi)
  radii_down <- 1.5 * (1 - (V(g)$fc_abs[down_idx] - min(V(g)$fc_abs[down_idx])) /
                         (max(V(g)$fc_abs[down_idx]) - min(V(g)$fc_abs[down_idx]) + 1e-6)) + 0.5
  layout_matrix[down_idx, 1] <- radii_down * cos(angles_down)
  layout_matrix[down_idx, 2] <- radii_down * sin(angles_down)
}

# Central node at origin
layout_matrix[central_idx, ] <- c(0, 0)

# Create layout for ggraph
g_layout <- create_layout(g, layout = "manual",
                          x = layout_matrix[,1], y = layout_matrix[,2])

# -------------------------------
# 5. Plot & Save
# -------------------------------
png("D:/NTU/Research/Projects/Microbiome_driven_metabolite_interaction_mapping/MDMIM/3_data_analysis/rplc_pos/ecoli/annotatedecoli_starnetwork_fc5-2.png",
    width = 6000, height = 6000, res = 400)

ggraph(g_layout) +
  geom_edge_link(arrow = arrow(length = unit(4, "mm")), alpha = 0.6) +
  geom_node_point(aes(size = size, color = direction)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "central" = "yellow")) +
  scale_size_continuous(range = c(5, 15), name = "|log2FC|") +
  theme_void()

dev.off()


# ===========================================
# Radial star network for E. coli metabolites
# Filter: fc > 2 or fc < 1/2
# Symmetry: fc and 1/fc have same distance & node size
# ===========================================

library(tidymass)
library(dplyr)
library(igraph)
library(ggraph)
library(ggrepel)
library(scales)

# -------------------------------
# 1. Prepare nodes
# -------------------------------

# Extract variable info
df <- as.data.frame(object_ecoli_sig_annotated@variable_info)
df$variable_id <- as.character(df$variable_id)

# Extract annotation table
anno <- as.data.frame(extract_annotation_table(object_ecoli_sig_annotated))
anno$variable_id <- as.character(anno$variable_id)

# Merge compound names
df <- df %>%
  left_join(anno[, c("variable_id", "Compound.name")], by = "variable_id") %>%
  mutate(Compound.name = ifelse(is.na(Compound.name) | Compound.name == "", variable_id, Compound.name))

# Collapse nodes by variable_id, compute log2 FC
nodes_collapsed <- df %>%
  group_by(variable_id) %>%
  summarize(
    label = paste(unique(Compound.name), collapse = ";"),
    mean_fc = mean(fc, na.rm = TRUE),
    log2_fc = mean(log2(fc), na.rm = TRUE),   # log2 transform
    direction = case_when(
      mean_fc > 2 ~ "up",
      mean_fc < 0.5 ~ "down",
      TRUE ~ "neutral"
    ),
    fc_abs = abs(log2_fc),   # symmetric magnitude
    .groups = "drop"
  ) %>%
  filter(direction != "neutral")

# Add central node
central_node <- tibble(
  variable_id = "Escherichia coli MG-1655 ATCC 7000926",
  label = "Escherichia coli MG-1655 ATCC 7000926",
  direction = "central",
  fc_abs = NA
)

nodes <- bind_rows(central_node, nodes_collapsed)

# Scale node sizes (based on |log2FC|)
size_range <- c(5, 15)
nodes <- nodes %>%
  mutate(size_scaled = ifelse(direction != "central",
                              scales::rescale(fc_abs, to = size_range),
                              20))

# -------------------------------
# 2. Create edges
# -------------------------------
edges <- nodes %>%
  filter(direction != "central") %>%
  transmute(from = variable_id, to = central_node$variable_id)

# -------------------------------
# 3. Build igraph object
# -------------------------------
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# Assign attributes to vertices
V(g)$label <- nodes$label
V(g)$direction <- nodes$direction
V(g)$fc_abs <- nodes$fc_abs
V(g)$size <- nodes$size_scaled

# -------------------------------
# 4. Radial layout
# -------------------------------
vcount_g <- vcount(g)
layout_matrix <- matrix(0, nrow = vcount_g, ncol = 2)

central_idx <- which(V(g)$direction == "central")
up_idx <- which(V(g)$direction == "up")
down_idx <- which(V(g)$direction == "down")

# Function to compute evenly spaced angles
compute_angles <- function(idx_vec, start_angle, end_angle){
  n <- length(idx_vec)
  if(n == 1) return((start_angle + end_angle)/2)
  seq(start_angle, end_angle, length.out = n + 1)[-1]
}

# Up nodes: top half
if(length(up_idx) > 0){
  angles_up <- compute_angles(up_idx, 0, pi)
  radii_up <- 1.5 * (1 - (V(g)$fc_abs[up_idx] - min(V(g)$fc_abs[up_idx])) /
                       (max(V(g)$fc_abs[up_idx]) - min(V(g)$fc_abs[up_idx]) + 1e-6)) + 0.5
  layout_matrix[up_idx, 1] <- radii_up * cos(angles_up)
  layout_matrix[up_idx, 2] <- radii_up * sin(angles_up)
}

# Down nodes: bottom half
if(length(down_idx) > 0){
  angles_down <- compute_angles(down_idx, pi, 2*pi)
  radii_down <- 1.5 * (1 - (V(g)$fc_abs[down_idx] - min(V(g)$fc_abs[down_idx])) /
                         (max(V(g)$fc_abs[down_idx]) - min(V(g)$fc_abs[down_idx]) + 1e-6)) + 0.5
  layout_matrix[down_idx, 1] <- radii_down * cos(angles_down)
  layout_matrix[down_idx, 2] <- radii_down * sin(angles_down)
}

# Central node at origin
layout_matrix[central_idx, ] <- c(0, 0)

# Create layout for ggraph
g_layout <- create_layout(g, layout = "manual",
                          x = layout_matrix[,1], y = layout_matrix[,2])

# -------------------------------
# 5. Plot & Save
# -------------------------------
png("D:/NTU/Research/Projects/Microbiome_driven_metabolite_interaction_mapping/MDMIM/3_data_analysis/rplc_pos/ecoli/annotatedecoli_starnetwork_fc2.png",
    width = 6000, height = 6000, res = 400)

ggraph(g_layout) +
  geom_edge_link(arrow = arrow(length = unit(4, "mm")), alpha = 0.6) +
  geom_node_point(aes(size = size, color = direction)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "central" = "yellow")) +
  scale_size_continuous(range = c(5, 15), name = "|log2FC|") +
  theme_void()

dev.off()
