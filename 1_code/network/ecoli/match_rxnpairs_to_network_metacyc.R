library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')
##Annotation of reaction pairs for ecoli

#1. add in the inchi column for the annotated object
#load database files and ecoli filtered object
load("~/compound_database_shenlab/metacyc/metacyc_ms1_edges_filtered_unique.rds")
load("~/compound_database_shenlab/metacyc/metacyc_ms1_nodes_all.rda")
object_ecoli_sig_annotated <- readRDS("~/MS2Library_ChiFang/In_Vitro/R_objects/object_ecoli_sig_annotated.rds")
#add column "inchi" back to the annotation tabler of the object
load("~/compound_database_shenlab/MS1/ms1_database.rda")
colnames(ms1_database@spectra.info)#  "Lab.ID"  "INCHI_ID"  "INCHI_ID_all"  "INCHIKEY_ID" "INCHIKEY_ID_all"  
annotation_table <- object_ecoli_sig_annotated@annotation_table
ms1_database_df <- ms1_database@spectra.info
inchi_info <- ms1_database_df[, c("INCHI_ID", "Lab.ID")]
annotation_table_with_inchi <- merge(
  annotation_table, 
  inchi_info, 
  by = "Lab.ID", 
  all.x = TRUE
)
object_ecoli_sig_annotated@annotation_table <- annotation_table_with_inchi
object_ecoli_sig_annotated_with_inchi <- object_ecoli_sig_annotated
save(object_ecoli_sig_annotated_with_inchi, file = "~/MS2Library_ChiFang/In_Vitro/R_objects/object_ecoli_sig_annotated_with_inchi.rds")

#2. annotate reaction pairs using inchi
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
df <- as.data.frame(object_ecoli_sig_annotated_with_inchi@variable_info)
df$variable_id <- as.character(df$variable_id)

# Extract annotation table
anno <- as.data.frame(extract_annotation_table(object_ecoli_sig_annotated_with_inchi))
anno$variable_id <- as.character(anno$variable_id)
colnames(anno) # 341 X 19 variables, "INCHI_ID" included

#remove replicates and all rows with missing INCHI_IDs
sum(duplicated(anno$INCHI_ID)) #33
sum(is.na(anno$INCHI_ID)) #32
  #remove rows with missing inchi_ids
anno_clean <- anno %>%
  filter(!is.na(INCHI_ID)) # 309 X 19 VARIABLES
sum(is.na(anno_clean$INCHI_ID)) #0
sum(duplicated(anno_clean$INCHI_ID)) #2
  #collapse duplicated inchi_ids

anno_collapsed <- anno_clean %>%
  group_by(INCHI_ID) %>%
  summarise(
    # Concatenate all variable_ids for this INCHI_ID
    variable_id = paste(unique(variable_id), collapse = "; "),
    
    # Concatenate all compound names
    Compound.name = paste(unique(Compound.name), collapse = "; "),
    
    .groups = "drop"
  )
sum(duplicated(anno_collapsed$INCHI_ID)) #0
#add in column of directions from variable_info by variable_id
    #extract variable_info for regulations
var_info <- as.data.frame(object_ecoli_sig_annotated_with_inchi@variable_info)
var_info$variable_id <- as.character(var_info$variable_id)
    # select variable_id and directions
var_info <- var_info %>%
  dplyr::select(variable_id, direction)
    #merge using collapse
anno_collapsed <- anno_clean %>%
  group_by(INCHI_ID) %>%
  summarise(
    variable_id = paste(unique(variable_id), collapse = "; "),
    Compound.name = paste(unique(Compound.name), collapse = "; "),
    mz = mean(mz, na.rm = TRUE),
    rt = mean(rt, na.rm = TRUE),
    CE = paste(unique(CE), collapse = "; "),
    Adduct = paste(unique(Adduct), collapse = "; "),
    Total.score = max(Total.score, na.rm = TRUE),
    SS = max(SS, na.rm = TRUE),
    Level = paste(unique(Level), collapse = "; "),
    Database = paste(unique(Database), collapse = "; "),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    vars = list(strsplit(variable_id, "; ")[[1]]),
    direction = {
      dir_vals <- object_ecoli_sig_annotated_with_inchi@variable_info$direction[
        match(vars, object_ecoli_sig_annotated_with_inchi@variable_info$variable_id)
      ]
      paste(unique(na.omit(dir_vals)), collapse = "; ")
    },
    fc = {
      fc_vals <- object_ecoli_sig_annotated_with_inchi@variable_info$fc[
        match(vars, object_ecoli_sig_annotated_with_inchi@variable_info$variable_id)
      ]
      paste(unique(na.omit(fc_vals)), collapse = "; ")
    },
    p_value = {
      pvals <- object_ecoli_sig_annotated_with_inchi@variable_info$p_value[
        match(vars, object_ecoli_sig_annotated_with_inchi@variable_info$variable_id)
      ]
      paste(unique(na.omit(pvals)), collapse = "; ")
    },
    p_value_adjust = {
      adjpvals <- object_ecoli_sig_annotated_with_inchi@variable_info$p_value_adjust[
        match(vars, object_ecoli_sig_annotated_with_inchi@variable_info$variable_id)
      ]
      paste(unique(na.omit(adjpvals)), collapse = "; ")
    }
  ) %>%
  ungroup() %>%
  dplyr::select(-vars)

colnames(anno_collapsed) #307 x 15 variables
ecoli_nodes <- anno_collapsed
# Add central node and combine with the nodes
central_node <- tibble(
  variable_id = "Escherichia coli MG-1655 ATCC 7000926",
  label = "Escherichia coli MG-1655 ATCC 7000926",
  direction = "central",
  fc = NA,        # no fold change for central node
  fc_abs = NA     # will be used for scaling
)
ecoli_nodes <- bind_rows(central_node, ecoli_nodes)

#make absolute fc, fc_abs
ecoli_nodes <- ecoli_nodes %>%
  mutate(
    fc_abs = ifelse(direction != "central",
                    abs(as.numeric(fc)),  # convert fc to numeric and take absolute
                    NA)
  )
#Scale node sizes based on absolute fold change
size_range <- c(5, 15)  # adjust to your preference
ecoli_nodes <- ecoli_nodes %>%
  mutate(
    size_scaled = ifelse(direction != "central",
                         scales::rescale(fc_abs, to = size_range),
                         20)   # fixed size for central node
  )
save(ecoli_nodes, file ="~/MS2Library_ChiFang/In_Vitro/R_objects/ecoli_nodes" ) #308 x 18 variables
# -------------------------------
# 2. Create edges
# -------------------------------
load("~/compound_database_shenlab/metacyc/metacyc_edges_filtered_unique")
   #separate nodes into ups and downs
up_nodes <- ecoli_nodes %>% filter(direction == "up")
down_nodes <- ecoli_nodes %>% filter(direction == "down")
  #connect with metacyc
metacyc_edges_clean <- metacyc_edges_filtered_unique %>%
  filter(
    grepl("^InChI=1S/", from_inchi) &
      grepl("^InChI=1S/", to_inchi)
  )
nrow(metacyc_edges_clean) #17538
sum(ecoli_nodes$INCHI_ID %in% metacyc_edges_clean$from_inchi |
      ecoli_nodes$INCHI_ID %in% metacyc_edges_clean$to_inchi) #15
# keep only edges where either side is in your ecoli_nodes list
metacyc_sub_edges <- metacyc_edges_clean %>%
  filter(
    from_inchi %in% ecoli_nodes$INCHI_ID |
      to_inchi %in% ecoli_nodes$INCHI_ID
  )
nrow(metacyc_sub_edges) #121
# collect unique InChIs involved in those reactions
sub_inchis <- unique(c(metacyc_sub_edges$from_inchi, metacyc_sub_edges$to_inchi))

length(sub_inchis) #95
# Join with your ecoli_nodes to add readable names
reaction_pairs_named <- metacyc_sub_edges %>%
  left_join(
    ecoli_nodes %>%
      dplyr::select(INCHI_ID, Compound.name) %>%
      dplyr::rename(from_name = Compound.name),
    by = c("from_inchi" = "INCHI_ID")
  ) %>%
  left_join(
    ecoli_nodes %>%
      dplyr::select(INCHI_ID, Compound.name) %>%
      dplyr::rename(to_name = Compound.name),
    by = c("to_inchi" = "INCHI_ID")
  ) %>%
  mutate(
    from_label = ifelse(!is.na(from_name), from_name, from_inchi),
    to_label   = ifelse(!is.na(to_name), to_name, to_inchi)
  ) %>%
  dplyr::select(from_label, to_label)


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