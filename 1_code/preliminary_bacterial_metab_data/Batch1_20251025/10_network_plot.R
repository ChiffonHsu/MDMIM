#Generate network plots for all bacteria
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/hybrid_collapsed_results")

plot_bacteria_network <- function(df, central_name, output_file = NULL) {
  
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(scales)
  library(ggrepel)
  
  # ---------------------------------
  # 1. Prepare nodes
  # ---------------------------------
  df <- df %>%
    mutate(
      variable_id   = as.character(variable_id),
      Compound.name = ifelse(is.na(Compound.name) | Compound.name == "",
                             variable_id, Compound.name)
    )
  
  nodes <- df %>%
    group_by(variable_id) %>%
    summarise(
      label     = paste(unique(Compound.name), collapse = ";"),
      mean_fc   = mean(fc_final, na.rm = TRUE),
      mean_pval = mean(pval_final, na.rm = TRUE),   # <-- FIXED COMMA !!!
      
      direction = case_when(
        mean_fc >  2.322 & mean_pval < 0.05 ~ "up",
        mean_fc < -2.322 & mean_pval < 0.05 ~ "down",
        TRUE ~ "neutral"
      ),
      
      fc_abs = abs(mean_fc),
      .groups = "drop"
    ) %>%
    filter(direction != "neutral")
  
  # central node
  central_node <- tibble(
    variable_id = central_name,
    label       = central_name,
    direction   = "central",
    fc_abs      = NA_real_
  )
  
  nodes_all <- bind_rows(central_node, nodes) %>%
    mutate(
      size_scaled = if_else(direction == "central",
                            20,
                            scales::rescale(fc_abs, to = c(5, 15))),
      name  = variable_id,
      size  = size_scaled,
      fontface = if_else(direction == "central", "bold", "plain")
    )
  
  # -------------------------
  # top 10 up/down labels
  # -------------------------
  top_up <- nodes %>% 
    filter(direction == "up") %>% 
    arrange(desc(fc_abs)) %>% 
    slice_head(n = 10) %>% 
    pull(variable_id)
  
  top_down <- nodes %>% 
    filter(direction == "down") %>% 
    arrange(desc(fc_abs)) %>% 
    slice_head(n = 10) %>% 
    pull(variable_id)
  
  nodes_all <- nodes_all %>%
    mutate(
      label_to_plot = case_when(
        variable_id == central_name ~ label,
        variable_id %in% top_up ~ label,
        variable_id %in% top_down ~ label,
        TRUE ~ ""
      )
    )
  
  # -------------------------
  # Edges
  # -------------------------
  edges <- nodes_all %>%
    filter(direction != "central") %>%
    transmute(from = variable_id, to = central_name)
  
  # -------------------------
  # Graph
  # -------------------------
  g <- graph_from_data_frame(edges, vertices = nodes_all, directed = TRUE)
  
  layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
  
  V(g)$label_to_plot <- nodes_all$label_to_plot
  V(g)$fontface      <- nodes_all$fontface
  
  central_idx <- which(V(g)$direction == "central")
  up_idx      <- which(V(g)$direction == "up")
  down_idx    <- which(V(g)$direction == "down")
  
  compute_angles <- function(idx_vec, start_angle, end_angle) {
    n <- length(idx_vec)
    if (n == 1) return((start_angle + end_angle) / 2)
    seq(start_angle, end_angle, length.out = n + 1)[-1]
  }
  
  if (length(up_idx) > 0) {
    angles <- compute_angles(up_idx, 0, pi)
    radii  <- 1.5 * (1 - scales::rescale(V(g)$size[up_idx], to = c(0, 1))) + 0.5
    layout_matrix[up_idx, ] <- cbind(radii * cos(angles), radii * sin(angles))
  }
  
  if (length(down_idx) > 0) {
    angles <- compute_angles(down_idx, pi, 2*pi)
    radii  <- 1.5 * (1 - scales::rescale(V(g)$size[down_idx], to = c(0, 1))) + 0.5
    layout_matrix[down_idx, ] <- cbind(radii * cos(angles), radii * sin(angles))
  }
  
  layout_matrix[central_idx, ] <- c(0, 0)
  
  g_layout <- create_layout(
    g,
    layout = "manual",
    x = layout_matrix[, 1],
    y = layout_matrix[, 2]
  )
  
  # -------------------------
  # PLOT
  # -------------------------
  p <- ggraph(g_layout) +
    geom_edge_link(arrow = arrow(length = unit(1, "mm")), alpha = 0.6) +
    geom_node_point(aes(size = size, color = direction)) +
    geom_node_text(
      aes(label = label_to_plot, fontface = fontface),
      repel = TRUE,
      max.overlaps = Inf,
      force = 2,
      point.padding = 0.3,
      box.padding = 0.5,
      size = 4
    ) +
    scale_color_manual(values = c(
      up = "red",
      down = "blue",
      central = "yellow"
    )) +
    theme_void()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 10, height = 10, dpi = 300)
  }
  
  return(p)
}



# -------------------------------
# LOOP OVER ALL BACTERIA
# -------------------------------
output_dir <- "~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/rplc_pos/Single_bacteria_map"
dir.create(output_dir, showWarnings = FALSE)

for (bac in names(collapsed_results)) {
  
  message("Processing: ", bac)
  
  df <- collapsed_results[[bac]]
  
  outfile <- file.path(
    output_dir,
    paste0(gsub(" ", "_", bac), "_network.png")
  )
  
  plot_bacteria_network(
    df = df,
    central_name = bac,
    output_file = outfile
  )
}
