#Build a bacteria-metabolite interaction network (bipartite network)
# Step One: Build a combined table of all bacteria and all significant metabolites
#A. Combine the entire tibble object “collapsd_results”
library(purrr)
library(dplyr)
combined <- imap_dfr(
  collapsed_results, 
  ~mutate(.x,  bacteria = .y)
)
#B. Filter out and keep only significant metabolties
sig_combined <- combined %>%
  mutate(variable_id = as.character(variable_id),
         Compound.name = ifelse(is.na(Compound.name) | Compound.name == "", 
                                variable_id, Compound.name)) %>%
  filter(abs(fc_final) > 2.322,
         pval_final < 0.05)
#C. Select only the top 20 metabolites for each bacteria
top20_metabs <- sig_combined %>%
  arrange(desc(abs(fc_final))) %>%
  distinct(variable_id, .keep_all = TRUE) %>%
  slice_head(n = 40) %>%
  pull(variable_id)

sig_top <- sig_combined %>%
  filter(variable_id %in% top20_metabs)

#Step Two: Build the bipartite edge list (bacteria->metabilite)
 #include directions so that up and down have different coloured edges
edges_bip <- sig_top %>%
  transmute(
    from = bacteria,
    to   = paste0("M:", variable_id),
    fc   = fc_final
  ) %>% 
  mutate(direction = ifelse(fc > 0, "up", "down"))

#Step Three: Build nodes list with node types, two nodes types -> nodes_bacteria and nodes_metabolite
# A. Bacteria nodes (type = "bacteria")
nodes_bateria <- tibble(
  name = unique(sig_combined$bacteria),
  type = "bacteria",
  label = name
)
# B. Metabolite nodes (type = "metabolite")
nodes_metabolite <- sig_combined %>%
  distinct(variable_id, Compound.name) %>%
  transmute(
    name = paste0("M:", variable_id),
    type = "metabolite",
    label = Compound.name
  )
# C. Combine bothe nodes_bacteria and nodes_metabolite
nodes_all <- bind_rows(nodes_bateria, nodes_metabolite)
# D. Sizing according to absolute fc
met_size <- sig_combined %>%
  mutate(node = paste0("M:", variable_id)) %>%
  group_by(node) %>%
  summarise(mean_fc_abs = mean(abs(fc_final), na.rm = TRUE))

#Step Four: Build igraph bipartite object
library(igraph)
nodes_all <- nodes_all %>%
  filter(name %in% c(edges_bip$from, edges_bip$to))
g <- graph_from_data_frame(
  d = edges_bip,
  vertices = nodes_all,
  directed = FALSE
)

 #for labeling nodes for the key
# Keep type as character
V(g)$type <- nodes_all$type

 # for sizing
V(g)$size <- 6  
met_idx <- which(V(g)$type == "metabolite")
V(g)$size[met_idx] <- scales::rescale(
  met_size$mean_fc_abs[ match(V(g)$name[met_idx], met_size$node) ],
  to = c(3, 6)
)

#Step Five: Plot bipartite plot using Fruchterman–Reingold
library(ggraph)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_color = fc), alpha = 0.3) +
  scale_edge_color_gradient2(
    low = "blue",
    high = "red",
    mid = "grey80",
    midpoint = 0) +
  geom_node_point(aes(color = type, size = size)) +
  
  #geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  scale_color_manual(values = c(
    "bacteria" = "orange",
    "metabolite" = "skyblue"
  )) +
   
  theme_void()

#Build interaction plot (will show labels when hover over to the edge or nodes)
#Color of edges according to fc 
library(igraph)
library(plotly)
library(dplyr)

set.seed(123)
layout_fr <- layout_with_fr(g)

layout_df <- as.data.frame(layout_fr)
colnames(layout_df) <- c("x", "y")

layout_df$name  <- V(g)$name
layout_df$type  <- V(g)$type
layout_df$size  <- V(g)$size
layout_df$label <- V(g)$label

node_colors <- ifelse(layout_df$type == "bacteria", "orange", "skyblue")
edges_df <- igraph::as_data_frame(g, "edges") %>%
  mutate(
    x0 = layout_df$x[match(from, layout_df$name)],
    y0 = layout_df$y[match(from, layout_df$name)],
    x1 = layout_df$x[match(to, layout_df$name)],
    y1 = layout_df$y[match(to, layout_df$name)],
    hover = paste0("<b>Bacteria:</b> ", from,
                   "<br><b>Metabolite:</b> ", to,
                   "<br><b>FC:</b> ", round(fc, 3))
  )

p <- plot_ly()

# Add each edge separately so gradient works
for (i in seq_len(nrow(edges_df))) {
  
  fc_i <- edges_df$fc[i]
  
  # map FC to color manually
  if (fc_i < 0) {
    col_i <- colorRampPalette(c("blue", "grey80"))(100)[
      cut(scales::rescale(fc_i, to=c(1, -1)), 100)
    ]
  } else {
    col_i <- colorRampPalette(c("grey80", "red"))(100)[
      cut(scales::rescale(fc_i, to=c(-1, 1)), 100)
    ]
  }
  
  p <- p %>% add_trace(
    x = c(edges_df$x0[i], edges_df$x1[i]),
    y = c(edges_df$y0[i], edges_df$y1[i]),
    type = "scatter",
    mode = "lines",
    line = list(color = col_i, width = 0.4),
    text = edges_df$hover[i],
    hoverinfo = "text",
    showlegend = FALSE
  )
}
p <- p %>% add_trace(
  x = layout_df$x,
  y = layout_df$y,
  type = "scatter",
  mode = "markers",
  text = layout_df$label,
  hoverinfo = "text",
  marker = list(
    size = layout_df$size * 1.8,
    color = node_colors,
    line = list(width = 1, color = "black")
  ),
  showlegend = TRUE,
  name = "Nodes",
  legendgroup = "nodes"
)
p <- p %>% add_trace(
  x = c(NA),
  y = c(NA),
  mode = "markers",
  marker = list(
    colorbar = list(title = "Fold Change"),
    colorscale = list(
      list(0, "blue"),
      list(0.5, "grey80"),
      list(1, "red")
    ),
    showscale = TRUE,
    color = 0
  ),
  showlegend = FALSE,
  hoverinfo = "none"
)
p <- p %>% layout(
  title = "Interactive Bacteria–Metabolite Network",
  xaxis = list(visible = FALSE),
  yaxis = list(visible = FALSE),
  plot_bgcolor = "#FFFFFF",
  paper_bgcolor = "#FFFFFF",
  legend = list(orientation = "v")
)

# --------------------------------------------------------------------
# ⭐ FIXED LEGEND (bacteria + metabolites separately)
# --------------------------------------------------------------------

# Bacteria legend entry
p <- p %>% add_trace(
  x = NA, y = NA,
  type = "scatter", mode = "markers",
  marker = list(color = "orange", size = 14),
  name = "Bacteria",
  legendgroup = "bacteria",
  inherit = FALSE,
  showlegend = TRUE,
  hoverinfo = "none"
)

# Metabolite legend entry
p <- p %>% add_trace(
  x = NA, y = NA,
  type = "scatter", mode = "markers",
  marker = list(color = "skyblue", size = 14),
  name = "Metabolites",
  legendgroup = "metabolites",
  inherit = FALSE,
  showlegend = TRUE,
  hoverinfo = "none"
)
  
p
p_fixed <- plotly::plotly_build(p)

htmlwidgets::saveWidget(
  widget = p_fixed,
  file = "3_data_analysis/preliminary_determination_bacterial_metabolomics/differential_analysis/best_annot/network_plots/bacteria_metabolite_interactive.html",
  selfcontained = TRUE
)