# #8_differential analysis hybrid
# load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/object_pos3_ms2")
# 
# colnames(object_pos3_ms2@sample_info)
# head(object_pos3_ms2@sample_info)
# # sample_id experiment sample_type Bacteria medium class batch group
# # 1  Blank_1_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
# # 2 Blank_10_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
# # 3 Blank_11_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
# # 4  Blank_2_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
# # 5  Blank_3_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
# # 6  Blank_4_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
# 
# #under object_pos3_ms2@sample_info$class, samples are grouped as Control and Subject, compare these two groups
# 
# #1. Get bacteria samples and their corresponding medium information
# bacteria_info <- object_pos3_ms2@sample_info %>%
#   filter(group == "bacteria_samples") %>%
#   dplyr::select(sample_id, Bacteria, medium, class)
# 
# #2. Create a vector of bacteria sample ids
# bacteria_classes <- unique(bacteria_info$Bacteria)
# #3. extra medium info based on medium values
# medium_6_samples <- object_pos3_ms2@sample_info %>%
#   filter(group == "medium_samples", medium == "MGAM-6" ) %>%
#   pull(sample_id)
# 
# #4. calculate p value for the differential analysis
# safe_fc <- function(object, control_ids, case_ids) {
#   expr <- object@expression_data
#   
#   # FIX: add pseudocount to avoid log2(0) issues
#   expr <- expr + 1e-6
#   
#   control_mean <- rowMeans(expr[, control_ids, drop = FALSE], na.rm = TRUE)
#   case_mean    <- rowMeans(expr[, case_ids, drop = FALSE], na.rm = TRUE)
#   
#   log2(case_mean / control_mean)
# }
# 
# ### Raw p-value (later adjust using BH)
# my_pvalue <- function(object, control_ids, case_ids) {
#   expr <- object@expression_data
#   
#   apply(expr, 1, function(x) {
#     x_control <- x[control_ids]
#     x_case    <- x[case_ids]
#     
#     x_control <- x_control[!is.na(x_control)]
#     x_case    <- x_case[!is.na(x_case)]
#     
#     if (length(x_control) < 2 | length(x_case) < 2) return(1)
#     if (var(x_control)==0 & var(x_case)==0) return(1)
#     
#     tryCatch(t.test(x_control, x_case)$p.value, error=function(e) 1)
#   })
# }
# #5. Loop through bacteria samples and perform comparisons
# ###############################################
# ### 5. Hybrid analysis based on sample size ###
# ###############################################
# comparison_results <- purrr::map(bacteria_classes, function(cls) {
#   
#   case_samples <- bacteria_info$sample_id[bacteria_info$Bacteria == cls]
#   bacteria_medium <- unique(bacteria_info$medium[bacteria_info$Bacteria == cls])[1]
#   
#   control_samples <- medium_6_samples
#   
#   n_case <- length(case_samples)
#   n_control <- length(control_samples)
#   
#   message("\nClass: ", cls, " | case=", n_case, " | control=", n_control)
#   
#   #################################################
#   ### CASE A: Full analysis (>=3 samples/group)
#   #################################################
#   if (n_case >= 3 & n_control >= 3) {
#     message(" → FULL analysis (mutate_fc + BH p-values)")
#     
#     # compute FC using pseudocount FC
#     FC_full <- safe_fc(object_pos3_ms2, control_samples, case_samples)
#     
#     # compute p-values using massdataset (raw)
#     obj_p <- object_pos3_ms2 %>%
#       mutate_p_value(control_sample_id = control_samples,
#                      case_sample_id = case_samples,
#                      method = "t.test",
#                      p_adjust_methods = "BH")
#     
#     df <- tibble(
#       variable_id = rownames(object_pos3_ms2@expression_data),
#       class = cls,
#       fc_final = FC_full,
#       pval_final = obj_p@variable_info$p_value_adjust
#     )
#     
#     
#     return(df)
#     
#   }
#   #################################################
#   ### CASE B: FALLBACK — SAFE hybrid version
#   ### Use mutate_fc for FC, manual t-test for p-value
#   #################################################
#   message(" → FALLBACK analysis (safe FC + BH p-values)")
#   
#   FC     <- safe_fc(object_pos3_ms2, control_samples, case_samples)
#   PVAL_r <- my_pvalue(object_pos3_ms2, control_samples, case_samples)
#   PVAL_a <- p.adjust(PVAL_r, method = "BH")   ### FIX ⭐ adjust p-values!
#   
#   tibble(
#     variable_id = rownames(object_pos3_ms2@expression_data),
#     class = cls,
#     fc_final = FC,
#     pval_final = PVAL_a
#   )
# })
# 
# names(comparison_results) <- bacteria_classes
# 
# ########################################
# ### 6.5 Add compound names (DO NOT LOSE NAMES)
# ########################################
# compound_lookup <- object_pos3_ms2@annotation_table %>%
#   dplyr::select(variable_id, Compound.name)
# 
# comparison_results <- purrr::map(comparison_results, function(df) {
#   dplyr::left_join(df, compound_lookup, by = "variable_id")
# })
# 
# ########################################
# ### Collapse duplicates
# ########################################
# collapsed_results <- lapply(comparison_results, function(df) {
#   df %>%
#     mutate(
#       ### FIX ⭐ Replace NA names with Unknown_varid
#       Compound.name = ifelse(
#         is.na(Compound.name),
#         paste0("Unknown_", variable_id),
#         Compound.name
#       )
#     ) %>%
#     group_by(Compound.name) %>%
#     arrange(pval_final, desc(abs(fc_final))) %>%
#     slice(1) %>%
#     ungroup()
# })
# 
# names(collapsed_results) <- bacteria_classes
# 
# #7. volcano plot
# #a. combine into one big tibble
# library(dplyr)
# library(purrr)
# all_results <- bind_rows(collapsed_results)
# 
# all_results <- all_results %>%
#   mutate(
#     logFC = fc_final,
#     neg_log10_p = -log10(pval_final)
#   )
# 
# volcano_prepare <- function(df) {
#   df %>%
#     mutate(
#       sig = pval_final < 0.05 & abs(fc_final) > 1,
#       regulation = case_when(
#         !sig ~ "none",
#         fc_final > 1 ~ "up",
#         fc_final < -1 ~ "down"
#       )
#     )
# }
# 
# volcano_plot2 <- function(df, class_name) {
#   
#   df <- volcano_prepare(df)
#   
#   label_df <- df %>%
#     filter(regulation != "none") %>%
#     arrange(pval_final) %>%
#     slice_head(n = 10)
#   
#   p <- ggplot(df, aes(x = fc_final, y = -log10(pval_final), color = regulation)) +
#     geom_point(alpha = 0.7) +
#     scale_color_manual(values = c("none"="grey70","up"="red","down"="blue")) +
#     theme_minimal() +
#     labs(
#       title = paste("Volcano Plot -", class_name),
#       x = "Log2 Fold Change",
#       y = "-log10(P-value)"
#     )
#   
#   if (nrow(label_df) > 0) {
#     p <- p + geom_text_repel(
#       data = label_df,
#       aes(label = Compound.name),
#       size = 3,
#       max.overlaps = 20
#     )
#   }
#   
#   return(p)
# }
# #e.save each plot automatically
# volcano_list <- purrr::map(
#   names(collapsed_results),
#   ~ volcano_plot2(collapsed_results[[.x]], .x)
# )
# 
# names(volcano_list) <- names(collapsed_results)
# 
# walk2(
#   volcano_list,
#   names(volcano_list),
#   ~ ggsave(
#     filename = paste0("volcano_", .y, ".pdf"),
#     path = "~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/batch2_rplc_pos_differential_analysis/hybrid_adj_fc_pval",
#     plot = .x,
#     width = 8,
#     height = 6
#   )
# )
# 
# ###############PART TWO: COLLATE SIGNIFICANT FEATURES IN EACH BACTERIA
# #(1) Define significance function
# sig_details <- lapply(collapsed_results, function(df) {
#   
#   df2 <- df %>%
#     mutate(
#       sig = pval_final < 0.05 & abs(fc_final) > 1,
#       regulation = case_when(
#         !sig ~ "none",
#         fc_final > 1 ~ "up",
#         fc_final < -1 ~ "down"
#       )
#     )
#   
#   tibble(
#     up = sum(df2$regulation == "up"),
#     down = sum(df2$regulation == "down"),
#     total_sig = sum(df2$regulation != "none"),
#     percent = round(100 * sum(df2$regulation != "none") / nrow(df2), 2)
#   )
# })
# 
# sig_details <- bind_rows(sig_details, .id = "bacterium")
# 
# sig_details
# # # A tibble: 9 × 5
# # bacterium                     up  down total_sig percent
# # <chr>                      <int> <int>     <int>   <dbl>
# #   1 Citrobacter youngae          313   497       810    2.01
# # 2 Citrobacter braakii         1699  3728      5427   13.4 
# # 3 Bacteroides fragilis        1738  5256      6994   17.3 
# # 4 Enterococcus casseliflavus  1361  3401      4762   11.8 
# # 5 Enterobacter ludwigii       1656  3153      4809   11.9 
# # 6 Pediococcus acidilactici    1225  4256      5481   13.6 
# # 7 Bacteroides finegoldii      1676  4139      5815   14.4 
# # 8 Enterococcus faecalis       1383  3586      4969   12.3 
# # 9 Enterococcus durans         1725  5095      6820   16.9 
# 
# hybrid_sig_summary <- sig_summary
# save(hybrid_sig_summary, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/hybrid_sig_summary")
# 
# hybrid_collapsed_results <- collapsed_results
# save(hybrid_collapsed_results, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/hybrid_collapsed_results")
# 
# hybrid_comparison_results <- comparison_results
# save(hybrid_comparison_results, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/hybrid_comparison_results")
# 
# #9_differential analysis hybrid
# load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/object_pos3_ms2")
# load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/hybrid_collapsed_results")
# #1. Create simple and cleaned tibble for each individual bacterium
# network_input <- lapply(hybrid_collapsed_results, function(df){
#   df %>%
#     mutate(
#       regulation = case_when(
#         pval_final < 0.05 & fc_final > 1 ~ "up",
#         pval_final < 0.05 & fc_final < -1 ~ "down",
#         TRUE ~ "none"
#       )
#     ) %>%
#     filter(regulation != "none") %>%
#     dplyr::select(class, Compound.name, regulation, fc_final, pval_final)
# })
# names(network_input) <- bacteria_classes
# 
# #2. build network edges from downregulated to upregulated
# #A. ADD in inchi_key and KEGG, CAS.ID
# #add column "inchi" back to the annotation tabler of the object
# load("~/compound_database_shenlab/MS1/ms1_database.rda")
# colnames(ms1_database@spectra.info)#  "Lab.ID"  "INCHI_ID"  "INCHI_ID_all"  "INCHIKEY_ID" "INCHIKEY_ID_all"  
# annotation_table <- object_pos3_ms2@annotation_table
# ms1_database_df <- ms1_database@spectra.info
# inchi_info <- ms1_database_df[, c("INCHIKEY_ID", "Lab.ID")]
# object_pos3_ms2_inchi <- merge(
#   annotation_table, 
#   inchi_info, 
#   by = "Lab.ID", 
#   all.x = TRUE
# )
# annotation_lookup <- object_pos3_ms2@annotation_table %>%
#   dplyr::select(
#     variable_id,
#     KEGG.ID,
#     HMDB.ID,
#     CAS.ID,
#     Lab.ID
#   ) %>%
#   dplyr::left_join(
#     ms1_database@spectra.info %>%
#       dplyr::select(Lab.ID, INCHIKEY_ID),
#     by = "Lab.ID"
#   )
# 
# hybrid_collapsed_results <- lapply(hybrid_collapsed_results, function(df) {
#   dplyr::left_join(df, annotation_lookup, by = "variable_id")
# })
# 
# 
# network_input <- lapply(hybrid_collapsed_results, function(df){
#   df %>%
#     mutate(
#       regulation = case_when(
#         pval_final < 0.05 & fc_final > 1 ~ "up",
#         pval_final < 0.05 & fc_final < -1 ~ "down",
#         TRUE ~ "none"
#       )
#     ) %>%
#     filter(regulation != "none") %>%
#     dplyr::select(class, Compound.name, regulation, fc_final, pval_final, KEGG.ID, HMDB.ID, CAS.ID, INCHIKEY_ID)
# })
# names(network_input) <- bacteria_classes
# 
# #B. Prepare edges
# network_nodes <- lapply(network_input, function(df){
#   df %>%
#     mutate(
#       direction = case_when(
#         regulation == "down" ~ "from",
#         regulation == "up"   ~ "to",
#         TRUE ~ NA_character_
#       )
#     )
# })
# 
# #C. Check for pairing with metacyc, annotate using kegg.id
# # per bacterium
# metacyc_edges_kegg <- readRDS("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/metcyc_database/metacyc_edges_kegg.rds")
# metacyc_nodes_all <- readRDS("~/compound_database_shenlab/metacyc/metacyc_nodes_all.rds")
# 
# metacyc_pairs_kegg <- lapply(names(network_nodes), function(bact){
#   
#   df <- hybrid_collapsed_results[[bact]]
#   
#   kegg_ids <- df$KEGG.ID %>% na.omit() %>% unique()
#   
#   # reaction pairs where either substrate OR product matches
#   matched <- metacyc_edges_kegg %>%
#     filter(
#       substrate_kegg %in% kegg_ids |
#         product_kegg   %in% kegg_ids
#     )
#   
#   matched
# })
# 
# names(metacyc_pairs_kegg) <- names(network_nodes)
# 
# 
# 
# library(tidymass)
# library(dplyr)
# library(igraph)
# library(ggraph)
# library(ggrepel)
# library(scales)
# 
# # -------------------------------
# # 1. Prepare nodes
# # -------------------------------
# 
# # Extract info
# df <- as.data.frame(hybrid_collapsed_results$`Citrobacter youngae`)
# 
# # Make sure info are correct
# df <- df %>%
#   mutate(
#     variable_id = as.character(variable_id),
#     Compound.name = ifelse(is.na(Compound.name) | Compound.name == "",
#                            variable_id, Compound.name)
#   )
# 
# # build nodes
# nodes <- df %>%
#   group_by(variable_id) %>%
#   summarise(
#     label     = paste(unique(Compound.name), collapse = ";"),
#     mean_fc   = mean(fc_final, na.rm = TRUE),
#     direction = case_when(
#       mean_fc >  4.644 ~ "up",    # tweak thresholds to taste
#       mean_fc < -4.644 ~ "down",
#       TRUE        ~ "neutral"
#     ),
#     fc_abs = abs(mean_fc),
#     .groups = "drop"
#   ) %>%
#   filter(direction != "neutral")   # keep only up/down
# 
# # Add central node
# central_name <- "Citrobacter youngae"
# central_node <- tibble(
#   variable_id = central_name,
#   label       = central_name,
#   direction   = "central",
#   fc_abs      = NA_real_
# )
# # Combine nodes
# nodes_all <- bind_rows(central_node, nodes) %>%
#   mutate(
#     size_scaled = if_else(direction == "central",
#                           20,
#                           scales::rescale(fc_abs, to = c(5, 15))),
#     # igraph expects a "name" column for vertex IDs
#     name = variable_id,
#     # create a 'size' column to use directly in igraph/ggraph
#     size = size_scaled
#   )
# nodes_all <- nodes_all %>%
#   mutate(
#     fontface = if_else(direction == "central", "bold", "plain")
#   )
# 
# # --------------------------------------------------------
# # Pick top 10 up and top 10 down regulated nodes for labels
# # --------------------------------------------------------
# top_up <- nodes_all %>%
#   filter(direction == "up") %>%
#   arrange(desc(fc_abs)) %>%
#   slice_head(n = 10) %>%
#   pull(variable_id)
# 
# top_down <- nodes_all %>%
#   filter(direction == "down") %>%
#   arrange(desc(fc_abs)) %>%
#   slice_head(n = 10) %>%
#   pull(variable_id)
# 
# # Label only these nodes + the central node
# nodes_all <- nodes_all %>%
#   mutate(
#     label_to_plot = case_when(
#       variable_id == central_name      ~ label,     # always keep bacterium label
#       variable_id %in% top_up          ~ label,     # top 10 up
#       variable_id %in% top_down        ~ label,     # top 10 down
#       TRUE                              ~ ""         # hide all others
#     )
#   )
# 
# # Scale node sizes
# size_range <- c(5, 15)
# nodes <- nodes %>%
#   mutate(size_scaled = ifelse(direction != "central",
#                               scales::rescale(fc_abs, to = size_range),
#                               20))
# 
# # -------------------------------
# # 2. Create edges
# # -------------------------------
# edges <- nodes_all %>%
#   filter(direction != "central") %>%
#   transmute(from = variable_id, to = central_name)
# 
# # -------------------------------
# # 3. Build igraph object
# # -------------------------------
# g <- graph_from_data_frame(edges, vertices = nodes_all, directed = TRUE)
# layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
# 
# V(g)$label_to_plot <- nodes_all$label_to_plot
# V(g)$fontface <- nodes_all$fontface
# 
# 
# central_idx <- which(V(g)$direction == "central")
# up_idx      <- which(V(g)$direction == "up")
# down_idx    <- which(V(g)$direction == "down")
# 
# compute_angles <- function(idx_vec, start_angle, end_angle) {
#   n <- length(idx_vec)
#   if (n == 1) return((start_angle + end_angle) / 2)
#   seq(start_angle, end_angle, length.out = n + 1)[-1]
# }
# 
# # Up-regulated nodes: top semicircle
# if (length(up_idx) > 0) {
#   angles_up <- compute_angles(up_idx, 0, pi)
#   radii_up  <- 1.5 * (1 - scales::rescale(V(g)$size[up_idx], to = c(0, 1))) + 0.5
#   layout_matrix[up_idx, ] <- cbind(radii_up * cos(angles_up),
#                                    radii_up * sin(angles_up))
# }
# 
# # Down-regulated nodes: bottom semicircle
# if (length(down_idx) > 0) {
#   angles_down <- compute_angles(down_idx, pi, 2 * pi)
#   radii_down  <- 1.5 * (1 - scales::rescale(V(g)$size[down_idx], to = c(0, 1))) + 0.5
#   layout_matrix[down_idx, ] <- cbind(radii_down * cos(angles_down),
#                                      radii_down * sin(angles_down))
# }
# 
# # Central node at origin
# layout_matrix[central_idx, ] <- c(0, 0)
# 
# g_layout <- create_layout(
#   g,
#   layout = "manual",
#   x = layout_matrix[, 1],
#   y = layout_matrix[, 2]
# )
# 
# #plot
# ggraph(g_layout) +
#   geom_edge_link(arrow = arrow(length = unit(1, "mm")), alpha = 0.6) +
#   geom_node_point(aes(size = size, color = direction)) +
#   geom_node_text(
#     aes(label = label_to_plot, fontface = fontface),
#     repel = TRUE,
#     size = 4
#   )+
#   scale_color_manual(values = c(up = "red",
#                                 down = "blue",
#                                 central = "yellow")) +
#   theme_void()
# 
# 
# png("~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/rplc_pos/Single_bacteria_map/Citrobacter youngae_starnetwork_fc2.png",
#     width = 6000, height = 6000, res = 300)
# 
# ggraph(g_layout) +
#   geom_edge_link(arrow = arrow(length = unit(4, "mm")), alpha = 0.6) +
#   geom_node_point(aes(size = size, color = direction)) +
#   geom_node_text(aes(label = label), repel = TRUE, size = 4) +
#   scale_color_manual(values = c("up" = "red", "down" = "blue", "central" = "yellow")) +
#   scale_size_continuous(range = c(5, 15), name = "|FC|") +
#   theme_void()
# 
# dev.off()
# 
# #10_network_plot
# #Generate network plots for all bacteria
# load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/hybrid_collapsed_results")
# 
# plot_bacteria_network <- function(df, central_name, output_file = NULL) {
#   
#   library(dplyr)
#   library(igraph)
#   library(ggraph)
#   library(scales)
#   library(ggrepel)
#   
#   # ---------------------------------
#   # 1. Prepare nodes
#   # ---------------------------------
#   df <- df %>%
#     mutate(
#       variable_id   = as.character(variable_id),
#       Compound.name = ifelse(is.na(Compound.name) | Compound.name == "",
#                              variable_id, Compound.name)
#     )
#   
#   nodes <- df %>%
#     group_by(variable_id) %>%
#     summarise(
#       label     = paste(unique(Compound.name), collapse = ";"),
#       mean_fc   = mean(fc_final, na.rm = TRUE),
#       mean_pval = mean(pval_final, na.rm = TRUE),   # <-- FIXED COMMA !!!
#       
#       direction = case_when(
#         mean_fc >  2.322 & mean_pval < 0.05 ~ "up",
#         mean_fc < -2.322 & mean_pval < 0.05 ~ "down",
#         TRUE ~ "neutral"
#       ),
#       
#       fc_abs = abs(mean_fc),
#       .groups = "drop"
#     ) %>%
#     filter(direction != "neutral")
#   
#   # central node
#   central_node <- tibble(
#     variable_id = central_name,
#     label       = central_name,
#     direction   = "central",
#     fc_abs      = NA_real_
#   )
#   
#   nodes_all <- bind_rows(central_node, nodes) %>%
#     mutate(
#       size_scaled = if_else(direction == "central",
#                             20,
#                             scales::rescale(fc_abs, to = c(5, 15))),
#       name  = variable_id,
#       size  = size_scaled,
#       fontface = if_else(direction == "central", "bold", "plain")
#     )
#   
#   # -------------------------
#   # top 10 up/down labels
#   # -------------------------
#   top_up <- nodes %>% 
#     filter(direction == "up") %>% 
#     arrange(desc(fc_abs)) %>% 
#     slice_head(n = 10) %>% 
#     pull(variable_id)
#   
#   top_down <- nodes %>% 
#     filter(direction == "down") %>% 
#     arrange(desc(fc_abs)) %>% 
#     slice_head(n = 10) %>% 
#     pull(variable_id)
#   
#   nodes_all <- nodes_all %>%
#     mutate(
#       label_to_plot = case_when(
#         variable_id == central_name ~ label,
#         variable_id %in% top_up ~ label,
#         variable_id %in% top_down ~ label,
#         TRUE ~ ""
#       )
#     )
#   
#   # -------------------------
#   # Edges
#   # -------------------------
#   edges <- nodes_all %>%
#     filter(direction != "central") %>%
#     transmute(from = variable_id, to = central_name)
#   
#   # -------------------------
#   # Graph
#   # -------------------------
#   g <- graph_from_data_frame(edges, vertices = nodes_all, directed = TRUE)
#   
#   layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
#   
#   V(g)$label_to_plot <- nodes_all$label_to_plot
#   V(g)$fontface      <- nodes_all$fontface
#   
#   central_idx <- which(V(g)$direction == "central")
#   up_idx      <- which(V(g)$direction == "up")
#   down_idx    <- which(V(g)$direction == "down")
#   
#   compute_angles <- function(idx_vec, start_angle, end_angle) {
#     n <- length(idx_vec)
#     if (n == 1) return((start_angle + end_angle) / 2)
#     seq(start_angle, end_angle, length.out = n + 1)[-1]
#   }
#   
#   if (length(up_idx) > 0) {
#     angles <- compute_angles(up_idx, 0, pi)
#     radii  <- 1.5 * (1 - scales::rescale(V(g)$size[up_idx], to = c(0, 1))) + 0.5
#     layout_matrix[up_idx, ] <- cbind(radii * cos(angles), radii * sin(angles))
#   }
#   
#   if (length(down_idx) > 0) {
#     angles <- compute_angles(down_idx, pi, 2*pi)
#     radii  <- 1.5 * (1 - scales::rescale(V(g)$size[down_idx], to = c(0, 1))) + 0.5
#     layout_matrix[down_idx, ] <- cbind(radii * cos(angles), radii * sin(angles))
#   }
#   
#   layout_matrix[central_idx, ] <- c(0, 0)
#   
#   g_layout <- create_layout(
#     g,
#     layout = "manual",
#     x = layout_matrix[, 1],
#     y = layout_matrix[, 2]
#   )
#   
#   # -------------------------
#   # PLOT
#   # -------------------------
#   p <- ggraph(g_layout) +
#     geom_edge_link(arrow = arrow(length = unit(1, "mm")), alpha = 0.6) +
#     geom_node_point(aes(size = size, color = direction)) +
#     geom_node_text(
#       aes(label = label_to_plot, fontface = fontface),
#       repel = TRUE,
#       max.overlaps = Inf,
#       force = 2,
#       point.padding = 0.3,
#       box.padding = 0.5,
#       size = 4
#     ) +
#     scale_color_manual(values = c(
#       up = "red",
#       down = "blue",
#       central = "yellow"
#     )) +
#     theme_void()
#   
#   if (!is.null(output_file)) {
#     ggsave(output_file, p, width = 10, height = 10, dpi = 300)
#   }
#   
#   return(p)
# }
# 
# 
# 
# # -------------------------------
# # LOOP OVER ALL BACTERIA
# # -------------------------------
# output_dir <- "~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/rplc_pos/Single_bacteria_map"
# dir.create(output_dir, showWarnings = FALSE)
# 
# for (bac in names(collapsed_results)) {
#   
#   message("Processing: ", bac)
#   
#   df <- collapsed_results[[bac]]
#   
#   outfile <- file.path(
#     output_dir,
#     paste0(gsub(" ", "_", bac), "_network.png")
#   )
#   
#   plot_bacteria_network(
#     df = df,
#     central_name = bac,
#     output_file = outfile
#   )
# }
# #11_build_interaction_network
# #Build a bacteria-metabolite interaction network (bipartite network)
# # Step One: Build a combined table of all bacteria and all significant metabolites
# #Includes merging of metabolites from the new batch
# #A. Combine the entire tibble object “collapsed_results”
# load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251010/MS2/differential_analysis_each_bacterium/collapsed_results")
# batch1_collapsed_results <- collapsed_results
# load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/hybrid_collapsed_results")
# batch2_collapsed_results <- hybrid_collapsed_results
# 
# library(purrr) 
# library(dplyr) 
# batch1_combined <- imap_dfr( 
#   batch1_collapsed_results, 
#   ~mutate(.x, bacteria = .y) ) 
# colnames(batch1_combined) #"variable_id" "class" "fc_final" "pval_final" "Compound.name" "bacteria" 
# nrow(batch1_combined) #32615
# 
# batch2_combined <- imap_dfr( 
#   batch2_collapsed_results, 
#   ~mutate(.x, bacteria = .y) ) 
# colnames(batch2_combined). # "variable_id" "class" "fc_final" "pval_final" "Compound.name" "bacteria" 
# nrow(batch2_combined) #363528
# 
# batches_combined <- rbind(batch1_combined, batch2_combined)
# nrow(batches_combined) #396143
# 
# #B. Filter out and keep only significant metabolties
# sig_combined <- batches_combined %>%
#   mutate(variable_id = as.character(variable_id),
#          Compound.name = ifelse(is.na(Compound.name) | Compound.name == "", 
#                                 variable_id, Compound.name)) %>%
#   filter(abs(fc_final) > 2.322,
#          pval_final < 0.05)
# #C. Select only the top 20 metabolites for each bacteria
# top20_metabs <- sig_combined %>%
#   arrange(desc(abs(fc_final))) %>%
#   distinct(variable_id, .keep_all = TRUE) %>%
#   slice_head(n = 40) %>%
#   pull(variable_id)
# 
# sig_top <- sig_combined %>%
#   filter(variable_id %in% top20_metabs)
# 
# #Step Two: Build the bipartite edge list (bacteria->metabilite)
# #include directions so that up and down have different coloured edges
# edges_bip <- sig_top %>%
#   transmute(
#     from = bacteria,
#     to   = paste0("M:", variable_id),
#     fc   = fc_final
#   ) %>% 
#   mutate(direction = ifelse(fc > 0, "up", "down"))
# 
# #the code below is for all metabolites, not just the top f20
# edges_bip <- sig_combined %>%
#   transmute(
#     from = bacteria,
#     to   = paste0("M:", variable_id),
#     fc   = fc_final
#   ) %>% 
#   mutate(direction = ifelse(fc > 0, "up", "down"))
# 
# #Step Three: Build nodes list with node types, two nodes types -> nodes_bacteria and nodes_metabolite
# # A. Bacteria nodes (type = "bacteria")
# nodes_bateria <- tibble(
#   name = unique(sig_combined$bacteria),
#   type = "bacteria",
#   label = name
# )
# # B. Metabolite nodes (type = "metabolite")
# nodes_metabolite <- sig_combined %>%
#   distinct(variable_id, Compound.name) %>%
#   transmute(
#     name = paste0("M:", variable_id),
#     type = "metabolite",
#     label = Compound.name
#   )
# # C. Combine both nodes_bacteria and nodes_metabolite and de-duplicate
# nodes_all <- bind_rows(nodes_bateria, nodes_metabolite)
# nodes_all <- nodes_all %>%
#   dplyr::distinct(name, .keep_all = TRUE)
# 
# # D. Sizing according to absolute fc
# met_size <- sig_combined %>%
#   mutate(node = paste0("M:", variable_id)) %>%
#   group_by(node) %>%
#   summarise(mean_fc_abs = mean(abs(fc_final), na.rm = TRUE))
# 
# #Step Four: Build igraph bipartite object
# library(igraph)
# nodes_all <- nodes_all %>%
#   filter(name %in% c(edges_bip$from, edges_bip$to))
# g <- graph_from_data_frame(
#   d = edges_bip,
#   vertices = nodes_all,
#   directed = FALSE
# )
# 
# #for labeling nodes for the key
# # Keep type as character
# V(g)$type <- nodes_all$type
# 
# # for sizing
# V(g)$size <- 6  
# met_idx <- which(V(g)$type == "metabolite")
# V(g)$size[met_idx] <- scales::rescale(
#   met_size$mean_fc_abs[ match(V(g)$name[met_idx], met_size$node) ],
#   to = c(3, 6)
# )
# 
# #Step Five: Plot bipartite plot using Fruchterman–Reingold
# library(ggraph)
# ggraph(g, layout = "fr") +
#   geom_edge_link(aes(edge_color = fc), alpha = 0.3) +
#   scale_edge_color_gradient2(
#     low = "blue",
#     high = "red",
#     mid = "grey80",
#     midpoint = 0) +
#   geom_node_point(aes(color = type, size = size)) +
#   
#   #geom_node_text(aes(label = label), repel = TRUE, size = 3) +
#   scale_color_manual(values = c(
#     "bacteria" = "orange",
#     "metabolite" = "skyblue"
#   )) +
#   
#   theme_void()
# 
# #Build interaction plot (will show labels when hover over to the edge or nodes)
# #Color of edges according to fc 
# library(igraph)
# library(plotly)
# library(dplyr)
# 
# set.seed(123)
# layout_fr <- layout_with_fr(g)
# 
# layout_df <- as.data.frame(layout_fr)
# colnames(layout_df) <- c("x", "y")
# 
# layout_df$name  <- V(g)$name
# layout_df$type  <- V(g)$type
# layout_df$size  <- V(g)$size
# layout_df$label <- V(g)$label
# 
# node_colors <- ifelse(layout_df$type == "bacteria", "orange", "skyblue")
# edges_df <- igraph::as_data_frame(g, "edges") %>%
#   mutate(
#     x0 = layout_df$x[match(from, layout_df$name)],
#     y0 = layout_df$y[match(from, layout_df$name)],
#     x1 = layout_df$x[match(to, layout_df$name)],
#     y1 = layout_df$y[match(to, layout_df$name)],
#     hover = paste0("<b>Bacteria:</b> ", from,
#                    "<br><b>Metabolite:</b> ", to,
#                    "<br><b>FC:</b> ", round(fc, 3))
#   )
# 
# p <- plot_ly()
# 
# # Add each edge separately so gradient works
# for (i in seq_len(nrow(edges_df))) {
#   
#   fc_i <- edges_df$fc[i]
#   
#   # map FC to color manually
#   if (fc_i < 0) {
#     col_i <- colorRampPalette(c("blue", "grey80"))(100)[
#       cut(scales::rescale(fc_i, to=c(1, -1)), 100)
#     ]
#   } else {
#     col_i <- colorRampPalette(c("grey80", "red"))(100)[
#       cut(scales::rescale(fc_i, to=c(-1, 1)), 100)
#     ]
#   }
#   
#   p <- p %>% add_trace(
#     x = c(edges_df$x0[i], edges_df$x1[i]),
#     y = c(edges_df$y0[i], edges_df$y1[i]),
#     type = "scatter",
#     mode = "lines",
#     line = list(color = col_i, width = 0.4),
#     text = edges_df$hover[i],
#     hoverinfo = "text",
#     showlegend = FALSE
#   )
# }
# p <- p %>% add_trace(
#   x = layout_df$x,
#   y = layout_df$y,
#   type = "scatter",
#   mode = "markers",
#   text = layout_df$label,
#   hoverinfo = "text",
#   marker = list(
#     size = layout_df$size * 1.8,
#     color = node_colors,
#     line = list(width = 1, color = "black")
#   ),
#   showlegend = TRUE,
#   name = "Nodes",
#   legendgroup = "nodes"
# )
# p <- p %>% add_trace(
#   x = c(NA),
#   y = c(NA),
#   mode = "markers",
#   marker = list(
#     colorbar = list(title = "Fold Change"),
#     colorscale = list(
#       list(0, "blue"),
#       list(0.5, "grey80"),
#       list(1, "red")
#     ),
#     showscale = TRUE,
#     color = 0
#   ),
#   showlegend = FALSE,
#   hoverinfo = "none"
# )
# p <- p %>% layout(
#   title = "Interactive Bacteria–Metabolite Network",
#   xaxis = list(visible = FALSE),
#   yaxis = list(visible = FALSE),
#   plot_bgcolor = "#FFFFFF",
#   paper_bgcolor = "#FFFFFF",
#   legend = list(orientation = "v")
# )
# 
# # --------------------------------------------------------------------
# # ⭐ FIXED LEGEND (bacteria + metabolites separately)
# # --------------------------------------------------------------------
# 
# # Bacteria legend entry
# p <- p %>% add_trace(
#   x = NA, y = NA,
#   type = "scatter", mode = "markers",
#   marker = list(color = "orange", size = 14),
#   name = "Bacteria",
#   legendgroup = "bacteria",
#   inherit = FALSE,
#   showlegend = TRUE,
#   hoverinfo = "none"
# )
# 
# # Metabolite legend entry
# p <- p %>% add_trace(
#   x = NA, y = NA,
#   type = "scatter", mode = "markers",
#   marker = list(color = "skyblue", size = 14),
#   name = "Metabolites",
#   legendgroup = "metabolites",
#   inherit = FALSE,
#   showlegend = TRUE,
#   hoverinfo = "none"
# )
# 
# p
# p_fixed <- plotly::plotly_build(p)
# 
# htmlwidgets::saveWidget(
#   widget = p_fixed,
#   file = "~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/rplc_pos/network_interaction_2batches/bacteria_metabolite_interactive.html",
#   selfcontained = TRUE
# )

#This function runs differential analysis, differential metabolite and bacteria map for single bacteria, then de-duplicating all
#repeating line, to finally build the the interaction network. Before moving on to this function, there must be a merged object 
#for the hilic and rplc neg pos
generate_volcano_plots <- function(
    object_path,
    out_root,
    medium_control_map = c("MGAM-6" = "MGAM-6"),
    p_cut_volcano = 0.05,
    fc_cut_volcano = 1,
    p_cut_network = 0.05,
    fc_cut_network = 2.322,
    use_tidymass = TRUE
) {
  # ---- packages ----
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
    library(ggplot2)
    library(ggrepel)
    if (use_tidymass) library(tidymass)
  })
  
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  # ---- load object ----
  load(object_path)  # expects object_pos3_ms2 in env; adjust if name differs
  if (!exists("object_pos3_ms2")) stop("object_pos3_ms2 not found after load()")
  
  obj <- object_pos3_ms2
  
  # ---- helper functions ----
  safe_fc <- function(object, control_ids, case_ids, pseudocount = 1e-6) {
    expr <- object@expression_data + pseudocount
    control_mean <- rowMeans(expr[, control_ids, drop = FALSE], na.rm = TRUE)
    case_mean    <- rowMeans(expr[, case_ids, drop = FALSE], na.rm = TRUE)
    log2(case_mean / control_mean)
  }
  
  my_pvalue <- function(object, control_ids, case_ids) {
    expr <- object@expression_data
    apply(expr, 1, function(x) {
      x_control <- x[control_ids]; x_case <- x[case_ids]
      x_control <- x_control[!is.na(x_control)]
      x_case    <- x_case[!is.na(x_case)]
      if (length(x_control) < 2 || length(x_case) < 2) return(1)
      if (var(x_control) == 0 && var(x_case) == 0) return(1)
      tryCatch(t.test(x_control, x_case)$p.value, error = function(e) 1)
    })
  }
  
  # ---- inputs from sample_info ----
  bacteria_info <- obj@sample_info %>%
    filter(group == "bacteria_samples") %>%
    select(sample_id, Bacteria, medium)
  
  bacteria_names <- unique(bacteria_info$Bacteria)
  
  get_control_samples <- function(case_medium) {
    ctrl_medium <- unname(medium_control_map[as.character(case_medium)])
    if (is.na(ctrl_medium)) stop("No control mapping for medium: ", case_medium)
    obj@sample_info %>%
      filter(group == "medium_samples", medium == ctrl_medium) %>%
      pull(sample_id)
  }
  
  # ==========================================================
  # STEP 1: differential per bacterium  (AUTOMATIC LOOP)
  # ==========================================================
  comparison_results <- purrr::map(bacteria_names, function(bac) {
    
    case_samples <- bacteria_info$sample_id[bacteria_info$Bacteria == bac]
    case_medium  <- unique(bacteria_info$medium[bacteria_info$Bacteria == bac])[1]
    control_samples <- get_control_samples(case_medium)
    
    message("Bacteria: ", bac, " | case=", length(case_samples), " | ctrl=", length(control_samples))
    
    # tidymass path if available
    can_tidymass <- use_tidymass &&
      exists("mutate_p_value", mode = "function") &&
      length(case_samples) >= 3 && length(control_samples) >= 3
    
    if (can_tidymass) {
      
      obj_p <- mutate_p_value(
        object = obj,
        control_sample_id = control_samples,
        case_sample_id = case_samples,
        method = "t.test",
        p_adjust_methods = "BH"
      )
      
      tibble(
        variable_id = rownames(obj@expression_data),
        bacteria    = bac,
        fc_final    = safe_fc(obj, control_samples, case_samples),
        pval_final  = obj_p@variable_info$p_value_adjust
      )
      
    } else {
      p_raw <- my_pvalue(obj, control_samples, case_samples)
      tibble(
        variable_id = rownames(obj@expression_data),
        bacteria    = bac,
        fc_final    = safe_fc(obj, control_samples, case_samples),
        pval_final  = p.adjust(p_raw, method = "BH")
      )
    }
  })
  names(comparison_results) <- bacteria_names
  
  # ==========================================================
  # STEP 2: add compound names + collapse duplicates
  # ==========================================================
  compound_lookup <- obj@annotation_table %>%
    select(variable_id, Compound.name)
  
  comparison_results <- purrr::map(comparison_results, ~ left_join(.x, compound_lookup, by = "variable_id"))
  
  collapsed_results <- purrr::map(comparison_results, function(df) {
    df %>%
      mutate(Compound.name = ifelse(is.na(Compound.name) | Compound.name == "",
                                    paste0("Unknown_", variable_id),
                                    Compound.name)) %>%
      group_by(Compound.name) %>%
      arrange(pval_final, desc(abs(fc_final))) %>%
      slice(1) %>%
      ungroup()
  })
  names(collapsed_results) <- bacteria_names
  
  # ==========================================================
  # STEP 3: volcano PDFs (loop + save)
  # ==========================================================
  volcano_dir <- file.path(out_root, "volcano")
  dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)
  
  volcano_plot <- function(df, name) {
    df <- df %>%
      mutate(regulation = case_when(
        pval_final < p_cut_volcano & fc_final >  fc_cut_volcano ~ "up",
        pval_final < p_cut_volcano & fc_final < -fc_cut_volcano ~ "down",
        TRUE ~ "none"
      ))
    
    label_df <- df %>% filter(regulation != "none") %>% arrange(pval_final) %>% slice_head(n = 10)
    
    p <- ggplot(df, aes(fc_final, -log10(pval_final), color = regulation)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c(none = "grey70", up = "red", down = "blue")) +
      theme_minimal() +
      labs(title = paste("Volcano -", name), x = "Log2FC", y = "-log10(adj p)")
    
    if (nrow(label_df) > 0) {
      p <- p + geom_text_repel(data = label_df, aes(label = Compound.name), size = 3, max.overlaps = 50)
    }
    p
  }
  
  volcano_list <- purrr::imap(collapsed_results, volcano_plot)
  
  purrr::iwalk(volcano_list, function(p, nm) {
    ggsave(
      filename = file.path(volcano_dir, paste0("volcano_", gsub(" ", "_", nm), ".pdf")),
      plot = p, width = 8, height = 6
    )
  })
  
  # ==========================================================
  # STEP 4: significance summary table
  # ==========================================================
  sig_details <- purrr::imap_dfr(collapsed_results, function(df, bac) {
    df2 <- df %>%
      mutate(sig = pval_final < p_cut_volcano & abs(fc_final) > fc_cut_volcano,
             regulation = case_when(!sig ~ "none", fc_final > 0 ~ "up", fc_final < 0 ~ "down"))
    tibble(
      bacterium = bac,
      up = sum(df2$regulation == "up"),
      down = sum(df2$regulation == "down"),
      total_sig = sum(df2$regulation != "none"),
      percent = round(100 * sum(df2$regulation != "none") / nrow(df2), 2)
    )
  })
  
  # ==========================================================
  # STEP 5: save outputs as RDS (one folder)
  # ==========================================================
  saveRDS(comparison_results, file.path(out_root, "comparison_results.rds"))
  saveRDS(collapsed_results,  file.path(out_root, "collapsed_results.rds"))
  saveRDS(sig_details,        file.path(out_root, "sig_details.rds"))
  
  invisible(list(
    comparison_results = comparison_results,
    collapsed_results  = collapsed_results,
    volcano_list       = volcano_list,
    sig_details        = sig_details
  ))
}

res <- generate_volcano_plots(
  object_path = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/object_pos3_ms2",
  out_root    = "~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/test/hybrid_adj_fc_pval",
  medium_control_map = c("MGAM-6" = "MGAM-6"),
  use_tidymass = TRUE
)





run_bacteria_metabolomics_pipeline <- function(
    object_path,
    out_root,
    medium_control_map = c("MGAM-6" = "MGAM-6"),
    use_tidymass = TRUE,
    
    # volcano thresholds (labels + coloring)
    p_cut_volcano = 0.05,
    fc_cut_volcano = 1,
    
    # network thresholds (for including edges)
    p_cut_network = 0.05,
    fc_cut_network = 2.322,
    
    # global network options
    topN_per_bacteria = NULL,   # e.g. 40; NULL = keep all sig
    make_interactive_html = TRUE
) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(purrr)
    library(ggrepel)
    library(scales)
    library(igraph)
    library(ggraph)
    if (make_interactive_html) {
      library(plotly)
      library(htmlwidgets)
    }
    if (use_tidymass) library(tidymass)
  })
  
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  volcano_dir <- file.path(out_root, "volcano")
  star_dir    <- file.path(out_root, "star_networks")
  global_dir  <- file.path(out_root, "global_network")
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(star_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(global_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---- load object ----
  load(object_path)
  if (!exists("object_pos3_ms2")) stop("object_pos3_ms2 not found after load(). Rename inside function if needed.")
  obj <- object_pos3_ms2
  
  # ---- helpers ----
  safe_fc <- function(object, control_ids, case_ids, pseudocount = 1e-6) {
    expr <- object@expression_data + pseudocount
    control_mean <- rowMeans(expr[, control_ids, drop = FALSE], na.rm = TRUE)
    case_mean    <- rowMeans(expr[, case_ids, drop = FALSE], na.rm = TRUE)
    log2(case_mean / control_mean)
  }
  
  my_pvalue <- function(object, control_ids, case_ids) {
    expr <- object@expression_data
    apply(expr, 1, function(x) {
      x_control <- x[control_ids]; x_case <- x[case_ids]
      x_control <- x_control[!is.na(x_control)]
      x_case    <- x_case[!is.na(x_case)]
      if (length(x_control) < 2 || length(x_case) < 2) return(1)
      if (var(x_control) == 0 && var(x_case) == 0) return(1)
      tryCatch(t.test(x_control, x_case)$p.value, error = function(e) 1)
    })
  }
  
  get_control_samples <- function(case_medium) {
    ctrl_medium <- unname(medium_control_map[as.character(case_medium)])
    if (is.na(ctrl_medium)) stop("No control mapping for medium: ", case_medium)
    obj@sample_info %>%
      dplyr::filter(group == "medium_samples", medium == ctrl_medium) %>%
      dplyr::pull(sample_id)
  }
  
  # ---- bacteria list ----
  bacteria_info <- obj@sample_info %>%
    dplyr::filter(group == "bacteria_samples") %>%
    dplyr::select(sample_id, Bacteria, medium)
  
  bacteria_names <- unique(bacteria_info$Bacteria)
  
  # ---- compound lookup ----
  compound_lookup <- obj@annotation_table %>%
    dplyr::select(variable_id, Compound.name)
  
  # ==========================================================
  # STEP 1: differential per bacterium
  # ==========================================================
  comparison_results <- purrr::map(bacteria_names, function(bac) {
    
    case_samples <- bacteria_info$sample_id[bacteria_info$Bacteria == bac]
    case_medium  <- unique(bacteria_info$medium[bacteria_info$Bacteria == bac])[1]
    control_samples <- get_control_samples(case_medium)
    
    message("Bacteria: ", bac, " | case=", length(case_samples), " | ctrl=", length(control_samples))
    
    can_tidymass <- use_tidymass &&
      exists("mutate_p_value", mode = "function") &&
      length(case_samples) >= 3 && length(control_samples) >= 3
    
    if (can_tidymass) {
      
      obj_p <- mutate_p_value(
        object = obj,
        control_sample_id = control_samples,
        case_sample_id = case_samples,
        method = "t.test",
        p_adjust_methods = "BH"
      )
      
      tibble(
        variable_id = rownames(obj@expression_data),
        bacteria    = bac,
        fc_final    = safe_fc(obj, control_samples, case_samples),
        pval_final  = obj_p@variable_info$p_value_adjust
      )
      
    } else {
      p_raw <- my_pvalue(obj, control_samples, case_samples)
      tibble(
        variable_id = rownames(obj@expression_data),
        bacteria    = bac,
        fc_final    = safe_fc(obj, control_samples, case_samples),
        pval_final  = p.adjust(p_raw, method = "BH")
      )
    }
  })
  names(comparison_results) <- bacteria_names
  
  # ==========================================================
  # STEP 2: add compound names + collapse duplicates
  # ==========================================================
  comparison_results <- purrr::map(comparison_results, ~ dplyr::left_join(.x, compound_lookup, by = "variable_id"))
  
  collapsed_results <- purrr::map(comparison_results, function(df) {
    df %>%
      dplyr::mutate(Compound.name = ifelse(is.na(Compound.name) | Compound.name == "",
                                           paste0("Unknown_", variable_id),
                                           Compound.name)) %>%
      dplyr::group_by(Compound.name) %>%
      dplyr::arrange(pval_final, dplyr::desc(abs(fc_final))) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  })
  names(collapsed_results) <- bacteria_names
  
  # ==========================================================
  # STEP 3: volcano plots (save PDFs)
  # ==========================================================
  volcano_plot <- function(df, name) {
    df <- df %>%
      dplyr::mutate(regulation = dplyr::case_when(
        pval_final < p_cut_volcano & fc_final >  fc_cut_volcano ~ "up",
        pval_final < p_cut_volcano & fc_final < -fc_cut_volcano ~ "down",
        TRUE ~ "none"
      ))
    
    label_df <- df %>%
      dplyr::filter(regulation != "none") %>%
      dplyr::arrange(pval_final) %>%
      dplyr::slice_head(n = 10)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(fc_final, -log10(pval_final), color = regulation)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_color_manual(values = c(none = "grey70", up = "red", down = "blue")) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Volcano -", name), x = "Log2FC", y = "-log10(adj p)")
    
    if (nrow(label_df) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(label = Compound.name),
        size = 3,
        max.overlaps = 50
      )
    }
    p
  }
  
  volcano_list <- purrr::imap(collapsed_results, volcano_plot)
  
  purrr::iwalk(volcano_list, function(p, nm) {
    ggplot2::ggsave(
      filename = file.path(volcano_dir, paste0("volcano_", gsub(" ", "_", nm), ".pdf")),
      plot = p, width = 8, height = 6
    )
  })
  
  # ==========================================================
  # STEP 4: per-bacterium "star" network PNGs
  # ==========================================================
  plot_bacteria_star <- function(df, central_name, outfile) {
    
    df <- df %>%
      dplyr::mutate(
        variable_id = as.character(variable_id),
        Compound.name = ifelse(is.na(Compound.name) | Compound.name == "", variable_id, Compound.name)
      )
    
    # up/down by network cutoffs (more stringent)
    nodes <- df %>%
      dplyr::mutate(direction = dplyr::case_when(
        pval_final < p_cut_network & fc_final >  fc_cut_network ~ "up",
        pval_final < p_cut_network & fc_final < -fc_cut_network ~ "down",
        TRUE ~ "neutral"
      )) %>%
      dplyr::filter(direction != "neutral") %>%
      dplyr::transmute(
        name      = paste0("M:", variable_id),   # metabolite IDs
        label     = Compound.name,
        direction = direction,
        fc_abs    = abs(fc_final)
      ) %>%
      dplyr::distinct(name, .keep_all = TRUE)    # ✅ critical
    
    if (nrow(nodes) == 0) return(invisible(NULL))
    
    # central node (namespaced)
    central_node <- tibble::tibble(
      name      = paste0("B:", central_name),
      label     = central_name,
      direction = "central",
      fc_abs    = NA_real_
    )
    
    nodes_all <- dplyr::bind_rows(central_node, nodes) %>%
      dplyr::distinct(name, .keep_all = TRUE) %>%
      dplyr::mutate(
        size = dplyr::if_else(direction == "central", 20, scales::rescale(fc_abs, to = c(5, 15))),
        fontface = dplyr::if_else(direction == "central", "bold", "plain")
      )
    
    # edges: metabolite -> bacteria (or bacteria -> metabolite; visually star either way)
    edges <- nodes_all %>%
      dplyr::filter(direction != "central") %>%
      dplyr::transmute(from = name, to = paste0("B:", central_name))

    g <- igraph::graph_from_data_frame(edges, vertices = nodes_all, directed = TRUE)
    
    # label only top 10 up/down
    top_up <- nodes_all %>% dplyr::filter(direction == "up") %>% dplyr::arrange(dplyr::desc(fc_abs)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(name)
    top_down <- nodes_all %>% dplyr::filter(direction == "down") %>% dplyr::arrange(dplyr::desc(fc_abs)) %>% dplyr::slice_head(n = 10) %>% dplyr::pull(name)
    
    central_id <- paste0("B:", central_name)
    
    nodes_all <- nodes_all %>%
      dplyr::mutate(label_to_plot = dplyr::case_when(
        name == central_id ~ label,
        name %in% top_up ~ label,
        name %in% top_down ~ label,
        TRUE ~ ""
      ))
    
    V(g)$direction <- nodes_all$direction[match(V(g)$name, nodes_all$name)]
    V(g)$size      <- nodes_all$size[match(V(g)$name, nodes_all$name)]
    V(g)$label_to_plot <- nodes_all$label_to_plot[match(V(g)$name, nodes_all$name)]
    V(g)$fontface  <- nodes_all$fontface[match(V(g)$name, nodes_all$name)]
    
    # manual semicircle layout
    layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
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
    
    g_layout <- ggraph::create_layout(g, layout = "manual", x = layout_matrix[,1], y = layout_matrix[,2])
    
    p <- ggraph::ggraph(g_layout) +
      ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(1, "mm")), alpha = 0.6) +
      ggraph::geom_node_point(ggplot2::aes(size = size, color = direction)) +
      ggraph::geom_node_text(
        ggplot2::aes(label = label_to_plot, fontface = fontface),
        repel = TRUE, max.overlaps = Inf, size = 4
      ) +
      ggplot2::scale_color_manual(values = c(up = "red", down = "blue", central = "gold")) +
      ggplot2::theme_void()
    
    ggplot2::ggsave(outfile, p, width = 10, height = 10, dpi = 300)
    invisible(p)
  }
  
  purrr::iwalk(collapsed_results, function(df, bac) {
    outfile <- file.path(star_dir, paste0(gsub(" ", "_", bac), "_star.png"))
    plot_bacteria_star(df, bac, outfile)
  })
  
  # ==========================================================
  # STEP 5: global bipartite network (edges + static + interactive)
  # ==========================================================
  combined <- purrr::imap_dfr(collapsed_results, ~ dplyr::mutate(.x, bacteria = .y))
  
  sig_combined <- combined %>%
    dplyr::mutate(
      variable_id = as.character(variable_id),
      Compound.name = ifelse(is.na(Compound.name) | Compound.name == "", variable_id, Compound.name)
    ) %>%
    dplyr::filter(abs(fc_final) > fc_cut_network, pval_final < p_cut_network)
  
  if (!is.null(topN_per_bacteria)) {
    sig_combined <- sig_combined %>%
      dplyr::group_by(bacteria) %>%
      dplyr::arrange(dplyr::desc(abs(fc_final))) %>%
      dplyr::slice_head(n = topN_per_bacteria) %>%
      dplyr::ungroup()
  }
  
  edges_bip <- sig_combined %>%
    dplyr::transmute(
      from = bacteria,
      to   = paste0("M:", variable_id),
      fc   = fc_final
    )
  
  nodes_bacteria <- tibble(
    name = unique(sig_combined$bacteria),
    type = "bacteria",
    label = name
  )
  
  nodes_metabolite <- sig_combined %>%
    dplyr::distinct(variable_id, Compound.name) %>%
    dplyr::transmute(
      name = paste0("M:", variable_id),
      type = "metabolite",
      label = Compound.name
    )
  
  nodes_all <- dplyr::bind_rows(nodes_bacteria, nodes_metabolite) %>%
    dplyr::distinct(name, .keep_all = TRUE)
  
  # limit vertices to those appearing in edges
  nodes_all <- nodes_all %>% dplyr::filter(name %in% c(edges_bip$from, edges_bip$to))
  
  g_bip <- igraph::graph_from_data_frame(edges_bip, vertices = nodes_all, directed = FALSE)
  V(g_bip)$type  <- nodes_all$type[match(V(g_bip)$name, nodes_all$name)]
  V(g_bip)$label <- nodes_all$label[match(V(g_bip)$name, nodes_all$name)]
  
  # node size by metabolite mean |fc|
  met_size <- sig_combined %>%
    dplyr::mutate(node = paste0("M:", variable_id)) %>%
    dplyr::group_by(node) %>%
    dplyr::summarise(mean_fc_abs = mean(abs(fc_final), na.rm = TRUE), .groups = "drop")
  
  V(g_bip)$size <- 6
  met_idx <- which(V(g_bip)$type == "metabolite")
  V(g_bip)$size[met_idx] <- scales::rescale(
    met_size$mean_fc_abs[match(V(g_bip)$name[met_idx], met_size$node)],
    to = c(3, 6)
  )
  
  # ---- static PNG ----
  static_png <- file.path(global_dir, "bacteria_metabolite_bipartite.png")
  p_static <- ggraph::ggraph(g_bip, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(edge_color = fc), alpha = 0.25) +
    ggraph::scale_edge_color_gradient2(low = "blue", mid = "grey80", high = "red", midpoint = 0) +
    ggraph::geom_node_point(ggplot2::aes(color = type, size = size)) +
    ggplot2::scale_color_manual(values = c(bacteria = "orange", metabolite = "skyblue")) +
    ggplot2::theme_void()
  
  ggplot2::ggsave(static_png, p_static, width = 12, height = 12, dpi = 300)
  
  # ---- interactive HTML (optional) ----
  html_file <- NULL
  if (make_interactive_html) {
    set.seed(123)
    layout_fr <- igraph::layout_with_fr(g_bip)
    layout_df <- as.data.frame(layout_fr)
    colnames(layout_df) <- c("x", "y")
    layout_df$name  <- V(g_bip)$name
    layout_df$type  <- V(g_bip)$type
    layout_df$size  <- V(g_bip)$size
    layout_df$label <- V(g_bip)$label
    
    node_colors <- ifelse(layout_df$type == "bacteria", "orange", "skyblue")
    
    edges_df <- igraph::as_data_frame(g_bip, "edges") %>%
      dplyr::mutate(
        x0 = layout_df$x[match(from, layout_df$name)],
        y0 = layout_df$y[match(from, layout_df$name)],
        x1 = layout_df$x[match(to, layout_df$name)],
        y1 = layout_df$y[match(to, layout_df$name)],
        hover = paste0("<b>Bacteria:</b> ", from,
                       "<br><b>Metabolite:</b> ", to,
                       "<br><b>FC:</b> ", round(fc, 3))
      )
    
    p <- plotly::plot_ly()
    
    # edges (per-edge color)
    for (i in seq_len(nrow(edges_df))) {
      fc_i <- edges_df$fc[i]
      col_i <- if (fc_i < 0) {
        colorRampPalette(c("blue", "grey80"))(100)[cut(scales::rescale(fc_i, to = c(1, -1)), 100)]
      } else {
        colorRampPalette(c("grey80", "red"))(100)[cut(scales::rescale(fc_i, to = c(-1, 1)), 100)]
      }
      
      p <- p %>% plotly::add_trace(
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
    
    # nodes
    p <- p %>% plotly::add_trace(
      x = layout_df$x, y = layout_df$y,
      type = "scatter", mode = "markers",
      text = layout_df$label,
      hoverinfo = "text",
      marker = list(
        size = layout_df$size * 2,
        color = node_colors,
        line = list(width = 1, color = "black")
      ),
      showlegend = FALSE
    )
    
    p <- p %>% plotly::layout(
      title = "Interactive Bacteria–Metabolite Network",
      xaxis = list(visible = FALSE),
      yaxis = list(visible = FALSE),
      plot_bgcolor = "#FFFFFF",
      paper_bgcolor = "#FFFFFF"
    )
    
    html_file <- file.path(global_dir, "bacteria_metabolite_interactive.html")
    htmlwidgets::saveWidget(plotly::plotly_build(p), file = html_file, selfcontained = TRUE)
  }
  
  # ==========================================================
  # save key R objects
  # ==========================================================
  saveRDS(comparison_results, file.path(out_root, "comparison_results.rds"))
  saveRDS(collapsed_results,  file.path(out_root, "collapsed_results.rds"))
  saveRDS(edges_bip,          file.path(out_root, "edges_bip.rds"))
  saveRDS(nodes_all,          file.path(out_root, "nodes_all.rds"))
  saveRDS(sig_combined,       file.path(out_root, "sig_combined.rds"))
  
  invisible(list(
    comparison_results = comparison_results,
    collapsed_results  = collapsed_results,
    volcano_list       = volcano_list,
    edges_bip          = edges_bip,
    nodes_all          = nodes_all,
    g_bip              = g_bip,
    static_png         = static_png,
    interactive_html   = html_file
  ))
}
res <- run_bacteria_metabolomics_pipeline(
  object_path = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/object_pos3_ms2",
  out_root    = "~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/test/full_pipeline",
  medium_control_map = c("MGAM-6" = "MGAM-6"),
  use_tidymass = TRUE,
  topN_per_bacteria = NULL,          # or 40
  make_interactive_html = TRUE
)
