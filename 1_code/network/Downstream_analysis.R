
load("MS2Library_ChiFang/In_Vitro/R_objects/annotated_object_merged.rda")

#The file with all significant features: sig_all, with differential analysis results
#The one with all annotated features: annotated_object_merged
#What now: Annotate features in "sig_all" using "annotated_object_merged"
str(annotated_object_merged) #check structure of annotated_object_merged
str(sig_all) #check structure of sig_all
annotation_info <- extract_annotation_table(annotated_object_merged)
head(annotation_info)
library(dplyr)
sig_all_annotated <- sig_all %>%
  left_join(annotation_info, by = c("feature_id" = "variable_id"))
head(sig_all_annotated) #check if it is joined properly
sum(!is.na(sig_all_annotated$Compound.name))

#Plot heatmap with top 5 most significant features for both upregulated and downregulates
library(dplyr)
sig_filtered <- sig_all_annotated %>%
  filter(adj_pval < 0.05, abs(log2FC) > 1, !is.na(Compound.name))
top_up <- sig_filtered %>%
  filter(log2FC > 0) %>%
  arrange(adj_pval, desc(log2FC)) %>%
  head(5)
top_down <- sig_filtered %>%
  filter(log2FC <0 ) %>%
  arrange(adj_pval, log2FC) %>%
  head(5)
top_hits <- bind_rows(top_up, top_down)

library(ggplot2)
log2FC_cutoff <- 1
pval_cutoff <- 0.05
results_df$significance <- "ns"
results_df$significance[results_df$log2FC > log2FC_cutoff & results_df$adj_pval < pval_cutoff] <- "up"
results_df$significance[results_df$log2FC < -log2FC_cutoff & results_df$adj_pval < pval_cutoff] <- "down"

library(ggrepel)
p <- ggplot(results_df, aes(x = log2FC, y = -log10(adj_pval), color = significance)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c( -log2FC_cutoff, log2FC_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed") +
  geom_text_repel(
    data = top_hits,
    aes(label = Compound.name),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf
  )
print(p)
ggsave(
  filename = "volcano_top5sig_annotated.png",
  plot = p,
  path = "R_objects",
  width = 16, height = 12, dpi = 500
)

#Plot heatmap with top 100 most significant features for both upregulated and downregulated features
library(dplyr)
library(pheatmap)

# Step 1: Select top 100 up and down within "mm" medium
sig_mm <- sig_all %>%
  filter(medium == "mm") %>%
  filter(!is.na(Compound.name), adj_pval < 0.05, abs(log2FC) > 1)

top_up <- sig_mm %>%
  filter(log2FC > 0) %>%
  arrange(adj_pval, desc(log2FC)) %>%
  head(100)

top_down <- sig_mm %>%
  filter(log2FC < 0) %>%
  arrange(adj_pval, log2FC) %>%
  head(100)

top_features_mm <- bind_rows(top_up, top_down)
feat_ids_mm <- unique(top_features_mm$feature_id)
length(feat_ids_mm)  # number of features to plot

# Step 2: Extract intensity matrix for mm medium
sig_exp <- exp_data[rownames(exp_data) %in% feat_ids_mm,
                    colnames(exp_data) %in% sample_info$sample_id[
                      sample_info$`Sample_Data:media` == "mm"
                    ],
                    drop = FALSE]

# Step 3: Clean NA values
row_na_pct <- rowMeans(is.na(sig_exp))
sig_exp <- sig_exp[row_na_pct <= 0.5, , drop = FALSE]
if(nrow(sig_exp) == 0) stop("No features left after NA filtering, Try lowering threshold or picking another medium")

# Step 4: Log-transform and scale
mat <- as.matrix(sig_exp)
mat_log <- log2(mat + 1)
mat_scaled <- t(scale(t(mat_log)))

# Remove remaining rows/columns with NA
mat_scaled <- mat_scaled[rowMeans(is.na(mat_scaled)) < 0.5, , drop = FALSE]
mat_scaled <- mat_scaled[, colMeans(is.na(mat_scaled)) < 0.5, drop = FALSE]
mat_scaled[!is.finite(mat_scaled)] <- 0

# Step 5: Prepare column annotations
annotation_col <- sample_info[match(colnames(mat_scaled), sample_info$sample_id),
                              c("sample_type", "Sample_Data:media"), drop = FALSE]
rownames(annotation_col) <- colnames(mat_scaled)

# Step 6: Draw heatmap
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "top100_sig_features_mm_heatmap.png",
         width = 20, height = 16)