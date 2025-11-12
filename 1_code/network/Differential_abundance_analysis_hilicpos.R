##Check for column 'media'in annotated_object_merge, annotated_object_merge
load("MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/annotated_object_merge")
library(massdataset)
head(extract_sample_info(annotated_object_merge))
#The column with media  "Sample_Data:media"

#Check for QC clustering ->check for reproducibility =>no qc samples - skipped this part
load("MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/annotated_object_merge")
object <- annotated_object_merge
head(extract_sample_info(object))
unique(extract_sample_info(object)$sample_type)
table(extract_sample_info(annotated_object_merge)$sample_type)
library(massdataset)
library(massqc)
#  "sample_type" marks QC vs supernatant/media_blank
qc_result<- massqc(
  object = object,
  qc_tag = "QC",
  qc_column = "sample_type"
)

##Differential analysis
##uSE TIDYMASS+ dplyr
library(dplyr)
exp_data <- extract_expression_data(annotated_object_merged)
sample_info <- extract_sample_info(annotated_object_merged)
results <- apply(exp_data,1 , function(x){
  t_res <- t.test (x ~ sample_info$sample_type)
  data.frame(p_value = t_res$p.value,
             mean_supernatant = mean(x[sample_info$sample_type=="supernatant"], na.rm=TRUE),
             mean_medium = mean(x[sample_info$sample_type=="media_blank"], na.rm=TRUE),
             log2FC = log2(mean(x[sample_info$sample_type=="supernatant"],na.rm=TRUE)/
                             mean(x[sample_info$sample_type=="media_blank"], na.rm=TRUE)))
})

results_df <- do.call(rbind, results)
results_df$adj_pval <- p.adjust (results_df$p_value, method="BH")
##Now I have p-values, adjusted p-values, means, and log2 fold-changes. 
##The next step is to add in up or downregulated features using log2FC
results_df$direction <- ifelse(results_df$log2FC > 0, "up", "down")
#Filter significant features by significant p values <0.05
sig_features <- subset(results_df, adj_pval <0.05)

#Visualize volcano plot
library(ggplot2)
p <- ggplot(results_df, aes(x = log2FC, y = -log10(adj_pval), color = direction)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c ("up" = "red", "down" = "blue"))

ggsave(
  filename = "volcano_plot_HILIC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/HILIC/POS",
  width = 8, height = 6, dpi = 300
)

#OR filter by a log2FC cutoff (e.g., |log2FC| > 1)
#Plot to find out the cutoff point to choose
library(ggplot2)
log2FC_cutoff <- 1
pval_cutoff <- 0.05
results_df$significance <- "ns"
results_df$significance[results_df$log2FC > log2FC_cutoff & results_df$adj_pval < pval_cutoff] <- "up"
results_df$significance[results_df$log2FC < -log2FC_cutoff & results_df$adj_pval < pval_cutoff] <- "down"
p <- ggplot(results_df, aes(x = log2FC, y = -log10(adj_pval), color = significance)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c( -log2FC_cutoff, log2FC_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed")
print(p)
ggsave(
  filename = "log2FC_plot_HILIC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/HILIC/POS",
  width = 8, height = 6, dpi = 300
)

#By comparing between mediums
library(dplyr)
#Make sure columns of exp_data match sample_info$sample_id
colnames(exp_data) <- sample_info$sample_id

results_list <- lapply(unique(sample_info$"Sample_Data:media"), function(med){
  # Identify samples in the current medium
  idx <- !is.na(sample_info$`Sample_Data:media`) & sample_info$`Sample_Data:media` == med
  #Assign first
  samples_in_med <- sample_info$sample_id[idx]     
  #Remove NA
  samples_in_med <- samples_in_med[!is.na(samples_in_med)] 
  #Keep only existing columns
  samples_in_med <- intersect(samples_in_med, colnames(exp_data)) 
  
  #Subset expression data and sample info
  exp_sub <- exp_data[, samples_in_med, drop = FALSE]
  sample_sub <- sample_info[sample_info$sample_id %in% samples_in_med, ]
  
  #Apply the above for each feature
  res <-apply(exp_sub, 1, function(x){
    super_vals <- x[sample_sub$sample_type=="supernatant"]
    media_vals <- x[sample_sub$sample_type=="media_blank"]
    #Only run the t test if both groups have observations
    if(length(super_vals[!is.na(super_vals)]) >= 2 & length(media_vals[!is.na(media_vals)]) >= 2){
      mean_sup <- mean(super_vals, na.rm=TRUE)
      mean_med <- mean(media_vals, na.rm=TRUE)
      log2fc <- log2(mean_sup / mean_med)
      t_res <- t.test(super_vals, media_vals)
      
      data.frame(
        medium = med,
        p_value = t_res$p.value,
        mean_supernatant = mean_sup,
        mean_medium = mean_med,
        log2FC = log2fc
      )
      
    } else {
      data.frame(
        medium = med,
        p_value = NA,
        mean_supernatant = NA,
        mean_medium = NA,
        log2FC = NA
      )
    }
  })
  do.call(rbind, res)
})

#Combine all media results into one dataframe
results_df <- do.call(rbind, results_list)

#Add feature_id column
results_df$feature_id <- rep(rownames(exp_data), times = length(unique(sample_info$`Sample_Data:media`)))

#Adjust pvalues
results_df$adj_pval <- p.adjust(results_df$p_value, method = "BH")

##Now I have p-values, adjusted p-values, means, and                   log2 fold-changes. 
##The next step is to add in up or downregulated features using log2FC
results_df$direction <- ifelse(results_df$log2FC > 0, "up", "down")
#Filter significant features by significant p values <0.05
sig_features <- subset(results_df, adj_pval <0.05)
#Visualize volcano plot
 om_point() +
  theme_minimal() +
  scale_color_manual(values = c ("up" = "red", "down" = "blue"))

ggsave(
  filename = "comparebtwmediums_log2fc_volcano_plot_HILIC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/HILIC/POS",
  width = 8, height = 6, dpi = 300
)

#According to log2fc cutoffs, put in the dashes
library(ggplot2)
log2FC_cutoff <- 1
pval_cutoff <- 0.05
results_df$significance <- "ns"
results_df$significance[results_df$log2FC > log2FC_cutoff & results_df$adj_pval < pval_cutoff] <- "up"
results_df$significance[results_df$log2FC < -log2FC_cutoff & results_df$adj_pval < pval_cutoff] <- "down"
p <- ggplot(results_df, aes(x = log2FC, y = -log10(adj_pval), color = significance)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +
  geom_vline(xintercept = c( -log2FC_cutoff, log2FC_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed")
print(p)
ggsave(
  filename = "comparebtwmediums_log2FCcutoffs_plot_HILIC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/HILIC/POS",
  width = 8, height = 6, dpi = 300
)

#Save
# Subset all significant metabolites
sig_all <- subset(results_df, significance != "ns")

# Save to CSV
write.csv(sig_all, 
          file = "MS2Library_ChiFang/In_Vitro/HILIC/POS/significant_metabolites.csv", 
                    row.names = FALSE)

#summary
#results_df - results after differential analysis, including NAs
#exp_data - with normalized abundance values
#sig_all - all significant features, no NAs


#Heatmaps
#Step 1: Check if sg_all is usable
nrow(sig_all) #Check for the number of rows
length(unique(sig_all$feature_id)) #check number of unique features
colSums(is.na(sig_all)[, c("feature_id", "log2FC", "adj_pval")]) #check for nas in these columns
sum(is.na(sig_all$log2FC)) #check for na values in log2FC
sum(is.na(sig_all$adj_pval)) #Check for na values in adjusted pvalues
  #Results for the above: all '0' values, means there are no nas in any of these three columns

#Step 2: Choose one medium - choose 'mm' first
feat_ids_mm <- unique(sig_all$feature_id[sig_all$medium == "mm"])
length(feat_ids_mm) #Result: 79

#Step 3: Extract intensity matrix from exp_data and transform
sig_exp <- exp_data[rownames(exp_data) %in% feat_ids_mm, ,  drop = FALSE]
#Clean na values
row_na_pct <- rowMeans(is.na(sig_exp))
sig_exp <- sig_exp[row_na_pct <= 0.5, , drop = FALSE]
if(nrow(sig_exp) == 0) stop("No features left after NA filtering, Try lowering threshold or picking another medium")


#Log-transform and convert sig_exp to matrix
mat <- as.matrix(sig_exp)
mat_log <- log2(mat + 1)
#Scale by row (z score per feature) for the heatmap
mat_scaled <- t(scale(t(mat_log)))
#Remove rows and columns with NA
mat_scaled <- mat_scaled[rowMeans(is.na(mat_scaled)) < 0.5, , drop = FALSE]
mat_scaled <- mat_scaled[, colMeans(is.na(mat_scaled)) < 0.5, drop = FALSE]
nrow(sig_all)

#Ensure no inf remains
mat_scaled[!is.finite(mat_scaled)] <- 0

#Step 4: Prepare column annotations
#Check sample_info, make sure it has sample_id that matches exp_data columns
#Columns: sample_type, Sample_Data:media
annotation_col <- sample_info[match(colnames(mat_scaled), sample_info$sample_id),
                              c("sample_type", "Sample_Data:media"), drop = FALSE]
rownames(annotation_col) <- colnames(mat_scaled)
#Check again
head(annotation_col)

#Draw heatmap mm vs others(to see the difference across media)
library(pheatmap)
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         scale = "none")
pheatmap::pheatmap(mat_scaled, 
                   annotation_col = annotation_col, 
                   filename = "sig_features_mmvsothers_heatmap.png",
                   width = 20, height = 16)


#Visualize heatmap of only mm samples (to visualize variation withinmm)
#Step 1: Check if sg_all is usable
nrow(sig_all) #Check for the number of rows
length(unique(sig_all$feature_id)) #check number of unique features
colSums(is.na(sig_all)[, c("feature_id", "log2FC", "adj_pval")]) #check for nas in these columns
sum(is.na(sig_all$log2FC)) #check for na values in log2FC
sum(is.na(sig_all$adj_pval)) #Check for na values in adjusted pvalues
#Results for the above: all '0' values, means there are no nas in any of these three columns

#Step 2: Choose one medium - choose 'mm' first
feat_ids_mm <- unique(sig_all$feature_id[sig_all$medium == "mm"])
length(feat_ids_mm) #Result: 79

#Step 3: Extract intensity matrix from exp_data and transform
##To visualize vartiation within mm medium
sig_exp <- exp_data[rownames(exp_data) %in% feat_ids_mm,
                    colnames(exp_data) %in% sample_info$sample_id[
                      sample_info$"Sample_Data:media" == "mm"
                      ],
                    drop = FALSE]

#Clean na values
row_na_pct <- rowMeans(is.na(sig_exp))
sig_exp <- sig_exp[row_na_pct <= 0.5, , drop = FALSE]
if(nrow(sig_exp) == 0) stop("No features left after NA filtering, Try lowering threshold or picking another medium")


#Log-transform and convert sig_exp to matrix
mat <- as.matrix(sig_exp)
mat_log <- log2(mat + 1)
#Scale by row (z score per feature) for the heatmap
mat_scaled <- t(scale(t(mat_log)))
#Remove rows and columns with NA
mat_scaled <- mat_scaled[rowMeans(is.na(mat_scaled)) < 0.5, , drop = FALSE]
mat_scaled <- mat_scaled[, colMeans(is.na(mat_scaled)) < 0.5, drop = FALSE]
nrow(sig_all)

#Ensure no inf remains
mat_scaled[!is.finite(mat_scaled)] <- 0

#Step 4: Prepare column annotations
#Check sample_info, make sure it has sample_id that matches exp_data columns
#Columns: sample_type, Sample_Data:media
annotation_col <- sample_info[match(colnames(mat_scaled), sample_info$sample_id),
                              c("sample_type", "Sample_Data:media"), drop = FALSE]
rownames(annotation_col) <- colnames(mat_scaled)
#Check again
head(annotation_col)

#Draw heatmap mm vs others(to see the difference across media)
library(pheatmap)
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         scale = "none")
pheatmap::pheatmap(mat_scaled, 
                   annotation_col = annotation_col, 
                   filename = "sig_featuresannotated_mm_heatmap.png",
                   width = 20, height = 16)
