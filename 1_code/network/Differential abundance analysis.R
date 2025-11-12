##Check for column 'media'in object_pos, object_neg
load("MS2Library_ChiFang/In_Vitro/RPLC/NEG/Result/data_cleaning/object_pos")
library(massdataset)
head(extract_sample_info(object_pos))
#The column with media  "Sample_Data:media"

#Check for QC clustering ->check for reproducibility =>no qc samples - skipped this part
load("MS2Library_ChiFang/In_Vitro/RPLC/POS/Result/data_cleaning/object_pos")
object <- object_pos
head(extract_sample_info(object))
unique(extract_sample_info(object)$sample_type)
table(extract_sample_info(object_pos)$sample_type)
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
exp_data <- extract_expression_data(object_pos)
sample_info <- extract_sample_info(object_pos)
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
 ##Now I have p-values, adjusted p-values, means, and                   log2 fold-changes. 
   ##The next step is to add in up or downregulated features using log2FC
results_df$direction <- ifelse(results_df$log2FC > 0, "up", "down")
   #Filter significant features by significant p values <0.05
sig_features <- subset(results_df, adj_pval <0.05)
  

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
  filename = "log2FC_plot_RPLC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/RPLC/POS",
  width = 8, height = 6, dpi = 300
)

#Visualize volcano plot
library(ggplot2)
p <- ggplot(results_df, aes(x = log2FC, y = -log10(adj_pval), color = direction)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c ("up" = "red", "down" = "blue"))

ggsave(
  filename = "volcano_plot_RPLC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/RPLC/POS",
  width = 8, height = 6, dpi = 300
  )

#Run t-test and differential analysis by comparison between the respective mediums
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

  #Adjust pvalues
results_df$adj_pval <- p.adjust(results_df$p_value, method = "BH")

##Now I have p-values, adjusted p-values, means, and                   log2 fold-changes. 
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
  filename = "comparebtwmediums_log2fc_volcano_plot_RPLC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/RPLC/POS",
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
  filename = "comparebtwmediums_log2FCcutoffs_plot_RPLC_pos.png",
  plot = p,
  path = "MS2Library_ChiFang/In_Vitro/RPLC/POS",
  width = 8, height = 6, dpi = 300
)






#Add in significance levels
log2FC_cutoff <- 1
pval_cutoff <- 0.05
results_df$significance <- "ns"
results_df$significance[results_df$log2FC > log2FC_cutoff & results_df$adj_pval < pval_cutoff] <- "up"
results_df$significance[results_df$log2FC < -log2FC_cutoff & results_df$adj_pval < pval_cutoff <- "down"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ]

# OR  Use MetaboAnalystR
library(MetaboAnalystR)