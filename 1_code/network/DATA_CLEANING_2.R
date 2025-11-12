library(massqc)
library(ggplot2)
# 1. Boxplot per sample
 p_box <- massqc_sample_boxplot(
      object_neg,
      color_by = "class",    
      fill_by = "class"
    )
 ggsave("sample_boxplot.png", p_box, width = 8, height = 6, dpi = 300)
 # 2. RSD plot for QC samples
 p_rsd <- massqc_rsd_plot(object_pos)
 ggsave("rsd_plot.png", p_rsd, width = 8, height = 6, dpi = 300)
 # 3. Cumulative RSD plot
 p_cum_rsd <- massqc_cumulative_rsd_plot(object_pos)
 ggsave("cumulative_rsd_plot.png", p_cum_rsd, width = 8, height = 6, dpi = 300)
 # 4. PCA plot of samples
 p_pca <- massqc_pca(object_pos)
 ggsave("pca_plot.png", p_pca, width = 8, height = 6, dpi = 300)
 # 5. QC sample correlation
 p_corr <- massqc_sample_correlation(object_pos)
 ggsave("sample_correlation.png", p_corr, width = 30, height = 30, dpi = 300)
 setwd()


 sample_info_df <- object_pos@sample_info
 
 # Filter for media_blank only
blank_ids <- sample_info_df$sample_id[sample_info_df$sample_type == "media_blank"]
 
 object_pos_blanks <- object_pos[, blank_ids]
 p_corr <- massqc_sample_correlation(object_pos_blanks)
 ggsave("sample_correlation_blanks.png", p_corr, width = 8, height = 6, dpi = 300)
 massqc_sample_correlation

library(massdataset)

# 1. Boxplot per sample
#> p_box <- massqc_sample_boxplot(
 # +   object_pos,
  #+   color_by = "class",     # or whichever metadata you prefer
  #+   fill_by = "class"
  #+ )
#> ggsave("sample_boxplot.png", p_box, width = 8, height = 6, dpi = 300)
#value change is not obvious , need to transform to log-intensity units ->modify first

library(massdataset)

# Extract expression data matrix
exprs <- object_pos@expression_data

# Apply log2(x + 1) transform
exprs_log <- log2(exprs + 1)

# Replace expression data with transformed data
object_pos_log <- object_pos
object_pos_log@expression_data <- exprs_log


library(massqc)

p_box <- massqc_sample_boxplot(
  object_neg_log,
  color_by = "class",
  fill_by  = "class"
)

ggplot2::ggsave("sample_boxplot_log2.png", p_box, width = 8, height = 6, dpi = 300)

massqc::massqc_report(object = object_neg, 
                      path = "data_cleaning/NEG/data_quality_before_data_cleaning")

