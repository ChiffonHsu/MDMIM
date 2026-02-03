library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(tidymass)
#1. process raw data
process_data(
  path = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/rawdata",
  polarity = "positive",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 2,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  group_for_figure = "QC"
)
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/rawdata/Result/object")
object #4813 x 79 data frame (variable_info)
#2. generate interactive plot
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/rawdata/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "bacteria_samples",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_case_bacteria_samples.html", selfcontained = TRUE)
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "Blank",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_case_medium_samples.html", selfcontained = TRUE)

p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "Blank",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_blank.html", selfcontained = TRUE)
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "QC",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_QC.html", selfcontained = TRUE)
######################################################################################################################################################
#1. process raw data
process_data(
  path = "2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/rawdata",
  polarity = "negative",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 5,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  group_for_figure = "QC"
)
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result/object")
object #1838 x 118  data frame (variable_info); 120 x 4 (sample_info)
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "bacteria_samples",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_case_neg_bacteria_samples.html", selfcontained = TRUE)
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "medium_samples",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_case_medium_samples.html", selfcontained = TRUE)

p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "Blank",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_blank.html", selfcontained = TRUE)
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "all",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_all.html", selfcontained = TRUE)
    #troubleshooting
output_path <- "2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result"
path = "."
f.in <- list.files(path = path, pattern = "\\.(mz[X|x]{0,1}[M|m][L|l]|cdf)", 
                   recursive = TRUE, full.names = TRUE)

# Example path (adjust to your directory)
# Set the path to the raw data files
# Extract retention times from the featureData slot
retention_times <- xdata2@featureData@data$retentionTime

# Check the range of the retention times
rt_range <- range(retention_times, na.rm = TRUE)
message("Retention time range: ", paste(rt_range, collapse = " to ")) # Retention time range: 0.252 to 1799.68999999998

path <- "2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/rawdata"

# List files with the desired extensions (.mzML, .mzXML, .cdf)
f.in <- list.files(path = path, pattern = "\\.(mz[X|x]{0,1}[M|m][L|l]|cdf)", recursive = TRUE, full.names = TRUE)

# Extract the sample group (the second-to-last part of the directory path)
sample_group <- unlist(lapply(stringr::str_split(f.in, "/"), function(x) {
  x[length(x) - 1]  # This grabs the group name from the directory structure
}))

# Create a data frame 'pd' with sample names (without file extensions) and their corresponding groups
pd <- data.frame(
  sample_name = sub(basename(f.in), pattern = "\\.(mz[X|x]{0,1}[M|m][L|l]|cdf)", replacement = "", fixed = TRUE), 
  sample_group = sample_group, 
  stringsAsFactors = FALSE, check.names = FALSE
)

# Verify pd structure
head(pd)

# Ensure you have valid file indices (idx) selected
figure_sample_name <- pd %>% dplyr::filter(sample_group %in% group_for_figure) %>% dplyr::pull(sample_name)
idx <- which(basename(xdata2@processingData@files) %in% figure_sample_name)

# Ensure you have BiocParallel set up correctly
library(BiocParallel)

# Define the number of threads (workers) you want to use
threads <- 4  # Set this to the desired number of parallel threads/workers

# Determine if you're on Windows or a Unix-like OS
if (masstools::get_os() == "windows") {
  # Use SnowParam for Windows
  bpparam <- BiocParallel::SnowParam(workers = threads, progressbar = TRUE)
} else {
  # Use MulticoreParam for non-Windows OS
  bpparam <- BiocParallel::MulticoreParam(workers = threads, progressbar = TRUE)
}

# Check if idx has valid files
if (length(idx) > 0) {
  # Generate the TIC plot
  tic.plot <- xcms::chromatogram(object = filterFile(object = xdata2, file = idx), aggregationFun = "sum", BPPARAM = bpparam)
  
  # Check if the plot was generated successfully
  if (!is.null(tic.plot)) {
    # Create the plot with the desired title
    plot_tic <- plot_chromatogram(object = tic.plot, title = "TIC", group_for_figure = group_for_figure)
    
    # Save the TIC plot manually as a PDF
    ggplot2::ggsave(filename = file.path(output_path, "TIC.pdf"), plot = plot_tic, width = 20, height = 7)
    message("TIC plot saved successfully.")
  } else {
    message("Error: TIC plot is NULL.")
  }
} else {
  message("Error: No valid indices for TIC plot.")
}




# Ensure you have the correct sample indices
figure_sample_name <- pd %>% dplyr::filter(sample_group %in% group_for_figure) %>% dplyr::pull(sample_name)
idx <- which(basename(xdata2@processingData@files) %in% figure_sample_name)

# Generate the TIC plot (Total Ion Chromatogram)
if (length(idx) > 0) {
  tic.plot <- tryCatch(xcms::chromatogram(object = filterFile(object = xdata2, file = idx), 
                                          aggregationFun = "sum", BPPARAM = bpparam), 
                       error = function(e) NULL)
  if (!is.null(tic.plot)) {
    plot_tic <- plot_chromatogram(object = tic.plot, title = "TIC", 
                                  group_for_figure = group_for_figure)
    ggplot2::ggsave(filename = file.path(output_path, "TIC.pdf"), plot = plot_tic, width = 20, height = 7)
    message("TIC plot saved successfully.")
  } else {
    message("Error: TIC plot is NULL.")
  }
} else {
  message("Error: No valid indices for TIC plot.")
}






# Generate the BPC plot (base peak chromatogram)
if (length(idx) > 0) {
  # Generate the BPC plot
  bpc.plot <- xcms::chromatogram(object = filterFile(object = xdata2, file = idx), aggregationFun = "max", BPPARAM = bpparam)
  
  # Check if the plot was generated successfully
  if (!is.null(bpc.plot)) {
    # Create the plot with the desired title
    plot_bpc <- plot_chromatogram(object = bpc.plot, title = "BPC", group_for_figure = group_for_figure)
    
    # Save the BPC plot manually as a PDF
    ggplot2::ggsave(filename = file.path(output_path, "BPC.pdf"), plot = plot_bpc, width = 20, height = 7)
    message("BPC plot saved successfully.")
  } else {
    message("Error: BPC plot is NULL.")
  }
} else {
  message("Error: No valid indices for BPC plot.")
}








#2. generate interactive plot
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata3,
  group_for_figure = "Case",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_case.html", selfcontained = TRUE)
load("2_data/preliminary_determination_bacterial_metabolomics/NEG/rawdata/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "QC",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_QC.html", selfcontained = TRUE)
load("2_data/preliminary_determination_bacterial_metabolomics/NEG/rawdata/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "Blank",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_blank.html", selfcontained = TRUE)
