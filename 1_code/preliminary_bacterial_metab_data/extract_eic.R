library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(tidymass)
library(tidyverse)

library(massprocesser)
library(massdataset)
targeted_table = 
  readr::read_csv("2_data/preliminary_determination_bacterial_metabolomics/POS/rawdata/Result/Peak_table_for_cleaning.csv")
nrow(targeted_table) #2723

load("2_data/preliminary_determination_bacterial_metabolomics/POS/rawdata/Result/intermediate_data/xdata3")
extract_eic(
  targeted_table = targeted_table,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = 15,
  rt_tolerance = 30,
  threads = 5,
  add_point = FALSE,
  path = "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics",
  group_for_figure = "Case", 
  feature_type = "png"
)

#check number of eic pngs generated for rplc_pos  
# 1️⃣ Folder with PNGs
folder_path <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics/feature_EIC"
num_png <- length(list.files(path = folder_path, pattern = "\\.png$", full.names = TRUE))
num_png #6313

# 2️⃣ List all PNG files (just names, no path)
png_files <- list.files(path = folder_path, pattern = "\\.png$", full.names = FALSE)

# 3️⃣ Remove the .png extension to get variable_id
png_ids <- sub("\\.png$", "", png_files)

# 4️⃣ Count how many times each PNG appears
png_counts <- table(png_ids)

# 5️⃣ Match targeted_rplc_pos variable_id to PNG counts
output_df <- targeted_table %>%
  dplyr::mutate(
    num_png = as.integer(png_counts[variable_id]),  # get counts from table
    num_png = ifelse(is.na(num_png), 0, num_png)   # replace NAs with 0
  )
# 6️⃣ Write to CSV
write.csv(output_df, file = "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics/pos_counts.csv",
          row.names = FALSE)

extra_pngs <- setdiff(png_ids, targeted_table$variable_id)
extra_pngs #0

#remove all extra pngs from the folder
#extra_png_paths <- file.path(folder_path, paste0(extra_pngs, ".png"))
#file.remove(extra_png_paths)
##check again 
folder_path <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics/feature_EIC"
num_png <- length(list.files(path = folder_path, pattern = "\\.png$", full.names = TRUE))
num_png #6313

######################################################################################################################################
library(tidymass)
library(tidyverse)

library(massprocesser)
library(massdataset)
targeted_table = 
  readr::read_csv("2_data/preliminary_determination_bacterial_metabolomics/NEG/rawdata/Result/Peak_table_for_cleaning.csv")
nrow(targeted_table) #1231

load("2_data/preliminary_determination_bacterial_metabolomics/NEG/rawdata/Result/intermediate_data/xdata3")
extract_eic(
  targeted_table = targeted_table,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = 15,
  rt_tolerance = 30,
  threads = 5,
  add_point = FALSE,
  path = "3_data_analysis/preliminary_determination_bacterial_metabolomics/NEG/extract_eics",
  group_for_figure = "Case", 
  feature_type = "png"
)

#check number of eic pngs generated for rplc_pos  
# 1️⃣ Folder with PNGs
folder_path <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics/feature_EIC"
num_png <- length(list.files(path = folder_path, pattern = "\\.png$", full.names = TRUE))
num_png #6313

# 2️⃣ List all PNG files (just names, no path)
png_files <- list.files(path = folder_path, pattern = "\\.png$", full.names = FALSE)

# 3️⃣ Remove the .png extension to get variable_id
png_ids <- sub("\\.png$", "", png_files)

# 4️⃣ Count how many times each PNG appears
png_counts <- table(png_ids)

# 5️⃣ Match targeted_rplc_pos variable_id to PNG counts
output_df <- targeted_table %>%
  dplyr::mutate(
    num_png = as.integer(png_counts[variable_id]),  # get counts from table
    num_png = ifelse(is.na(num_png), 0, num_png)   # replace NAs with 0
  )
# 6️⃣ Write to CSV
write.csv(output_df, file = "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics/pos_counts.csv",
          row.names = FALSE)

extra_pngs <- setdiff(png_ids, targeted_table$variable_id)
extra_pngs #0

#remove all extra pngs from the folder
#extra_png_paths <- file.path(folder_path, paste0(extra_pngs, ".png"))
#file.remove(extra_png_paths)
##check again 
folder_path <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics/feature_EIC"
num_png <- length(list.files(path = folder_path, pattern = "\\.png$", full.names = TRUE))
num_png #6313