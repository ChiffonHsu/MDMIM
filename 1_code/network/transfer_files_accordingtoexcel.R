library(readr)
library(dplyr)
library(stringr)
library(fs)

# Path to your Excel file
csv_path <- "MS2Library_ChiFang/In_vivo/in_vivo_data.csv"

# Read the sheet that contains your data
df <- read_csv(csv_path)

# Define source and destination folders
source_folder <- "MS2Library_ChiFang/In_vivo/ST001683_c18negative/reverse_phase_c18tive"
destination_folder <- file.path(source_folder, "samples")

# Create destination folder if it doesn't exist
if (!dir_exists(destination_folder)) {
  dir_create(destination_folder)
}

# Filter for rows where ms_dial_sample_name starts with "" and pull Sample_name
con_files <- df %>%
  filter(!str_starts(Colonization, "germ-free")) %>%
  pull(Sample_name)          

# Copy corresponding .mzml files to destination folder
for (file_name in con_files) {
  source_file <- file.path(source_folder, paste0(file_name, ".mzml"))
  dest_file <- file.path(destination_folder, paste0(file_name, ".mzml"))
  
  if (file_exists(source_file)) {
    file_copy(source_file, dest_file, overwrite = TRUE)
    message("Copied: ", file_name)
  } else {
    warning("File not found: ", source_file)
  }
}