library(dplyr)
library(massdataset)
library(masstools)
library(ggplot2)

load("2_data/preliminary_determination_bacterial_metabolomics/MS2/POS/object_pos4_ms2")

# Step 3: Create output folder for plots
output_dir <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/ms2_matching_plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Step 4: Load all MS2 database objects
load("~/compound_database_shenlab/MS2/gnps_ms2.rda")
load("~/compound_database_shenlab/MS2/hmdb_ms2.rda")
load("~/compound_database_shenlab/MS2/massbank_ms2.rda")
load("~/compound_database_shenlab/MS2/metlin_ms2.rda")
load("~/compound_database_shenlab/MS2/mona_ms2.rda")
load("~/compound_database_shenlab/MS2/mpsnyder_rplc_ms2.rda")
load("~/compound_database_shenlab/MS2/nist_ms2.rda")

library(dplyr)
object_pos4_ms2@annotation_table %>%
  count(Database, name = "num_matches") %>%
  arrange(desc(num_matches))
#Database                     num_matches
#1 NIST_20220425                          57
#2 METLIN_20220425                        43
#3 MoNA_2022-04-27                        15
#4 HMDB_2022-04-11                        12
#5 GNPS_2022-04-29                         9
#6 MassBank_2022-04-27                     5
#7 Michael_Snyder_rplc_20220424           3

annotation_table <- object_pos4_ms2@annotation_table

database <-
  unique(annotation_table$Database)
database
#[1] "Shen-Lab_0.0.1"               "NIST_20220425"                "METLIN_20220425"             
#[4] "GNPS_2022-04-29"              "MoNA_2022-04-27"              "Michael_Snyder_RPLC_20220424"
#[7] "MassBank_2022-04-27"          "HMDB_2022-04-11"     
variable_id <-
  annotation_table$variable_id %>%
  unique()

database_list <-
  list(
    "gnps" = gnps_ms2,
    "hmdb" = hmdb_ms2,
    "massbank" = massbank_ms2,
    "metlin" = metlin_ms2,
    "mona" = mona_ms2,
    "mpsnyder_rplc" = mpsnyder_rplc_ms2,
    "nist" = nist_ms2
  )

names(database_list) <-
  c(
    paste(
      gnps_ms2@database.info$Source,
      gnps_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      hmdb_ms2@database.info$Source,
      hmdb_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      massbank_ms2@database.info$Source,
      massbank_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      metlin_ms2@database.info$Source,
      metlin_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      mona_ms2@database.info$Source,
      mona_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      mpsnyder_rplc_ms2@database.info$Source,
      mpsnyder_rplc_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      nist_ms2@database.info$Source,
      nist_ms2@database.info$Version,
      sep = "_"
    )
  )


# --- Define tolerances
ms1.match.ppm <- 25
ms2.match.ppm <- 30

# --- Pick database
database_list <- list(
  "NIST_20220425" = nist_ms2,
  "METLIN_20220425" = metlin_ms2,
  "GNPS_2022-04-29" = gnps_ms2,
  "HMDB_2022-04-11" = hmdb_ms2,
  "MoNA_2022-04-27" = mona_ms2,
  "Michael_Snyder_RPLC_20220424" = mpsnyder_rplc_ms2,
  "MassBank_2022-04-27" = massbank_ms2
)

#create main folder
main_outdir <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/ms2_matching_plots"
dir.create(main_outdir, showWarnings = FALSE)

# --- Extract MS2 data
annotation_table <- object_pos4_ms2@annotation_table
ms2_data <- object_pos4_ms2@ms2_data
dataset_names <- names(ms2_data)

# --- Find dataset containing the target variable
for (target_id in unique(annotation_table$variable_id)) {
  # Find the MS2 object containing this variable
  ms2_obj <- NULL
  ms2_files_id <- NULL   ########## ADD THIS LINE
  
  for (dataset_name in dataset_names) {
    if (target_id %in% ms2_data[[dataset_name]]@variable_id) {
      ms2_obj <- ms2_data[[dataset_name]]
      ms2_files_id <- dataset_name   ########## ADD THIS LINE
      break
    }
  }
  
  if (is.null(ms2_obj)) {
    message("No MS2 data found for ", target_id)
    next
  }
  
  idx <- which(ms2_obj@variable_id == target_id)
  spectrum1 <- ms2_obj@ms2_spectra[[idx]]
  
  # create one folder per variable_id
  var_dir <- file.path(main_outdir, target_id)
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get all annotation matches for this variable
  annots <- annotation_table %>% filter(variable_id == target_id)
  
  for (i in seq_len(nrow(annots))) {
    db_name <- annots$Database[i]
    lab_id <- annots$Lab.ID[i]
    ce_value <- annots$CE[i]
    
    if (!db_name %in% names(database_list)) {
      message("Database ", db_name, " not found in list, skipping.")
      next
    }
    
    db_object <- database_list[[db_name]]
    
    # Get library MS2 spectrum
    spectrum2 <- tryCatch({
      get_ms2_spectrum(
        lab.id = lab_id,
        polarity = "positive",
        database = db_object,
        ce = ce_value
      )
    }, error = function(e) NULL)
    
    if (is.null(spectrum1) & is.null(spectrum2)) {
      message("No spectra available for ", target_id, " in ", db_name)
      next
    }
    
    # Generate mirror plot
    p <- masstools::ms2_plot(
      spectrum1 = spectrum1,
      spectrum2 = spectrum2,
      spectrum1_name = paste0(target_id),
      spectrum2_name = paste0(annots$Compound.name[i]),
      ppm.tol = ms1.match.ppm,
      mz.ppm.thr = ms2.match.ppm,
      interactive_plot = FALSE
    )
    
    ############## ADD THIS SECTION #################
    # If you want to include a simple ms2_spectrum_id label
    # you can construct it using the target_id (since that’s unique)
    ms2_spectrum_id <- paste0("mz", target_id)
    
    # Build the full annotation text
    info_text <- paste0(
      "variable_id: ", target_id, "\n",
      "ms2_files_id: ", ms2_files_id, "\n",
      "ms2_spectrum_id: ", ms2_spectrum_id, "\n",
      "Compound.name: ", annots$Compound.name[i], "\n",
      "CAS.ID: ", annots$CAS.ID[i], "\n",
      "HMDB.ID: ", annots$HMDB.ID[i], "\n",
      "KEGG.ID: ", annots$KEGG.ID[i], "\n",
      "Lab.ID: ", lab_id, "\n",
      "Adduct: ", annots$Adduct[i], "\n",
      "mz.error: ", annots$mz.error[i], "\n",
      "mz.match.score: ", annots$mz.match.score[i], "\n",
      "RT.error: ", annots$RT.error[i], "\n",
      "RT.match.score: ", annots$RT.match.score[i], "\n",
      "CE: ", ce_value, "\n",
      "SS: ", annots$SS[i], "\n",
      "Total.score: ", annots$Total.score[i], "\n",
      "Database: ", db_name, "\n",
      "Level: ", annots$Level[i]
    )
    ###############################################
    
    p <- p + annotate("text", x = -Inf, y = Inf, label = info_text, hjust = 0, vjust = 1, size = 4)
    
    ggsave(
      filename = file.path(var_dir, paste0(target_id, "_", lab_id, "_", db_name, "_ms2_plot_", ".pdf")),
      plot = p,
      width = 8,
      height = 6
    )
    
    message("✅ Saved: ", target_id, " (", lab_id, ") from ", db_name)
  }
}

#check and delete all empty folders
# List all folders in the directory
folder_list <- list.dirs("3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/ms2_matching_plots", recursive = FALSE)

# Function to check if the folder is empty
is_empty_folder <- function(folder) {
  # Check if the folder has files or subfolders
  files_in_folder <- list.files(folder)
  return(length(files_in_folder) == 0)
}
total_folders <-length(folder_list)
empty_folders <- 0

for (folder in folder_list) {
  if (is_empty_folder(folder)) {
    empty_folders <- empty_folders + 1  # Increment empty folder count
  }
}

# Print results
cat("Total number of folders:", total_folders, "\n") #1681
cat("Number of empty folders:", empty_folders, "\n") #1534

delete_empty <- readline(prompt = "Do you want to delete empty folders? (yes/no): ")

if (tolower(delete_empty) == "yes") {
  # Loop through each folder and delete if it's empty
  for (folder in folder_list) {
    if (is_empty_folder(folder)) {
      cat("Deleting empty folder:", folder, "\n")
      unlink(folder, recursive = TRUE)  # Remove the empty folder
    }
  }
}
############################################################################################################

load("2_data/preliminary_determination_bacterial_metabolomics/MS2/NEG/object_neg4_ms2")

# Step 3: Create output folder for plots
output_dir <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/NEG/ms2_matching_plots"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Step 4: Load all MS2 database objects
load("~/compound_database_shenlab/MS2/gnps_ms2.rda")
load("~/compound_database_shenlab/MS2/hmdb_ms2.rda")
load("~/compound_database_shenlab/MS2/massbank_ms2.rda")
load("~/compound_database_shenlab/MS2/metlin_ms2.rda")
load("~/compound_database_shenlab/MS2/mona_ms2.rda")
load("~/compound_database_shenlab/MS2/mpsnyder_rplc_ms2.rda")
load("~/compound_database_shenlab/MS2/nist_ms2.rda")

library(dplyr)
object_neg4_ms2@annotation_table %>%
  count(Database, name = "num_matches") %>%
  arrange(desc(num_matches))
#Database                     num_matches
#1 Shen-Lab_0.0.1                      2832
#2 NIST_20220425                          8
#3 Michael_Snyder_RPLC_20220424           5
#4 METLIN_20220425                        4
#5 MoNA_2022-04-27                        4
#6 HMDB_2022-04-11                        1

annotation_table <- object_neg4_ms2@annotation_table

database <-
  unique(annotation_table$Database)
database
#[1] "Shen-Lab_0.0.1"               "METLIN_20220425"              "NIST_20220425"               
#[4] "MoNA_2022-04-27"              "Michael_Snyder_RPLC_20220424" "HMDB_2022-04-11"     
variable_id <-
  annotation_table$variable_id %>%
  unique()

database_list <-
  list(
    "gnps" = gnps_ms2,
    "hmdb" = hmdb_ms2,
    "massbank" = massbank_ms2,
    "metlin" = metlin_ms2,
    "mona" = mona_ms2,
    "mpsnyder_rplc" = mpsnyder_rplc_ms2,
    "nist" = nist_ms2
  )

names(database_list) <-
  c(
    paste(
      gnps_ms2@database.info$Source,
      gnps_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      hmdb_ms2@database.info$Source,
      hmdb_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      massbank_ms2@database.info$Source,
      massbank_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      metlin_ms2@database.info$Source,
      metlin_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      mona_ms2@database.info$Source,
      mona_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      mpsnyder_rplc_ms2@database.info$Source,
      mpsnyder_rplc_ms2@database.info$Version,
      sep = "_"
    ),
    paste(
      nist_ms2@database.info$Source,
      nist_ms2@database.info$Version,
      sep = "_"
    )
  )


# --- Define tolerances
ms1.match.ppm <- 25
ms2.match.ppm <- 30

# --- Pick database
database_list <- list(
  "NIST_20220425" = nist_ms2,
  "METLIN_20220425" = metlin_ms2,
  "GNPS_2022-04-29" = gnps_ms2,
  "HMDB_2022-04-11" = hmdb_ms2,
  "MoNA_2022-04-27" = mona_ms2,
  "Michael_Snyder_RPLC_20220424" = mpsnyder_rplc_ms2,
  "MassBank_2022-04-27" = massbank_ms2
)

#create main folder
main_outdir <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/NEG/ms2_matching_plots"
dir.create(main_outdir, showWarnings = FALSE)

# --- Extract MS2 data
annotation_table <- object_neg4_ms2@annotation_table
ms2_data <- object_neg4_ms2@ms2_data
dataset_names <- names(ms2_data)

# --- Find dataset containing the target variable
for (target_id in unique(annotation_table$variable_id)) {
  # Find the MS2 object containing this variable
  ms2_obj <- NULL
  ms2_files_id <- NULL   ########## ADD THIS LINE
  
  for (dataset_name in dataset_names) {
    if (target_id %in% ms2_data[[dataset_name]]@variable_id) {
      ms2_obj <- ms2_data[[dataset_name]]
      ms2_files_id <- dataset_name   ########## ADD THIS LINE
      break
    }
  }
  
  if (is.null(ms2_obj)) {
    message("No MS2 data found for ", target_id)
    next
  }
  
  idx <- which(ms2_obj@variable_id == target_id)
  spectrum1 <- ms2_obj@ms2_spectra[[idx]]
  
  # create one folder per variable_id
  var_dir <- file.path(main_outdir, target_id)
  dir.create(var_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get all annotation matches for this variable
  annots <- annotation_table %>% filter(variable_id == target_id)
  
  for (i in seq_len(nrow(annots))) {
    db_name <- annots$Database[i]
    lab_id <- annots$Lab.ID[i]
    ce_value <- annots$CE[i]
    
    if (!db_name %in% names(database_list)) {
      message("Database ", db_name, " not found in list, skipping.")
      next
    }
    
    db_object <- database_list[[db_name]]
    
    # Get library MS2 spectrum
    spectrum2 <- tryCatch({
      get_ms2_spectrum(
        lab.id = lab_id,
        polarity = "negative",
        database = db_object,
        ce = ce_value
      )
    }, error = function(e) NULL)
    
    if (is.null(spectrum1) & is.null(spectrum2)) {
      message("No spectra available for ", target_id, " in ", db_name)
      next
    }
    
    # Generate mirror plot
    p <- masstools::ms2_plot(
      spectrum1 = spectrum1,
      spectrum2 = spectrum2,
      spectrum1_name = paste0(target_id),
      spectrum2_name = paste0(annots$Compound.name[i]),
      ppm.tol = ms1.match.ppm,
      mz.ppm.thr = ms2.match.ppm,
      interactive_plot = FALSE
    )
    
    ############## ADD THIS SECTION #################
    # If you want to include a simple ms2_spectrum_id label
    # you can construct it using the target_id (since that’s unique)
    ms2_spectrum_id <- paste0("mz", target_id)
    
    # Build the full annotation text
    info_text <- paste0(
      "variable_id: ", target_id, "\n",
      "ms2_files_id: ", ms2_files_id, "\n",
      "ms2_spectrum_id: ", ms2_spectrum_id, "\n",
      "Compound.name: ", annots$Compound.name[i], "\n",
      "CAS.ID: ", annots$CAS.ID[i], "\n",
      "HMDB.ID: ", annots$HMDB.ID[i], "\n",
      "KEGG.ID: ", annots$KEGG.ID[i], "\n",
      "Lab.ID: ", lab_id, "\n",
      "Adduct: ", annots$Adduct[i], "\n",
      "mz.error: ", annots$mz.error[i], "\n",
      "mz.match.score: ", annots$mz.match.score[i], "\n",
      "RT.error: ", annots$RT.error[i], "\n",
      "RT.match.score: ", annots$RT.match.score[i], "\n",
      "CE: ", ce_value, "\n",
      "SS: ", annots$SS[i], "\n",
      "Total.score: ", annots$Total.score[i], "\n",
      "Database: ", db_name, "\n",
      "Level: ", annots$Level[i]
    )
    ###############################################
    
    p <- p + annotate("text", x = -Inf, y = Inf, label = info_text, hjust = 0, vjust = 1, size = 4)
    
    ggsave(
      filename = file.path(var_dir, paste0(target_id, "_", lab_id, "_", db_name, "_ms2_plot_", ".pdf")),
      plot = p,
      width = 8,
      height = 6
    )
    
    message("✅ Saved: ", target_id, " (", lab_id, ") from ", db_name)
  }
}

#check and delete all empty folders
# List all folders in the directory
folder_list <- list.dirs("3_data_analysis/preliminary_determination_bacterial_metabolomics/NEG/ms2_matching_plots", recursive = FALSE)

# Function to check if the folder is empty
is_empty_folder <- function(folder) {
  # Check if the folder has files or subfolders
  files_in_folder <- list.files(folder)
  return(length(files_in_folder) == 0)
}
total_folders <-length(folder_list)
empty_folders <- 0

for (folder in folder_list) {
  if (is_empty_folder(folder)) {
    empty_folders <- empty_folders + 1  # Increment empty folder count
  }
}

# Print results
cat("Total number of folders:", total_folders, "\n") #683
cat("Number of empty folders:", empty_folders, "\n") #669

delete_empty <- readline(prompt = "Do you want to delete empty folders? (yes/no): ")

if (tolower(delete_empty) == "yes") {
  # Loop through each folder and delete if it's empty
  for (folder in folder_list) {
    if (is_empty_folder(folder)) {
      cat("Deleting empty folder:", folder, "\n")
      unlink(folder, recursive = TRUE)  # Remove the empty folder
    }
  }
}
