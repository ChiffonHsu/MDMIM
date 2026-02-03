library(r4projects)
setwd(get_project_wd())
rm(list = ls())

load("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/data_cleaning/object_pos")
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/data_cleaning/object_neg")
object_pos
object_neg
#1. Data quality assessment before data cleaning:
object_pos <- 
  object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::mutate(experiment = as.character(experiment))
object_neg <- 
  object_neg %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::mutate(experiment = as.character(experiment))

# #Get html format quality assessment report (saved in folder “data_cleaning/POS/data_quality_before_data_cleaning/Report.”:
# object_pos@sample_info$group <- as.factor(object_pos@sample_info$group)
# object_pos@sample_info$sample_type <- as.factor(object_pos@sample_info$sample_type)

# # Directly modify the 'group' column in sample_info
# object_pos@sample_info$group <- gsub("bacterium_samples", "bacteria_samples", object_pos@sample_info$group)
# 
# # Verify the changes
# unique(object_pos@sample_info$group)

massqc::massqc_report(object = object_pos,
    path = "2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/data_quality_before_data_cleaning")

 #2. Remove noisy metabolic features: Remove variables which have MVs in more than 20% QC samples and in at least 50% samples in control group or case group.
   object_pos %>% 
   activate_mass_dataset(what = "sample_info") %>% 
   dplyr::count(group)
#3. MV percentage in QC samples: 
   show_variable_missing_values(object = object_pos %>% 
                                  activate_mass_dataset(what = "sample_info") %>% 
                                  filter(class == "QC"), 
                                percentage = TRUE) +
     scale_size_continuous(range = c(0.01, 2))
   qc_id =
     object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "QC") %>%
     pull(sample_id)
   blank_id =
     object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "Blank") %>%
     pull(sample_id)
   bacteria_id <- object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group %in% c("bacteria_samples")) %>%
     pull(sample_id)
   
   medium_id <- object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group %in% c("medium_samples")) %>%
     pull(sample_id)
   object_pos =
     object_pos %>%
     mutate_variable_na_freq(according_to_samples = qc_id) %>%
     mutate_variable_na_freq(according_to_samples = bacteria_id) %>%
     mutate_variable_na_freq(according_to_samples = medium_id) 
   object_pos =
     object_pos %>%
     mutate_variable_na_freq(according_to_samples = medium_id)
   head(extract_variable_info(object_pos))
   #do a preliminary check on t he diestribution of missing values in each groups
   # Check the missing value frequencies for Bacteria and Medium groups
   # Load necessary library
   library(ggplot2)
   
   # Create a function to calculate missing values for each variable in a specific group
   missing_values_histogram <- function(data, group_column, group_name) {
     missing_data <- data %>%
       filter(group_column == group_name) %>%  # Filter data by group
       summarise(across(everything(), ~ sum(is.na(.)))) %>%  # Count NAs for each variable
       gather(key = "variable", value = "missing_count")  # Reshape for plotting
     
     # Plot the histogram for missing values
     ggplot(missing_data, aes(x = missing_count)) +
       geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
       labs(title = paste("Missing Values Distribution in", group_name),
            x = "Number of Missing Values",
            y = "Frequency") +
       theme_minimal()
   }
   
   # Call the function to create histograms for Bacteria and Medium groups
   missing_values_histogram(object_pos, object_pos@sample_info$group, "bacteria_samples")
   missing_values_histogram(object_pos, object_pos@sample_info$group, "medium_samples")
   missing_values_histogram(object_pos, object_pos@sample_info$group, "QC")
   
   #CONTINUE CHECKING
   # Check missing values for specific variables in Medium samples
   medium_missing_details <- object_pos %>%
     filter(group == "medium_samples") %>%
     summarise(across(everything(), ~ sum(is.na(.))))  # Count missing values for each variable
   
   # Print the missing details for each variable in the Medium group
   print(medium_missing_details)
   bacteria_missing_details <- object_pos %>%
     filter(group == "bacteria_samples") %>%
     summarise(across(everything(), ~ sum(is.na(.))))  # Count missing values for each variable
   
   # Print the missing details for each variable in the Bacteria group
   print(bacteria_missing_details)
   qc_missing_details <- object_pos %>%
     filter(group == "qc_samples") %>%
     summarise(across(everything(), ~ sum(is.na(.))))  # Count missing values for each variable
   
   # Print the missing details for each variable in the QC group
   print(qc_missing_details)
   
   #4. Remove variables:
   object_pos <- object_pos %>%
     activate_mass_dataset(what = "variable_info") %>%
     filter(
       na_freq.1 < 0.2 |           # QC good
         na_freq.2 < 0.5       # detected in medium → stable and real
     )
    # object_pos <- 
  #   object_pos %>% 
  #   activate_mass_dataset(what = "variable_info") %>%
  #   filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5))
   object_pos
   #Check number of variables left
#5. Filter outlier sampels:
   massdataset::show_sample_missing_values(object = object_pos,
                                           color_by = "group",
                                           order_by = "injection.order",
                                           percentage = TRUE) +
     theme(axis.text.x = element_text(size = 2)) +
     scale_size_continuous(range = c(0.1, 2)) +
     ggsci::scale_color_aaas()                                            
                                             
#6. Detect outlier samples:
   outlier_samples =
     object_pos %>%
     `+`(1) %>% 
     log(2) %>%
     scale() %>%
     detect_outlier()
   outlier_samples
   outlier_table <-
     extract_outlier_table(outlier_samples)
   outlier_table %>% 
     head()
   outlier_table %>% 
     apply(1, function(x){
       sum(x)
     }) %>% 
     `>`(0) %>% 
     which()
# Check if there is any outlier samples left (should be 0)
##Troubleshooting: there is one QC sample with 100% missing value
   # Calculate the percentage of missing values per sample
   missing_percent_per_sample <- apply(object_pos@expression_data, 2, function(x) mean(is.na(x)) * 100)
   
   # Identify the sample with 100% missing values
   qc_sample_with_missing_values <- names(missing_percent_per_sample)[missing_percent_per_sample == 100]
   
   # Print the QC sample with 100% missing values
   print(qc_sample_with_missing_values) # "POS blk 1 (2)."
   
   # Remove the sample with 100% missing values from the object
   object_pos <- object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(!sample_id %in% qc_sample_with_missing_values)
   
   # Verify that the sample has been removed
   object_pos
   
   #7. missing value imputation
   get_mv_number(object_pos)
   #148779
   #object_pos <- 
  #   impute_mv(object = object_pos, method = "knn")
   
   #get_mv_number(object_pos) #0

   #8. data normalization
#   object_pos <- 
#     normalize_data(object_pos, method = "median")
      #imputation method: half minimum (manually)
   #(1)extract matrix
   expr <- object_pos@expression_data
   #(2)compute half-min values per feature (row)
   half_min_values <- apply(expr, 1, function(x) {
     min(x, na.rm = TRUE) / 2
   })
   #(3)replace NA with half-min of that feature
   for (i in 1:nrow(expr)) {
     expr[i, is.na(expr[i, ])] <- half_min_values[i]
   }
   #(4)put matrix back into massdataset
   object_pos@expression_data <- expr
     
   object_pos2 <- 
     integrate_data(object_pos, method = "subject_median")
   object_pos2 %>% 
     `+`(1) %>% 
     log(2) %>% 
     massqc::massqc_pca(color_by = "experiment", line = FALSE)
   save(object_pos2, file = "2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/data_cleaning/object_pos2")
   
################################################################################################################
   #1. Data quality assessment before data cleaning:
   #Get html format quality assessment report (saved in folder “data_cleaning/POS/data_quality_before_data_cleaning/ Report.”:
   load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/data_cleaning/object_neg")
   massqc::massqc_report(object = object_neg,
                         path = "3_data_analysis/preliminary_determination_bacterial_metabolomics/NEG/data_quality_before_data_cleaning")
   
   #2. Remove noisy metabolic features: Remove variables which have MVs in more than 20% QC samples and in at least 50% samples in control group or case group.
   object_neg %>% 
     activate_mass_dataset(what = "sample_info") %>% 
     dplyr::count(group)
   #3. MV percentage in QC samples: 
   show_variable_missing_values(object = object_neg %>% 
                                  activate_mass_dataset(what = "sample_info") %>% 
                                  filter(class == "QC"), 
                                percentage = TRUE) +
     scale_size_continuous(range = c(0.01, 2))
   qc_id =
     object_neg %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "QC") %>%
     pull(sample_id)
   blank_id =
     object_neg %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "Blank") %>%
     pull(sample_id)
   bacteria_id <- object_neg %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group %in% c("bacteria_samples")) %>%
     pull(sample_id)
   
   medium_id <- object_neg %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group %in% c("medium_samples")) %>%
     pull(sample_id)
   object_neg =
     object_neg %>%
     mutate_variable_na_freq(according_to_samples = qc_id) %>%
     mutate_variable_na_freq(according_to_samples = bacteria_id) %>%
     mutate_variable_na_freq(according_to_samples = medium_id) 
   object_neg =
     object_neg %>%
     mutate_variable_na_freq(according_to_samples = medium_id)
   head(extract_variable_info(object_neg))
   #do a preliminary check on t he diestribution of missing values in each groups
   # Check the missing value frequencies for Bacteria and Medium groups
   # Load necessary library
   library(ggplot2)
   
   # Create a function to calculate missing values for each variable in a specific group
   missing_values_histogram <- function(data, group_column, group_name) {
     missing_data <- data %>%
       filter(group_column == group_name) %>%  # Filter data by group
       summarise(across(everything(), ~ sum(is.na(.)))) %>%  # Count NAs for each variable
       gather(key = "variable", value = "missing_count")  # Reshape for plotting
     
     # Plot the histogram for missing values
     ggplot(missing_data, aes(x = missing_count)) +
       geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
       labs(title = paste("Missing Values Distribution in", group_name),
            x = "Number of Missing Values",
            y = "Frequency") +
       theme_minimal()
   }
   
   # Call the function to create histograms for Bacteria and Medium groups
   missing_values_histogram(object_neg, object_neg@sample_info$group, "bacteria_samples")
   missing_values_histogram(object_neg, object_neg@sample_info$group, "medium_samples")
   missing_values_histogram(object_neg, object_neg@sample_info$group, "QC")
   
   #CONTINUE CHECKING
   # Check missing values for specific variables in Medium samples
   medium_missing_details <- object_neg %>%
     filter(group == "medium_samples") %>%
     summarise(across(everything(), ~ sum(is.na(.))))  # Count missing values for each variable
   
   # Print the missing details for each variable in the Medium group
   print(medium_missing_details)
   bacteria_missing_details <- object_neg %>%
     filter(group == "bacteria_samples") %>%
     summarise(across(everything(), ~ sum(is.na(.))))  # Count missing values for each variable
   
   # Print the missing details for each variable in the Bacteria group
   print(bacteria_missing_details)
   qc_missing_details <- object_neg %>%
     filter(group == "qc_samples") %>%
     summarise(across(everything(), ~ sum(is.na(.))))  # Count missing values for each variable
   
   # Print the missing details for each variable in the QC group
   print(qc_missing_details)
   
   #4. Remove variables:
   object_neg <- object_neg %>%
     activate_mass_dataset(what = "variable_info") %>%
     filter(
       na_freq.1 < 0.2 |           # QC good
         na_freq.2 < 0.5       # detected in medium → stable and real
     )
   # object_neg <- 
   #   object_neg %>% 
   #   activate_mass_dataset(what = "variable_info") %>%
   #   filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5))
   object_neg
   #Check number of variables left
   #5. Filter outlier sampels:
   massdataset::show_sample_missing_values(object = object_neg,
                                           color_by = "group",
                                           order_by = "injection.order",
                                           percentage = TRUE) +
     theme(axis.text.x = element_text(size = 2)) +
     scale_size_continuous(range = c(0.1, 2)) +
     ggsci::scale_color_aaas()                                            
   
   #6. Detect outlier samples:
   outlier_samples =
     object_neg %>%
     `+`(1) %>% 
     log(2) %>%
     scale() %>%
     detect_outlier()
   outlier_samples
   outlier_table <-
     extract_outlier_table(outlier_samples)
   outlier_table %>% 
     head()
   outlier_table %>% 
     apply(1, function(x){
       sum(x)
     }) %>% 
     `>`(0) %>% 
     which()
   # Check if there is any outlier samples left (should be 0)
   ##Troubleshooting: there is one QC sample with 100% missing value
   # Calculate the percentage of missing values per sample
   missing_percent_per_sample <- apply(object_neg@expression_data, 2, function(x) mean(is.na(x)) * 100)
   
   # Identify the sample with 100% missing values
   qc_sample_with_missing_values <- names(missing_percent_per_sample)[missing_percent_per_sample == 100]
   
   # Print the QC sample with 100% missing values
   print(qc_sample_with_missing_values) # "POS blk 1 (2)."
   
   # Remove the sample with 100% missing values from the object
   object_neg <- object_neg %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(!sample_id %in% qc_sample_with_missing_values)
   
   # Verify that the sample has been removed
   object_neg
   
   #7. missing value imputation
   get_mv_number(object_neg)
   #54540
   #object_neg <- 
   #   impute_mv(object = object_neg, method = "knn")
   
   #get_mv_number(object_neg) #0
   
   #8. data normalization
   #   object_neg <- 
   #     normalize_data(object_neg, method = "median")
   #imputation method: half minimum (manually)
   #(1)extract matrix
   expr <- object_neg@expression_data
   #(2)compute half-min values per feature (row)
   half_min_values <- apply(expr, 1, function(x) {
     min(x, na.rm = TRUE) / 2
   })
   #(3)replace NA with half-min of that feature
   for (i in 1:nrow(expr)) {
     expr[i, is.na(expr[i, ])] <- half_min_values[i]
   }
   #(4)put matrix back into massdataset
   object_neg@expression_data <- expr
   
   object_neg2 <- 
     integrate_data(object_neg, method = "subject_median")
   object_neg2@sample_info$experiment <- 
     factor(object_neg2@sample_info$experiment)
   
     object_neg2 %>% 
     `+`(1) %>% 
     log(2) %>% 
     massqc::massqc_pca(
       color_by = "experiment",
       line = FALSE,
       ellipse = FALSE,
       point_size = 3,
       label = FALSE,
       order = FALSE
     )
   
   save(object_neg2, file = "2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/data_cleaning/object_neg2")
   