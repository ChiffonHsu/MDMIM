library(r4projects)
setwd(get_project_wd())
rm(list = ls())

load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/data_cleaning/object_pos")
  #load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/data_cleaning/object_neg")
object_pos
  #object_neg
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
                      path = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/data_quality_before_data_cleaning")

#2. Remove noisy metabolic features: Remove variables which have MVs in more than 20% QC samples and in at least 50% samples in control group or case group.
object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count(group)
    # group  n
    # 1            Blank 11
    # 2               QC 15
    # 3 bacteria_samples 40
    # 4   medium_samples  4
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
bacteria_samples_id = 
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "bacteria_samples") %>%
  pull(sample_id)
Citrobacter_youngae_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Citrobacter youngae")) %>%
  pull(sample_id)
Citrobacter_braakii_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Citrobacter braakii")) %>%
  pull(sample_id)
Bacteroides_fragilis_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Bacteroides fragilis")) %>%
  pull(sample_id)
Enterococcus_casseliflavus_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Enterococcus casseliflavus")) %>%
  pull(sample_id)
Enterobacter_ludwigii_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Enterobacter ludwigii")) %>%
  pull(sample_id)
Pediococcus_acidilactici <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Pediococcus acidilactici")) %>%
  pull(sample_id)
Bacteroides_finegoldii_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Bacteroides finegoldii")) %>%
  pull(sample_id)
Enterococcus_faecalis_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Enterococcus faecalis")) %>%
  pull(sample_id)
Enterococcus_durans_id <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(Bacteria %in% c("Enterococcus durans")) %>%
  pull(sample_id)
MGAM_6h <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class %in% c("Control")) %>%
  pull(sample_id)
    
object_pos =
  object_pos %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>%
  mutate_variable_na_freq(according_to_samples = blank_id) %>%
  mutate_variable_na_freq(according_to_samples = bacteria_samples_id) %>%
  mutate_variable_na_freq(according_to_samples = Citrobacter_youngae_id) %>%
  mutate_variable_na_freq(according_to_samples = Citrobacter_braakii_id) %>%
  mutate_variable_na_freq(according_to_samples = Bacteroides_fragilis_id) %>%
  mutate_variable_na_freq(according_to_samples = Enterococcus_casseliflavus_id) %>%
  mutate_variable_na_freq(according_to_samples = Enterobacter_ludwigii_id) %>%
  mutate_variable_na_freq(according_to_samples = Pediococcus_acidilactici) %>%
  mutate_variable_na_freq(according_to_samples = Bacteroides_finegoldii_id) %>%
  mutate_variable_na_freq(according_to_samples = Enterococcus_faecalis_id) %>%
  mutate_variable_na_freq(according_to_samples = Enterococcus_durans_id) %>%
  mutate_variable_na_freq(according_to_samples = MGAM_6h) 
head(extract_variable_info(object_pos))
    # variable_id       mz        rt   na_freq na_freq.1 na_freq.2 na_freq.3 na_freq.4
    # 1 M50T125_POS 50.01472 125.31911 0.4666667         1     0.075         0      0.00
    # 2 M51T125_POS 51.02256 125.32843 0.4666667         1     0.075         0      0.00
    # 3  M51T35_POS 51.05700  35.48802 0.4000000         1     0.000         0      0.00
    # 4  M52T34_POS 52.06481  34.49569 0.6000000         1     0.475         0      0.25
    # 5  M53T68_POS 53.00180  68.34187 0.4000000         1     0.025         0      0.00
    # 6 M53T126_POS 53.00178 125.53596 0.4666667         1     0.075         0      0.00
    # na_freq.5 na_freq.6 na_freq.7 na_freq.8 na_freq.9 na_freq.10 na_freq.11 na_freq.12
    # 1         0      0.00      0.00      0.00      0.00       0.25       0.50          0
    # 2         0      0.00      0.00      0.00      0.25       0.00       0.50          0
    # 3         0      0.00      0.00      0.00      0.00       0.00       0.00          0
    # 4         1      0.75      0.50      0.75      0.25       0.50       0.25          0
    # 5         0      0.25      0.00      0.00      0.00       0.00       0.00          0
    # 6         0      0.25      0.25      0.00      0.00       0.00       0.00          0

# Assuming you have already applied mutate_variable_na_freq for each sample
# and stored them in object_pos

# Select variables (columns) related to the 'na_freq' values
# Here, we are selecting `variable_id` and all columns starting with 'na_freq'
na_freq_data <- object_pos@variable_info %>% 
  dplyr::select(variable_id, starts_with("na_freq"))

# Check the first few rows to ensure correct extraction
head(na_freq_data)

# Reshape the data into a long format, where 'na_freq' values are grouped by sample
na_freq_long <- na_freq_data %>%
  gather(key = "sample_group", value = "na_freq", -variable_id)

na_freq_long$sample_group <- factor(na_freq_long$sample_group, 
                                    levels = colnames(na_freq_data)[-1])  # Exclude 'variable_id' column

# Check the structure to verify the changes
str(na_freq_long)
# Visualize the missing values across different sample groups
library(ggplot2)

ggplot(na_freq_long, aes(x = sample_group, y = na_freq)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for readability
  labs(title = "Distribution of Missing Values (NA Frequency) Across Sample Groups", 
       x = "Sample Group", 
       y = "Missing Value Proportion (na_freq)") +
  theme_minimal()




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
  # sample_id experiment sample_type Bacteria medium class batch group
  # 1         0          0           0        0      0     0     0     0
# Print the missing details for each variable in the Bacteria group
print(bacteria_missing_details)
qc_missing_details <- object_pos %>%
  filter(group == "qc_samples") %>%
  summarise(across(everything(), ~ sum(is.na(.))))  # Count missing values for each variable

# Print the missing details for each variable in the QC group
print(qc_missing_details)
  # sample_id experiment sample_type Bacteria medium class batch group
  # 1         0          0           0        0      0     0     0     0

#4. Remove variables:
object_pos <- object_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(
    na_freq < 0.2 |           # QC good
      na_freq.2 < 0.8  # detected in medium → stable and real
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
#791120
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
save(object_pos2, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/data_cleaning/object_pos2")

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

