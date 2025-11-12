load("MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/object_pos")

object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count(sample_type)

show_variable_missing_values(object = object_pos %>% 
                               activate_mass_dataset(what = "sample_info") %>% 
                               filter(class == "media_blank"), 
                             percentage = TRUE) +
  scale_size_continuous(range = c(0.01, 2))

mix_id =
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(sample_type == "mix") %>%
  pull(sample_id)

media_blank_id =
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(sample_type == "media_blank") %>%
  pull(sample_id)

supernatant_id =
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(sample_type == "supernatant") %>%
  pull(sample_id)

object_pos =
  object_pos %>%
  #no mix samples, in diff folder, mutate_variable_na_freq(according_to_samples = mix_id) %>%
  mutate_variable_na_freq(according_to_samples = media_blank_id) %>%
  mutate_variable_na_freq(according_to_samples = supernatant_id)

head(extract_variable_info(object_pos))
                         
object_pos <- 
  object_pos %>% 
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2 & (na_freq < 0.5 | na_freq.1 < 0.5))
object_pos


#Detect outlier samples
load("MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/object_pos")
massdataset::show_sample_missing_values(object = object_pos,
                                        color_by = "sample_type",
                                        order_by = "injection.order",
                                        percentage = TRUE) +
  theme(axis.text.x = element_text(size = 2)) +
  scale_size_continuous(range = c(0.1, 2)) +
  ggsci::scale_color_aaas()

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

#Missing values imputation
load("MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/object_pos")
library(massqc)

get_mv_number(object_pos)
#> [1] 148965
library(tidymass)

library(massstat)
object_pos <- 
  impute_mv(object = object_pos, method = "knn")

#If there is any column with too many NAs, remove this column before imputating
# Calculate missing value frequencies first
object_pos <- object_pos %>%
  mutate_variable_na_freq()

# Filter variables with high missingness and impute with minimum
library(tidymass)
library(dplyr)

# 1. Compute missing value frequency
object_pos <- object_pos %>% mutate_variable_na_freq()

# 2. Subset and impute sparse variables (>= 80% missing)
object_pos_sparse <- object_pos %>%
  activate_mass_dataset("variable_info") %>%
  filter(na_freq >= 0.8) %>%
  impute_mv(method = "minimum")

# 3. Subset and impute dense variables (< 80% missing)
object_pos_dense <- object_pos %>%
  activate_mass_dataset("variable_info") %>%
  filter(na_freq < 0.8) %>%
  impute_mv(method = "knn")

object_pos_imputed <- combine_mass_dataset(object_pos_sparse, object_pos_dense)

# Combine back
object_pos <- combine_mass_dataset(object_pos_sparse, object_pos_dense)

#OR
# Impute with two steps:
  #1. For features with >80% missing, replace with min/median.
  #2. For the rest, use KNN.
# Step 1: fill very sparse vars with min
object_neg_sparse <- object_neg %>%
  activate_mass_dataset("variable_info") %>%
  filter(na_freq >= 0.8) %>%
  extract_mass_dataset() %>%
  impute_mv(method = "minimum")

object_neg_sparse <- object_neg %>%
  activate_mass_dataset("variable_info") %>%
  filter(na_freq >= 0.8) %>%
  impute_mv(method = "minimum")


# Step 2: use knn for the rest
object_neg_dense <- object_neg %>%
  activate_mass_dataset("variable_info") %>%
  filter(na_freq < 0.8) %>%
  extract_mass_dataset() %>%
  impute_mv(method = "knn")

# Merge back
object_neg <- combine_mass_dataset(object_neg_sparse, object_neg_dense)

remove.packages(c("tidymass", "masstools", "massdataset", "massqc", "massstat"))
remove.packages(
  c("tidymass", "masstools", "massdataset", "massqc", "massstat"),
  lib = "/usr/lib64/R/library"
)

##########WATCH OUT##############TOO MANY MISSING VALUES, NONE OF THE ABOVE WORKED!!!!!
#######ALL THE ABOVE DIDNT WORK, THIS ALSO DIDN WORK FOR DATA NORMALISATION SO PLS DUN USE THIS:#######
#load("MS2Library_ChiFang/In_Vitro/RPLC/POS/Result/data_cleaning/object_pos")
# Load necessary libraries
#library(tidymass)
#library(dplyr)
#library(impute)  # for KNN imputation
# 1. Compute missing value frequency for all variables
#object_pos <- object_pos %>% mutate_variable_na_freq()
# 2. Extract expression data as a data frame
#expr <- object_pos@expression_data  # directly access the matrix
#expr_df <- as.data.frame(expr)
# 3. Get variable IDs and their missing value frequencies
#na_freq <- colSums(is.na(expr_df)) / nrow(expr_df)
#summary(na_freq)
# 4. Separate sparse vs dense variables
 # Sparse variables (≥ 80% missing)
#parse_vars <- names(na_freq[na_freq >= 0.8])
# Dense variables (< 80% missing)
#dense_vars <- names(na_freq[na_freq < 0.8])
# 5. Impute sparse variables with minimum value
# expr_df[sparse_vars] <- lapply(expr_df[sparse_vars], function(x) {
  #  +     x[is.na(x)] <- min(x, na.rm = TRUE)
 #   +     return(x)
#    + })
 # 6. Impute dense variables using KNN
#library(impute)
#dense_matrix <- as.matrix(expr_df[dense_vars])
#dense_imputed <- impute.knn(dense_matrix)$data
# 7. Result: fully imputed expression matrix
 #expr_df[dense_vars] <- dense_imputed #STILL WONT WORK, USE THE ONE BELOW:

#THE ABOVE DOESNTY WORK FOR DATA NORMALISATION, SO HERE IS THE NEW ONE:
library(tidymass)
library(dplyr)
library(impute)

# 1. Load object
load("MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/object_pos")

# 2. Compute NA frequency
object_pos <- object_pos %>% mutate_variable_na_freq()

# 3. Extract expression data as a data.frame
expr_df <- as.data.frame(object_pos@expression_data)

# 4. Identify sparse and dense variables
na_freq <- colSums(is.na(expr_df)) / nrow(expr_df)
sparse_vars <- names(na_freq[na_freq >= 0.8])
dense_vars  <- names(na_freq[na_freq < 0.8])

# 5. Impute sparse variables (minimum)
expr_df[sparse_vars] <- lapply(expr_df[sparse_vars], function(x) {
  x[is.na(x)] <- min(x, na.rm = TRUE)
  x
})

# 6. Impute dense variables (KNN)
dense_matrix <- as.matrix(expr_df[dense_vars])
dense_imputed <- impute.knn(dense_matrix)$data
expr_df[dense_vars] <- dense_imputed

# 7. Assign back to object (data.frame!)
object_pos@expression_data <- expr_df

# 8. Normalize
object_pos <- normalize_data(object_pos, method = "median")

# 9. Integrate by subject
object_pos2 <- integrate_data(object_pos, method = "subject_median")


#Data normalisation and integration
load("MS2Library_ChiFang/In_Vitro/RPLC/NEG/Result/data_cleaning/object_neg")
object_pos <- 
  normalize_data(object_pos, method = "median")

object_pos2 <- 
  integrate_data(object_pos, method = "subject_median")

#The following didnt work, check whats in the sample_info first
#object_neg2 %>% 
#  `+`(1) %>% 
#  log(2) %>% 
#  massqc::massqc_pca(color_by = "batch", line = FALSE)

#There's no columns with "batch", so take the column called "experiment" as batch
   # Copy 'experiment' to a new 'batch' column
object_pos2@sample_info$batch <- object_pos2@sample_info$experiment

# Now you can run PCA colored by batch
library(massqc)
object_pos2@sample_info$batch <- as.factor(object_pos2@sample_info$batch)
object_pos2 %>%
  log(2) %>%  # log-transform if needed
  massqc_pca(color_by = "batch", line = FALSE)
#save(object_neg2, file = "data_cleaning/POS/object_pos2")

#FOR NEG DATASET
library(tidymass)
library(dplyr)
library(massqc)

# Use the already integrated object_neg2

# 1. Extract expression data
expr <- object_neg2@expression_data

# 2. Replace zero or negative values with a small positive number
expr[] <- lapply(expr, function(x) {
  x[x <= 0 | is.infinite(x) | is.na(x)] <- min(x[x > 0], na.rm = TRUE)
  x
})

# 3. Assign back to object
object_neg2@expression_data <- as.data.frame(expr)

# 4. Now do log2 and PCA safely
object_neg2 %>%
  log(2) %>%
  massqc_pca(color_by = "batch", line = FALSE)

