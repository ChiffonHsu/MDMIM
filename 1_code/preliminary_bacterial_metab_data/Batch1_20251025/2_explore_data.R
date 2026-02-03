library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(tidymass)
library(tidyverse)
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/Result/object")
object_pos <- object
object_pos
sample_info_pos <- readxl::read_xlsx("2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/sample_info_lpos.xlsx")
head(sample_info_pos)

#1. Add sample_info_pos to object_pos:
  object_pos %>% 
  extract_sample_info() %>% 
  head()
object_pos <- 
  object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::select(-c("group", "class", "injection.order"))
object_pos =
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  left_join(sample_info_pos,
            by = "sample_id")
object_pos %>% 
  extract_sample_info() %>% 
  head()
#2. Save the object_pos in a new folder named “data_cleaning”
dir.create("2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/data_cleaning/", showWarnings = FALSE, recursive = TRUE)
save(object_pos, file = "2_data/preliminary_determination_bacterial_metabolomics/20251025/MS1/POS/Result/data_cleaning/object_pos")
object_pos
dim(object_pos)
# variables   samples 
# 32011        70 
object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count()
object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count(group) #70
object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count(experiment)
# experiment  n
# 1   20251025 70
#3. Generate peak distribution plots:
  object_pos %>%
  `+`(1) %>% 
  log(10) %>%
  show_mz_rt_plot() +
  scale_size_continuous(range = c(0.01, 2))

#4. Explore missing values (mvs):
  get_mv_number(object = object_pos)  #1206278
#5. Get missing number in each samples:
  get_mv_number(object = object_pos, by = "sample") %>% 
  head()
  # Blank_1_POS (2). Blank_10_POS (2). Blank_11_POS (2).  Blank_2_POS (2). 
  # 27002             26175             25789             25449 
  # Blank_3_POS (2).  Blank_4_POS (2). 
  # 25103             25882 
#6. Get missing number in each variables:
  get_mv_number(object = object_pos, by = "variable") %>% 
  head()
  # M50T125_POS M51T125_POS  M51T35_POS  M52T34_POS  M53T68_POS M53T126_POS 
  # 21          21          17          39          18          21
#7. Show missing value information:
  show_missing_values(object = object_pos, show_column_names = FALSE, percentage = TRUE)

#8. Show missing value in samples:
  show_sample_missing_values(object = object_pos, percentage = TRUE)

#9. Show missing value in variables:
  show_variable_missing_values(
    object = object_pos,
    percentage = TRUE,
    show_x_text = FALSE,
    show_x_ticks = FALSE
  ) +
  scale_size_continuous(range = c(0.01, 1))
##############################################################################################################################
  load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result/object")
  object_neg <- object
  object_neg
  sample_info_neg <- readxl::read_xlsx("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/sample_info_lneg.xlsx")
  head(sample_info_neg)
  
  #1. Add sample_info_pos to object_neg:
  object_neg %>% 
    extract_sample_info() %>% 
    head()
  object_neg <- 
    object_neg %>% 
    activate_mass_dataset(what = "sample_info") %>%
    dplyr::select(-c("group", "class", "injection.order"))
  object_neg =
    object_neg %>%
    activate_mass_dataset(what = "sample_info") %>%
    left_join(sample_info_neg,
              by = "sample_id")
  object_neg %>% 
    extract_sample_info() %>% 
    head()
  #2. Save the object_neg in a new folder named “data_cleaning”
  dir.create("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/data_cleaning", showWarnings = FALSE, recursive = TRUE)
  save(object_neg, file = "2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/data_cleaning/object_neg")
  object_neg
  dim(object_neg)
  object_neg %>% 
    activate_mass_dataset(what = "sample_info") %>% 
    dplyr::count()
  object_neg %>% 
    activate_mass_dataset(what = "sample_info") %>% 
    dplyr::count(group)
  object_neg %>% 
    activate_mass_dataset(what = "sample_info") %>% 
    dplyr::count(experiment)
  #3. Generate peak distribution plots:
  object_neg %>%
    `+`(1) %>% 
    log(10) %>%
    show_mz_rt_plot() +
    scale_size_continuous(range = c(0.01, 2))
  
  #4. Explore missing values (mvs):
  get_mv_number(object = object_neg)  #31166
  #5. Get missing number in each samples:
  get_mv_number(object = object_neg, by = "sample") %>% 
    head()
  #6. Get missing number in each variables:
  get_mv_number(object = object_neg, by = "variable") %>% 
    head()
  #7. Show missing value information:
  show_missing_values(object = object_neg, show_column_names = FALSE, percentage = TRUE)
  
  #8. Show missing value in samples:
  show_sample_missing_values(object = object_neg, percentage = TRUE)
  
  #9. Show missing value in variables:
  show_variable_missing_values(
    object = object_neg,
    percentage = TRUE,
    show_x_text = FALSE,
    show_x_ticks = FALSE
  ) +
    scale_size_continuous(range = c(0.01, 1))
  