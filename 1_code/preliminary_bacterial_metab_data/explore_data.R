library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(tidymass)
library(tidyverse)
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/Result/object")
object_pos <- object
object_pos
sample_info_pos <- readxl::read_xlsx("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/sample_info_lpos.xlsx")
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
dir.create("data_cleaning/POS", showWarnings = FALSE, recursive = TRUE)
save(object_pos, file = "data_cleaning/POS/object_pos")
object_pos
dim(object_pos)
object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count()
object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count(group)
object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::count(experiment)
#3. Generate peak distribution plots:
  object_pos %>%
  `+`(1) %>% 
  log(10) %>%
  show_mz_rt_plot() +
  scale_size_continuous(range = c(0.01, 2))

#4. Explore missing values (mvs):
  get_mv_number(object = object_pos)  #72124
#5. Get missing number in each samples:
  get_mv_number(object = object_pos, by = "sample") %>% 
  head()
#6. Get missing number in each variables:
  get_mv_number(object = object_pos, by = "variable") %>% 
  head()
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
  dir.create("data_cleaning/NEG", showWarnings = FALSE, recursive = TRUE)
  save(object_neg, file = "data_cleaning/NEG/object_neg")
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
  