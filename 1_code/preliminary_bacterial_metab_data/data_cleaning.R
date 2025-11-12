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

#Get html format quality assessment report (saved in folder “data_cleaning/POS/data_quality_before_data_cleaning/Report.”:
object_pos@sample_info$group <- as.factor(object_pos@sample_info$group)
object_pos@sample_info$sample_type <- as.factor(object_pos@sample_info$sample_type)

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
     filter(class == "QC") %>%
     pull(sample_id)
   control_id =
     object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "Control") %>%
     pull(sample_id)
   case_id =
     object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "Case") %>%
     pull(sample_id)
   object_pos =
     object_pos %>%
     mutate_variable_na_freq(according_to_samples = qc_id) %>%
     mutate_variable_na_freq(according_to_samples = control_id) %>%
     mutate_variable_na_freq(according_to_samples = case_id)
   head(extract_variable_info(object_pos))  
   
#4. Remove variables:
   object_pos <- 
     object_pos %>% 
     activate_mass_dataset(what = "variable_info") %>%
     filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5))
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
   
################################################################################################################
   #1. Data quality assessment before data cleaning:
   #Get html format quality assessment report (saved in folder “data_cleaning/POS/data_quality_before_data_cleaning/ Report.”:
   massqc::massqc_report(object = object_pos,
                         path = "data_cleaning/POS/data_quality_before_data_cleaning")
   
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
     filter(class == "QC") %>%
     pull(sample_id)
   control_id =
     object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "Control") %>%
     pull(sample_id)
   case_id =
     object_pos %>%
     activate_mass_dataset(what = "sample_info") %>%
     filter(group == "Case") %>%
     pull(sample_id)
   object_pos =
     object_pos %>%
     mutate_variable_na_freq(according_to_samples = qc_id) %>%
     mutate_variable_na_freq(according_to_samples = control_id) %>%
     mutate_variable_na_freq(according_to_samples = case_id)
   head(extract_variable_info(object_pos))  
   
   #4. Remove variables:
   object_pos <- 
     object_pos %>% 
     activate_mass_dataset(what = "variable_info") %>%
     filter(na_freq < 0.2 & (na_freq.1 < 0.5 | na_freq.2 < 0.5))
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