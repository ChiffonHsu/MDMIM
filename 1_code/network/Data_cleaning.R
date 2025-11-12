library(tidymass)
load("MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/object_pos")
load("MS2Library_ChiFang/In_Vitro/RPLC/POS/Result/data_cleaning/object_pos")

#Change batch to character
object_pos
object_pos <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(batch = as.character(batch))

object_pos <- object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::mutate(batch = as.character(batch))

#Check quality of the data
massqc::massqc_report(object = object_pos,
                      path = "MS2Library_ChiFang/In_Vitro/HILIC/POS/Result/data_cleaning/data_quality_before_data_cleaing")

