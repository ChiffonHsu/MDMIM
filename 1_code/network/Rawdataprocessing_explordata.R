library(tidymass)
library(Biobase)

message("Starting processing of raw data...")
start_time <- Sys.time()

result <- process_data(
  path = "MS2Library_ChiFang/In_vivo/RPLC/POS",
  polarity = "positive",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 10,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = FALSE,
  min_fraction = 0.5,
  group_for_figure = "mix"
)



end_time <- Sys.time()
message("Done. Total time: ", round(end_time - start_time, 2))

load("MS2Library_ChiFang/In_vivo/RPLC/POS/Result/object")
object


load("MS2Library_ChiFang/In_vivo/RPLC/POS/Result/intermediate_data/xdata3")
load("MS2Library_ChiFang/In_vivo/RPLC/POS/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
plot = 
  massprocesser::plot_adjusted_rt(object = xdata2, 
                                  group_for_figure = "mix", 
                                  interactive = TRUE)                                             
plot

png("adjusted_rt_plot_rplc_pos_xdata2.png", width = 1200, height = 800)
massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "mix",
  interactive = FALSE
)
dev.off()


#Extract EIC
targeted_table = 
  readr::read_csv("MS2Library_ChiFang/In_vivo/rplc/POS/Result/Peak_table_for_cleaning.csv")

mean_int = targeted_table %>% 
  dplyr::select(-c(variable_id:rt)) %>% 
  apply(1, function(x){
    mean(x, na.rm = TRUE)
  })

targeted_table = 
  targeted_table %>% 
  dplyr::select(variable_id:rt) %>% 
  dplyr::mutate(mean_int = mean_int) %>% 
  dplyr::arrange(dplyr::desc(mean_int)) %>% 
  head(10) %>% 
  dplyr::select(-mean_int)

targeted_table

#Load EIC data
load("MS2Library_ChiFang/In_vivo/RPLC/POS/Result/intermediate_data/xdata3")

extract_eic(
  targeted_table = targeted_table,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = 15,
  rt_tolerance = 30,
  threads = 5,
  add_point = FALSE,
  path = "MS2Library_ChiFang/In_vivo/RPLC/POS/Result",
  group_for_figure = "mix", 
  feature_type = "png"
)

library(tidymass)
library(tidyverse)

#load object
load("MS2Library_ChiFang/In_vivo/RPLC/NEG/Result/object")
object_pos <- object
object_pos

#Read sample information
sample_info_pos <- readr::read_csv("MS2Library_ChiFang/In_vivo/in_vivo_data.csv")
head(sample_info_pos)

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



dir.create("MS2Library_ChiFang/In_vivo/HILIC/POS/Result/data_cleaning", showWarnings = FALSE, recursive = TRUE)
save(object_pos, file = "MS2Library_ChiFang/In_vivo/HILIC/POS/Result/data_cleaning/object_pos")
object_pos
dim(object_pos)
object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(class)
object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(group)
object_pos %>%
  activate_mass_dataset(what = "sample_info")
  dplyr::count(batch)
  
######if the above doesnt work, use this:
dir.create("MS2Library_ChiFang/In_vivo/HILIC/POS/Result/data_cleaning", showWarnings = FALSE, recursive = TRUE)
save(object_pos, file = "MS2Library_ChiFang/In_vivo/HILIC/POS/Result/data_cleaning/object_pos")
  
object_pos
dim(object_pos)
  
# Use extract_sample_info() to safely count
library(dplyr)
  
si <- extract_sample_info(object_pos)
  
si %>% dplyr::count(class)
si %>% dplyr::count(group)
si %>% dplyr::count(batch)
  
  
#Peak distribution plot of posative mode
p <- object_pos %>%
  '+'(1) %>%
  log(10) %>%
  show_mz_rt_plot() +
  scale_size_continuous(range = c(0.01,2))
ggsave(
  "MS2Library_ChiFang/In_vivo/HILIC/POS/Result/mz_rt_plot_HILIC_pos.png", 
  plot = p,
  width = 8, 
  height = 6, 
  dpi = 300
  )

#Explore missing values (mvs) in posative mode data
get_mv_number(object = object_pos)
#Get missing value number in each samples
get_mv_number(object = object_pos, by = "sample") %>%
  head()
#Get missing value number in each sample
get_mv_number(object = object_pos, by = "variable") %>%
  head()
#Show missing value information
p <- show_missing_values(object = object_pos, show_column_names = FALSE, percentage = TRUE)
print(p)
ggsave(
  "MS2Library_ChiFang/In_vivo/HILIC/POS/Result/missing_values.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
  )

#Show the missing values in samples
show_sample_missing_values(object = object_pos, percentage = TRUE)

#Show the missing values in variables
show_variable_missing_values(
  object = object_pos,
  percentage = TRUE,
  show_x_text = FALSE,
  show_x_ticks = FALSE
) +
  scale_size_continuous((range = c(0.1, 1)))

show_variable_missing_values(object = object_pos, 
                             percentage = TRUE, 
                             show_x_text = FALSE, 
                             show_x_ticks = FALSE) +
  scale_size_continuous(range = c(0.01, 1))