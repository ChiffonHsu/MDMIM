library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(tidymass)
#1. process raw data
process_data(
  path = "2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/rawdata",
  polarity = "positive",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 4,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  group_for_figure = "QC"
)
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/rawdata/Result/object")
object #4813 x 79 data frame (variable_info)
#2. generate interactive plot
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "Case",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_case.html", selfcontained = TRUE)
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "Blank",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_blank.html", selfcontained = TRUE)
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "QC",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_QC.html", selfcontained = TRUE)
######################################################################################################################################################
#1. process raw data
process_data(
  path = "2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/rawdata",
  polarity = "negative",
  ppm = 10,
  peakwidth = c(10, 60),
  threads = 4,
  output_tic = TRUE,
  output_bpc = TRUE,
  output_rt_correction_plot = TRUE,
  min_fraction = 0.5,
  group_for_figure = "QC"
)
load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result/object")
object #1206 x 118  data frame (variable_info); 120 x 4 (sample_info)
#2. generate interactive plot
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata3,
  group_for_figure = "Case",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_case.html", selfcontained = TRUE)
load("2_data/preliminary_determination_bacterial_metabolomics/NEG/rawdata/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "QC",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_QC.html", selfcontained = TRUE)
load("2_data/preliminary_determination_bacterial_metabolomics/NEG/rawdata/Result/intermediate_data/xdata2")
##set the group_for_figure if you want to show specific groups. And set it as "all" if you want to show all samples.
p <- massprocesser::plot_adjusted_rt(
  object = xdata2,
  group_for_figure = "Blank",
  interactive = TRUE
)
print(p)
htmlwidgets::saveWidget(p, "rt_adjusted_plot_blank.html", selfcontained = TRUE)
