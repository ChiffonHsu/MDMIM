load("MS2Library_ChiFang/In_Vitro/RPLC/POS/Result/data_cleaning/object_pos")
#Generate table including metabolite id, log2FC and p values
top_up <- subset(results_df, significance == "up") %>% arrange(desc(log2FC))
top_down <- subset(results_df, significance == "down") %>% arrange(log2FC)
#Generate heatmaps
library(massdataset)

#Expression matrix=intensity values, look at intensity per metabolite/sample first
load("MS2Library_ChiFang/In_Vitro/RPLC/POS/Result/data_cleaning/object_pos")
exp_data <- extract_expression_data(object_pos)
head(exp_data[, 11:5])
#Rows=metabolites+ids, columns=samples
#Link metabolite info
var_info <- extract_variable_info(object_pos)
head(var_info)
#Subset to significant metabolites only -> create sig_exp, a matrix for only significant metabolites 
sig <- read.csv("MS2Library_ChiFang/In_Vitro/RPLC/POS/significant_metabolites.csv")
sig_exp <- exp_data[rownames(exp_data) %in% sig$feature_id, ]

#generate heatmap of only the significant metabolites