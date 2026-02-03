#merge data
load("2_data/preliminary_determination_bacterial_metabolomics/MS2/object_merged")
colnames(object_merged@sample_info)
#[1] "sample_id"        "experiment"       "sample_type"      "Bacteria"         "medium"           "class"            "injectionorder"  
#[8] "group"            "sample_id_merge"  "sample_id_2"      "experiment_2"     "sample_type_2"    "Bacteria_2"       "medium_2"        
#[15] "class_2"          "injectionorder_2" "group_2" 

#keep to one annotation evey feature, by taking the highest similarity score for each feature
# Select best annotation per variable_id
annot1 <- object_merged@annotation_table %>%
  group_by(variable_id) %>%
  arrange(desc(Total.score)) %>%
  slice(1) %>%
  ungroup()
object_annot_best <- object_merged
object_annot_best@annotation_table <- annot1 #2965 variables left
save(object_annot_best, file = "2_data/preliminary_determination_bacterial_metabolomics/MS2/object_annot_best")

#1. Get bacteria samples and their corresponding medium information
bacteria_info <- object_annot_best@sample_info %>%
  filter(group == "bacteria_samples") %>%
  dplyr::select(sample_id, Bacteria, medium, class)

#2. Create a vector of bacteria sample ids
bacteria_classes <- unique(bacteria_info$class)
#3. extra medium info based on medium values
medium_6_samples <- object_annot_best@sample_info %>%
  filter(group == "medium_samples", medium == 6) %>%
  pull(sample_id)
medium_10_samples <- object_annot_best@sample_info %>%
  filter(group == "medium_samples", medium == 10) %>%
  pull(sample_id)
#4. calculate p value for the differential analysis
safe_fc <- function(object, control_ids, case_ids) {
  expr <- object@expression_data
  
  # FIX: add pseudocount to avoid log2(0) issues
  expr <- expr + 1e-6
  
  control_mean <- rowMeans(expr[, control_ids, drop = FALSE], na.rm = TRUE)
  case_mean    <- rowMeans(expr[, case_ids, drop = FALSE], na.rm = TRUE)
  
  log2(case_mean / control_mean)
}

### Raw p-value (later adjust using BH)
my_pvalue <- function(object, control_ids, case_ids) {
  expr <- object@expression_data
  
  apply(expr, 1, function(x) {
    x_control <- x[control_ids]
    x_case    <- x[case_ids]
    
    x_control <- x_control[!is.na(x_control)]
    x_case    <- x_case[!is.na(x_case)]
    
    if (length(x_control) < 2 | length(x_case) < 2) return(1)
    if (var(x_control)==0 & var(x_case)==0) return(1)
    
    tryCatch(t.test(x_control, x_case)$p.value, error=function(e) 1)
  })
}
#5. Loop through bacteria samples and perform comparisons
###############################################
### 5. Hybrid analysis based on sample size ###
###############################################
comparison_results <- purrr::map(bacteria_classes, function(cls) {
  
  case_samples <- bacteria_info$sample_id[bacteria_info$class == cls]
  bacteria_medium <- unique(bacteria_info$medium[bacteria_info$class == cls])[1]
  
  control_samples <- if (bacteria_medium == 6) medium_6_samples else medium_10_samples
  
  n_case <- length(case_samples)
  n_control <- length(control_samples)
  
  message("\nClass: ", cls, " | case=", n_case, " | control=", n_control)
  
  #################################################
  ### CASE A: Full analysis (>=3 samples/group)
  #################################################
  if (n_case >= 3 & n_control >= 3) {
    message(" → FULL analysis (mutate_fc + BH p-values)")
    
    # compute FC using pseudocount FC
    FC_full <- safe_fc(object_annot_best, control_samples, case_samples)
    
    # compute p-values using massdataset (raw)
    obj_p <- object_annot_best %>%
      mutate_p_value(control_sample_id = control_samples,
                     case_sample_id = case_samples,
                     method = "t.test",
                     p_adjust_methods = "BH")
    
    df <- tibble(
      variable_id = rownames(object_annot_best@expression_data),
      class = cls,
      fc_final = FC_full,
      pval_final = obj_p@variable_info$p_value_adjust
    )
    
    
    return(df)
    
  }
  #################################################
  ### CASE B: FALLBACK — SAFE hybrid version
  ### Use mutate_fc for FC, manual t-test for p-value
  #################################################
  message(" → FALLBACK analysis (safe FC + BH p-values)")
  
  FC     <- safe_fc(object_annot_best, control_samples, case_samples)
  PVAL_r <- my_pvalue(object_annot_best, control_samples, case_samples)
  PVAL_a <- p.adjust(PVAL_r, method = "BH")   ### FIX ⭐ adjust p-values!
  
  tibble(
    variable_id = rownames(object_annot_best@expression_data),
    class = cls,
    fc_final = FC,
    pval_final = PVAL_a
  )
})

names(comparison_results) <- bacteria_classes

########################################
### 6.5 Add compound names (DO NOT LOSE NAMES)
########################################
compound_lookup <- object_annot_best@annotation_table %>%
  dplyr::select(variable_id, Compound.name)

comparison_results <- map(comparison_results, function(df) {
  left_join(df, compound_lookup, by = "variable_id")
})

########################################
### Collapse duplicates
########################################
collapsed_results <- lapply(comparison_results, function(df) {
  df %>%
    mutate(
      ### FIX ⭐ Replace NA names with Unknown_varid
      Compound.name = ifelse(
        is.na(Compound.name),
        paste0("Unknown_", variable_id),
        Compound.name
      )
    ) %>%
    group_by(Compound.name) %>%
    arrange(pval_final, desc(abs(fc_final))) %>%
    slice(1) %>%
    ungroup()
})

names(collapsed_results) <- bacteria_classes

#7. volcano plot
#a. combine into one big tibble
library(dplyr)
library(purrr)
all_results <- bind_rows(collapsed_results)

all_results <- all_results %>%
  mutate(
    logFC = fc_final,
    neg_log10_p = -log10(pval_final)
  )

volcano_prepare <- function(df) {
  df %>%
    mutate(
      sig = pval_final < 0.05 & abs(fc_final) > 1,
      regulation = case_when(
        !sig ~ "none",
        fc_final > 1 ~ "up",
        fc_final < -1 ~ "down"
      )
    )
}

volcano_plot2 <- function(df, class_name) {
  
  df <- volcano_prepare(df)
  
  label_df <- df %>%
    filter(regulation != "none") %>%
    arrange(pval_final) %>%
    slice_head(n = 10)
  
  p <- ggplot(df, aes(x = fc_final, y = -log10(pval_final), color = regulation)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("none"="grey70","up"="red","down"="blue")) +
    theme_minimal() +
    labs(
      title = paste("Volcano Plot -", class_name),
      x = "Log2 Fold Change",
      y = "-log10(P-value)"
    )
  
  if (nrow(label_df) > 0) {
    p <- p + geom_text_repel(
      data = label_df,
      aes(label = Compound.name),
      size = 3,
      max.overlaps = 20
    )
  }
  
  return(p)
}
#e.save each plot automatically
volcano_list <- map(
  names(collapsed_results),
  ~ volcano_plot2(collapsed_results[[.x]], .x)
)

names(volcano_list) <- names(collapsed_results)

walk2(
  volcano_list,
  names(volcano_list),
  ~ ggsave(
    filename = paste0("volcano_", .y, ".pdf"),
    path = "3_data_analysis/preliminary_determination_bacterial_metabolomics/differential_analysis/best_annot/hybrid_adj_fc_pval",
    plot = .x,
    width = 8,
    height = 6
  )
)

###############PART TWO: COLLATE SIGNIFICANT FEATURES IN EACH BACTERIA
#(1) Define significance function
sig_details <- lapply(collapsed_results, function(df) {
  
  df2 <- df %>%
    mutate(
      sig = pval_final < 0.05 & abs(fc_final) > 1,
      regulation = case_when(
        !sig ~ "none",
        fc_final > 1 ~ "up",
        fc_final < -1 ~ "down"
      )
    )
  
  tibble(
    up = sum(df2$regulation == "up"),
    down = sum(df2$regulation == "down"),
    total_sig = sum(df2$regulation != "none"),
    percent = round(100 * sum(df2$regulation != "none") / nrow(df2), 2)
  )
})

sig_details <- bind_rows(sig_details, .id = "bacterium")

sig_details
