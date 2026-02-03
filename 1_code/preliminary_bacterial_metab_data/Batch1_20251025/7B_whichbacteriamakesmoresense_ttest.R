#This script uses just neg object to find out which bacteria makes more sense (releases or consumes more metabolites)
load("~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/object_pos3_ms2")

    colnames(object_pos3_ms2@sample_info)
    head(object_pos3_ms2@sample_info)
    # sample_id experiment sample_type Bacteria medium class batch group
    # 1  Blank_1_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
    # 2 Blank_10_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
    # 3 Blank_11_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
    # 4  Blank_2_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
    # 5  Blank_3_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank
    # 6  Blank_4_POS (2).   20251025        lpos  Subject   <NA> Blank     2 Blank

  #under object_pos3_ms2@sample_info$class, samples are grouped as Control and Subject, compare these two groups
    
    #1. Get bacteria samples and their corresponding medium information
    bacteria_info <- object_pos3_ms2@sample_info %>%
      filter(group == "bacteria_samples") %>%
      dplyr::select(sample_id, Bacteria, medium, class)
    
    #2. Create a vector of bacteria sample ids
    bacteria_classes <- unique(bacteria_info$Bacteria)
    #3. extra medium info based on medium values
    medium_6_samples <- object_pos3_ms2@sample_info %>%
      filter(group == "medium_samples", medium == "MGAM-6" ) %>%
      pull(sample_id)
    # medium_10_samples <- object_pos3_ms2@sample_info %>%
    #   filter(group == "medium_samples", medium == 10) %>%
    #   pull(sample_id)
    #4. calculate p value for the differential analysis
    my_fc <- function(object, control_ids, case_ids) {
      expr <- object@expression_data
      
      control_mean <- rowMeans(expr[, control_ids, drop = FALSE], na.rm = TRUE)
      case_mean    <- rowMeans(expr[, case_ids, drop = FALSE], na.rm = TRUE)
      
      log2(case_mean / control_mean)
    }
    
    my_pvalue <- function(object, control_ids, case_ids) {
      
      expr <- object@expression_data
      
      apply(expr, 1, function(x) {
        x_control <- x[control_ids]
        x_case    <- x[case_ids]
        
        # Remove NA values
        x_control <- x_control[!is.na(x_control)]
        x_case    <- x_case[!is.na(x_case)]
        
        # Not enough valid values → return p=1
        if (length(x_control) < 2 | length(x_case) < 2) {
          return(1)
        }
        
        # No variance → t-test cannot run → return p=1
        if (var(x_control) == 0 & var(x_case) == 0) {
          return(1)
        }
        
        # Try t-test safely
        tryCatch(
          t.test(x_control, x_case)$p.value,
          error = function(e) 1
        )
      })
    }
    #5. Loop through bacteria samples and perform comparisons
    comparison_results <- purrr::map(bacteria_classes, function(cls) {
      
      case_samples <- bacteria_info$sample_id[bacteria_info$Bacteria == cls]
      bacteria_medium <- unique(bacteria_info$medium[bacteria_info$Bacteria == cls])
      
      control_samples <- medium_6_samples
      # if (bacteria_medium == 6) {
      #   control_samples <- medium_6_samples
      # } else {
      #   control_samples <- medium_10_samples
      # }
      
      
      # FC and p-value
      FC <- my_fc(object_pos3_ms2, control_samples, case_samples)
      PVAL <- my_pvalue(object_pos3_ms2, control_samples, case_samples)
      
      tibble(
        variable_id = rownames(object_pos3_ms2@expression_data),
        class = cls,
        fold_change = FC,
        p_value = PVAL
      )
    })
    
    names(comparison_results) <- bacteria_classes
    
    #6. add compound names
    compound_lookup <- object_pos3_ms2@annotation_table %>%
      dplyr::select(variable_id, Compound.name)
    
    comparison_results <- map(comparison_results, function(df) {
      df %>% left_join(compound_lookup, by = "variable_id")
    })
    
    #7. volcano plot
    #a. combine into one big tibble
    library(dplyr)
    library(purrr)
    
    ########################################
    ### 6.5 Collapse duplicated Compound.name per bacterium
    ########################################
    
    collapsed_results <- lapply(comparison_results, function(df){
      
      df %>%
        group_by(Compound.name) %>%
        arrange(p_value, desc(abs(fold_change))) %>%   # best based on pval then FC
        slice(1) %>%
        ungroup()
    })
    
    names(collapsed_results) <- names(comparison_results)
    
    all_results <- bind_rows(collapsed_results)
    
    #b. prepare columns
    all_results <- all_results %>%
      mutate(
        logFC = fold_change,
        neg_log10_p = -log10(p_value)
      )
    #c. volcano plot
    volcano_prepare <- function(df) {
      df %>%
        mutate(
          sig = ifelse(!is.na(p_value) & p_value < 0.05 & abs(fold_change) > 1,
                       "yes", "no"),
          regulation = case_when(
            sig == "no"                      ~ "none",
            sig == "yes" & fold_change > 1  ~ "up",
            sig == "yes" & fold_change < -1 ~ "down",
            TRUE ~ "none"
          )
        )
    }
    
    volcano_plot2 <- function(df, class_name) {
      
      df <- volcano_prepare(df)
      
      ###############################  CHANGED
      # Select TOP-10 by smallest p-value
      label_df <- df %>%
        filter(regulation != "none", !is.na(p_value)) %>%
        arrange(p_value) %>%
        slice_head(n = 10)
      
      ggplot(df, aes(x = fold_change, y = -log10(p_value), color = regulation)) +
        geom_point(alpha = 0.7) +
        scale_color_manual(
          values = c("none" = "grey70", "up" = "red", "down" = "blue")
        ) +
        geom_text_repel(
          ###############################  CHANGED
          data = label_df,
          aes(label = Compound.name),   # labels = compound names
          ###############################
          size = 3,
          max.overlaps = 20
        ) +
        labs(
          title = paste("Volcano Plot -", class_name),
          x = "Log2 Fold Change",
          y = "-log10(P-value)"
        ) +
        theme_minimal()
    }
    
    #d. loop to generate one volcano plot per class
    volcano_list <- map(
      names(collapsed_results),
      ~ volcano_plot2(collapsed_results[[.x]], .x)
    )
    names(volcano_list) <- names(collapsed_results)
    
    #e.save each plot automatically
    walk2(
      volcano_list,
      names(volcano_list),
      ~ ggsave(
        filename = paste0("volcano_", .y, ".pdf"),
        path = "~/MDMIM/3_data_analysis/preliminary_determination_bacterial_metabolomics/20251025/batch2_rplc_pos_differential_analysis/all_ttest_fc",
        plot = .x,
        width = 8,
        height = 6
      )
    )
    
    ###############PART TWO: COLLATE SIGNIFICANT FEATURES IN EACH BACTERIA
    #(1) Define significance function
    sig_summary <- lapply(collapsed_results, function(df) {
      
      df2 <- df %>%
        mutate(
          sig = p_value < 0.05 & abs(fold_change) > 1,
          regulation = case_when(
            !sig ~ "none",
            fold_change > 1  ~ "up",
            fold_change < -1 ~ "down"
          )
        )
      
      up_count   <- sum(df2$regulation == "up")
      down_count <- sum(df2$regulation == "down")
      total_sig  <- up_count + down_count
      total_feat <- nrow(df2)
      
      tibble(
        up = up_count,
        down = down_count,
        total_sig = total_sig,
        percent = round(100 * total_sig / total_feat, 2)
      )
    })
    
    # Add bacterium names
    sig_summary <- bind_rows(sig_summary, .id = "bacterium")
    
    sig_summary
    # bacterium                     up  down total_sig percent
    # <chr>                      <int> <int>     <int>   <dbl>
    #   1 Citrobacter youngae         1591  2051      3642    9.09
    # 2 Citrobacter braakii         2731  5559      8290   20.7 
    # 3 Bacteroides fragilis        2630  6970      9600   24.0 
    # 4 Enterococcus casseliflavus  2368  5275      7643   19.1 
    # 5 Enterobacter ludwigii       3390  4262      7652   19.1 
    # 6 Pediococcus acidilactici    2503  6187      8690   21.7 
    # 7 Bacteroides finegoldii      2628  5885      8513   21.2 
    # 8 Enterococcus faecalis       2483  5482      7965   19.9 
    # 9 Enterococcus durans         2566  6904      9470   23.6 
    allttest_sig_summary <- sig_summary
    save(allttest_sig_summary, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/allttest_sig_summary")
    
    allttest_collapsed_results <- collapsed_results
    save(allttest_collapsed_results, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/allttest_collapsed_results")
   
    allttest_comparison_results <- comparison_results
    save(allttest_comparison_results, file = "~/MDMIM/2_data/preliminary_determination_bacterial_metabolomics/20251025/MS2/ms2_annotation/differential_analysis/allttest_comparison_results")
    