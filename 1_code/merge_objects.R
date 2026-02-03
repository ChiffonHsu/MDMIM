#merge objects
load("2_data/preliminary_determination_bacterial_metabolomics/MS2/POS/object_pos4_ms2")
load("2_data/preliminary_determination_bacterial_metabolomics/MS2/NEG/object_neg4_ms2")
#REMOVE FEATURES AIWHOUT ANNOTATION
object_pos4_ms2 <- 
  object_pos4_ms2 %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2 | Level == 3)
object_pos4_ms2
object_neg4_ms2<- 
  object_neg4_ms2%>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  filter(!is.na(Level)) %>% 
  filter(Level == 1 | Level == 2 | Level == 3)
object_neg4_ms2
#merge objects
head(colnames(object_pos4_ms2))
head(colnames(object_neg4_ms2))
#normalize naming before merging
clean_sample_id <- function(x) {
  # Add this line to remove POS or NEG from the start of the string
  # followed by zero or more whitespace characters
  x <- gsub("^(POS|NEG)\\s*", "", x, ignore.case = TRUE)
  
  # remove (2) or (3) or any (...) patterns
  x <- gsub("\\s*\\(.*\\)", "", x)
  
  # remove trailing "."
  x <- gsub("\\.$", "", x)
  
  # remove polarity suffixes
  x <- gsub("-POS$|-NEG$|POS-$|NEG-$", "", x, ignore.case = TRUE)
  
  # trim whitespace
  x <- trimws(x)
  
  # return cleaned IDs
  return(x)
}

# Clean sample_info
object_pos4_ms2 <- object_pos4_ms2 %>%
  activate_mass_dataset("sample_info") %>%
  mutate(sample_id_merge = clean_sample_id(sample_id))

# Apply the same cleaned IDs to expression_data
new_names_pos <- clean_sample_id(colnames(object_pos4_ms2@expression_data))
colnames(object_pos4_ms2@expression_data) <- new_names_pos

# Also update sample_id to match the expression matrix (IMPORTANT)
object_pos4_ms2 <- object_pos4_ms2 %>%
  activate_mass_dataset("sample_info") %>%
  mutate(sample_id = sample_id_merge)

# Clean sample_info
object_neg4_ms2<- object_neg4_ms2%>%
  activate_mass_dataset("sample_info") %>%
  mutate(sample_id_merge = clean_sample_id(sample_id))

# Apply to expression_data
new_names_neg <- clean_sample_id(colnames(object_neg4_ms2@expression_data))
colnames(object_neg4_ms2@expression_data) <- new_names_neg

# Update sample_id
object_neg4_ms2<- object_neg4_ms2%>%
  activate_mass_dataset("sample_info") %>%
  mutate(sample_id = sample_id_merge)

object_merged <- 
  merge_mass_dataset(x = object_pos4_ms2, 
                     y = object_neg4_ms2, 
                     sample_direction = "inner",
                     variable_direction = "full", 
                     sample_by = "sample_id_merge", 
                     variable_by = c("variable_id", "mz", "rt"))
object_merged
save(object_merged, file = "2_data/preliminary_determination_bacterial_metabolomics/MS2/object_merged")
#remove redundant metabolites
object_merged <-
  object_merged %>% 
  activate_mass_dataset(what = "annotation_table") %>% 
  group_by(Compound.name) %>% 
  filter(Level == min(Level)) %>% 
  filter(SS == max(SS)) %>% 
  slice_head(n = 1)

