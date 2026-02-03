library(r4projects)
setwd(get_project_wd())
rm(list = ls())
library(tidymass)
library(tidyverse)
internal_standards <- readxl::read_excel("3_data_analysis/preliminary_determination_bacterial_metabolomics/internal_standards.xlsx")
library(dplyr)
library(stringr)

# 1. calculate Monoisotopic atomic masses
atom_mass <- c(
  C = 12.0000,        # Carbon
  H = 1.007825,       # Hydrogen
  D = 2.014102,       # Deuterium (Hydrogen isotope)
  O = 15.994915,      # Oxygen
  N = 14.003074,      # Nitrogen
  P = 30.973762,      # Phosphorus
  S = 31.972071,      # Sulfur
  Na = 22.989769,     # Sodium (Na)
  Cl = 35.453,        # Chlorine (Cl)
  K = 39.0983,        # Potassium (K)
  Li = 6.941         # Lithium (Li)
)

H_mass <- 1.007276  # proton mass for adducts
Na_mass <- 22.989769
K_mass <- 39.0983
Li_mass <- 6.941
NH4_mass <- 18.806204
Cl_mass <- 35.453
CH3OH_mass <- 32.039    # Methanol (CH₃OH) - for methanol adducts

# Negative ions
OH_mass <- 17.00734    # Hydroxide ion (OH⁻)
CH3COO_mass <- 59.044  # Acetate ion (CH₃COO⁻)
HCOO_mass <- 45.017    # Formate ion (HCOO⁻)
F_mass <- 18.998403    # Fluoride ion (F⁻)
NO3_mass <- 62.0049    # Nitrate ion (NO₃⁻)
SO4_mass <- 96.063     # Sulfate ion (SO₄²⁻)
H2PO4_mass <- 97.976   # Dihydrogen phosphate ion (H₂PO₄⁻)
PO4_mass <- 94.9714
    
   # Function to calculate neutral mass from molecular formula
calc_neutral_mass <- function(formula) {
  # extract each element and its count using regex
  elems <- str_match_all(formula, "([A-Z][a-z]?)(\\d*)")[[1]]
  mass <- 0
  for (i in seq_len(nrow(elems))) {
    elem <- elems[i, 2]
    count <- as.numeric(elems[i, 3])
    if (is.na(count) || count == 0) count <- 1
    mass <- mass + atom_mass[[elem]] * count
  }
  return(mass)
}

  # Example: your internal standards table
  # internal standard list should have at least: name, molecular_formula
internal_standards_list <- internal_standards %>%
  rowwise() %>%
  mutate(
    neutral_mass = calc_neutral_mass(molecular_formula),
    positive_actual_mz_H = neutral_mass + H_mass,   # [M+H]+
    negative_actual_mz_H = neutral_mass - H_mass,   # [M−H]−
    positive_actual_mz_Na = neutral_mass + Na_mass, # [M+Na]+
    negative_actual_mz_Cl = neutral_mass - Cl_mass, # [M−Cl]−
    positive_actual_mz_NH4 = neutral_mass + NH4_mass, # [M+NH4]+
    positive_actual_mz_K = neutral_mass + K_mass,   # [M+K]+
    positive_actual_mz_Li = neutral_mass + Li_mass, # [M+Li]+
    negative_actual_mz_PO4 = neutral_mass - PO4_mass, # [M−PO4]−
    negative_actual_mz_SO4 = neutral_mass - SO4_mass, # [M−SO4]−
    negative_actual_mz_F = neutral_mass - F_mass,   # [M−F]−
    negative_actual_mz_CH3COO = neutral_mass - CH3COO_mass, # [M−CH3COO]−
    negative_actual_mz_HCOO = neutral_mass - HCOO_mass, # [M−HCOO]−
    negative_actual_mz_NO3 = neutral_mass - NO3_mass 
    
  ) %>%
  ungroup()

  # View the table
internal_standards_list
colnames(internal_standards_list)

colnames(internal_standards_list)#2. extract eics
library(massprocesser)
library(dplyr)

mz_tol <- 25
rt_tol <- 60
threads <- 5
path <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/POS/extract_eics/internal_standards"

load("2_data/preliminary_determination_bacterial_metabolomics/MS1/POS/Result/intermediate_data/xdata3")

# 1. Safely extract from the NAnnotatedDataFrame structure
npheno <- xdata3@phenoData

# check slots
slotNames(npheno)

# 2. Pull out the data and varMetadata manually
pdat <- as.data.frame(npheno@data)
vmeta <- as.data.frame(npheno@varMetadata)

# 3. Rebuild a normal AnnotatedDataFrame
pd_fixed <- Biobase::AnnotatedDataFrame(data = pdat, varMetadata = vmeta)

# 4. Assign back to xdata3
phenoData(xdata3) <- pd_fixed

# 5. Verify
class(phenoData(xdata3))

class(xdata3)
range(mz(xdata3))
range(rtime(xdata3))
unique(polarity(xdata3))

# Create targeted_table for positive-mode standards
# targeted_table_pos <- internal_standards_list %>%
#   # filter(positive_ion_response != "\\") %>%
  # transmute(
  #   variable_id = name,             # must be first column
  #   mz = positive_mode_actual_mz,        # second column
  #   rt = rt *60                         # third column (or rt_s if in seconds)
  #   # optional: add other columns if needed
  # )
library(dplyr)
library(tidyr)

# View the table for negative-mode standards
# Replace blank (empty) values with NA in the specific column
internal_standards$positive_mode_actual_mz[internal_standards$positive_mode_actual_mz == ""] <- NA
internal_standards$positive_mode_actual_mz <- as.numeric(internal_standards$positive_mode_actual_mz)
head(internal_standards$positive_mode_actual_mz)


  #1. H ion  
targeted_table_pos <- internal_standards%>%
  transmute(
    variable_id = name,
    mz = positive_mode_actual_mz,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+H]+"
  )
targeted_table_pos <- targeted_table_pos[-3, ]
# Extract EICs 
extract_eic(
  targeted_table = targeted_table_pos,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
   #2. Na ion
targeted_table_pos_Na <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = positive_actual_mz_Na,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+Na]+"
  )
extract_eic(
  targeted_table = targeted_table_pos_Na,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
     #3. NH4 ion
targeted_table_pos_NH4 <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = positive_actual_mz_NH4,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+NH4]+"
  )
extract_eic(
  targeted_table = targeted_table_pos_NH4,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
   #4. K ion
targeted_table_pos_K <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = positive_actual_mz_K,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+K]+"
  )
extract_eic(
  targeted_table = targeted_table_pos_K,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
  #5. Li ion
targeted_table_pos_Li <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = positive_actual_mz_Li,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+Li]+"
  )
extract_eic(
  targeted_table = targeted_table_pos_Li,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
####################################################################################################################
library(massprocesser)
library(dplyr)

mz_tol <- 25
rt_tol <- 60
threads <- 5
path <- "3_data_analysis/preliminary_determination_bacterial_metabolomics/NEG/extract_eics/internal_standards"

load("2_data/preliminary_determination_bacterial_metabolomics/MS1/NEG/Result/intermediate_data/xdata3")

# 1. Safely extract from the NAnnotatedDataFrame structure
npheno <- xdata3@phenoData

# check slots
slotNames(npheno)

# 2. Pull out the data and varMetadata manually
pdat <- as.data.frame(npheno@data)
vmeta <- as.data.frame(npheno@varMetadata)

# 3. Rebuild a normal AnnotatedDataFrame
pd_fixed <- Biobase::AnnotatedDataFrame(data = pdat, varMetadata = vmeta)

# 4. Assign back to xdata3
phenoData(xdata3) <- pd_fixed

# 5. Verify
class(phenoData(xdata3))

class(xdata3)
range(mz(xdata3))
range(rtime(xdata3))
unique(polarity(xdata3))


# Create targeted_table for positive-mode standards
#1. H ion
targeted_table_neg <- internal_standards %>%
  transmute(
    variable_id = name,
    mz = negative_mode_actual_mz,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+H]-"
  )
targeted_table_neg <- targeted_table_neg[-5, ]
targeted_table_neg$mz <- as.numeric(targeted_table_neg$mz)
extract_eic(
  targeted_table = targeted_table_neg,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
#2. Cl ion
targeted_table_neg_Cl <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = negative_actual_mz_Cl,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+Cl]+"
  )
extract_eic(
  targeted_table = targeted_table_neg_Cl,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
#3. PO4 ion
targeted_table_neg_PO4 <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = negative_actual_mz_PO4,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+PO4]+"
  )
extract_eic(
  targeted_table = targeted_table_neg_PO4,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
#4. SO4 ion
targeted_table_neg_SO4 <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = negative_actual_mz_SO4,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+SO4]+"
  )
extract_eic(
  targeted_table = targeted_table_neg_SO4,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
#5. F ion
targeted_table_neg_F <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = negative_actual_mz_F,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+F]+"
  )
extract_eic(
  targeted_table = targeted_table_neg_F,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
#6. CH3COO ion
targeted_table_neg_CH3COO <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = negative_actual_mz_CH3COO,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+CH3COO]+"
  )
extract_eic(
  targeted_table = targeted_table_neg_CH3COO,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
#7. HCOO ion
targeted_table_neg_HCOO <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = negative_actual_mz_HCOO,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+HCOO]+"
  )
extract_eic(
  targeted_table = targeted_table_neg_HCOO,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
#8. NO3 ion
targeted_table_neg_NO3 <- internal_standards_list %>%
  transmute(
    variable_id = name,
    mz = negative_actual_mz_NO3,
    rt = rt * 60,  # Convert RT to seconds
    ion_type = "[M+NO3]+"
  )
extract_eic(
  targeted_table = targeted_table_neg_NO3,
  object = xdata3,
  polarity = "negative",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "Blank",
  feature_type = "png"
)
 