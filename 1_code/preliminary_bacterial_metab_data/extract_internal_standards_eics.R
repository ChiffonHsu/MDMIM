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
  C = 12.0000,
  H = 1.007825,
  D = 2.014102,   # Deuterium
  O = 15.994915,
  N = 14.003074,
  P = 30.973762,
  S = 31.972071
)

H_mass <- 1.007276  # proton mass for adducts

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
    positive_actual_mz = neutral_mass + H_mass,   # [M+H]+
    negative_actual_mz = neutral_mass - H_mass    # [M−H]−
  ) %>%
  ungroup()

  # View the table
internal_standards_list

#2. extract eics
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
targeted_table_pos <- internal_standards_list %>%
  filter(positive_ion_response != "\\") %>%
  transmute(
    variable_id = name,             # must be first column
    mz = positive_mode_actual_mz,        # second column
    rt = rt *60                         # third column (or rt_s if in seconds)
    # optional: add other columns if needed
  )


# Extract EICs — no loop needed
extract_eic(
  targeted_table = targeted_table_pos,
  object = xdata3,
  polarity = "positive",
  mz_tolerance = mz_tol,
  rt_tolerance = rt_tol,
  threads = threads,
  add_point = FALSE,
  path = path,
  group_for_figure = "QC",
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
targeted_table_neg <- internal_standards_list %>%
  filter(negative_ion_response != "\\") %>%
  transmute(
    variable_id = name,             # must be first column
    mz = negative_mode_actual_mz,        # second column
    rt = rt *60                         # third column (or rt_s if in seconds)
    # optional: add other columns if needed
  )


# Extract EICs — no loop needed
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


