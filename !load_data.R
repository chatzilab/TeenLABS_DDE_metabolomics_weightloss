# Load covariates outcomes dataset wide format
data_wide <- readRDS(fs::path(dir_data,
                              "tl_covariates_outcomes_environmentalchemicals_w.RDS"))

# Load covariates outcomes dataset long format
data_long <- readRDS(fs::path(dir_data,
                              "tl_covariates_outcomes_environmentalchemicals_l.RDS"))

# Load plasma metabolome dataset
data_plasma_metabo <- readRDS(fs::path(dir_data,
                                       "Metabolomics",
                                       "plasma_metabolomics_fts.rds"))

# Load metabolic feature annotatoin
feat_annot_l <- readRDS(fs::path(dir_home,
                                 "4_Projects",
                                 "TL_POPs AT_plasma metabolome",
                                 "0_data",
                                 "confirmed_unique_annotation.rds"))

# Load AT metabolome dataset
data_AT_metabo <- readRDS(fs::path(dir_data,
                                       "Metabolomics",
                                       "adipose_tissue_metabolomics_fts.rds"))

feat_annot_l_AT <- readRDS(fs::path(
  dir_home,
  "1_Data Cleaning",
  "5_Omics",
  "1_Metabolomics",
  "3_Adipose tissue",
  "0_Data",
  "feature_metadata_with_summaries.RDS"
))
