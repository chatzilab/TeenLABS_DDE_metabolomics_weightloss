# load AT transcriptomics data
RNAseq <- readRDS(
  fs::path(dir_trans,
           "TL_baseline_AT_transcriptomics.rds")
)

# Load TL cohort data
data_long <- readRDS(fs::path(dir_data,
                              "tl_covariates_outcomes_environmentalchemicals_l.RDS"))

data_wide <- readRDS(fs::path(dir_data,
                              "tl_covariates_outcomes_environmentalchemicals_w.RDS"))

# Load plasma metabolome dataset
data_plasma_metabo <- readRDS(fs::path(dir_data,
                                       "Metabolomics",
                                       "plasma_metabolomics_fts.rds"))
