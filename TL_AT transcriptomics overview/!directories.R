# Home
dir_home <- here::here() %>% dirname %>% dirname %>% dirname %>% dirname %>%
  fs::path(.)

# Cleaned data
dir_data <- fs::path(
  dir_home,
  "2_Cleaned data"
)

# Transcriptomics data
dir_trans <- fs::path(
  dir_home,
  "1_Data Cleaning",
  "5_Omics",
  "3_RNAseq",
  "1_human",
  "0_data"
)

# Project
dir_project <- fs::path(
  here::here() %>% dirname %>% dirname
)
