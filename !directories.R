# Directory

# home directory for project
dir_home <- here::here() %>% 
  dirname() %>% dirname() %>% dirname() %>% fs::path()

# project folder
dir_project <- here::here() %>%
  dirname() %>% fs::path()

# data folder
dir_data <- fs::path(dir_home,"2_Cleaned Data")

# report folder
dir_report <- fs::path(dir_project, "2_reports")

# figure folder
dir_figure <- fs::path(dir_project, "3_figures")

# Temporary datasets
dir_temp_data <- fs::path(here::here(), "0_temporary_data")
