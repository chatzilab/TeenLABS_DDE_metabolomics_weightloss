# Linear-mixed effect model for vAT metabolites and plasma metabolites
# Created by Zhenjiang Li
# 03202024
# Run the code chunks in the !analytic_workflow_AT metabolomics.qmd first

colnames_var <- colnames(data_AT_metabo_select)[-1]
colnames_omics <- colnames(data_plasma_metabo_endo_l1)[-c(1:3)]

df <- left_join(
  data,
  data_AT_metabo_select,
  by = "key"
)
df <- left_join(
  df,
  data_plasma_metabo_endo_l1,
  by = c("key", "visit", "visit_new")
)

lmm_AT_plasma_covars_output <- lapply(
  colnames_var,
  function(x) summarize_AT_plasma(
    df = df,
    var = x,
    colnames_omics = colnames_omics,
    covariates
  )
)

saveRDS(
  lmm_AT_plasma_covars_output,
  fs::path(
    dir_project,
    "1_scripts",
    "0_temporary_data",
    "model_statistics_vAT_plasma_selected.rds"
  )
)
