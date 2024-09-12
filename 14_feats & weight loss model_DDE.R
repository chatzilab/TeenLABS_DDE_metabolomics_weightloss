# Linear-mixed effect model for metabolic features and bmi/uwaist
# Adjusted for numeric DDE
# Created by Zhenjiang Li
# 01252024
# Run the codes of setup-select chunk in the !analytic_workflow.Rmd first

lmm_random_intercept_knot1_covars_bmi_DDE_output <- {}
lmm_random_intercept_knot1_covars_uwaist__DDE_output <- {}
lmm_random_intercept_knot1_covars_outcome_DDE_output <- list(
  lmm_random_intercept_knot1_covars_bmi_DDE_output,
  lmm_random_intercept_knot1_covars_uwaist__DDE_output
)
for(i in 1:length(platforms))
{
  data_abundance_df <- data_abundance[[i]]
  model_df <- data_abundance_df
  
  covariates_df <- data %>%
    ungroup() %>%
    select(all_of(covariates), dde_cat_num)
  key <- data[,"key"]
  visit <- data$visit
  
  for(k in 1:length(outcome))
  {
    y <- data[,outcome[k]]
    cl <- makeCluster(detectCores()-2)
    registerDoParallel(cl)
    model_output <- foreach(j = 1:dim(model_df)[2], .combine = rbind, .packages = c("nlme", "splines2")) %dopar%
      fit_lmm_random_intercept_knot1_outcome_covars_ddenum(x = model_df[,j], y = y, key = key, visit = visit, covariates_df)
    stopCluster(cl)
    rownames(model_output) <- colnames(model_df)
    
    model_output <- model_output %>%
      as.data.frame %>%
      dplyr::mutate(p_adjust = p.adjust(pvalue_anova, method = "BH"))
    
    model_output <- cbind(colnames(model_df), model_output)
    colnames(model_output)[1] <- "row_name"
    
    model_output <- model_output %>%
      as.data.frame %>%
      mutate_at(vars(coef, se, statistic, pvalue, pvalue_anova, L_ratio, p_adjust), as.numeric)
    
    lmm_random_intercept_knot1_covars_outcome_DDE_output[[k]][[i]] <- as.data.frame(model_output)
    rm(model_output) 
  }
  rm(model_df)
}
names(lmm_random_intercept_knot1_covars_outcome_DDE_output[[1]]) <- platforms
names(lmm_random_intercept_knot1_covars_outcome_DDE_output[[2]]) <- platforms
names(lmm_random_intercept_knot1_covars_outcome_DDE_output) <- outcome

save(lmm_random_intercept_knot1_covars_outcome_DDE_output,
     file = fs::path(dir_temp_data, "model_statistics_bmi_uwaist_DDE.RData"))