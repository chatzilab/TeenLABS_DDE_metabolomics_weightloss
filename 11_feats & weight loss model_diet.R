# Linear-mixed effect model for metabolic features and bmi/uwaist
# Created by Zhenjiang Li
# 12142023
# Run the codes of setup-select chunk in the !analytic_workflow.Rmd first

# Create data for the analysis
key_include <- data %>%
  dplyr::filter(visit_new %in% c(0, 6, 12)) %>%
  dplyr::select(key, visit, kcal) %>%
  na.omit %>%
  dplyr::select(key) %>% unlist %>% unique
length(key_include)

data_sub <- data %>%
  dplyr::filter(key %in% intersect(key_include, subject_include)) %>%
  dplyr::filter(visit_new %in% c(0, 6, 12))
dim(data_sub)

data_sub$row_name <- with(data_sub, paste(key, "_", visit_new))

# EDA
table1(as.formula(paste("~", paste0(c(covariates, "kcal"), collapse = "+"))),
                  data = data_sub)

# Run model
lmm_random_intercept_covars_diet_bmi_output <- {}
lmm_random_intercept_covars_diet_uwaist_output <- {}
lmm_random_intercept_covars_diet_outcome_output <- list(
  lmm_random_intercept_covars_diet_bmi_output,
  lmm_random_intercept_covars_diet_uwaist_output
)
for(i in 1:length(platforms))
{
  data_abundance_df <- data_abundance[[i]]
  model_df <- data_abundance_df[match(data_sub$row_name, rownames(data_abundance_df)),]
  
  covariates_df <- data_sub %>%
    ungroup() %>%
    select(all_of(covariates), kcal)
  key <- data_sub[,"key"]
  visit <- data_sub$visit
  
  for(k in 1:length(outcome))
  {
    y <- data_sub[,outcome[k]]
    cl <- makeCluster(detectCores()-2)
    registerDoParallel(cl)
    model_output <- foreach(j = 1:dim(model_df)[2], .combine = rbind, .packages = c("nlme", "dplyr")) %dopar%
      fit_lmm_random_intercept_diet_outcome_covars(x = model_df[,j], y = y, key = key, visit = visit, covariates_df)
    stopCluster(cl)
    rownames(model_output) <- colnames(model_df)
    
    model_output <- model_output %>%
      as.data.frame %>%
      dplyr::mutate(p_adjust = p.adjust(pvalue_anova, method = "BH"))
    
    model_output <- model_output %>%
      as.data.frame %>%
      dplyr::mutate(p_adjust_nodiet = p.adjust(pvalue_anova_nodiet, method = "BH"))
    
    model_output <- cbind(colnames(model_df), model_output)
    colnames(model_output)[1] <- "row_name"
    
    model_output <- model_output %>%
      as.data.frame %>%
      mutate_at(vars(coef, se, statistic, pvalue, pvalue_anova, L_ratio, p_adjust, pvalue_anova_nodiet, L_ratio_nodiet, p_adjust_nodiet), as.numeric)
    
    lmm_random_intercept_covars_diet_outcome_output[[k]][[i]] <- as.data.frame(model_output)
    rm(model_output) 
  }
  rm(model_df)
}
names(lmm_random_intercept_covars_diet_outcome_output[[1]]) <- platforms
names(lmm_random_intercept_covars_diet_outcome_output[[2]]) <- platforms
names(lmm_random_intercept_covars_diet_outcome_output) <- outcome

save(lmm_random_intercept_covars_diet_outcome_output,
     file = fs::path(dir_temp_data, "model_statistics_bmi_uwaist_diet.RData"))
