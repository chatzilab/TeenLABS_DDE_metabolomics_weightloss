# Linear-mixed effect model for vAT DDE (numeric 3-level) and metabolic features
# Created by Zhenjiang Li
# 01252024
# Run the codes of setup-select chunk in the !analytic_workflow.Rmd first

lmm_random_intercept_knot1_covars_DDE_num_output <- {}
for(i in 1:length(platforms))
{
  data_abundance_df <- data_abundance[[i]]
  model_df <- data_abundance_df
  
  covariates_df <- data %>%
    ungroup() %>%
    select(all_of(covariates))
  x <- data[,exposure[4]]
  key <- data[,"key"]
  visit <- data$visit
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  model_output <- foreach(j = 1:dim(model_df)[2], .combine = rbind, .packages = c("nlme", "splines2")) %dopar%
    fit_lmm_random_intercept_knot1_covars_ddenum(x = x, y = model_df[,j], key = key, visit = visit, covariates_df)
  stopCluster(cl)
  rownames(model_output) <- rep(colnames(model_df), each = 3)
  model_output <- cbind(model_output, rep(c(0,1,2), each = 1, time = dim(model_df)[2]))
  model_output <- cbind(model_output, rep(c("continuous 3-level"), time = dim(model_df)[2]*3))
  colnames(model_output)[7:8] <- c("period", "exposure")
  
  model_output <- model_output %>% 
    as.data.frame %>%
    dplyr::group_by(period) %>%
    dplyr::mutate(p_adjust = p.adjust(pvalue, method = "BH"))
  
  model_output <- model_output %>%
    as.data.frame %>%
    dplyr::mutate(p_adjust_interaction = p.adjust(pvalue_interaction, method = "BH"))
  
  model_output <- cbind(rep(colnames(model_df), each = 3), model_output)
  colnames(model_output)[1] <- "row_name"
  
  model_output <- model_output %>%
    as.data.frame %>%
    mutate_at(vars(coef, se, statistic, pvalue, pvalue_interaction, L_ratio, p_adjust, p_adjust_interaction), as.numeric)
  
  lmm_random_intercept_knot1_covars_DDE_num_output[[i]] <- as.data.frame(model_output)
  rm(model_output, model_df)
}
names(lmm_random_intercept_knot1_covars_DDE_num_output) <- platforms

save(lmm_random_intercept_knot1_covars_DDE_num_output,
     file = fs::path(dir_temp_data, "model_statistics_DDE_num.RData"))