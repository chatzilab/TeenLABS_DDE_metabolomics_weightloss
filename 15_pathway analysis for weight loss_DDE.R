# Pathway enrichment analysis via mummichog for metabolic features and bmi/uwaist
# Adjusted for numeric DDE
# Created by Zhenjiang Li
# 01252024
# Run the codes of setup-select chunk in the !analytic_workflow.Rmd first

for(k in 1:length(mum_wkpath_outcome_DDE))
{
  if(file.exists(mum_wkpath_outcome_DDE[k])){
    do.call(function(x) unlink(x, recursive = T), list(list.files(mum_wkpath_outcome_DDE[k], full.names = T)))
  }else{
    dir.create(mum_wkpath_outcome_DDE[k])
  }
}


for(k in 1:length(lmm_random_intercept_knot1_covars_outcome_DDE_output))
{
  lmm_random_intercept_knot1_covars_output_temp <- lmm_random_intercept_knot1_covars_outcome_DDE_output[[k]]
  for(i in 1:length(lmm_random_intercept_knot1_covars_output_temp))
  {
    feat_table_name <- platforms[i]
    model_stats <- lmm_random_intercept_knot1_covars_output_temp[[i]]
    input <- model_stats %>%
      na.omit %>%
      dplyr::select(row_name, p_adjust, L_ratio)
    input$mz <- as.numeric(sub('_.*', '', input$row_name))
    input$rt <- as.numeric(sub('.*_', '', input$row_name))
    input <- input %>%
      dplyr::select(mz, rt, p_adjust, L_ratio)
    colnames(input) <- c("mz", "rtime", "p-value", "t-score")
    
    write.table(
      input,
      paste0(mum_wkpath_outcome_DDE[k],"/sum_table_", feat_table_name,".txt"),
      quote = F, row.names = F, sep = "\t"
    )
    
    rm(feat_table_name, model_stats, input)
  }
  rm(lmm_random_intercept_knot1_covars_output_temp)
}

# cmd code for mummichog
for(k in 1:length(outcome))
{
  for(i in 1:length(platforms))
  {
    feat_table_name <- platforms[i]
    if(grepl("neg", feat_table_name))
    {
      mode <- "negative"
    }else{
      mode <- "positive"
    }
    cmd_code <- paste("python -m mummichog.main -f",
                      paste0("sum_table_", feat_table_name,".txt"),
                      "-o", paste0("mummichog_", feat_table_name),
                      "-p 1000 -c ", q_cutoff_sig_feat, " -m ", mode, " -z TRUE")
    print(cmd_code)
    setwd(mum_wkpath_outcome_DDE[k])
    system(cmd_code)
    
    rm(feat_table_name, mode, cmd_code)
  }
}

