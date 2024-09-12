create_shortcut <- function(path, short_cut, drive)
{
  return(
    sub(gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", short_cut)), paste0(drive, ":"), path))
}

check_visit <- function(df){
  return(with(df, table(visit, visit_new, exclude = F)))
}

check_sd <- function(X){
  # var.coef
  sd <- unlist(lapply(as.data.frame(X), 
                      function(x) abs(sd(x))))
  return(sd)
}

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

fit_lmm_random_intercept_knot1_covars <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  model <- as.formula(paste("y ~ x*bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  reduced_model <- as.formula(paste("y ~ x + bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  fit <- lme(model,
             random = ~1|key,
             data = df, na.action=na.omit)
  coef_table <- summary(fit)$tTable
  start_index <- grep("xModerate", rownames(coef_table))[2]
  end_index <- start_index+3
  output <- coef_table[c(2,3,start_index:end_index), c(1,2,4,5)]
  
  fit_reduced <- lme(reduced_model,
                     random = ~1|key,
                     data = df, na.action = na.omit)
  overall_sig <- anova(update(fit, . ~ ., method = "ML"),
                       update(fit_reduced, . ~ ., method = "ML"))
  output <- cbind(output,
                  c(overall_sig$`p-value`[2], rep(NA, time = 5)))
  output <- cbind(output,
                  c(overall_sig$L.Ratio[2], rep(NA, time = 5)))
  colnames(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_interaction", "L_ratio")
  return(output)
}

fit_lmm_random_intercept_slope_knot1_covars <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  model <- as.formula(paste("y ~ x*bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  reduced_model <- as.formula(paste("y ~ x + bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  fit <- lme(model,
             random = ~1 + x|key,
             data = df, na.action=na.omit,
             control = lmeControl(opt='optim'))
  coef_table <- summary(fit)$tTable
  start_index <- grep("xModerate", rownames(coef_table))[2]
  end_index <- start_index+3
  output <- coef_table[c(2,3,start_index:end_index), c(1,2,4,5)]
  
  fit_reduced <- lme(reduced_model,
                     random = ~1 + x|key,
                     data = df, na.action = na.omit,
                     control = lmeControl(opt='optim'))
  overall_sig <- anova(update(fit, . ~ ., method = "ML"),
                       update(fit_reduced, . ~ ., method = "ML"))
  output <- cbind(output,
                  c(overall_sig$`p-value`[2], rep(NA, time = 5)))
  output <- cbind(output,
                  c(overall_sig$L.Ratio[2], rep(NA, time = 5)))
  colnames(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_interaction", "L_ratio")
  return(output)
}

fit_lmm_random_intercept_knot1_covars_coefequality <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  model <- as.formula(paste("y ~ x*bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  fit <- lme(model,
             random = ~1|key,
             data = df, na.action=na.omit)
  coef_table <- summary(fit)$tTable
  start_index <- grep("xModerate", rownames(coef_table))[2]
  end_index <- start_index+3
  coefequality_t1 <- linearHypothesis(fit, paste(rownames(coef_table)[start_index], "=", rownames(coef_table)[start_index+1]))
  coefequality_t2 <- linearHypothesis(fit, paste(rownames(coef_table)[end_index-1], "=", rownames(coef_table)[end_index]))
  output <- matrix(nrow = 2, ncol = 2)
  output[1,1] <- "1"
  output[2,1] <- "2"
  output[1,2] <- coefequality_t1$`Pr(>Chisq)`[2]
  output[2,2] <- coefequality_t2$`Pr(>Chisq)`[2]
  colnames(output) <- c("period", "p_coefequality")
  return(output)
}

fit_lmm_random_intercept_knot1_outcome_covars <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  df$x <- as.numeric(x)
  model <- as.formula(paste("y ~ x + bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  reduced_model <- as.formula(paste("y ~ bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  fit <- lme(model,
             random = ~1|key,
             data = df, na.action=na.omit)
  coef_table <- summary(fit)$tTable
  # start_index <- grep("x", rownames(coef_table))[1]
  # end_index <- start_index+1
  output <- coef_table[2, c(1,2,4,5)]
  
  fit_reduced <- lme(reduced_model,
                     random = ~1|key,
                     data = df, na.action = na.omit)
  overall_sig <- anova(update(fit, . ~ ., method = "ML"),
                       update(fit_reduced, . ~ ., method = "ML"))
  output <- c(output,
              overall_sig$`p-value`[2])
  output <- c(output,
              overall_sig$L.Ratio[2])
  names(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_anova", "L_ratio")
  return(output)
}

fit_lmm_random_intercept_knot1_outcome_covars_interaction <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  df$x <- as.numeric(x)
  model <- as.formula(paste("y ~ x*bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + "), " + (1|key)"))
  reduced_model <- as.formula(paste("y ~ x + bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + "), " + (1|key)"))
  fit <- lmer(model,
              data = df, na.action=na.omit)
  coef_table <- summary(fit)$coefficients
  start_index <- grep("x:", rownames(coef_table))[1]
  end_index <- start_index+1
  output <- coef_table[c(2,start_index:end_index), c(1,2,4,5)]
  
  fit_reduced <- lmer(reduced_model,
                      data = df, na.action = na.omit)
  overall_sig <- anova(fit, fit_reduced)
  output <- cbind(output,
                  c(overall_sig$`Pr(>Chisq)`[2], rep(NA, time = 2)))
  output <- cbind(output,
                  c(overall_sig$Chisq[2], rep(NA, time = 2)))
  colnames(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_interaction", "Chisq")
  return(output)
}

fit_lmm_random_intercept_diet_outcome_covars <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  df$x <- as.numeric(x)
  model <- as.formula(paste("y ~ x + visit", "+", paste(colnames(covariates_df)[! colnames(covariates_df) %in% c("key", "site")], collapse = " + ")))
  reduced_model <- as.formula(paste("y ~ visit", "+", paste(colnames(covariates_df)[! colnames(covariates_df) %in% c("key", "site")], collapse = " + ")))
  model_nodiet <- as.formula(paste("y ~ x + visit", "+", paste(colnames(covariates_df)[!colnames(covariates_df) %in% c("key", "kcal", "site")], collapse = " + ")))
  reduced_model_nodiet <- as.formula(paste("y ~ visit", "+", paste(colnames(covariates_df)[!colnames(covariates_df) %in% c("key", "kcal", "site")], collapse = " + ")))
  fit <- lme(model,
             random = ~1|key,
             data = df, na.action=na.omit)
  df_index <- as.numeric(rownames(fit$fitted))
  coef_table <- summary(fit)$tTable
  # start_index <- grep("x", rownames(coef_table))[1]
  # end_index <- start_index+1
  output <- coef_table[2, c(1,2,4,5)]
  
  fit_reduced <- lme(reduced_model,
                     random = ~1|key,
                     data = df,
                     na.action = na.omit)
  overall_sig <- anova(update(fit, . ~ ., method = "ML"),
                       update(fit_reduced, . ~ ., method = "ML"))
  output <- c(output,
              overall_sig$`p-value`[2])
  output <- c(output,
              overall_sig$L.Ratio[2])
  names(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_anova", "L_ratio")
  rm(overall_sig)
  
  fit_nodiet <- lme(model_nodiet,
                    random = ~1|key,
                    data = df %>% slice(df_index),
                    na.action = na.omit)
  fit_reduced_nodiet <- lme(reduced_model_nodiet,
                            random = ~1|key,
                            data = df %>% slice(df_index),
                            na.action = na.omit)
  coef_table_nodiet <- summary(fit_nodiet)$tTable
  output_nodiet <- coef_table_nodiet[2, c(1,2,4,5)]
  names(output_nodiet) <- c("coef_nodiet", "se_nodiet", "statistic_nodiet", "pvalue_nodiet")
  
  overall_sig <- anova(update(fit_nodiet, . ~ ., method = "ML"),
                       update(fit_reduced_nodiet, . ~ ., method = "ML"))
  output <- c(output,
              output_nodiet)
  output <- c(output,
              overall_sig$`p-value`[2])
  output <- c(output,
              overall_sig$L.Ratio[2])
  names(output)[length(output)-1] <- "pvalue_anova_nodiet"
  names(output)[length(output)] <- "L_ratio_nodiet"
  
  return(output)
}

fit_lmm_random_intercept_knot1_covars_ddenum <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  model <- as.formula(paste("y ~ x*bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  reduced_model <- as.formula(paste("y ~ x + bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  fit <- lme(model,
             random = ~1|key,
             data = df, na.action=na.omit)
  coef_table <- summary(fit)$tTable
  start_index <- grep("x:", rownames(coef_table))[1]
  end_index <- start_index+1
  output <- coef_table[c(2,start_index:end_index), c(1,2,4,5)]
  
  fit_reduced <- lme(reduced_model,
                     random = ~1|key,
                     data = df, na.action = na.omit)
  overall_sig <- anova(update(fit, . ~ ., method = "ML"),
                       update(fit_reduced, . ~ ., method = "ML"))
  output <- cbind(output,
                  c(overall_sig$`p-value`[2], rep(NA, time = 2)))
  output <- cbind(output,
                  c(overall_sig$L.Ratio[2], rep(NA, time = 2)))
  colnames(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_interaction", "L_ratio")
  return(output)
}

fit_lmm_random_intercept_knot1_outcome_covars_ddenum <- function(x, y, key, visit, covariates_df){
  df <- as.data.frame(cbind(key, x, y, visit))
  names(df) <- c("key", "x", "y", "visit")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  df$x <- as.numeric(x)
  model <- as.formula(paste("y ~ x + dde_cat_num*bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) %in% c("key", "dde_cat_num")], collapse = " + ")))
  reduced_model <- as.formula(paste("y ~ dde_cat_num*bSpline(visit, knots = 1, degree = 1)", "+", paste(colnames(covariates_df)[! colnames(covariates_df) %in% c("key", "dde_cat_num")], collapse = " + ")))
  fit <- lme(model,
             random = ~1|key,
             data = df, na.action=na.omit)
  coef_table <- summary(fit)$tTable
  # start_index <- grep("x", rownames(coef_table))[1]
  # end_index <- start_index+1
  output <- coef_table[2, c(1,2,4,5)]
  
  fit_reduced <- lme(reduced_model,
                     random = ~1|key,
                     data = df, na.action = na.omit)
  overall_sig <- anova(update(fit, . ~ ., method = "ML"),
                       update(fit_reduced, . ~ ., method = "ML"))
  output <- c(output,
              overall_sig$`p-value`[2])
  output <- c(output,
              overall_sig$L.Ratio[2])
  names(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_anova", "L_ratio")
  return(output)
}


fit_lm_covars <- function(x, y, covariates_df){
  df <- as.data.frame(cbind(x, y))
  names(df) <- c("x", "y")
  df <- as.data.frame(cbind(df, covariates_df))
  df$y <- as.numeric(df$y)
  model <- as.formula(paste("y ~ x +", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  reduced_model <- as.formula(paste("y ~ ", paste(colnames(covariates_df)[! colnames(covariates_df) == "key"], collapse = " + ")))
  fit <- lm(model,
            data = df, na.action=na.omit)
  coef_table <- summary(fit)$coefficient

  output <- as.data.frame(coef_table)[2,]
  
  fit_reduced <- lm(reduced_model,
                    data = df, na.action = na.omit)
  overall_sig <- anova(fit, fit_reduced)
  
  output <- cbind(output,
                  overall_sig$`Pr(>F)`[2])
  output <- cbind(output,
                  overall_sig$`F`[2])
  colnames(output) <- c("coef", "se", "statistic", "pvalue", "pvalue_anova", "F_stats")
  return(output)
}


fit_lmm_random_intercept_covars_AT_plasma <- function(df, var, feature, covars)
{
  model <- as.formula(paste(paste0("X", feature), " ~ ", paste0("X", var), "*bSpline(visit, knots = 1, degree = 1) + dde_cat_num", "+", paste(covars, collapse = " + "), "+ (1|key)"))
  reduced_model <- as.formula(paste(paste0("X", feature), " ~ ", paste0("X", var), "+bSpline(visit, knots = 1, degree = 1) + dde_cat_num", "+", paste(covars, collapse = " + "), "+ (1|key)"))
  fit <- lmer(model,
              data = df, na.action=na.omit)
  coef_table <- summary(fit)$coefficients
  conf <- confint(fit)
  start_index <- grep("X", rownames(coef_table))[2]
  end_index <- start_index+1
  output <- cbind(
    coef_table[c(start_index:end_index), c(1,2,4,5)],
    conf[c((start_index+2):(end_index+2)),]
  )
  
  fit_reduced <- lmer(reduced_model,
                      data = df, na.action = na.omit)
  overall_sig <- anova(fit,
                       fit_reduced)
  output <- cbind(output,
                  rep(overall_sig$`Pr(>Chisq)`[2], time = 2))
  output <- cbind(output,
                  rep(overall_sig$Chisq[2], time = 2))
  colnames(output) <- c("estimate", "se", "statistic", "pvalue", "lower", "upper", "pvalue_interaction", "chisq")
  
  output <- as.data.frame(output)
  
  return(output)
}

summarize_AT_plasma <- function(df, var, colnames_omics, covariates)
{
  var <- as.character(var)
  df <- data.frame(df)
  
  cl <- makeCluster(detectCores()-2)
  registerDoParallel(cl)
  model_output <- foreach(
    i = 1:length(colnames_omics),
    .combine = rbind,
    .packages = c("lme4", "lmerTest", "splines2"),
    .export=c("fit_lmm_random_intercept_covars_AT_plasma")) %dopar%
    fit_lmm_random_intercept_covars_AT_plasma(
      df = df,
      var = var,
      feature = colnames_omics[i],
      covars = covariates)
  stopCluster(cl)
  
  model_output <- cbind(model_output, rep(c(1,2), time = length(colnames_omics)))
  model_output <- cbind(model_output, rep(colnames_omics, each = 2))
  model_output <- cbind(model_output, var)
  colnames(model_output)[9:10] <- c("period", "feature_name")
  
  model_output <- as_tibble(model_output) %>% 
    select(var, feature_name, period, everything())
  
  return(model_output)
}

add_CI <- function(df)
{
  if(any(grepl("estimate", colnames(df))))
  {
    colnames(df)[colnames(df) == "estimate"] <- "coef"
  }
  t_score <- qt(p = p_cutoff_sig_feat/2, df = 60-1, lower.tail = F)
  margin_error <- t_score*as.numeric(df$se)
  lower_bound <- as.numeric(df$coef) - margin_error
  upper_bound <- as.numeric(df$coef) + margin_error
  
  df$coef_num <- as.numeric(df$coef)
  df <- Reduce(
    cbind,
    list(df, lower_bound, upper_bound))
  colnames(df)[length(colnames(df))-1] <- "lower_bound"
  colnames(df)[length(colnames(df))] <- "upper_bound"
  return(df)
}

create_two_digits <- function(x)
{
  x <- format(round(x, 2), nsmall = 2)
}
