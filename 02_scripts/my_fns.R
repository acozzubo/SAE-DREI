# ---
# title: "Functions Script for Simulation for ACS, DREI SAE 2023 Project"
# subtitle: "DREI2023_SAE"
# author: "[Angelo Cozzubo] (https://sites.google.com/pucp.pe/acozz)"
# date creation: "13/April/2023"
# date last edit: "`r Sys.Date()`"
# ---

# 1. Draw samples function ####
draw_sts_samples <- function(n_samples, strata_vars, strata_sizes, sampling_frame){
  #Draw n stratified samples from frame and store in list   
  #input: number of samples (n_samples), 
  #       vector of strata variables as strings (strata_vars),
  #       vector of sample size for each stratum (strata_sizes)
  #       sampling frame (sampling_frame) SHOULD BE ORDERED BY STRATA VARS!
  #return: list of dataframes with drawn samples
  
  #Order df, required by strata
  sampling_frame <- sampling_frame %>% 
    arrange_at(strata_vars)
    
  #Replications 
  ## different seed for each sample
  reps_df <- future_replicate(n_samples, 
                              strata(sampling_frame, 
                                     strata_vars,  
                                     size = strata_sizes,
                                     method = "srswor"))
  
  ## list of sampled df from replications 
  samples <- list()
  for (k in 1:n_samples) {
    samples[[k]] <- getdata(sampling_frame, data.frame(reps_df[,k]))
  }
  
  return(samples)
}


# 2. Pre-process list of strata_size  ####
fix_list <- function(orig_list, popsize_list){
  #Fix list of strata sample sizes: round, replace 0s, set to min pop size    
  #input: original list in the required order by strata (orig_list), 
  #       vector of strata pop size in same order than orig_list (popsize_list),
  #return: fixed list 
  
  rnd_list <- round(orig_list)
  nozero_list <- ifelse(rnd_list %in% c(0, 1), 2, rnd_list) # replace samples = 0 or 1 with 2 (to compute variances) 
  fixed_list <- pmin(nozero_list, popsize_list)   # if pop. size < sample, replace with pop value
  ## harder fix https://stackoverflow.com/questions/14735411/strata-from-sampling-returns-an-error-arguments-imply-differing-number-of-r
  
  return(fixed_list)
}


# 3. Obtain direct estimates
direct_estim <- function(list_dfs, outcome_finalvar, grouping_var,
                         ID_unit=ID_unit, Prob=Prob, Stratum=Stratum){
  #Obtain survey direct estimates, s.e., CV for a list of samples     
  #input: list_dfs (list of samples as dataframes)
  #       outcome_finalvar: outcome variable to obtain mean
  #       grouping_var: grouping var, areas 
  #       ID_unit: unique identifier 
  #       Prob: probability of selection
  #       Stratum: design strata
  #return: df with estimates, s.e., CV for each df in list_dfs  

  means_df <- list()
  for (i in 1:length(list_dfs)) {

    svy <- svydesign(id=list_dfs[[i]]$ID_unit, 
                     strata=list_dfs[[i]]$Stratum, 
                     nest=TRUE,
                     weights = (list_dfs[[i]]$Prob)^-1, #inverse selec. probability
                     data=list_dfs[[i]])
    
    svy_means <- svyby(formula=as.formula(paste("~", outcome_finalvar)), 
                       by=as.formula(paste("~", paste(grouping_var, collapse= "+"))),  
                       design=svy, 
                       FUN=svymean, 
                       na.rm=T) 
    svy_means["CV"] <- (svy_means$se/svy_means[outcome_finalvar])*100
    svy_means["var"] <- (svy_means$se)^2
    
    if (i>1){svy_means <- svy_means %>% select(-grouping_var)} #not repeat IDs
    means_df <- c(means_df, svy_means) #populate list 
  }
  return(as.data.frame(means_df))
}

#3.5.1 Fix zero variances in direct estaimtes (due to small samplesize)
fix_zero_var_simple <- function(estimates_df, samplesize_list) {
  #Replace the direct estimates zero variances with p(1-p)/n, taking p=.5 for max variance      
  #input: estimates_df: dataframe with direct estimates
  #       samplesize_list: list with sample size for each domain (fix_list fn)
  #return: df with replaced zero variances
  
  
  estimates_df <- estimates_df[order(estimates_df$race, estimates_df$ST), ]
  estimates_df <- cbind(estimates_df, samplesize = samplesize_list)
  
  estimates_df <- estimates_df %>%
    mutate(across(starts_with("var"), ~ ifelse(. %in% c(0), 0.25/samplesize, .)))
  
  estimates_df <- estimates_df %>%
    mutate(across(starts_with("se"), ~ ifelse(. %in% c(0), sqrt(0.25/samplesize), .)))
  
  estimates_df <- subset(estimates_df, select = -samplesize)
  
  return(estimates_df)
  
}

#3.5.2 Fix zero variances in direct estimates using deff
fix_zero_var <- function(estimates_df, samplesize_list, niter=499) {
  #Replace the direct estimates zero variances with p(1-p)/n, taking p=.5 for max variance      
  #input: estimates_df: dataframe with direct estimates
  #       samplesize_list: list with sample size for each domain (fix_list fn)
  #return: df with replaced zero variances
  
  #sort and merge sample size 
  estimates_df <- estimates_df[order(estimates_df$race, estimates_df$ST), ]
  estimates_df <- cbind(estimates_df, samplesize = samplesize_list)
  
  for (i in 1:niter) {
    var_column_name <- paste("var.", i, sep = "")
    house_ten_column_name <- paste("house_ten.", i, sep = "")
    deff_column_name <- paste("deff.", i, sep = "")

    estimates_df[[deff_column_name]] <- (estimates_df[[var_column_name]] * estimates_df[["samplesize"]]) / 
      ((1-estimates_df[[house_ten_column_name]])*estimates_df[[house_ten_column_name]])
    # var(p)*n/(1-p)(p)
    
  }
  
  # Compute the average of columns starting with "deff." Add the average as new column
  estimates_df[sapply(estimates_df, is.infinite)] <- NA # replace all Infs with NA
  estimates_df$average_deff <- rowMeans(estimates_df[, grep("^deff\\.", names(estimates_df), value = TRUE)], na.rm = TRUE)
  
  # Drop the columns starting with "deff."
  estimates_df <- estimates_df[, !names(estimates_df) %in% grep("^deff\\.", names(estimates_df))]
  
  # Calculate the average_deff by "race" and merge
  deff_race <- aggregate(average_deff ~ race, data = subset(estimates_df, samplesize > 15), FUN = mean)
  estimates_df <- estimates_df[, !names(estimates_df) %in% "average_deff"]
  estimates_df <- merge(estimates_df, deff_race, by = "race", all.x = TRUE)
  
  #Calculate the average_estimate by "race" and merge (for replacing 0, 1)
  estimates_df$average_estim <- rowMeans(estimates_df[, grep("^house_ten\\.", names(estimates_df), value = TRUE)]) #avg across iterations
  estim_race <- aggregate(average_estim ~ race, data = subset(estimates_df, samplesize > 15), FUN = mean) #average by race
  estimates_df <- estimates_df[, !names(estimates_df) %in% "average_estim"]
  estimates_df <- merge(estimates_df, estim_race, by = "race", all.x = TRUE)
  
  for (i in 1:niter) {
    CV_column_name <- paste("CV.", i, sep = "")
    se_column_name <- paste("se.", i, sep = "")
    var_column_name <- paste("var.", i, sep = "")
    house_ten_column_name <- paste("house_ten.", i, sep = "")
    house_ten_column_name_aux <- paste("house_ten_aux.", i, sep = "")
    
    # drop variance, std error and CVs
    estimates_df <- estimates_df[, !names(estimates_df) %in% c(var_column_name, se_column_name, CV_column_name)]
    
    # generate auxiliary columns with estimate, as they can be replaced if 0 or 1
    estimates_df[[house_ten_column_name_aux]] <- estimates_df[[house_ten_column_name]] 
    
    # replace estimates 0,1 for average race estim
    rows_to_replace <- estimates_df[[house_ten_column_name_aux]] %in% c(0, 1)
    estimates_df[rows_to_replace, house_ten_column_name_aux] <- estimates_df[rows_to_replace, "average_estim"]
    
    # Compute variance, std error and CVs using the auxiliary column
    estimates_df[[var_column_name]] <- ((1-estimates_df[[house_ten_column_name_aux]])*estimates_df[[house_ten_column_name_aux]]) * 
                                       (estimates_df[["average_deff"]] / estimates_df[["samplesize"]])
                                       # (1-p)(p)*deff/n 
    estimates_df[[se_column_name]] <- sqrt(estimates_df[[var_column_name]])
    estimates_df[[CV_column_name]] <- estimates_df[[se_column_name]] / estimates_df[[house_ten_column_name]]
    
    # drop auxiliary column 
    estimates_df <- estimates_df[, !names(estimates_df) %in% c(house_ten_column_name_aux)]
    
  }
  
  estimates_df <- subset(estimates_df, select = -samplesize)
  return(estimates_df)
}

# 4. Obtain Fay-Herriot estimates
FH_estim <- function(dir_estimates_df, covars_df, 
                     covars_names, grouping_var, outcome_finalvar, FE=F){
  #Obtain FH estimates, synth estimates, mse and CV for list of direct estimates     
  #input: dir_estimates_df: dataframe with direct estimates from direct_estim fn
  #       covars_df: dataframe that includes covariates and area IDs
  #       covars_names: list of covariates name to subset covars_df
  #       outcome_finalvar: outcome variable to obtain mean
  #       grouping_var: grouping var, areas 
  #return: df with FH estimates, synth estimates, mse and CV   
  
  # merge and order covariates to the direct estimates
  dir_estimates_df <- inner_join(dir_estimates_df, covars_df, 
                                 by = grouping_var, na_matches = "never") 
  
  # select the covariates columns + constant for the FH
  covars_filtered <- data.matrix(cbind(const=rep(1, dim(dir_estimates_df)[1]), 
                                       dir_estimates_df[,(names(dir_estimates_df) %in% covars_names)]))
  
  if(FE){
    design_matrix <- model.matrix(~ as.factor(race) - 1, data=dir_estimates_df)
    covars_filtered <- cbind(covars_filtered, design_matrix[, -ncol(design_matrix)]) #drop last column for multicollineraity
  }
  
  # estimates and variances dfs 
  estimates <- dir_estimates_df[, grepl(paste0("^", outcome_finalvar), names(dir_estimates_df))]
  variances <- dir_estimates_df[, grepl(paste0("^", "var"), names(dir_estimates_df))]
  
  # loop thorugh every estimate 
  FH_df <- data.frame(matrix(nrow = dim(estimates)[1], ncol = 0)) #empty df
  for (col_index in seq_along(estimates)) {

    # compute FH 
    area.FH.res <- eblupFH(estimates[[col_index]]~(covars_filtered)-1,
                           vardir = (variances[[col_index]]),
                           method = "REML", MAXITER = 100, PRECISION = 0.0001)
    ## FH
    area.FH <- area.FH.res$eblup
    
    # ## synth estimator 
    # area.rsyn1 <- covars_filtered%*%area.FH.res$fit$estcoef[,1]
    
    ## MSE and CV 
    area.FH.mse.res <- mseFH(estimates[[col_index]]~(covars_filtered)-1,
                             vardir = (variances[[col_index]]),
                             method = "REML", MAXITER = 100, PRECISION = 0.0001)
    area.FH.mse <- area.FH.mse.res$mse
    area.FH.cv <- 100*sqrt(area.FH.mse)/area.FH
    
    aux_df <- data.frame("FH.estimates" = area.FH,
                         # "FH.synth" = area.rsyn1,
                         "FH.mse" = area.FH.mse,
                         "FH.cv" = area.FH.cv)
    
    #fix names with suffix from the rep
    new_colnames <- paste0(names(aux_df), ".", col_index)
    aux_df <- setNames(aux_df, new_colnames)
    
    FH_df <- cbind(FH_df, aux_df)
  }
  
  # paste grouping vars
  FH_df <- cbind(dir_estimates_df[, grouping_var], FH_df)
  return(FH_df)
}


#5. Compute MSE from replicates estimates and population value
MSE <- function(estimates_df, outcome_df, grouping_var, outcome_finalvar, FH=F, diff=F) {
  #Obtain MSE estimate for each domain      
  #input: estimates_df: dataframe with estimates from direct_estim or FH_estim fns
  #       outcome_df: dataframe with population value of parameter 
  #       outcome_finalvar: outcome variable to obtain mean
  #       grouping_var: grouping var, areas
  #       FH: T if computing the MSE for FH results (renaming)
  #       diff: T if computing the MSE for race difference (sort only 1 var)
  #return: df with MSE and domain IDs    
  
  if(FH){
    estimates_df <- estimates_df %>%
                    select(starts_with("FH.estimates"), all_of(grouping_var)) %>%
                    rename_with(~paste0("house_ten_", .), starts_with("FH.estimates"))
  }
  
  # Create a new dataframe with selected replications, and grouping vars 
  mse_df <- data.frame(estimates_df[, grouping_var],
                       estimates_df[, grep(paste0("^", outcome_finalvar), names(estimates_df))])
  
  # Order dfs
  if (diff) { #MSE for difference between races (only sort by ST)
    mse_df <- mse_df[order(mse_df$ST), ]
    outcome_df <- outcome_df[order(outcome_df$ST), ]
    names(outcome_df)[names(outcome_df) == "difference"] <- outcome_finalvar
  } else {
    mse_df <- mse_df[order(mse_df$race, mse_df$ST), ]
    outcome_df <- outcome_df[order(outcome_df$race, outcome_df$ST), ]
  }
  
  # Subtract population value of outcome from selected columns in new_dataframe
  mse_df[, -(1:length(grouping_var))] <- (mse_df[, -(1:length(grouping_var))] - 
                                          outcome_df[[outcome_finalvar]])^2
  
  # Average the difference and keep the final MSE only 
  mse_df <- mse_df %>%
    mutate(mse = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)) %>%
    select(all_of(grouping_var), mse)
  
  if(FH){mse_df <- mse_df %>% rename(mse_fh = mse)} #rename to compare with direct estim
  
  return(mse_df) 
}

#5.0 export MSE (for differences estimates)
MSE_direct_export <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(MSE_eq_1 = (mse_eq_1$mse),
                          MSE_eq_gss = (mse_eq_gss$mse),
                          MSE_pr_1 = (mse_pr_1$mse), 
                          MSE_pr_05 = (mse_pr_05$mse), 
                          MSE_pr_gss = (mse_pr_gss$mse))  
  
  overall_summ <- aggregate(cbind(MSE_eq_1, MSE_eq_gss, MSE_pr_1, MSE_pr_05, MSE_pr_gss) ~ NULL, 
                            data = summ_MSEs, 
                            FUN = function(x) c(Median = median(x)))

  return(overall_summ)
}

MSE_export <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(MSE_eq_1 = (mse_eq_1$mse_fh),
                          MSE_eq_gss = (mse_eq_gss$mse_fh),
                          MSE_pr_1 = (mse_pr_1$mse_fh), 
                          MSE_pr_05 = (mse_pr_05$mse_fh), 
                          MSE_pr_gss = (mse_pr_gss$mse_fh))  
  
  overall_summ <- aggregate(cbind(MSE_eq_1, MSE_eq_gss, MSE_pr_1, MSE_pr_05, MSE_pr_gss) ~ NULL, 
                            data = summ_MSEs, 
                            FUN = function(x) c(Median = median(x)))
  
  return(overall_summ)
}

#5.1 export MSE
MSE_direct_race <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$race, mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$race, mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$race, mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$race, mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$race, mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(race = mse_eq_1$race, 
                          MSE_eq_1 = (mse_eq_1$mse),
                          MSE_eq_gss = (mse_eq_gss$mse),
                          MSE_pr_1 = (mse_pr_1$mse), 
                          MSE_pr_05 = (mse_pr_05$mse), 
                          MSE_pr_gss = (mse_pr_gss$mse))  
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(MSE_eq_1, MSE_eq_gss, MSE_pr_1, MSE_pr_05, MSE_pr_gss) ~ race, 
                         data = summ_MSEs, 
                         FUN = function(x) c(Median = median(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(MSE_eq_1, MSE_eq_gss, MSE_pr_1, MSE_pr_05, MSE_pr_gss) ~ NULL, 
                            data = summ_MSEs, 
                            FUN = function(x) c(Median = median(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}

MSE_export_race <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$race, mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$race, mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$race, mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$race, mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$race, mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(race = mse_eq_1$race, 
                          MSE_eq_1 = (mse_eq_1$mse_fh),
                          MSE_eq_gss = (mse_eq_gss$mse_fh),
                          MSE_pr_1 = (mse_pr_1$mse_fh), 
                          MSE_pr_05 = (mse_pr_05$mse_fh), 
                          MSE_pr_gss = (mse_pr_gss$mse_fh))  
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(MSE_eq_1, MSE_eq_gss, MSE_pr_1, MSE_pr_05, MSE_pr_gss) ~ race, 
                         data = summ_MSEs, 
                         FUN = function(x) c(Median = median(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(MSE_eq_1, MSE_eq_gss, MSE_pr_1, MSE_pr_05, MSE_pr_gss) ~ NULL, 
                            data = summ_MSEs, 
                            FUN = function(x) c(Median = median(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}

#5.2 MSE improvements by race and overall
MSE_improve <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(state = mse_eq_1$ST, 
                          improve_eq_1 = (mse_eq_1$mse_fh - mse_eq_1$mse)/mse_eq_1$mse,
                          improve_eq_gss = (mse_eq_gss$mse_fh - mse_eq_gss$mse)/mse_eq_gss$mse,
                          improve_pr_1 = (mse_pr_1$mse_fh - mse_pr_1$mse)/mse_pr_1$mse,
                          improve_pr_05 = (mse_pr_05$mse_fh - mse_pr_05$mse)/mse_pr_05$mse,
                          improve_pr_gss = (mse_pr_gss$mse_fh - mse_pr_gss$mse)/mse_pr_gss$mse) 
  
  overall_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ NULL, 
                            data = summ_MSEs,
                            FUN = function(x) c(Median = median(x)))
                            # FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(overall_summ)
}


#5.3 MSE improvements by race and overall
MSE_improve_race <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$race, mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$race, mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$race, mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$race, mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$race, mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(race = mse_eq_1$race, 
                          improve_eq_1 = (mse_eq_1$mse_fh - mse_eq_1$mse)/mse_eq_1$mse,
                          improve_eq_gss = (mse_eq_gss$mse_fh - mse_eq_gss$mse)/mse_eq_gss$mse,
                          improve_pr_1 = (mse_pr_1$mse_fh - mse_pr_1$mse)/mse_pr_1$mse,
                          improve_pr_05 = (mse_pr_05$mse_fh - mse_pr_05$mse)/mse_pr_05$mse,
                          improve_pr_gss = (mse_pr_gss$mse_fh - mse_pr_gss$mse)/mse_pr_gss$mse) 
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ race, 
                         data = summ_MSEs, 
                         FUN = function(x) c(Median = median(x)))
                          # FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ NULL, 
                            data = summ_MSEs, 
                            FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}


#6. Scatterplot for MSE Direct vs FH, by race and sample size 
plotMSE <- function(mse_df, samplesize_list, subtitle, covars=T, FE=F) {
  #Scatterplot for MSE Direct vs FH, by race and sample size      
  #input: mse_df: dataframe with direct and FH mse (MSE fn)
  #       samplesize_list: list with sample size for each domain (fix_list fn) 
  #       subtitle: string with specific subtitle
  #return: df with MSE and domain IDs    
  
  #order 
  mse_df <- mse_df[order(mse_df$race, mse_df$ST), ]
  
  mse_df <- cbind(mse_df, samplesize = samplesize_list)
  title1 <- "MSE: Direct vs FH estimates,\n"
  
  covar_lbl <- "covariates"
  FE_lbl <- "race FE,"
  
  title2 <- paste0("w/ ", covar_lbl, " & no ", FE_lbl, "\n") #default
  if(!covars && !FE){ title2 <- paste0("no ", covar_lbl, " & no ", FE_lbl, "\n")}
  if(!covars && FE){ title2 <- paste0("no ", covar_lbl, " & ", FE_lbl, "\n")}
  if(covars && FE){ title2 <- paste0("w/ ", covar_lbl, " & ", FE_lbl, "\n")}
  
  labels_race = c("Nonhispanic White", "Nonhispanic black", 
                  "Hispanic", "Other race or multirace")
  
  scatter <- ggplot(mse_df, 
                    aes(x = mse_fh, y = mse, 
                        color = as.factor(race),
                        size = samplesize)) + # , size = n for size of bubble
    geom_point(alpha=0.7) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(title = paste0(title1, title2, subtitle), # title,  axis labels
         x = "FH MSE", y = "Direct MSE") + 
    scale_color_manual(values = cbPalette, labels = labels_race) +
    labs(color = "Race")
  
  return(scatter) 
}

#6.1 Scatterplot for MSE Direct vs FH, by race and sample size 
plot_diffMSE <- function(mse_df, subtitle, covars=T, FE=F) {
  #Scatterplot for MSE Direct vs FH, by race and sample size      
  #input: mse_df: dataframe with direct and FH mse (MSE fn)
  #       samplesize_list: list with sample size for each domain (fix_list fn) 
  #       subtitle: string with specific subtitle
  #return: df with MSE and domain IDs    
  
  #order 
  mse_df <- mse_df[order(mse_df$ST), ]
  
  title1 <- "MSE: Direct vs FH estimates (diff. in HH. tenancy),\n"
  
  covar_lbl <- "covariates"
  FE_lbl <- "race FE,"
  
  title2 <- paste0("w/ ", covar_lbl, " & no ", FE_lbl, "\n") #default
  if(!covars && !FE){ title2 <- paste0("no ", covar_lbl, " & no ", FE_lbl, "\n")}
  if(!covars && FE){ title2 <- paste0("no ", covar_lbl, " & ", FE_lbl, "\n")}
  if(covars && FE){ title2 <- paste0("w/ ", covar_lbl, " & ", FE_lbl, "\n")}
  
  scatter <- ggplot(mse_df, 
                    aes(x = mse_fh, y = mse)) + # , size = n for size of bubble
    geom_point(alpha=0.7) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(title = paste0(title1, title2, subtitle), # title,  axis labels
         x = "FH MSE", y = "Direct MSE") + 
    scale_color_manual(values = cbPalette) 
  
  return(scatter) 
}

#7. Scatterplot for MSE ratio (Direct vs FH) vs sample size, by race  
plotMSE2 <- function(mse_df, samplesize_list, subtitle, covars=T, FE=F){
  #Scatterplot for MSE Direct vs FH, by race and sample size      
  #input: mse_df: dataframe with direct and FH mse (MSE fn)
  #       samplesize_list: list with sample size for each domain (fix_list fn) 
  #       subtitle: string with specific subtitle
  #return: df with MSE and domain IDs
  
  #order 
  mse_df <- mse_df[order(mse_df$race, mse_df$ST), ]
  
  mse_df <- cbind(mse_df, samplesize = samplesize_list)
  mse_df$ln_samplesize <- log(mse_df$samplesize) 
  mse_df$ratio <- mse_df$mse/mse_df$mse_fh 
  title1 <- "MSE ratios: Direct/FH estimates,\n"
  
  covar_lbl <- "covariates"
  FE_lbl <- "race FE,"
  
  title2 <- paste0("w/ ", covar_lbl, " & no ", FE_lbl, "\n") #default
  if(!covars && !FE){ title2 <- paste0("no ", covar_lbl, " & no ", FE_lbl, "\n")}
  if(!covars && FE){ title2 <- paste0("no ", covar_lbl, " & ", FE_lbl, "\n")}
  if(covars && FE){ title2 <- paste0("w/ ", covar_lbl, " & ", FE_lbl, "\n")}
  
  labels_race = c("Nonhispanic White", "Nonhispanic black", 
                  "Hispanic", "Other race or multirace")
  
  scatter <- ggplot(mse_df, 
                    aes(x = ln_samplesize, y = ratio, 
                        color = as.factor(race))) + # , size = n for size of bubble
    geom_point(size=2, alpha=0.7) + 
    geom_abline(intercept = 1, slope = 0, linetype = "dashed", color = "black") +
    labs(title = paste0(title1, title2, subtitle), # title,  axis labels
         x = "log sample size", y = "MSE ratio (direct vs FH)") + 
    scale_color_manual(values = cbPalette, labels = labels_race) +
    labs(color = "Race")
  
  return(scatter) 
}

# 8. Estimates convergence plot  ####
converplot <- function(popsize, fh_estimates, direct_estimates, grouping_var, title, xlabel, ylabel, note) {
  #Draw a convergence plot for ratio of direct/model estimate or std errors   
  #input: domains pop size, estimates ratio column, title,labels and notes string
  #return: none, use jpg command header to save plot 
  
  par(mar = c(7, 4, 3, 3))
  
  fh_df <- fh_estimates %>%
    select(starts_with("FH.estimates"), all_of(grouping_var)) %>%
    rename_with(~paste0("house_ten_", .), starts_with("FH.estimates"))
  
  fh_df <- fh_df %>%
    mutate(avg_estim = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)) %>%
    select(all_of(grouping_var), avg_estim)
  
  direct_df <- direct_estimates %>%
    select(starts_with("house_ten"), all_of(grouping_var)) 
  
  direct_df <- direct_df %>%
    mutate(avg_direct = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)) %>%
    select(all_of(grouping_var), avg_direct)
  
  df <- inner_join(fh_df, direct_df, by = grouping_var, na_matches = "never")
  
  df$estim_ratio_col <- (df$avg_direct/df$avg_estim)*100 
  
  df <- df[order(df$race, df$ST), ]
  df$race <- factor(df$race)
  
  race_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
  labels_race = c("Nonhispanic White", "Nonhispanic black", 
                  "Hispanic", "Other race or multirace")
  
  plot(popsize, df$estim_ratio_col, 
       log="x", xaxt = 'n', 
       main=title,  
       xlab= xlabel,  
       ylab= ylabel, 
       pch=19,
       col=race_colors[df$race])
  abline(h=100, col="red")
  legend("topright", legend=labels_race, fill=race_colors)
  
  myTicks = axTicks(1)
  axis(1, at = myTicks, labels = formatC(myTicks, format = 'd'))
  mtext(note, side = 1, line = 5, cex = 0.7, adj = 0) 

}

# 8. Estimates convergence plot  ####
converplot_std <- function(popsize, fh_estimates, direct_estimates, grouping_var, title, xlabel, ylabel, note) {
  #Draw a convergence plot for ratio of direct/model estimate or std errors   
  #input: domains pop size, estimates ratio column, title,labels and notes string
  #return: none, use jpg command header to save plot 
  
  par(mar = c(7, 4, 3, 3))
  
  fh_df <- fh_estimates %>%
    select(starts_with("FH.mse"), all_of(grouping_var)) %>%
    rename_with(~paste0("house_ten_", .), starts_with("FH.mse"))
  
  fh_df <- fh_df %>%
    mutate(avg_mse = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)) %>%
    select(all_of(grouping_var), avg_mse)
  
  direct_df <- direct_estimates %>%
    select(starts_with("se"), all_of(grouping_var)) 
  
  direct_df <- direct_df %>%
    mutate(avg_se = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)) %>%
    select(all_of(grouping_var), avg_se)
  
  df <- inner_join(fh_df, direct_df, by = grouping_var, na_matches = "never")
  
  df$estim_ratio_col <- (df$avg_se/sqrt(df$avg_mse))*100 #
  
  df <- df[order(df$race, df$ST), ]
  df$race <- factor(df$race)
  
  race_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
  labels_race = c("Nonhispanic White", "Nonhispanic black", 
                  "Hispanic", "Other race or multirace")
  
  plot(popsize, df$estim_ratio_col, 
       log="x", xaxt = 'n', 
       main=title,  
       xlab= xlabel,  
       ylab= ylabel, 
       pch=19,
       col=race_colors[df$race])
  abline(h=100, col="red")
  legend("topright", legend=labels_race, fill=race_colors)
  
  myTicks = axTicks(1)
  axis(1, at = myTicks, labels = formatC(myTicks, format = 'd'))
  mtext(note, side = 1, line = 5, cex = 0.7, adj = 0) 
}


# ----- 2024 Additions - cv ----- #

#9. Compute Absolute from  estimates and population value
AbsBias <- function(estimates_df, outcome_df, grouping_var, outcome_finalvar, FH=F, diff=F) {
  #Obtain AbsBias estimate for each domain      
  #input: estimates_df: dataframe with estimates from direct_estim or FH_estim fns
  #       outcome_df: dataframe with population value of parameter 
  #       outcome_finalvar: outcome variable to obtain mean
  #       grouping_var: grouping var, areas
  #       FH: T if computing the for FH results (renaming)
  #       diff: T if computing for race difference (sort only 1 var)
  #return: df with AbsBias and domain IDs    
  
  if(FH){
    estimates_df <- estimates_df %>%
      select(starts_with("FH.estimates"), all_of(grouping_var)) %>%
      rename_with(~paste0("house_ten_", .), starts_with("FH.estimates"))
  }
  
  # Create a new dataframe with selected replications, and grouping vars 
  AB_df <- data.frame(estimates_df[, grouping_var],
                       estimates_df[, grep(paste0("^", outcome_finalvar), names(estimates_df))])
  
  # Order dfs
  if (diff) { #Abs Bias for difference between races (only sort by ST)
    AB_df <- AB_df[order(AB_df$ST), ]
    outcome_df <- outcome_df[order(outcome_df$ST), ]
    names(outcome_df)[names(outcome_df) == "difference"] <- outcome_finalvar
  } else {
    AB_df <- AB_df[order(AB_df$race, AB_df$ST), ]
    outcome_df <- outcome_df[order(outcome_df$race, outcome_df$ST), ]
  }
  
  # Subtract population value of outcome from selected columns in new_dataframe and take absolute value 
  AB_df[, -(1:length(grouping_var))] <- abs(AB_df[, -(1:length(grouping_var))] - 
                                            outcome_df[[outcome_finalvar]])
  
  # Average the difference and keep the final AbsBias only 
  AB_df <- AB_df %>%
    mutate(Abs_Bias = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)) %>%
    select(all_of(grouping_var),Abs_Bias)
  
  if(FH){AB_df <- AB_df %>% rename(Abs_Bias_FH = Abs_Bias)} #rename to compare with direct estim
  
  return(AB_df) 
}


#9.0 export Absolute Bias (for differences estimates)
AbsBias_direct_export <- function(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss){
  #Obtain csv file with summary of the Absolute Bias improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with Absolute Bias improvements by race and overall     
  
  AbsBias_eq_1 <- AbsBias_eq_1[order(AbsBias_eq_1$ST), ]
  AbsBias_eq_gss <- AbsBias_eq_gss[order(AbsBias_eq_gss$ST), ]
  AbsBias_pr_1 <- AbsBias_pr_1[order(AbsBias_pr_1$ST), ]
  AbsBias_pr_05 <- AbsBias_pr_05[order(AbsBias_pr_05$ST), ]
  AbsBias_pr_gss <- AbsBias_pr_gss[order(AbsBias_pr_gss$ST), ]
  
  summ_AbsBiases <- data.frame(AbsBias_eq_1 = (AbsBias_eq_1$Abs_Bias),
                          AbsBias_eq_gss = (AbsBias_eq_gss$Abs_Bias),
                          AbsBias_pr_1 = (AbsBias_pr_1$Abs_Bias), 
                          AbsBias_pr_05 = (AbsBias_pr_05$Abs_Bias), 
                          AbsBias_pr_gss = (AbsBias_pr_gss$Abs_Bias))  
  
  overall_summ <- aggregate(cbind(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss) ~ NULL, 
                            data = summ_AbsBiases, 
                            FUN = function(x) c(Median = median(x), Average = mean(x), Max = max(x), Min = min(x)) )
  
  return(overall_summ)
}

AbsBias_export <- function(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss){
  #Obtain csv file with summary of the AbsBias improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with AbsBias improvements by race and overall     
  
  AbsBias_eq_1 <- AbsBias_eq_1[order(AbsBias_eq_1$ST), ]
  AbsBias_eq_gss <- AbsBias_eq_gss[order(AbsBias_eq_gss$ST), ]
  AbsBias_pr_1 <- AbsBias_pr_1[order(AbsBias_pr_1$ST), ]
  AbsBias_pr_05 <- AbsBias_pr_05[order(AbsBias_pr_05$ST), ]
  AbsBias_pr_gss <- AbsBias_pr_gss[order(AbsBias_pr_gss$ST), ]
  
  summ_AbsBiases <- data.frame(AbsBias_eq_1 = (AbsBias_eq_1$AbsBias_FH),
                          AbsBias_eq_gss = (AbsBias_eq_gss$AbsBias_FH),
                          AbsBias_pr_1 = (AbsBias_pr_1$AbsBias_FH), 
                          AbsBias_pr_05 = (AbsBias_pr_05$AbsBias_FH), 
                          AbsBias_pr_gss = (AbsBias_pr_gss$AbsBias_FH))  
  
  overall_summ <- aggregate(cbind(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss) ~ NULL, 
                            data = summ_AbsBiases, 
                            FUN = function(x) c(Median = median(x), Average = mean(x), Max = max(x), Min = min(x)))
  
  return(overall_summ)
}

#9.1 export AbsBias
AbsBias_direct_race <- function(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the AbsBias improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with AbsBias improvements by race and overall     
  
  AbsBias_eq_1 <- AbsBias_eq_1[order(AbsBias_eq_1$race, AbsBias_eq_1$ST), ]
  AbsBias_eq_gss <- AbsBias_eq_gss[order(AbsBias_eq_gss$race, AbsBias_eq_gss$ST), ]
  AbsBias_pr_1 <- AbsBias_pr_1[order(AbsBias_pr_1$race, AbsBias_pr_1$ST), ]
  AbsBias_pr_05 <- AbsBias_pr_05[order(AbsBias_pr_05$race, AbsBias_pr_05$ST), ]
  AbsBias_pr_gss <- AbsBias_pr_gss[order(AbsBias_pr_gss$race, AbsBias_pr_gss$ST), ]
  
  summ_AbsBiases <- data.frame(race = AbsBias_eq_1$race, 
                          AbsBias_eq_1 = (AbsBias_eq_1$Abs_Bias),
                          AbsBias_eq_gss = (AbsBias_eq_gss$Abs_Bias),
                          AbsBias_pr_1 = (AbsBias_pr_1$Abs_Bias), 
                          AbsBias_pr_05 = (AbsBias_pr_05$Abs_Bias), 
                          AbsBias_pr_gss = (AbsBias_pr_gss$Abs_Bias))  
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate (cbind(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss) ~ race, 
                         data = summ_AbsBiases, 
                         FUN = function(x) c(Median = median(x), Average = mean(x), Max = max(x), Min = min(x) ))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss) ~ NULL, 
                            data = summ_AbsBiases, 
                            FUN = function(x) c(Median = median(x), Average = mean(x), Max = max(x), Min = min(x))) 
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}



AbsBias_export_race <- function(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the AbsBias improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with AbsBias improvements by race and overall     
  
  AbsBias_eq_1 <- AbsBias_eq_1[order(AbsBias_eq_1$race, AbsBias_eq_1$ST), ]
  AbsBias_eq_gss <- AbsBias_eq_gss[order(AbsBias_eq_gss$race, AbsBias_eq_gss$ST), ]
  AbsBias_pr_1 <- AbsBias_pr_1[order(AbsBias_pr_1$race, AbsBias_pr_1$ST), ]
  AbsBias_pr_05 <- AbsBias_pr_05[order(AbsBias_pr_05$race, AbsBias_pr_05$ST), ]
  AbsBias_pr_gss <- AbsBias_pr_gss[order(AbsBias_pr_gss$race, AbsBias_pr_gss$ST), ]
  
  summ_AbsBiases <- data.frame(race = AbsBias_eq_1$race, 
                          AbsBias_eq_1 = (AbsBias_eq_1$Abs_Bias_FH),
                          AbsBias_eq_gss = (AbsBias_eq_gss$Abs_Bias_FH),
                          AbsBias_pr_1 = (AbsBias_pr_1$Abs_Bias_FH), 
                          AbsBias_pr_05 = (AbsBias_pr_05$Abs_Bias_FH), 
                          AbsBias_pr_gss = (AbsBias_pr_gss$Abs_Bias_FH))  
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss) ~ race, 
                         data = summ_AbsBiases, 
                         FUN = function(x) c(Median = median(x), Average = mean(x), Max = max(x), Min = min(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss) ~ NULL, 
                            data = summ_AbsBiases, 
                            FUN = function(x) c(Median = median(x), Average = mean(x), Max = max(x), Min = min(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}

#9.2 AbsBias improvements by race and overall
AbsBias_improve <- function(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss){
  #Obtain csv file with summary of the AbsBias improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with AbsBias improvements by race and overall     
  
  AbsBias_eq_1 <- AbsBias_eq_1[order(AbsBias_eq_1$ST), ]
  AbsBias_eq_gss <- AbsBias_eq_gss[order(AbsBias_eq_gss$ST), ]
  AbsBias_pr_1 <- AbsBias_pr_1[order(AbsBias_pr_1$ST), ]
  AbsBias_pr_05 <- AbsBias_pr_05[order(AbsBias_pr_05$ST), ]
  AbsBias_pr_gss <- AbsBias_pr_gss[order(AbsBias_pr_gss$ST), ]
  
  summ_AbsBiases <- data.frame(state = AbsBias_eq_1$ST, 
                          improve_eq_1 = (AbsBias_eq_1$Abs_Bias_FH - AbsBias_eq_1$Abs_Bias)/AbsBias_eq_1$Abs_Bias,
                          improve_eq_gss = (AbsBias_eq_gss$Abs_Bias_FH - AbsBias_eq_gss$Abs_Bias)/AbsBias_eq_gss$Abs_Bias,
                          improve_pr_1 = (AbsBias_pr_1$Abs_Bias_FH - AbsBias_pr_1$Abs_Bias)/AbsBias_pr_1$Abs_Bias,
                          improve_pr_05 = (AbsBias_pr_05$Abs_Bias_FH - AbsBias_pr_05$Abs_Bias)/AbsBias_pr_05$Abs_Bias,
                          improve_pr_gss = (AbsBias_pr_gss$Abs_Bias_FH - AbsBias_pr_gss$Abs_Bias)/AbsBias_pr_gss$Abs_Bias) 
  
  overall_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ NULL, 
                            data = summ_AbsBiases,
                            FUN = function(x) c(Median = median(x), Average = mean(x), Max = max(x), Min = min(x)))
  # FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(overall_summ)
}


#9.3 AbsBias improvements by race and overall
AbsBias_improve_race <- function(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the AbsBias improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with AbsBias improvements by race and overall     
  
  AbsBias_eq_1 <- AbsBias_eq_1[order(AbsBias_eq_1$race, AbsBias_eq_1$ST), ]
  AbsBias_eq_gss <- AbsBias_eq_gss[order(AbsBias_eq_gss$race, AbsBias_eq_gss$ST), ]
  AbsBias_pr_1 <- AbsBias_pr_1[order(AbsBias_pr_1$race, AbsBias_pr_1$ST), ]
  AbsBias_pr_05 <- AbsBias_pr_05[order(AbsBias_pr_05$race, AbsBias_pr_05$ST), ]
  AbsBias_pr_gss <- AbsBias_pr_gss[order(AbsBias_pr_gss$race, AbsBias_pr_gss$ST), ]
  
  summ_AbsBiases <- data.frame(race = AbsBias_eq_1$race, 
                          improve_eq_1 = (AbsBias_eq_1$Abs_Bias_FH - AbsBias_eq_1$Abs_Bias)/AbsBias_eq_1$Abs_Bias,
                          improve_eq_gss = (AbsBias_eq_gss$Abs_Bias_FH - AbsBias_eq_gss$Abs_Bias)/AbsBias_eq_gss$Abs_Bias,
                          improve_pr_1 = (AbsBias_pr_1$Abs_Bias_FH - AbsBias_pr_1$Abs_Bias)/AbsBias_pr_1$Abs_Bias,
                          improve_pr_05 = (AbsBias_pr_05$Abs_Bias_FH - AbsBias_pr_05$Abs_Bias)/AbsBias_pr_05$Abs_Bias,
                          improve_pr_gss = (AbsBias_pr_gss$Abs_Bias_FH - AbsBias_pr_gss$Abs_Bias)/AbsBias_pr_gss$Abs_Bias) 
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ race, 
                         data = summ_AbsBiases, 
                         FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  # FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ NULL, 
                            data = summ_AbsBiases, 
                            FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}


#10. Scatterplot for AbsBias Direct vs FH, by race and sample size 
plotAbsBias <- function(AbsBias_df, samplesize_list, subtitle, covars=T, FE=F) {
  #Scatterplot for AbsBias Direct vs FH, by race and sample size      
  #input: AbsBias_df: dataframe with direct and FH AbsBias (AbsBias fn)
  #       samplesize_list: list with sample size for each domain (fix_list fn) 
  #       subtitle: string with specific subtitle
  #return: df with AbsBias and domain IDs    
  
  #order 
  AbsBias_df <- AbsBias_df[order(AbsBias_df$race, AbsBias_df$ST), ]
  
  AbsBias_df <- cbind(AbsBias_df, samplesize = samplesize_list)
  title1 <- "AbsBias: Direct vs FH estimates,\n"
  
  covar_lbl <- "covariates"
  FE_lbl <- "race FE,"
  
  title2 <- paste0("w/ ", covar_lbl, " & no ", FE_lbl, "\n") #default
  if(!covars && !FE){ title2 <- paste0("no ", covar_lbl, " & no ", FE_lbl, "\n")}
  if(!covars && FE){ title2 <- paste0("no ", covar_lbl, " & ", FE_lbl, "\n")}
  if(covars && FE){ title2 <- paste0("w/ ", covar_lbl, " & ", FE_lbl, "\n")}
  
  labels_race = c("Nonhispanic White", "Nonhispanic black", 
                  "Hispanic", "Other race or multirace")
  
  scatter <- ggplot(AbsBias_df, 
                    aes(x = Abs_Bias_FH, y = Abs_Bias, 
                        color = as.factor(race),
                        size = samplesize)) + # , size = n for size of bubble
    geom_point(alpha=0.7) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(title = paste0(title1, title2, subtitle), # title,  axis labels
         x = "FH Absolute Bias", y = "Direct Absolute Bias") + 
    scale_color_manual(values = cbPalette, labels = labels_race) +
    labs(color = "Race")
  
  return(scatter) 
}

#10.1 Scatterplot for AbsBias Direct vs FH, by race and sample size 
plot_diffAbsBias <- function(AbsBias_df, subtitle, covars=T, FE=F) {
  #Scatterplot for AbsBias Direct vs FH, by race and sample size      
  #input: AbsBias_df: dataframe with direct and FH AbsBias (AbsBias fn)
  #       samplesize_list: list with sample size for each domain (fix_list fn) 
  #       subtitle: string with specific subtitle
  #return: df with AbsBias and domain IDs    
  
  #order 
  AbsBias_df <- AbsBias_df[order(AbsBias_df$ST), ]
  
  title1 <- "Absolute Bias: Direct vs FH estimates (diff. in HH. tenancy),\n"
  
  covar_lbl <- "covariates"
  FE_lbl <- "race FE,"
  
  title2 <- paste0("w/ ", covar_lbl, " & no ", FE_lbl, "\n") #default
  if(!covars && !FE){ title2 <- paste0("no ", covar_lbl, " & no ", FE_lbl, "\n")}
  if(!covars && FE){ title2 <- paste0("no ", covar_lbl, " & ", FE_lbl, "\n")}
  if(covars && FE){ title2 <- paste0("w/ ", covar_lbl, " & ", FE_lbl, "\n")}
  
  scatter <- ggplot(AbsBias_df, 
                    aes(x = Abs_Bias_FH, y = Abs_Bias)) + # , size = n for size of bubble
    geom_point(alpha=0.7) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(title = paste0(title1, title2, subtitle), # title,  axis labels
         x = "FH AbsBias", y = "Direct AbsBias") + 
    scale_color_manual(values = cbPalette) 
  
  return(scatter) 
}

#11. Scatterplot for MSE ratio (Direct vs FH) vs sample size, by race  
plotMSE2 <- function(mse_df, samplesize_list, subtitle, covars=T, FE=F){
  #Scatterplot for MSE Direct vs FH, by race and sample size      
  #input: mse_df: dataframe with direct and FH mse (MSE fn)
  #       samplesize_list: list with sample size for each domain (fix_list fn) 
  #       subtitle: string with specific subtitle
  #return: df with MSE and domain IDs
  
  #order 
  mse_df <- mse_df[order(mse_df$race, mse_df$ST), ]
  
  mse_df <- cbind(mse_df, samplesize = samplesize_list)
  mse_df$ln_samplesize <- log(mse_df$samplesize) 
  mse_df$ratio <- mse_df$mse/mse_df$mse_fh 
  title1 <- "MSE ratios: Direct/FH estimates,\n"
  
  covar_lbl <- "covariates"
  FE_lbl <- "race FE,"
  
  title2 <- paste0("w/ ", covar_lbl, " & no ", FE_lbl, "\n") #default
  if(!covars && !FE){ title2 <- paste0("no ", covar_lbl, " & no ", FE_lbl, "\n")}
  if(!covars && FE){ title2 <- paste0("no ", covar_lbl, " & ", FE_lbl, "\n")}
  if(covars && FE){ title2 <- paste0("w/ ", covar_lbl, " & ", FE_lbl, "\n")}
  
  labels_race = c("Nonhispanic White", "Nonhispanic black", 
                  "Hispanic", "Other race or multirace")
  
  scatter <- ggplot(mse_df, 
                    aes(x = ln_samplesize, y = ratio, 
                        color = as.factor(race))) + # , size = n for size of bubble
    geom_point(size=2, alpha=0.7) + 
    geom_abline(intercept = 1, slope = 0, linetype = "dashed", color = "black") +
    labs(title = paste0(title1, title2, subtitle), # title,  axis labels
         x = "log sample size", y = "MSE ratio (direct vs FH)") + 
    scale_color_manual(values = cbPalette, labels = labels_race) +
    labs(color = "Race")
  
  return(scatter) 
}


#Relative MSE - normalizing by dividing by mean of the true population
#12. Compute MSE from replicates estimates and population value
mse_r <- function(estimates_df, outcome_df, grouping_var, outcome_finalvar, FH=F, diff=F) {
  #Obtain MSE estimate for each domain      
  #input: estimates_df: dataframe with estimates from direct_estim or FH_estim fns
  #       outcome_df: dataframe with population value of parameter 
  #       outcome_finalvar: outcome variable to obtain mean
  #       grouping_var: grouping var, areas
  #       FH: T if computing the MSE for FH results (renaming)
  #       diff: T if computing the MSE for race difference (sort only 1 var)
  #return: df with MSE and domain IDs    
  
  if(FH){
    estimates_df <- estimates_df %>%
      select(starts_with("FH.estimates"), all_of(grouping_var)) %>%
      rename_with(~paste0("house_ten_", .), starts_with("FH.estimates"))
  }
  
  # Create a new dataframe with selected replications, and grouping vars 
  mse_df <- data.frame(estimates_df[, grouping_var],
                       estimates_df[, grep(paste0("^", outcome_finalvar), names(estimates_df))])
  
  # Order dfs
  if (diff) { #MSE for difference between races (only sort by ST)
    mse_df <- mse_df[order(mse_df$ST), ]
    outcome_df <- outcome_df[order(outcome_df$ST), ]
    names(outcome_df)[names(outcome_df) == "difference"] <- outcome_finalvar
  } else {
    mse_df <- mse_df[order(mse_df$race, mse_df$ST), ]
    outcome_df <- outcome_df[order(outcome_df$race, outcome_df$ST), ]
  }
  
  # Subtract population value of outcome from selected columns in new_dataframe
  mse_df[, -(1:length(grouping_var))] <- (mse_df[, -(1:length(grouping_var))] - 
                                            outcome_df[[outcome_finalvar]])^2
  
  # Average the difference and keep the final MSE only 
  mse_df <- mse_df %>%
    mutate(mse_r = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)/outcome_df$house_ten) %>%
    select(all_of(grouping_var), mse_r)
  
  if(FH){mse_df <- mse_df %>% rename(mse_r_fh = mse_r)} #rename to compare with direct estim
  
  return(mse_df) 
}

#12.0 export MSE (for differences estimates)
mse_r_direct_export <- function(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss){
  #Obtain csv file with summary of the relative MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_r_eq_1 <- mse_r_eq_1[order(mse_r_eq_1$ST), ]
  mse_r_eq_gss <- mse_r_eq_gss[order(mse_r_eq_gss$ST), ]
  mse_r_pr_1 <- mse_r_pr_1[order(mse_r_pr_1$ST), ]
  mse_r_pr_05 <- mse_r_pr_05[order(mse_r_pr_05$ST), ]
  mse_r_pr_gss <- mse_r_pr_gss[order(mse_r_pr_gss$ST), ]
  
  summ_mse_rs <- data.frame(mse_r_eq_1 = (mse_r_eq_1$mse_r),
                          mse_r_eq_gss = (mse_r_eq_gss$mse_r),
                          mse_r_pr_1 = (mse_r_pr_1$mse_r), 
                          mse_r_pr_05 = (mse_r_pr_05$mse_r), 
                          mse_r_pr_gss = (mse_r_pr_gss$mse_r))  
  
  overall_summ <- aggregate(cbind(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss) ~ NULL, 
                            data = summ_mse_rs, 
                            FUN = function(x) c(Median = median(x)))
  
  return(overall_summ)
}

mse_r_export <- function(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_r_eq_1 <- mse_r_eq_1[order(mse_r_eq_1$ST), ]
  mse_r_eq_gss <- mse_r_eq_gss[order(mse_r_eq_gss$ST), ]
  mse_r_pr_1 <- mse_r_pr_1[order(mse_r_pr_1$ST), ]
  mse_r_pr_05 <- mse_r_pr_05[order(mse_r_pr_05$ST), ]
  mse_r_pr_gss <- mse_r_pr_gss[order(mse_r_pr_gss$ST), ]
  
  summ_mse_rs <- data.frame(mse_r_eq_1 = (mse_r_eq_1$mse_r_fh),
                          mse_r_eq_gss = (mse_r_eq_gss$mse_r_fh),
                          mse_r_pr_1 = (mse_r_pr_1$mse_r_fh), 
                          mse_r_pr_05 = (mse_r_pr_05$mse_r_fh), 
                          mse_r_pr_gss = (mse_r_pr_gss$mse_r_fh))  
  
  overall_summ <- aggregate(cbind(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss) ~ NULL, 
                            data = summ_mse_rs, 
                            FUN = function(x) c(Median = median(x)))
  
  return(overall_summ)
}

#12.1 export MSE -- 
mse_r_direct_race <- function(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the mse_r improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with mse_r improvements by race and overall     
  
  mse_r_eq_1 <- mse_r_eq_1[order(mse_r_eq_1$race, mse_r_eq_1$ST), ]
  mse_r_eq_gss <- mse_r_eq_gss[order(mse_r_eq_gss$race, mse_r_eq_gss$ST), ]
  mse_r_pr_1 <- mse_r_pr_1[order(mse_r_pr_1$race, mse_r_pr_1$ST), ]
  mse_r_pr_05 <- mse_r_pr_05[order(mse_r_pr_05$race, mse_r_pr_05$ST), ]
  mse_r_pr_gss <- mse_r_pr_gss[order(mse_r_pr_gss$race, mse_r_pr_gss$ST), ]
  
  summ_mse_rs <- data.frame(race = mse_r_eq_1$race, 
                          mse_r_eq_1 = (mse_r_eq_1$mse_r),
                          mse_r_eq_gss = (mse_r_eq_gss$mse_r),
                          mse_r_pr_1 = (mse_r_pr_1$mse_r), 
                          mse_r_pr_05 = (mse_r_pr_05$mse_r), 
                          mse_r_pr_gss = (mse_r_pr_gss$mse_r))  
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss) ~ race, 
                         data = summ_mse_rs, 
                         FUN = function(x) c(Median = median(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss) ~ NULL, 
                            data = summ_mse_rs, 
                            FUN = function(x) c(Median = median(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}

mse_r_export_race <- function(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the mse_r improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_r_eq_1 <- mse_r_eq_1[order(mse_r_eq_1$race, mse_r_eq_1$ST), ]
  mse_r_eq_gss <- mse_r_eq_gss[order(mse_r_eq_gss$race, mse_r_eq_gss$ST), ]
  mse_r_pr_1 <- mse_r_pr_1[order(mse_r_pr_1$race, mse_r_pr_1$ST), ]
  mse_r_pr_05 <- mse_r_pr_05[order(mse_r_pr_05$race, mse_r_pr_05$ST), ]
  mse_r_pr_gss <- mse_r_pr_gss[order(mse_r_pr_gss$race, mse_r_pr_gss$ST), ]
  
  summ_mse_rs <- data.frame(race = mse_r_eq_1$race, 
                          mse_r_eq_1 = (mse_r_eq_1$mse_r_fh),
                          mse_r_eq_gss = (mse_r_eq_gss$mse_r_fh),
                          mse_r_pr_1 = (mse_r_pr_1$mse_r_fh), 
                          mse_r_pr_05 = (mse_r_pr_05$mse_r_fh), 
                          mse_r_pr_gss = (mse_r_pr_gss$mse_r_fh))  
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss) ~ race, 
                         data = summ_mse_rs, 
                         FUN = function(x) c(Median = median(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss) ~ NULL, 
                            data = summ_mse_rs, 
                            FUN = function(x) c(Median = median(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}

#12.2 MSE improvements by race and overall -- not yet done 
MSE_improve <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(state = mse_eq_1$ST, 
                          improve_eq_1 = (mse_eq_1$mse_fh - mse_eq_1$mse)/mse_eq_1$mse,
                          improve_eq_gss = (mse_eq_gss$mse_fh - mse_eq_gss$mse)/mse_eq_gss$mse,
                          improve_pr_1 = (mse_pr_1$mse_fh - mse_pr_1$mse)/mse_pr_1$mse,
                          improve_pr_05 = (mse_pr_05$mse_fh - mse_pr_05$mse)/mse_pr_05$mse,
                          improve_pr_gss = (mse_pr_gss$mse_fh - mse_pr_gss$mse)/mse_pr_gss$mse) 
  
  overall_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ NULL, 
                            data = summ_MSEs,
                            FUN = function(x) c(Median = median(x)))
  # FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(overall_summ)
}


#12.3 MSE improvements by race and overall -- not yet done 
MSE_improve_race <- function(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, racevar, statevar){
  #Obtain csv file with summary of the MSE improvement direct vs FH estimates for all designs       
  #input: all designs and race+state variable to sort
  #return: df with MSE improvements by race and overall     
  
  mse_eq_1 <- mse_eq_1[order(mse_eq_1$race, mse_eq_1$ST), ]
  mse_eq_gss <- mse_eq_gss[order(mse_eq_gss$race, mse_eq_gss$ST), ]
  mse_pr_1 <- mse_pr_1[order(mse_pr_1$race, mse_pr_1$ST), ]
  mse_pr_05 <- mse_pr_05[order(mse_pr_05$race, mse_pr_05$ST), ]
  mse_pr_gss <- mse_pr_gss[order(mse_pr_gss$race, mse_pr_gss$ST), ]
  
  summ_MSEs <- data.frame(race = mse_eq_1$race, 
                          improve_eq_1 = (mse_eq_1$mse_fh - mse_eq_1$mse)/mse_eq_1$mse,
                          improve_eq_gss = (mse_eq_gss$mse_fh - mse_eq_gss$mse)/mse_eq_gss$mse,
                          improve_pr_1 = (mse_pr_1$mse_fh - mse_pr_1$mse)/mse_pr_1$mse,
                          improve_pr_05 = (mse_pr_05$mse_fh - mse_pr_05$mse)/mse_pr_05$mse,
                          improve_pr_gss = (mse_pr_gss$mse_fh - mse_pr_gss$mse)/mse_pr_gss$mse) 
  
  # Use aggregate to obtain Min, Median, Mean, and Max by 'race'
  race_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ race, 
                         data = summ_MSEs, 
                         FUN = function(x) c(Median = median(x)))
  # FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  race_summ$race <- as.character(race_summ$race)
  race_summ <- race_summ %>%
    mutate(race = case_when(race == 1 ~ "Nonhispanic White", race == 2 ~ "Nonhispanic black",
                            race == 3 ~ "Hispanic", race == 4 ~ "Other race or multirace",
                            TRUE ~ as.character(race))) 
  
  overall_summ <- aggregate(cbind(improve_eq_1, improve_eq_gss, improve_pr_1, improve_pr_05, improve_pr_gss) ~ NULL, 
                            data = summ_MSEs, 
                            FUN = function(x) c(Min = min(x), Median = median(x), Mean = mean(x), Max = max(x)))
  overall_summ <- data.frame(race = "Overall", overall_summ)
  
  return(rbind(race_summ, overall_summ))
}


