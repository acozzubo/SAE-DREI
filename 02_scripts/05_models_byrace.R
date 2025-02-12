#  
# title: "Run Simulation for ACS: FH estimates, separate models by race"
# subtitle: "DREI2023_SAE"
# author: "[Angelo Cozzubo] \n(https://sites.google.com/pucp.pe/acozz)"
# date creation: "28/June/2023"
#  

#**************************************
#0. Preliminaries ---- 
#**************************************
rm(list=ls())
#install.packages("pacman")
pacman::p_load(here, haven, tidyverse, labelled, dplyr, stringr, tableone, tidysmd, tictoc,
               beepr, ggplot2, writexl, survey, sampling, sae, janitor, epiDisplay, future.apply)
##fix priorities with dplyr
select <- dplyr::select
filter <- dplyr::filter

#routes 
here::i_am("02_scripts/05_models_byrace.R")
source(here("02_scripts","my_fns.R")) # load my functions 

#theme
theme <- theme_classic() #ggplot theme 
theme$plot.title$hjust <- 0.5

#seed
set.seed(123123) # for stochastic part 
#parallel 
plan(multisession)

#**************************************
#1. Read data, obtain population differences ---- 
#**************************************
load(here("04_results","direct_FH_500reps.RData")) #load direct estimates

outcome_finalvar <- "house_ten"
grouping_var <- c("race", "ST") 

#****************************************************************************
#2 .Compute FH estimates for each race separately, all desigsn and model w/o FE ---- 
#****************************************************************************

#2.1 No covariates no FE ----
folder <-  "nocovars_noFE"

# Create empty data frames to store the results
fh_eq_1 <- data.frame()  
fh_eq_gss <- data.frame() 
fh_pr_1 <- data.frame()
fh_pr_05 <- data.frame()
fh_pr_gss <- data.frame()

# Loop through each race and do FH, append rows 
for (i in c(1, 2, 3, 4)) {
  fh_eq_1 <- rbind(FH_estim(direct_eq_1[direct_eq_1$race == i,], 
                           collapsed_df[collapsed_df$race == i,],
                           nocovars, grouping_var, outcome_finalvar, FE=F),
                   fh_eq_1)

  fh_eq_gss <- rbind(FH_estim(direct_eq_gss[direct_eq_gss$race == i,], 
                            collapsed_df[collapsed_df$race == i,],
                            nocovars, grouping_var, outcome_finalvar, FE=F),
                     fh_eq_gss)
  
  fh_pr_1 <- rbind(FH_estim(direct_pr_1[direct_pr_1$race == i,], 
                            collapsed_df[collapsed_df$race == i,],
                            nocovars, grouping_var, outcome_finalvar, FE=F),
                   fh_pr_1)
  
  fh_pr_05 <- rbind(FH_estim(direct_pr_05[direct_pr_05$race == i,], 
                            collapsed_df[collapsed_df$race == i,],
                            nocovars, grouping_var, outcome_finalvar, FE=F),
                    fh_pr_05)
  
  fh_pr_gss <- rbind(FH_estim(direct_pr_gss[direct_pr_gss$race == i,], 
                            collapsed_df[collapsed_df$race == i,],
                            nocovars, grouping_var, outcome_finalvar, FE=F),
                     fh_pr_gss)
}

#MSE computations 
mse_eq_1 <- inner_join(MSE(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                       MSE(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                       by = grouping_var, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                         MSE(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                         by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                       MSE(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                       by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                        MSE(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                        by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                         MSE(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                         by = grouping_var, na_matches = "never") #join both MSE

#MSE export 
write.csv(MSE_export_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", "_models_by_race", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", "_models_by_race", folder, "summary_MSE_improve.csv"))

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  "_models_by_race", folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  "_models_by_race", folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

#2.2 Simple covariates & no FE ---- 
folder <- "covars_simple_noFE"

# Create empty data frames to store the results
fh_eq_1 <- data.frame()  
fh_eq_gss <- data.frame() 
fh_pr_1 <- data.frame()
fh_pr_05 <- data.frame()
fh_pr_gss <- data.frame()


# Loop through each race and do FH, append rows 
for (i in c(1, 2, 3, 4)) {
  fh_eq_1 <- rbind(FH_estim(direct_eq_1[direct_eq_1$race == i,],
                            collapsed_df[collapsed_df$race == i,],
                            covars, grouping_var, outcome_finalvar, FE=F),
                   fh_eq_1)
  
  fh_eq_gss <- rbind(FH_estim(direct_eq_gss[direct_eq_gss$race == i,],
                              collapsed_df[collapsed_df$race == i,],
                              c("AGEP", "MAR", "LAPTOP"), grouping_var, outcome_finalvar, FE=F), ## problem with covar "BDSP" and race=1
                     fh_eq_gss)                                                                          ## standardization does not solve
  
  fh_pr_1 <- rbind(FH_estim(direct_pr_1[direct_pr_1$race == i,],
                            collapsed_df[collapsed_df$race == i,],
                            covars, grouping_var, outcome_finalvar, FE=F),
                   fh_pr_1)

  fh_pr_05 <- rbind(FH_estim(direct_pr_05[direct_pr_05$race == i,],
                             collapsed_df[collapsed_df$race == i,],
                             covars, grouping_var, outcome_finalvar, FE=F),
                    fh_pr_05)

  fh_pr_gss <- rbind(FH_estim(direct_pr_gss[direct_pr_gss$race == i,],
                              collapsed_df[collapsed_df$race == i,],
                              covars, grouping_var, outcome_finalvar, FE=F),
                     fh_pr_gss)
}

#MSE computations 
mse_eq_1 <- inner_join(MSE(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                       MSE(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                       by = grouping_var, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar),
                         MSE(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T),
                         by = grouping_var, na_matches = "never") #join both MSE

mse_pr_1 <- inner_join(MSE(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                       MSE(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                       by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                        MSE(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                        by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                         MSE(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                         by = grouping_var, na_matches = "never") #join both MSE

#MSE export 
write.csv(MSE_export_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", "_models_by_race", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", "_models_by_race", folder, "summary_MSE_improve.csv"))

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  "_models_by_race", folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  "_models_by_race", folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()


#2.3 Stepwise covariates & no FE ----
folder <- "covars_step_noFE"

# Create empty data frames to store the results
fh_eq_1 <- data.frame()  
fh_eq_gss <- data.frame() 
fh_pr_1 <- data.frame()
fh_pr_05 <- data.frame()
fh_pr_gss <- data.frame()

# Loop through each race and do FH, append rows 
for (i in c(1, 2, 3, 4)) {
  fh_eq_1 <- rbind(FH_estim(direct_eq_1[direct_eq_1$race == i,], 
                            collapsed_df[collapsed_df$race == i,],
                            stepcovars, grouping_var, outcome_finalvar, FE=F),
                   fh_eq_1)
  
  fh_eq_gss <- rbind(FH_estim(direct_eq_gss[direct_eq_gss$race == i,], 
                              collapsed_df[collapsed_df$race == i,],
                              stepcovars, grouping_var, outcome_finalvar, FE=F),
                     fh_eq_gss)
  
  fh_pr_1 <- rbind(FH_estim(direct_pr_1[direct_pr_1$race == i,], 
                            collapsed_df[collapsed_df$race == i,],
                            stepcovars, grouping_var, outcome_finalvar, FE=F),
                   fh_pr_1)
  
  fh_pr_05 <- rbind(FH_estim(direct_pr_05[direct_pr_05$race == i,], 
                             collapsed_df[collapsed_df$race == i,],
                             stepcovars, grouping_var, outcome_finalvar, FE=F),
                    fh_pr_05)
  
  fh_pr_gss <- rbind(FH_estim(direct_pr_gss[direct_pr_gss$race == i,], 
                              collapsed_df[collapsed_df$race == i,],
                              stepcovars, grouping_var, outcome_finalvar, FE=F),
                     fh_pr_gss)
}

#MSE computations 
mse_eq_1 <- inner_join(MSE(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                       MSE(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                       by = grouping_var, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                         MSE(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                         by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                       MSE(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                       by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                        MSE(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                        by = grouping_var, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                         MSE(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                         by = grouping_var, na_matches = "never") #join both MSE

#MSE export 
write.csv(MSE_export_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", "_models_by_race", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", "_models_by_race", folder, "summary_MSE_improve.csv"))

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  "_models_by_race", folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  "_models_by_race", folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_models_by_race", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", "_models_by_race", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

##### END ####
beep(8)