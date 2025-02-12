#  
# title: "Run Simulation for ACS: direct and FH estimates"
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
here::i_am("02_scripts/03_models_simulation.R")
source(here("02_scripts","my_fns.R")) # load my functions 

#theme
theme <- theme_classic() #ggplot theme 
theme$plot.title$hjust <- 0.5

#seed
set.seed(123123) # for stochastic part 
#parallel 
plan(multisession)

#**************************************
#1. Read data, obtain direct estimates ---- 
#**************************************
load(here("04_results", "230818_500iter.RData"))
source(here("02_scripts","my_fns.R")) # load updated functions
options(survey.lonely.psu="adjust") # adjust for lonely PSU in survey design 

collapsed_stepw_df <- read.csv(here("03_data", "230926_additional_covariates.csv")) # Taylor's covariate. generated in other script
# collapsed_stepw_df <- readRDS(here("04_results", "stepwise_vars.rds"))
collapsed_stepw_df$ST <- sprintf("%02d", collapsed_stepw_df$ST) # make ST character and add leading 0s
collapsed_stepw_df$FINCP <- (collapsed_stepw_df$FINCP - mean(collapsed_stepw_df$FINCP))/sd(collapsed_stepw_df$FINCP) 
stepcovars <- colnames(collapsed_stepw_df)[!grepl("race|ST", colnames(collapsed_stepw_df), ignore.case = TRUE)] #no constant, race or ST 

outcome_finalvar <- "house_ten"
grouping_var <- c("race", "ST") 

# direct estimates (~100min for 500 reps)
direct_eq_1 <- direct_estim(df_eq_1, outcome_finalvar, grouping_var)
direct_eq_gss <- direct_estim(df_eq_gss, outcome_finalvar, grouping_var)
direct_pr_1 <- direct_estim(df_pr_1, outcome_finalvar, grouping_var)
direct_pr_05 <- direct_estim(df_pr_05, outcome_finalvar, grouping_var)
direct_pr_gss <- direct_estim(df_pr_gss, outcome_finalvar, grouping_var)

# replace 0.00 variances with p(1-p)/n assuming p=.5 por max variance
direct_eq_1 <- fix_zero_var(direct_eq_1, st_eq_1)
direct_eq_gss <- fix_zero_var(direct_eq_gss, st_eq_gss)
direct_pr_1 <- fix_zero_var(direct_pr_1, st_pr_1)
direct_pr_05 <- fix_zero_var(direct_pr_05, st_pr_05)
direct_pr_gss <- fix_zero_var(direct_pr_gss, st_pr_gss)


#**************************************
#2. FH & Plotting: Compute MSE for different cases ---- 
#**************************************

#labels for race 
labels_race = c("Nonhispanic White", "Nonhispanic black", 
                "Hispanic", "Other race or multirace")
# Colorblind friendly palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

nocovars <- c() 

#3.1 No covariates no FE ----
folder <-  "nocovars_noFE"

#FH estimates (~4min for 500 reps)
fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F)  
fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F)  
fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F)
fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F)  
fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F)

write.csv(fh_eq_1,  here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"))
write.csv(fh_eq_gss, here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"))
write.csv(fh_pr_05, here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"))
write.csv(fh_pr_1, here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"))
write.csv(fh_pr_gss, here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"))


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
##directs
write.csv(MSE_direct_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_MSE_direct.csv"))

## models 
write.csv(MSE_export_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_MSE_improve.csv"))

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()


#3.2 No covariates & FE ---- 
folder <- "nocovars_FE"

fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T)  
fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T)  
fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T)
fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T)  
fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T)  

write.csv(fh_eq_1,  here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"))
write.csv(fh_eq_gss, here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"))
write.csv(fh_pr_05, here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"))
write.csv(fh_pr_1, here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"))
write.csv(fh_pr_gss, here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"))

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
          here("04_results", "plots", folder, "summary_MSE.csv"))
#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", folder, "summary_MSE_improve.csv"))

#MSE plots 
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=F, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=F, FE=T),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()



#3.3 Simple covariates & no FE ---- 
folder <- "covars_simple_noFE"

fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F)  
fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F)  
fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F)
fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F)  
fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F)  

write.csv(fh_eq_1,  here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"))
write.csv(fh_eq_gss, here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"))
write.csv(fh_pr_05, here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"))
write.csv(fh_pr_1, here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"))
write.csv(fh_pr_gss, here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"))

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
          here("04_results", "plots", folder, "summary_MSE.csv"))
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_MSE_improve.csv"))

#MSE plots 
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

#3.4 Simple covariates & FE ---- 
folder <- "covars_simple_FE"

fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T)  
fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T)  
fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T)
fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T)  
fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T)

write.csv(fh_eq_1,  here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"))
write.csv(fh_eq_gss, here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"))
write.csv(fh_pr_05, here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"))
write.csv(fh_pr_1, here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"))
write.csv(fh_pr_gss, here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"))


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
          here("04_results", "plots", folder, "summary_MSE.csv"))
#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_MSE_improve.csv"))

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()


#3.5 Stepwise covariates & no FE ----
folder <- "covars_step_noFE"

fh_eq_1 <- FH_estim(direct_eq_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F)  
fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F)
fh_pr_1 <- FH_estim(direct_pr_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F)
fh_pr_05 <- FH_estim(direct_pr_05, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F)  
fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F)

write.csv(fh_eq_1,  here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"))
write.csv(fh_eq_gss, here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"))
write.csv(fh_pr_05, here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"))
write.csv(fh_pr_1, here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"))
write.csv(fh_pr_gss, here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"))


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
          here("04_results", "plots", folder, "summary_MSE.csv"))
#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_MSE_improve.csv"))

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()


#3.6 Stepwise covariates & FE ---- 
folder <- "covars_step_FE"

fh_eq_1 <- FH_estim(direct_eq_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T)  
fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T)
fh_pr_1 <- FH_estim(direct_pr_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T)
fh_pr_05 <- FH_estim(direct_pr_05, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T)  
fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T)

write.csv(fh_eq_1,  here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"))
write.csv(fh_eq_gss, here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"))
write.csv(fh_pr_05, here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"))
write.csv(fh_pr_1, here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"))
write.csv(fh_pr_gss, here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"))


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
          here("04_results", "plots", folder, "summary_MSE.csv"))
#MSE improvements 
write.csv(MSE_improve_race(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_MSE_improve.csv"))

#MSE plots 
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eq1.png"), 
       plot = plotMSE(mse_eq_1, st_eq_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_eqgss.png"), 
       plot = plotMSE(mse_eq_gss, st_eq_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr1.png"), 
       plot = plotMSE(mse_pr_1, st_pr_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_MSE_pr05.png"), 
       plot = plotMSE(mse_pr_05, st_pr_05, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_MSE_prgss.png"), 
       plot = plotMSE(mse_pr_gss, st_pr_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)


# Convergence plots and save 
xlabel <- "Sample Size"
ylabel <- "(Direct / FH) * 100"
note <- "Note: FH = Fay-Harriot model. X-axis in logarithmic scale. Compiled by authors."

## Design: eq 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eq1.jpg"), width = 500, height = 500)
converplot(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eq1.jpg"), width = 500, height = 500)
converplot_std(st_eq_1, fh_eq_1, direct_eq_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: eq GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_eqgss.jpg"), width = 500, height = 500)
converplot(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Equal sample size for all races (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_eqgss.jpg"), width = 500, height = 500)
converplot_std(st_eq_gss, fh_eq_gss, direct_eq_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 1% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr1.jpg"), width = 500, height = 500)
converplot(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (1% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr1.jpg"), width = 500, height = 500)
converplot_std(st_pr_1, fh_pr_1, direct_pr_1, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr 0.5% 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (0.5% ACS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_05, fh_pr_05, direct_pr_05, grouping_var, title, xlabel, ylabel, note)
dev.off()

## Design: pr GSS 
title <- "House Tenancy by Race. Ratio of Estimates \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_estim_pr05.jpg"), width = 500, height = 500)
converplot(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()

title <- "House Tenancy by Race. Ratio Std Errors \n Proportional sample size for each race (GSS size)"
jpeg(here("04_results", "plots", folder, "conver_std_pr05.jpg"), width = 500, height = 500)
converplot_std(st_pr_gss, fh_pr_gss, direct_pr_gss, grouping_var, title, xlabel, ylabel, note)
dev.off()


##### END ####
#save 
save.image(here("04_results","direct_FH_500reps.RData"))
beep(2)