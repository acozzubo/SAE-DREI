#  
# title: "Run Simulation for ACS: direct and FH estimates for difference between races"
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
here::i_am("02_scripts/04_models_difference.R")
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

popdiff_df <- outcome_df %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(difference = diff(range(house_ten)))

outcome_finalvar <- "house_ten"
grouping_var <- c("race", "ST") 
grouping_var_diff <- "ST" #dont consider race here as we compute race differences

#**************************************
#2. Compute difference of direct estimates for each design ---- 
#**************************************
folder <-  "direct"

diff_dir_eq_1 <- direct_eq_1 %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("house_ten"), ~ diff(range(.))))

diff_dir_eq_gss <- direct_eq_gss %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("house_ten"), ~ diff(range(.))))

diff_dir_pr_1 <- direct_pr_1 %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("house_ten"), ~ diff(range(.))))

diff_dir_pr_05 <- direct_pr_05 %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("house_ten"), ~ diff(range(.))))

diff_dir_pr_gss <- direct_pr_gss %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("house_ten"), ~ diff(range(.))))



write.csv(diff_dir_eq_1,  here("04_results", "estimates", folder, "Estimates_diff_dir_eq_1.csv"))
write.csv(diff_dir_eq_gss, here("04_results", "estimates", folder,"Estimates_diff_dir_eq_gss.csv"))
write.csv(diff_dir_pr_05, here("04_results", "estimates", folder, "Estimates_diff_dir_pr_05.csv"))
write.csv(diff_dir_pr_1, here("04_results", "estimates", folder, "Estimates_diff_dir_pr_1.csv"))
write.csv(diff_dir_pr_gss, here("04_results", "estimates", folder, "Estimates_diff_dir_pr_gss.csv"))

collapsed_df$ST <- as.numeric(collapsed_df$ST)

#**************************************
#3. Compute difference of FH estimates for each design and model ---- 
#**************************************

nocovars <- c() 

#3.1 No covariates no FE ----
folder <-  "nocovars_noFE"

#Difference FH estimates (~4min for 500 reps)
diff_fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F)  %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

write.csv(diff_fh_eq_1,  here("04_results", "estimates", folder, "Estimates_diff_fh_eq_1.csv"))
write.csv(diff_fh_eq_gss, here("04_results", "estimates", folder,"Estimates_diff_fh_eq_gss.csv"))
write.csv(diff_fh_pr_05, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_05.csv"))
write.csv(diff_fh_pr_1, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_1.csv"))
write.csv(diff_fh_pr_gss, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_gss.csv"))


#MSE computations 
mse_eq_1 <- inner_join(MSE(diff_dir_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), #direct 
                       MSE(diff_fh_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), #FH
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(diff_dir_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(diff_dir_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                       MSE(diff_fh_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(diff_dir_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                        MSE(diff_fh_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                        by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(diff_dir_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE

#MSE export 
##directs
write.csv(MSE_direct_export(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", "summary_MSE_direct.csv"))

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eq1.png"), 
       plot = plot_diffMSE(mse_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eqgss.png"), 
       plot = plot_diffMSE(mse_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr1.png"), 
       plot = plot_diffMSE(mse_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr05.png"), 
       plot = plot_diffMSE(mse_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_prgss.png"), 
       plot = plot_diffMSE(mse_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

#MSE export 
write.csv(MSE_export(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE_improve.csv"))


#3.2 No covariates & FE ---- 
folder <- "nocovars_FE"

#Difference FH estimates (~4min for 500 reps)
diff_fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, nocovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

write.csv(diff_fh_eq_1,  here("04_results", "estimates", folder, "Estimates_diff_fh_eq_1.csv"))
write.csv(diff_fh_eq_gss, here("04_results", "estimates", folder,"Estimates_diff_fh_eq_gss.csv"))
write.csv(diff_fh_pr_05, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_05.csv"))
write.csv(diff_fh_pr_1, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_1.csv"))
write.csv(diff_fh_pr_gss, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_gss.csv"))

#MSE computations 
mse_eq_1 <- inner_join(MSE(diff_dir_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), #direct 
                       MSE(diff_fh_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), #FH
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(diff_dir_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(diff_dir_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                       MSE(diff_fh_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(diff_dir_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                        MSE(diff_fh_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                        by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(diff_dir_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eq1.png"), 
       plot = plot_diffMSE(mse_eq_1, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eqgss.png"), 
       plot = plot_diffMSE(mse_eq_gss, subtitle, covars=F, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr1.png"), 
       plot = plot_diffMSE(mse_pr_1, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr05.png"), 
       plot = plot_diffMSE(mse_pr_05, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_prgss.png"), 
       plot = plot_diffMSE(mse_pr_gss, subtitle, covars=F, FE=T),  
       width = 8, height = 4, dpi = 300)

#MSE export 
write.csv(MSE_export(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE_improve.csv"))


#3.3 Simple covariates & no FE ---- 
folder <- "covars_simple_noFE"

#Difference FH estimates (~4min for 500 reps)
diff_fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

write.csv(diff_fh_eq_1,  here("04_results", "estimates", folder, "Estimates_diff_fh_eq_1.csv"))
write.csv(diff_fh_eq_gss, here("04_results", "estimates", folder,"Estimates_diff_fh_eq_gss.csv"))
write.csv(diff_fh_pr_05, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_05.csv"))
write.csv(diff_fh_pr_1, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_1.csv"))
write.csv(diff_fh_pr_gss, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_gss.csv"))

#MSE computations 
mse_eq_1 <- inner_join(MSE(diff_dir_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), #direct 
                       MSE(diff_fh_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), #FH
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(diff_dir_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(diff_dir_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                       MSE(diff_fh_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(diff_dir_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                        MSE(diff_fh_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                        by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(diff_dir_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE


#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eq1.png"), 
       plot = plot_diffMSE(mse_eq_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eqgss.png"), 
       plot = plot_diffMSE(mse_eq_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr1.png"), 
       plot = plot_diffMSE(mse_pr_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr05.png"), 
       plot = plot_diffMSE(mse_pr_05, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_prgss.png"), 
       plot = plot_diffMSE(mse_pr_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

#MSE export 
write.csv(MSE_export(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE_improve.csv"))


#3.4 Simple covariates & FE ---- 
folder <- "covars_simple_FE"

#Difference FH estimates (~4min for 500 reps)
diff_fh_eq_1 <- FH_estim(direct_eq_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T)  %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_1 <- FH_estim(direct_pr_1, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_05 <- FH_estim(direct_pr_05, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_df, covars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

write.csv(diff_fh_eq_1,  here("04_results", "estimates", folder, "Estimates_diff_fh_eq_1.csv"))
write.csv(diff_fh_eq_gss, here("04_results", "estimates", folder,"Estimates_diff_fh_eq_gss.csv"))
write.csv(diff_fh_pr_05, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_05.csv"))
write.csv(diff_fh_pr_1, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_1.csv"))
write.csv(diff_fh_pr_gss, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_gss.csv"))

#MSE computations 
mse_eq_1 <- inner_join(MSE(diff_dir_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), #direct 
                       MSE(diff_fh_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), #FH
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(diff_dir_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(diff_dir_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                       MSE(diff_fh_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(diff_dir_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                        MSE(diff_fh_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                        by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(diff_dir_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eq1.png"), 
       plot = plot_diffMSE(mse_eq_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eqgss.png"), 
       plot = plot_diffMSE(mse_eq_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr1.png"), 
       plot = plot_diffMSE(mse_pr_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr05.png"), 
       plot = plot_diffMSE(mse_pr_05, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_prgss.png"), 
       plot = plot_diffMSE(mse_pr_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

#MSE export 
write.csv(MSE_export(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE_improve.csv"))


#3.5 Stepwise covariates & no FE ----
collapsed_stepw_df$ST <- as.numeric(collapsed_stepw_df$ST)

folder <- "covars_step_noFE"

#Difference FH estimates (~4min for 500 reps)
diff_fh_eq_1 <- FH_estim(direct_eq_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_1 <- FH_estim(direct_pr_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_05 <- FH_estim(direct_pr_05, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=F) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

write.csv(diff_fh_eq_1,  here("04_results", "estimates", folder, "Estimates_diff_fh_eq_1.csv"))
write.csv(diff_fh_eq_gss, here("04_results", "estimates", folder,"Estimates_diff_fh_eq_gss.csv"))
write.csv(diff_fh_pr_05, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_05.csv"))
write.csv(diff_fh_pr_1, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_1.csv"))
write.csv(diff_fh_pr_gss, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_gss.csv"))

#MSE computations 
mse_eq_1 <- inner_join(MSE(diff_dir_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), #direct 
                       MSE(diff_fh_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), #FH
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(diff_dir_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(diff_dir_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                       MSE(diff_fh_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(diff_dir_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                        MSE(diff_fh_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                        by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(diff_dir_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eq1.png"), 
       plot = plot_diffMSE(mse_eq_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eqgss.png"), 
       plot = plot_diffMSE(mse_eq_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr1.png"), 
       plot = plot_diffMSE(mse_pr_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr05.png"), 
       plot = plot_diffMSE(mse_pr_05, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_prgss.png"), 
       plot = plot_diffMSE(mse_pr_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

#MSE export 
write.csv(MSE_export(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE.csv"))
#MSE improvements 
write.csv(MSE_improve(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE_improve.csv"))


#3.6 Stepwise covariates & FE
folder <- "covars_step_FE"

#Difference FH estimates (~4min for 500 reps)
diff_fh_eq_1 <- FH_estim(direct_eq_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_eq_gss <- FH_estim(direct_eq_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_1 <- FH_estim(direct_pr_1, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_05 <- FH_estim(direct_pr_05, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

diff_fh_pr_gss <- FH_estim(direct_pr_gss, collapsed_stepw_df, stepcovars, grouping_var, outcome_finalvar, FE=T) %>%
  filter(race %in% c(1, 2)) %>%
  group_by(ST) %>%
  summarise(across(starts_with("FH.estimates."), ~ diff(range(.))))

write.csv(diff_fh_eq_1,  here("04_results", "estimates", folder, "Estimates_diff_fh_eq_1.csv"))
write.csv(diff_fh_eq_gss, here("04_results", "estimates", folder,"Estimates_diff_fh_eq_gss.csv"))
write.csv(diff_fh_pr_05, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_05.csv"))
write.csv(diff_fh_pr_1, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_1.csv"))
write.csv(diff_fh_pr_gss, here("04_results", "estimates", folder, "Estimates_diff_fh_pr_gss.csv"))

#MSE computations 
mse_eq_1 <- inner_join(MSE(diff_dir_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), #direct 
                       MSE(diff_fh_eq_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), #FH
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_eq_gss <- inner_join(MSE(diff_dir_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_eq_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_1 <- inner_join(MSE(diff_dir_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                       MSE(diff_fh_pr_1, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                       by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_05 <- inner_join(MSE(diff_dir_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                        MSE(diff_fh_pr_05, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                        by = grouping_var_diff, na_matches = "never") #join both MSE 

mse_pr_gss <- inner_join(MSE(diff_dir_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, diff=T), 
                         MSE(diff_fh_pr_gss, popdiff_df, grouping_var_diff, outcome_finalvar, FH=T, diff=T), 
                         by = grouping_var_diff, na_matches = "never") #join both MSE

#MSE plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eq1.png"), 
       plot = plot_diffMSE(mse_eq_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_eqgss.png"), 
       plot = plot_diffMSE(mse_eq_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr1.png"), 
       plot = plot_diffMSE(mse_pr_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots", "_differences",  folder, "scatter_MSE_pr05.png"), 
       plot = plot_diffMSE(mse_pr_05, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", "_differences", folder, "scatter_MSE_prgss.png"), 
       plot = plot_diffMSE(mse_pr_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

#MSE export 
write.csv(MSE_export(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE.csv"))

#MSE improvements 
write.csv(MSE_improve(mse_eq_1, mse_eq_gss, mse_pr_1, mse_pr_05, mse_pr_gss), 
          here("04_results", "plots", "_differences", folder, "summary_MSE_improve.csv"))

##### END ####
beep(8)