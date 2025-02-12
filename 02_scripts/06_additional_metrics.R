#  
# title: "Additional Metrics"
# subtitle: "DREI2023_SAE"
# author: "[Carissa Villanueva]"
# date creation: "06/March/2024"
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

folder <- "direct"
diff_dir_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_diff_dir_eq_1.csv"), header = T)
diff_dir_eq_1 = diff_dir_eq_1[,-1]
diff_dir_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_diff_dir_eq_gss.csv"), header = T)
diff_dir_eq_gss = diff_dir_eq_gss[,-1]
diff_dir_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_diff_dir_pr_05.csv"), header = T)
diff_dir_pr_05 = diff_dir_pr_05[,-1]
diff_dir_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_diff_dir_pr_1.csv"), header = T)
diff_dir_pr_1 = diff_dir_pr_1[,-1]
diff_dir_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_diff_dir_pr_gss.csv"), header = T)
diff_dir_pr_gss = diff_dir_pr_gss[,-1]

#numbers converted in csv fixing for direct ests
direct_eq_1$ST <-as.numeric(direct_eq_1$ST)

direct_eq_gss$ST <-as.numeric(direct_eq_gss$ST)

direct_pr_1$ST <-as.numeric(direct_pr_1$ST)

direct_pr_05$ST <-as.numeric(direct_pr_05$ST)

direct_pr_gss$ST <-as.numeric(direct_pr_gss$ST)



#****************************************************************************
#2 .Compute Absolute Bias for FH and Direct Estimates ---- 
#****************************************************************************

#2.1 No covariates no FE ----
folder <-  "nocovars_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias


#AbsBias export 
##directs
write.csv(AbsBias_direct_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_AbsBias_direct.csv"))

## models 
write.csv(AbsBias_export_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias.csv"))

#AbsBias improvements 
write.csv(AbsBias_improve_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias_improve.csv"))


#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plot_diffAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plot_diffAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plot_diffAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plot_diffAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plot_diffAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

#2.2  No covariates & FE ----

folder <-  "nocovars_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias

#AbsBias export 
##directs
write.csv(AbsBias_direct_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_AbsBias_direct.csv"))

## models 
write.csv(AbsBias_export_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias.csv"))

#AbsBias improvements 
write.csv(AbsBias_improve_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias_improve.csv"))


#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=T),  
       width = 8, height = 4, dpi = 300)

#2.3 Simple covariates & no FE ---- 

folder <-  "covars_simple_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias

#AbsBias export 
##directs
write.csv(AbsBias_direct_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_AbsBias_direct.csv"))

## models 
write.csv(AbsBias_export_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias.csv"))

#AbsBias improvements 
write.csv(AbsBias_improve_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias_improve.csv"))


#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

#Additional Plots
#pr05
AbsBias_pr_05$siz = st_pr_05
ggplot(data = AbsBias_pr_05, aes(y=Abs_Bias_FH, x=siz)) +geom_point() + ylim(0,.42) + scale_x_log10() +
  geom_jitter() + labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and no Fixed Effects", 
                       subtitle = "Where the race sample is proportional to size in ACS pop (.5% size of ACS)", 
                       x = "Sample Size - scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_pr05.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_05, aes(y=Abs_Bias, x=(siz))) +geom_point() +
  geom_jitter() + ylim(0,.42) + scale_x_log10() +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is proportional to size in ACS pop (.5% size of ACS)",
       x = "Sample Size- scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_pr05.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)
#pr1
AbsBias_pr_1$size = st_pr_1
ggplot(data = AbsBias_pr_1, aes(y=Abs_Bias_FH, x=size)) +geom_point() + geom_jitter() +
  ylim(0,.45)  + scale_x_log10() +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is proportional to size in ACS pop (1% size of ACS)",
       x = "Sample Size - scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_pr1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_1, aes(y=Abs_Bias, x=size)) +geom_point() + geom_jitter() +
  ylim(0,.45)  + scale_x_log10() +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is proportional to size in ACS pop (1% size of ACS)",
       x = "Sample Size - scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_pr1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)
#prgss
AbsBias_pr_gss$size = st_pr_gss
ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias_FH, x=size)) +geom_point() + geom_jitter() +
  ylim(0,.45)  + scale_x_log10() +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is proportional to size in GSS pop",
       x = "Sample Size - scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_prgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias, x=size)) +geom_point() + geom_jitter() +
  ylim(0,.45)  + scale_x_log10() +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is proportional to size in GSS population",
       x = "Sample Size - scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_prgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

#prgss - race
ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias_FH, x=size)) +geom_point() +
  geom_jitter() + facet_wrap(vars(race)) + ylim(0,.42)  + scale_x_log10() +
  labs(title ="Bias of Fay-Harriot Estimates by Race", 
       subtitle = "Model with Simple Covariates and no Fixed Effects", 
       caption ="With race sample proportional to size in GSS population",
       x = "Sample Size - scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_prgss_Race.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias, x=size)) +geom_point() +
  geom_jitter() + facet_wrap(vars(race))+ ylim(0,.42)  + scale_x_log10() +
  labs(title ="Bias of Direct Estimates by Race", 
       subtitle = "Model with Simple Covariates and no Fixed Effects", 
       caption ="With race sample proportional to size in GSS population",
       x = "Sample Size - scaled to Log10",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_prgss_race.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)
#eqgss
AbsBias_eq_gss$size = st_eq_gss
ggplot(data = AbsBias_eq_gss, aes(x=Abs_Bias_FH)) +geom_histogram() + xlim(0.02,.12) + ylim(0,50) +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is equal from the GSS population (N=20)",
       x = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_eqgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_eq_gss, aes(x=Abs_Bias)) +geom_histogram() + xlim(0,.13) + 
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is equal from the GSS population (N=20)",
       x = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_eqgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

#eq1
AbsBias_eq_1$size = st_eq_1
ggplot(data = AbsBias_eq_1, aes(x=Abs_Bias_FH)) +geom_histogram() + xlim(-0.01,.10) + ylim(0,101) +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is equal from the ACS population (N=152)",
       x = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_eq1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_eq_1, aes(x=Abs_Bias)) +geom_histogram()   + xlim(-0.01,.10)  + ylim(0,101)+
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and no Fixed Effects", 
       subtitle = "Where the race sample is equal from the ACS population (N =152)",
       x  = "Absolute Bias of the Estimate" )

ggsave(here("04_results", "plots", folder, "size_AbsBias_eq1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)


#2.4 Simple covariates & FE ----

folder <-  "covars_simple_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias

#AbsBias export 
##directs
write.csv(AbsBias_direct_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_AbsBias_direct.csv"))

## models 
write.csv(AbsBias_export_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias.csv"))

#AbsBias improvements 
write.csv(AbsBias_improve_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias_improve.csv"))


#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=T, FE=), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars= T, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

#pr05
AbsBias_pr_05$siz = st_pr_05
ggplot(data = AbsBias_pr_05, aes(y=Abs_Bias_FH, x=log(siz))) +geom_point() + 
  geom_jitter() + labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and Fixed Effects for Race", 
                       subtitle = "Where the race sample is proportional to size in ACS pop (.5% size of ACS)", x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_pr05.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_05, aes(y=Abs_Bias, x=log(siz))) +geom_point() +
  geom_jitter() +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is proportional to size in ACS pop (.5% size of ACS)",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_pr05.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)
#pr1
AbsBias_pr_1$size = st_pr_1
ggplot(data = AbsBias_pr_1, aes(y=Abs_Bias_FH, x=log(size))) +geom_point() + geom_jitter() +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is proportional to size in ACS pop (1% size of ACS)",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_pr1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_1, aes(y=Abs_Bias, x=log(size))) +geom_point() + geom_jitter() +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is proportional to size in ACS pop (1% size of ACS)",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_pr1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)
#prgss
AbsBias_pr_gss$size = st_pr_gss
ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias_FH, x=log(size))) +geom_point() + geom_jitter() +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is proportional to size in GSS pop",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_prgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias, x=log(size))) +geom_point() + geom_jitter() +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is proportional to size in GSS population",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_prgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

#prgss - race   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias_FH, x=size)) +geom_point() +
  geom_jitter() + facet_wrap(vars(race)) + scale_x_log10() + ylim(0,.3) +
  labs(title ="Bias of Fay-Harriot Estimates by Race", 
       subtitle = "Model with Simple Covariates and Fixed Effects for Race", 
       caption ="With race sample proportional to size in GSS population",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_prgss_Race.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias, x=size)) +geom_point() +
  geom_jitter() + facet_wrap(vars(race))+  scale_x_log10() + ylim(0,.3) +
  labs(title ="Bias of Direct Estimates by Race", 
       subtitle = "Model with Simple Covariates and Fixed Effects for Race", 
       caption ="With race sample proportional to size in GSS population",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_prgss_race.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)
#eqgss
AbsBias_eq_gss$size = st_eq_gss
ggplot(data = AbsBias_eq_gss, aes(x=Abs_Bias_FH)) +geom_histogram() + xlim(0.0,.15) +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the GSS population (N=20)",
       x  = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_eqgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_eq_gss, aes(x=Abs_Bias)) +geom_histogram() +  xlim(0.0,.15) +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the GSS population (N=20)",
       x =  "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_eqgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

#eq1
AbsBias_eq_1$size = st_eq_1
ggplot(data = AbsBias_eq_1, aes(x=Abs_Bias_FH)) +geom_histogram() +  xlim(-0.01,.15) +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the ACS population (N=152)",
       x  = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_eq1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_eq_1, aes(x=Abs_Bias)) +geom_histogram() + xlim(-0.01,.15) +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the ACS population (N=152)",
       x = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_eq1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

#2.5 Stepwise covariates & no FE ----

folder <-  "covars_step_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias

#AbsBias export 
##directs
write.csv(AbsBias_direct_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_AbsBias_direct.csv"))

## models 
write.csv(AbsBias_export_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias.csv"))

#AbsBias improvements 
write.csv(AbsBias_improve_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias_improve.csv"))

simp_improve_pr1 = AbsBias_pr_05
simp_improve_pr1$biasdiff =  simp_improve_pr1$Abs_Bias - simp_improve_pr1$Abs_Bias_FH 

ggplot(data = simp_improve_pr1, aes(x=biasdiff)) + geom_histogram(bins = 35) +  facet_wrap(vars(race)) + 
  labs(title ="Reducation in Bias (Direct vs FH) for Model with Stepwise Covariates and no FE", 
       subtitle = "Where the race sample is proportional to size in ACS pop (.05% size of ACS)",
       x = "Absolute Bias of the Estimate") 


#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=T, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=T, FE=F),  
       width = 8, height = 4, dpi = 300)


#2.6 Stepwise covariates & FE ---- 
folder <-  "covars_step_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias



#AbsBias export 
##directs
write.csv(AbsBias_direct_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_AbsBias_direct.csv"))

## models 
write.csv(AbsBias_export_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias.csv"))

#AbsBias improvements 
write.csv(AbsBias_improve_race(AbsBias_eq_1, AbsBias_eq_gss, AbsBias_pr_1, AbsBias_pr_05, AbsBias_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_AbsBias_improve.csv"))


#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=T, FE=T), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=T, FE=T),  
       width = 8, height = 4, dpi = 300)


###NOT YET COMPLETE
#****************************************************************************
#3 .Compute difference in Absolute Bias by race for FH and Direct Estimates ----
#****************************************************************************

#2.1 No covariates no FE ----
folder <-  "nocovars_noFE"

diff_fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_diff_fh_eq_1.csv"), header = T)
diff_fh_eq_1 = diff_fh_eq_1[,-1]
diff_fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_diff_fh_eq_gss.csv"), header = T)
diff_fh_eq_gss = diff_fh_eq_gss[,-1]
diff_fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_diff_fh_pr_05.csv"), header = T)
diff_fh_pr_05 = diff_fh_pr_05[,-1]
diff_fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_diff_fh_pr_1.csv"), header = T)
diff_fh_pr_1 = diff_fh_pr_1[,-1]
diff_fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_diff_fh_pr_gss.csv"), header = T)
diff_fh_pr_gss = diff_fh_pr_gss[,-1]


## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias




#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

#2.2  No covariates & FE ----

folder <-  "nocovars_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias
#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

#2.3 Simple covariates & no FE ---- 

folder <-  "covars_simple_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias

#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)
#2.4 Simple covariates & FE ----

folder <-  "covars_simple_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias
#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

#2.5 Stepwise covariates & no FE ----

folder <-  "covars_step_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias
#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

#2.6 Stepwise covariates & FE ---- 
folder <-  "covars_step_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
AbsBias_eq_1 <- inner_join(AbsBias(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           AbsBias(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_eq_gss <- inner_join(AbsBias(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_1 <- inner_join(AbsBias(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           AbsBias(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_05 <- inner_join(AbsBias(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            AbsBias(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both AbsBias 

AbsBias_pr_gss <- inner_join(AbsBias(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             AbsBias(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both AbsBias





#AbsBias plots
subtitle <- "Equal sample size for all races (1% ACS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eq1.png"), 
       plot = plotAbsBias(AbsBias_eq_1, st_eq_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Equal sample size for all races (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_eqgss.png"), 
       plot = plotAbsBias(AbsBias_eq_gss, st_eq_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (1% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr1.png"), 
       plot = plotAbsBias(AbsBias_pr_1, st_pr_1, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (0.5% ACS size)"
ggsave(here("04_results", "plots",  folder, "scatter_AbsBias_pr05.png"), 
       plot = plotAbsBias(AbsBias_pr_05, st_pr_05, subtitle, covars=F, FE=F), 
       width = 8, height = 4, dpi = 300)

subtitle <- "Proportional sample size for each race (GSS size)"
ggsave(here("04_results", "plots", folder, "scatter_AbsBias_prgss.png"), 
       plot = plotAbsBias(AbsBias_pr_gss, st_pr_gss, subtitle, covars=F, FE=F),  
       width = 8, height = 4, dpi = 300)




#****************************************************************************
#4.Compute Relative MSE for FH and Direct Estimates from models ----
#****************************************************************************
#*Note: there are several ways to normalize MSE here we have used division by mean of the truth (house_ten)
#*

#4.1 No covariates no FE ----
folder <-  "nocovars_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
mse_r_eq_1 <- inner_join(mse_r(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           mse_r(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_eq_gss <- inner_join(mse_r(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_1 <- inner_join(mse_r(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           mse_r(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_05 <- inner_join(mse_r(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            mse_r(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_gss <- inner_join(mse_r(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r


#mse_r export 
##directs
write.csv(mse_r_direct_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_mse_r_direct.csv"))


#2.2  No covariates & FE ----

folder <-  "nocovars_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
mse_r_eq_1 <- inner_join(mse_r(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           mse_r(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_eq_gss <- inner_join(mse_r(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_1 <- inner_join(mse_r(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           mse_r(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_05 <- inner_join(mse_r(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            mse_r(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_gss <- inner_join(mse_r(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r

#mse_r export 
##directs
write.csv(mse_r_direct_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_mse_r_direct.csv"))

## models 
write.csv(mse_r_export_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_mse_r.csv"))



#2.3 Simple covariates & no FE ---- 

folder <-  "covars_simple_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
mse_r_eq_1 <- inner_join(mse_r(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           mse_r(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_eq_gss <- inner_join(mse_r(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_1 <- inner_join(mse_r(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           mse_r(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_05 <- inner_join(mse_r(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            mse_r(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_gss <- inner_join(mse_r(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r

#mse_r export 
##directs
write.csv(mse_r_direct_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_mse_r_direct.csv"))

## models 
write.csv(mse_r_export_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_mse_r.csv"))



#2.4 Simple covariates & FE ----

folder <-  "covars_simple_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
mse_r_eq_1 <- inner_join(mse_r(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           mse_r(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_eq_gss <- inner_join(mse_r(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_1 <- inner_join(mse_r(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           mse_r(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_05 <- inner_join(mse_r(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            mse_r(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_gss <- inner_join(mse_r(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r

#mse_r export 
##directs
write.csv(mse_r_direct_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_mse_r_direct.csv"))

## models 
write.csv(mse_r_export_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_mse_r.csv"))



#2.5 Stepwise covariates & no FE ----

folder <-  "covars_step_noFE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
mse_r_eq_1 <- inner_join(mse_r(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           mse_r(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_eq_gss <- inner_join(mse_r(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_1 <- inner_join(mse_r(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           mse_r(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_05 <- inner_join(mse_r(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            mse_r(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_gss <- inner_join(mse_r(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r

#mse_r export 
##directs
write.csv(mse_r_direct_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_mse_r_direct.csv"))

## models 
write.csv(mse_r_export_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_mse_r.csv"))





#2.6 Stepwise covariates & FE ---- 
folder <-  "covars_step_FE"

fh_eq_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_eq_1.csv"), header = T)
fh_eq_1 = fh_eq_1[,-1]
fh_eq_gss = read.csv( here("04_results", "estimates", folder,"Estimates_fh_eq_gss.csv"), header = T)
fh_eq_gss = fh_eq_gss[,-1]
fh_pr_05 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_05.csv"), header = T)
fh_pr_05 = fh_pr_05[,-1]
fh_pr_1 = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_1.csv"), header = T)
fh_pr_1 = fh_pr_1[,-1]
fh_pr_gss = read.csv( here("04_results", "estimates", folder, "Estimates_fh_pr_gss.csv"), header = T)
fh_pr_gss = fh_pr_gss[,-1]

## Calculating absolute bias from the file 04_results
mse_r_eq_1 <- inner_join(mse_r(direct_eq_1, outcome_df, grouping_var, outcome_finalvar), #direct 
                           mse_r(fh_eq_1, outcome_df, grouping_var, outcome_finalvar, FH=T), #FH
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_eq_gss <- inner_join(mse_r(direct_eq_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_eq_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_1 <- inner_join(mse_r(direct_pr_1, outcome_df, grouping_var, outcome_finalvar), 
                           mse_r(fh_pr_1, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                           by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_05 <- inner_join(mse_r(direct_pr_05, outcome_df, grouping_var, outcome_finalvar), 
                            mse_r(fh_pr_05, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                            by = grouping_var, na_matches = "never") #join both mse_r 

mse_r_pr_gss <- inner_join(mse_r(direct_pr_gss, outcome_df, grouping_var, outcome_finalvar), 
                             mse_r(fh_pr_gss, outcome_df, grouping_var, outcome_finalvar, FH=T), 
                             by = grouping_var, na_matches = "never") #join both mse_r



#mse_r export 
##directs
write.csv(mse_r_direct_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", "summary_mse_r_direct.csv"))

## models 
write.csv(mse_r_export_race(mse_r_eq_1, mse_r_eq_gss, mse_r_pr_1, mse_r_pr_05, mse_r_pr_gss, "race", "ST"), 
          here("04_results", "plots", folder, "summary_mse_r.csv"))




