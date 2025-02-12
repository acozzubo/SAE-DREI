





#AbsBias <- function(estimates_df, outcome_df, grouping_var, outcome_finalvar, FH=F, diff=F) {
  

#Obtain AbsBias estimate for each domain      
  #input: estimates_df: dataframe with estimates from direct_estim or FH_estim fns
  #       outcome_df: dataframe with population value of parameter 
  #       outcome_finalvar: outcome variable to obtain mean
  #       grouping_var: grouping var, areas
  #       FH: T if computing the for FH results (renaming)
  #       diff: T if computing for race difference (sort only 1 var)
  #return: df with AbsBias and domain IDs    
  
estimates_df = fh_eq_1 
outcome_df
outcome_finalvar <- "house_ten"
grouping_var <- c("race", "ST") 


    estimates_df <- estimates_df %>%
      select(starts_with("FH.estimates"), all_of(grouping_var)) %>%
      rename_with(~paste0("house_ten_", .), starts_with("FH.estimates"))

  
  # Create a new dataframe with selected replications, and grouping vars 
  AB_df <- data.frame(estimates_df[, grouping_var],
                      estimates_df[, grep(paste0("^", outcome_finalvar), names(estimates_df))])
  

  
  # Subtract population value of outcome from selected columns in new_dataframe and take absolute value 
  AB_df[, -(1:length(grouping_var))] <- abs(AB_df[, -(1:length(grouping_var))] - 
                                              outcome_df[[outcome_finalvar]])
  
  # Average the difference and keep the final AbsBias only 
  AB_df <- AB_df %>%
    mutate(Abs_Bias = rowMeans(select(., -one_of(grouping_var)), na.rm=TRUE)) %>%
    select(all_of(grouping_var),Abs_Bias)
  
  if(FH){AB_df <- AB_df %>% rename(Abs_Bias_FH = Abs_Bias)} #rename to compare with direct estim
  
  return(AB_df) 



library(ggplot2)
library(see)

outcome_df$race2 = sort(rep(c(1,2,3,4), 51))  
outcome_df$race3 = c(rep(c("Nonhispanic White"), 51), rep(c("Nonhispanic Black"), 51), rep(c("Hispanic"), 51), rep(c("Other Race or Multiracial"), 51))
  
hist(outcome_df$house_ten, breaks = 10) #qqplot as well, maybe a few other normality checks
ggplot(data = outcome_df, aes(x= house_ten)) +geom_histogram()
 
ggplot(data = outcome_df, aes(x= house_ten)) +geom_histogram() + facet_wrap(vars(as.character(race2)))

ggplot(data = outcome_df, aes(y= house_ten, x = race2)) +geom_point()
 
ggplot(data = outcome_df, aes(y= house_ten, x = as.character(race3))) +geom_violinhalf() + 
  labs(title ="Distribution of House Tenancy (Outcome Variable) in the US by Race", 
      subtitle = "House Tenancy is Calculated at the State Level from ACS data",
      x = "Racial Category",  y = "House Tenancy" )


qqplot(y = outcome_df$house_ten, x= outcome_df$race2)
   
  summary(outcome_df$house_ten)
  
  hist(direct_eq_1$house_ten)
  hist(direct_eq_1$house_ten.1)
  hist(direct_eq_1$house_ten.134)
  
  hist(direct_eq_1$house_ten, breaks = 20)
  hist(direct_eq_gss$house_ten, breaks = 20)
  hist(direct_pr_05$house_ten, breaks = 20)
  hist(direct_pr_1$house_ten, breaks = 20)
  hist(direct_pr_gss$house_ten, breaks = 20)

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

#prgss - race
ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias_FH, x=log(size))) +geom_point() +
  geom_jitter() + facet_wrap(vars(race)) +
  labs(title ="Bias of Fay-Harriot Estimates by Race", 
       subtitle = "Model with Simple Covariates and Fixed Effects for Race", 
       caption ="With race sample proportional to size in GSS population",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_prgss_Race.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_pr_gss, aes(y=Abs_Bias, x=log(size))) +geom_point() +
  geom_jitter() + facet_wrap(vars(race))+
  labs(title ="Bias of Direct Estimates by Race", 
       subtitle = "Model with Simple Covariates and Fixed Effects for Race", 
       caption ="With race sample proportional to size in GSS population",
       x = "Log of the Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_prgss_race.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)
#eqgss
AbsBias_eq_gss$size = st_eq_gss
ggplot(data = AbsBias_eq_gss, aes(y=Abs_Bias_FH, x=(size))) +geom_violin() + ylim(0.02,.12) +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the GSS population",
       x = "Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_eqgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_eq_gss, aes(y=Abs_Bias, x=(size))) +geom_violin() + ylim(0.02,.12) +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the GSS population",
       x = "Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_eqgss.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

#eq1
AbsBias_eq_1$size = st_eq_1
ggplot(data = AbsBias_eq_1, aes(y=Abs_Bias_FH, x=size)) +geom_violin(0,.1) +
  labs(title ="Bias of Fay-Harriot Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the ACS population",
       x = "Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBiasFH_eq1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

ggplot(data = AbsBias_eq_1, aes(y=Abs_Bias, x=size)) +geom_violin(0,.1) +
  labs(title ="Bias of Direct Estimates from Model with Simple Covariates and Fixed Effects for Race", 
       subtitle = "Where the race sample is equal from the ACS population",
       x = "Sample Size",  y = "Absolute Bias of the Estimate" )
ggsave(here("04_results", "plots", folder, "size_AbsBias_eq1.png"), 
       plot = last_plot(), 
       width = 8, height = 4, dpi = 300)

