#  
# title: "Draw repeated samples from ACS"
# subtitle: "DREI2023_SAE"
# author: "[Angelo Cozzubo] \n(https://sites.google.com/pucp.pe/acozz)"
# date creation: "07/April/2023"
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
here::i_am("02_scripts/01_loadACS.R")
source(here("02_scripts","my_fns.R")) # load my functions 

#theme
theme <- theme_classic() #ggplot theme 
theme$plot.title$hjust <- 0.5

#seed
set.seed(123123) # for stochastic part 
#parallel 
plan(multisession)

#**************************************
#1. Read data, labels, filter and params ---- 
#**************************************

ids <- c("ST", "SERIALNO", "DIVISION", "PUMA", "REGION") #SPORDER for individuals
grouping_var0 <- c("RAC1P", "HISP") #inputs to construct race
covars_p <- c("AGEP", "MAR") #ENG no variation
covars_h <- c("BDSP", "LAPTOP")
outcome_var <- c("TEN")

# Unzip file first, if not col_select doesnt work
ACS_p <- rbind(read_sas(here("03_data", "ACS", "psam_pusa.sas7bdat"), 
                        col_select = all_of(c(ids, covars_p, grouping_var0))), 
               read_sas(here("03_data", "ACS", "psam_pusb.sas7bdat"), 
                        col_select = all_of(c(ids, covars_p, grouping_var0))))
ACS_p_labels <- lapply(var_label(ACS_p), str_conv, encoding="UTF-8")

ACS_h <- rbind(read_sas(here("03_data", "ACS", "psam_husa.sas7bdat"),
                 col_select = all_of(c(ids, outcome_var, covars_h))),
               read_sas(here("03_data", "ACS", "psam_husb.sas7bdat"),
                        col_select = all_of(c(ids, outcome_var, covars_h))))
ACS_h_labels <- lapply(var_label(ACS_h), str_conv, encoding="UTF-8")

# Join 
df <- inner_join(ACS_p, ACS_h, by = ids, multiple = "all", na_matches = "never") 
labels <- c(ACS_p_labels, ACS_h_labels)
rm(list = ls()[startsWith(ls(), "ACS")])

# destring 
covars <- c(covars_p, covars_h)
df <- df %>% 
  mutate(across(all_of(c(outcome_var, covars, grouping_var0)), as.integer)) %>% 
  select(-ids[-1]) # keep ST only 


#**************************************
#2. Grouping variables and recodes ---- 
#**************************************

# 2.1 Grouping vars ----
grouping_var <- c("race", "ST") 

# Create the 'race' variable
## Nonhispanic white = white alone + not hispanic
## Nonhispanic black = black alone + not hispanic
## hispanic = hispanic 
## other = else 

df <- df %>%
  mutate(race = NA) %>%
  mutate(race = ifelse(RAC1P == 1 & HISP == 1, 1, race)) %>%
  mutate(race = ifelse(RAC1P == 2 & HISP == 1, 2, race)) %>%
  mutate(race = ifelse(HISP != 1 & !is.na(HISP), 3, race)) %>%
  mutate(race = ifelse(is.na(race), 4, race))

var_label(df$race) <- "Race"
df$race <- labelled(df$race, labels = c("Nonhispanic White" = 1, 
                                        "Nonhispanic black" = 2, 
                                        "Hispanic" = 3, 
                                        "Other race or multirace" = 4))
    # RC1P  (ACS documentation) 
    # "1": "White alone",
    # "2": "Black or African American alone",
    # "3": "American Indian alone",
    # "4": "Alaska Native alone",
    # "5": "American Indian and Alaska Native tribes specified",
    # "6": "Asian alone",
    # "7": "Native Hawaiian and Other Pacific Islander alone",
    # "8": "Some other race alone"
    # "9": "Two or More Races",
    
    # HISP (ACS documentation)
    # "1": "Not Spanish/Hispanic/Latino"
    # else (2-24): "Spanish/Hispanic/Latino"

# 2.2 Outcome var ----
## Recode outcome variable HH Tenancy (1=owned or mortgage, 0=otherwise)
df <- df %>%
  mutate(house_ten = case_when(
    TEN == 2 ~ 1,
    TEN == 3 ~ 0,
    TEN == 4 ~ 0,
    TRUE ~ TEN)) %>% 
  select(-TEN)

# Drop missings for outcome variable
df <- df %>% filter(!is.na(house_ten))

outcome_finalvar <- "house_ten"

# 2.3 Covariates var ----
df <- df %>% 
  mutate(LAPTOP = ifelse(LAPTOP == 2, 0, LAPTOP),      # dummy has a laptop  
         MAR = ifelse(MAR %in% c(2, 3, 4, 5), 0, MAR)) %>%  # dummy of married
  select(-grouping_var0) # keep ST only 

# 2.4 Covariates and pop counts/freqs ----
collapsed_df <- df %>%
  group_by(.dots = lapply(grouping_var, as.symbol)) %>%   # grouping vars as symbols 
  summarise(n = n(),                                      # total population by area
            across(all_of(covars), 
                   \(x) mean(x, na.rm = TRUE))) %>%       # mean considering NAs
  mutate(freq = n/dim(df)[1])                             # freq: %race1 between all ST
                                                          # sum across ST = 1 for e/ race

outcome_df <- df %>%
  group_by(.dots = lapply(grouping_var, as.symbol)) %>%   # grouping vars as symbols 
  summarise(across(all_of(outcome_finalvar), 
                   \(x) mean(x, na.rm = TRUE)))           # mean considering NAs

pdf(here("04_results", "outcome_histogram.pdf"))
hist(outcome_df$house_ten, 
     main = paste("Histogram of", colnames(outcome_df)[3]), 
     xlab = "Values", 
     ylab = "Frequency",
     col = "grey",
     border = "black",
     breaks = 5,
     xlim = c(0, 1),
     ylim = c(0, 60),
     freq = TRUE)
dev.off()

##!! check all states has all races! ---- 
##!! order matter for the sum,freqs! ---- 

# 2.5 Population sizes and number of areas 
total_num_area <- dim(collapsed_df)[1]
total_pop_size <- dim(df)[1]

#**************************************
#4. Draw stratified samples ---- 
#**************************************

# 4.0 Set params and filter data 

# GSS sample size 
gss_size <- 4032             #gss.norc.org/Documents/codebook/GSS%202021%20Codebook%20R1.pdf
# Number of simulations 
n_sim <- 500
# Drop unneeded covars 
df <- df %>% select(-covars) #dont need covars for sampling

# 4.1 List of strata sample sizes ----
## order of vector: first races then ST. Example
## [race1_ST1, ..., race1_STk, race2_ST1,...race2_STk,...,raceJ_STk]

## Race sample proportional to size in ACS pop (1% size of ACS)
st_pr_1 <- collapsed_df$freq*(total_pop_size/100)
st_pr_1 <- fix_list(st_pr_1, collapsed_df$n)

## Race sample proportional to size in ACS pop (0.5% size of ACS)
st_pr_05 <- collapsed_df$freq*(total_pop_size/200)
st_pr_05 <- fix_list(st_pr_05, collapsed_df$n)

## Race sample proportional to size in GSS pop 
st_pr_gss <- collapsed_df$freq*(gss_size) #/length(unique(df$race)))
st_pr_gss <- fix_list(st_pr_gss, collapsed_df$n)

## Equal sample size for each race (1% size of size)
st_eq_1 <- rep(length(collapsed_df$freq)^-1, length(collapsed_df$freq))*(total_pop_size/100)
st_eq_1 <- fix_list(st_eq_1, collapsed_df$n)

## Equal sample size for each race (GSS size)
st_eq_gss <- rep(length(collapsed_df$freq)^-1, length(collapsed_df$freq))*(gss_size)
st_eq_gss <- fix_list(st_eq_gss, collapsed_df$n) 

# 4.2 Draw n_reps stratified samples ---- 
## Watch out! This takes a lot of time to run
gc()
Sys.setenv(R_MAX_MEM_SIZE = "5Gb") # allocate 5gb 
tic()
df_eq_gss <- draw_sts_samples(n_sim, grouping_var, st_eq_gss, df) 
df_eq_1 <- draw_sts_samples(n_sim, grouping_var, st_eq_1, df) 
df_pr_gss <- draw_sts_samples(n_sim, grouping_var, st_pr_gss, df) 
df_pr_1 <- draw_sts_samples(n_sim, grouping_var, st_pr_1, df)
df_pr_05 <- draw_sts_samples(n_sim, grouping_var, st_pr_05, df)
toc()

##### END ####
#save 
save.image(here("04_results","230818_500iter.RData"))
beep(8)