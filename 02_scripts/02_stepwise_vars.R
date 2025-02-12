#  
# title: "Variable selection by Stepwise"
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
               StepReg, caret, beepr, ggplot2, writexl, survey, sampling, sae, janitor, epiDisplay, future.apply)
##fix priorities with dplyr
select <- dplyr::select
filter <- dplyr::filter

#routes 
here::i_am("02_scripts/02_stepwise_vars.R")
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

ids <- c("SERIALNO", "DIVISION", "PUMA", "REGION") #SPORDER for individuals
grouping_var0 <- c("RAC1P", "HISP") #inputs to construct race
covars_p <- c("AGEP", "MAR") #ENG no variation
covars_h <- c("BDSP", "LAPTOP")
outcome_var <- c("TEN")

# Unzip file first, if not col_select doesnt work
ACS_p <- rbind(read_sas(here("03_data", "ACS", "psam_pusa.sas7bdat")), 
               read_sas(here("03_data", "ACS", "psam_pusb.sas7bdat")))
ACS_p_labels <- lapply(var_label(ACS_p), str_conv, encoding="UTF-8")

ACS_h <- rbind(read_sas(here("03_data", "ACS", "psam_husa.sas7bdat")),
               read_sas(here("03_data", "ACS", "psam_husb.sas7bdat")))
ACS_h_labels <- lapply(var_label(ACS_h), str_conv, encoding="UTF-8")

# Join 
df <- inner_join(ACS_p, ACS_h, by = ids, multiple = "all", na_matches = "never") 
labels <- c(ACS_p_labels, ACS_h_labels)
rm(list = ls()[startsWith(ls(), "ACS")])

df <- df %>% 
  select(-ids) %>% # keep ST only 
  mutate_all(funs(as.numeric(as.factor(.))))  

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

outcome_finalvar <- "house_ten"


#**************************************
#4. Stepwise with population for covariates ---- 
#**************************************

# sample 1% for stepwise
sample_df <- sample_n(df, dim(df)[1]/100)

#Formula 
y <- "house_ten"
not_candidates <- c("race", grouping_var0, y,
                    c("WGTP", "NP", "PWGTP", "NAICSP", 
                      "NOP", "POBP", "POWSP", "RAC2P", 
                      "RAC3P", "RACASN", "RACBLK", "RACNUM", 
                      "RACSOR", "RACWHT", "WAOB"))

# Drop vars with 0 variance
sample_df <- sample_df[, -nearZeroVar(sample_df)]
candidates_df <- sample_df %>%
  select(-all_of(not_candidates), -matches("^PWGTP|^WGTP")) %>%
  select(where(~!any(is.na(.x))))
step_formula <- as.formula(paste(y, " ~ ", 
                                 paste(colnames(candidates_df), collapse= "+")))

#standardize
#candidates_df <- as.data.frame(scale(candidates_df))

#Stepwise
stp.sel <- stepwise(formula=step_formula, 
                    data=cbind(candidates_df, house_ten=sample_df[[y]]), #paste outcome
                    selection="bidirection", #bidirection, score
                    sle=0.01, # significance is irrelevant 
                    sls=0.01, # significance is irrelevant 
                    select="AICc") #2°order AIC, ↑penalization if ↓n

#Selected vars (no constant)
stepw.vars <- stp.sel[["Selected Varaibles"]][!stp.sel[["Selected Varaibles"]] %in% c('1')]

#Selected vars dataframe
stepwise_df <- cbind(const=rep(1,dim(df)[1]),
                     df[,names(df) %in% c(stepw.vars, grouping_var)])

stepwise_df <- cbind(stepwise_df, ST=df$ST.x)

# 2.4 Covariates and pop counts/freqs ----
nogrouping_vars <- colnames(stepwise_df)[1:(dim(stepwise_df)[2]-2)] 
collapsed_stepwise_df <- stepwise_df %>%
  group_by(.dots = lapply(grouping_var, as.symbol)) %>%   # grouping vars as symbols
  summarise(n = n(),                                      # total population by area
            across(all_of(nogrouping_vars),
                   \(x) mean(x, na.rm = TRUE))) %>%       # mean considering NAs
  mutate(freq = n/sum(n))                                 # freq: %race1 between all ST


##### END ####
#save 
saveRDS(stepwise_df, here("04_results", "stepwise_vars.rds"))
beep(8)