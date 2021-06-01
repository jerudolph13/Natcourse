
###############################################################################################
#
# Project: Natural course paper
#
# Purpose: Read in EAGeR data and prep for use in analyses
#
# Author: Jacqueline Rudolph
#
# Last Update: 11 Nov 2020
#
###############################################################################################


packages <- c("tidyverse")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

eager <- read.csv(file="../data/2020_11_06-eager_weekly.csv", header=TRUE) %>% 
  select(id, week, GA, treatment, age, BMI, smoke, eligibility, studymed_b_imputed, bleed_d_imputed, 
         nausea_d_imputed, outcome) %>% 
  #Remove ID 417199 due to missing information
  filter(id != 417199)


# Manipulate data ---------------------------------------------------------

#Define exposure and covariates
eager <- eager %>% 
    mutate(
      #baseline variables
      eligibility = as.numeric(eligibility == "new"),
           
      #t-varying compliance exposure
      compliance = as.numeric(studymed_b_imputed >= 5/7), 
      compliance_last = if_else(week == 1, 0, as.double(lag(compliance))),
           
      #t-varying confounders
      bleed = as.numeric(bleed_d_imputed > 0), 
      bleed_last = if_else(week == 1, 0, as.double(lag(bleed))),
      nausea = as.numeric(nausea_d_imputed > 0),
      nausea_last = if_else(week == 1, 0, as.double(lag(nausea))),
           
      #flag last record
      last = as.numeric(!duplicated(id, fromLast=T)))

#Define the outcome 
  #This treats women with unknown GA but a post-pregnancy outcome as censored when we last see them
  #We could instead impute GA for these women
eager_last <- eager %>% 
  filter(last==1) %>% 
  mutate(c_week = week - GA) %>% 
  select(id, c_week)
eager <- left_join(eager, eager_last, by="id") %>% 
  mutate(conception = ifelse(GA==99, 0, as.numeric((week - c_week) >= 0)))

#Define drop out
eager$drop <- ifelse(eager$last==1 & eager$outcome=="withdrawal", 1, 
              ifelse(eager$last==1 & eager$conception==0 & eager$outcome %in% c("live birth", "pregnancy loss"), 1, 0))

#Subset to records through event or censoring
eager <- filter(eager, week <= 26)

#Remove extra records after conception
eager <- eager %>% 
  group_by(id) %>% 
  mutate(cum_concep = cumsum(conception)) %>% 
  filter(cum_concep <= 1) %>% 
  ungroup(id) %>% 
  select(-c(cum_concep, last))

#For t-fixed analysis, obtain each woman's average compliance through 26 weeks of follow-up
  #Compare average compliance to the 5/7 rule
  comp <- eager %>% 
    group_by(id) %>% 
    summarize(avg.comp=mean(studymed_b_imputed)) %>% 
    mutate(compliance=as.numeric(avg.comp>=5/7)) %>% 
    select(-avg.comp)


# Select final variables --------------------------------------------------

#For time-fixed analysis, keep relevant baseline variables, outcome, and drop out
eager_base <- eager %>% 
    select(id, week, treatment, age, BMI, smoke, eligibility, conception, drop) %>% 
    inner_join(comp, by="id")  #Merge average compliance

#For time-varying analysis, keep relevant baseline variables, compliance, bleeding, nausea, outcome, and drop out
eager_tvar <- select(eager, id, week, treatment, age, BMI, smoke, eligibility, compliance, compliance_last, 
                     bleed, bleed_last, nausea, nausea_last, conception, drop)


###################################################################################################################################  
# Output data

write.table(eager_base, file="../data/eager_base.txt", sep="\t", row.names=FALSE)
write.table(eager_tvar, file="../data/eager_tvar.txt", sep="\t", row.names=FALSE)

