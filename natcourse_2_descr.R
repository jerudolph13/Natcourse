
############################################################################################  
#
# Project: g-formula natural course paper
#
# Purpose: Summarize and visualize the data used in time-fixed and time-varying analyses
#
# Author: Jacqueline Rudolph
#
# Last Update: 10 Nov 2020
#
############################################################################################  


library("tidyverse")


# Time-fixed analysis -----------------------------------------------------

eager_base <- read.table(file="../data/eager_base.txt", sep="\t", header=TRUE) %>% 
  filter(week == 1)

  #Summarize baseline variables
  summ_base <- eager_base %>% 
    summarize(n.comp = sum(compliance),
              avg.comp = mean(compliance),
              n.trt = sum(treatment),
              avg.trt = mean(treatment),
              n.eligi = sum(eligibility),
              avg.eligi = mean(eligibility),
              n.smoke = sum(smoke),
              avg.smoke = mean(smoke),
              avg.age = mean(age),
              sd.age = sd(age),
              avg.bmi = mean(BMI),
              sd.bmi = sd(BMI))
  
  #Summarize by compliance
  summ_bycomp <- eager_base %>% 
    group_by(compliance) %>% 
    summarize(n.trt = sum(treatment),
              avg.trt = mean(treatment),
              n.eligi = sum(eligibility),
              avg.eligi = mean(eligibility),
              n.smoke = sum(smoke),
              avg.smoke = mean(smoke),
              avg.age = mean(age),
              sd.age = sd(age),
              avg.bmi = mean(BMI),
              sd.bmi = sd(BMI))  


# Time-varying analysis ---------------------------------------------------

eager_tvar <- read.table(file="../data/eager_tvar.txt", sep="\t", header=TRUE)

  #Average across all person-weeks
  summ_tvar <- eager_tvar %>% 
    summarize(avg.comp = mean(compliance),
              avg.bleed = mean(bleed),
              avg.nausea = mean(nausea))
  
  summ_trt <- eager_tvar %>% 
    group_by(treatment) %>% 
    summarize(avg.comp = mean(compliance))
  
  #Average by week
  weekly <- eager_tvar %>% 
    group_by(week) %>% 
    summarize(comp.wk = mean(compliance),
              bleed.wk = mean(bleed),
              nausea.wk = mean(nausea))
  
  weekly_trt <- eager_tvar %>% 
    group_by(treatment, week) %>% 
    summarize(comp.wk = mean(compliance))
  
  #Visualize
  thm <- theme_classic() +
    theme(
      legend.position = "top",
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA)
    )
  
  #Compliance
  ggplot(data=weekly, aes(x=week, y=comp.wk)) + thm +
    geom_line() +
    geom_point() +
    xlab("Week") + scale_x_continuous(expand=c(0, 0), limits=c(0, 27)) +
    ylab("Average compliance") + scale_y_continuous(expand=c(0, 0), limits=c(0.5, 1.0))
  
  #Bleeding
  ggplot(data=weekly, aes(x=week, y=bleed.wk)) + thm +
    geom_line() +
    geom_point() +
    xlab("Week") + scale_x_continuous(expand=c(0, 0), limits=c(0, 27)) +
    ylab("Average bleeding") + scale_y_continuous(expand=c(0, 0), limits=c(0, 1))
  
  #Nausea
  ggplot(data=weekly, aes(x=week, y=nausea.wk)) + thm +
    geom_line() +
    geom_point() +
    xlab("Week") + scale_x_continuous(expand=c(0, 0), limits=c(0, 27)) +
    ylab("Average nausea") + scale_y_continuous(expand=c(0, 0), limits=c(0.15, 0.5))

