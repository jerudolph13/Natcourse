
###################################################################################################################################  
#
# Project: g-formula natural course paper
#
# Purpose: Estimate natural course in EAGeR -- Is compliance associated with becoming preg within 6 months? (time-varying)
#
# Author: Jacqueline Rudolph
#
# Last Update: 19 May 2020
#
###################################################################################################################################  

setwd("~/Documents/Pitt Projects/Natural Course/data")

packages <- c("survival", "tidyverse", "survminer", "splines")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

#What size Monte Carlo resample?
montecarlo <- 5000

#When do we admin censor?
K <- 26

#Set seed for predictions
set.seed(123)


###################################################################################################################################  
# Read in data

eager <- read.table(file="eager_tvar.txt", sep="\t", header=TRUE) %>% 
  mutate(last = as.numeric(!duplicated(id, fromLast=T)),
         age = scale(age),
         BMI = scale(BMI))

#Create a version of data with only the last record
  eager_last <- filter(eager, last == 1)

  
###################################################################################################################################  
# What was the observed survival?
  
surv_obs <- survfit(Surv(week, delta) ~ 1, data=eager_last)
  #Survival by assigned treatment
  surv_exp <- survfit(Surv(week, delta) ~ treatment, data=eager_last)


###################################################################################################################################  
# Run parametric g-formula
  
#Model nausea (N), bleeding (B), exposure (A), censoring (C), and outcome (T) <- in that temporal order
  #Run all models stratified by treatment arm
pgf.mod <- function(k) {
  mod.n <- glm(nausea ~ compliance_last + bleed_last + nausea_last + 
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=eager, subset=(treatment==k))
  
  mod.b <- glm(bleed ~ compliance_last + bleed_last + nausea + nausea_last + 
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=eager, subset=(treatment==k))
  
  mod.a <- glm(compliance ~ compliance_last + bleed + bleed_last + nausea + nausea_last +
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=eager, subset=(treatment==k))
  
  mod.c <- glm(drop ~ compliance + compliance_last + bleed + bleed_last + nausea + nausea_last +
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=eager, subset=(treatment==k))
  
  mod.t <- glm(delta ~ compliance + compliance_last + bleed + bleed_last + nausea + nausea_last +
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=eager, subset=(treatment==k))
  
  return(list(mod.n, mod.b, mod.a, mod.c, mod.t))
}
#Store model objects in a list we will use in the pgf function below
mods <- lapply(0:1, function(x){pgf.mod(x)})
  
#Monte Carlo (MC) implementation
  #Resample from all women
  MC0 <- eager %>% 
    filter(week == 1) %>%  #Just want the first record
    select(-c(id, week, last, delta, drop))  #Keep only baseline variables
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  #Sample with replacement
  MC<-MC0[index,]  #Grab the sampled records
  MC$id<-1:montecarlo  #Assign new IDs
    
  #Resample from those assigned to aspirin
  aspirin <- filter(eager, treatment==1)
  MC0 <- aspirin %>% 
    filter(week == 1) %>%  
    select(-c(id, week, last, delta, drop))  
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  
  MC_aspirin <- MC0[index,]  
  MC_aspirin$id <- 1:montecarlo 
    
  #Resample from those assigned to placebo
  placebo <- filter(eager, treatment==0)
  MC0 <- placebo %>% 
    filter(week == 1) %>%  
    select(-c(id, week, last, delta, drop))  
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  
  MC_placebo <- MC0[index,]  
  MC_placebo$id <- 1:montecarlo 
  
#Use parametric g-formula (pgf) function to reconstruct follow-up for each person
  pgf <- function(i, data, maxint){
    #Grab record for individual i
    g_data <- data[data$id==i, ]  
    #Initialize variables
    treatment<-eligibility<-age<-BMI<-smoke<-nausea<-nausea_last<-bleed<-bleed_last <- numeric()
    compliance<-compliance_last<-week<-drop<-delta <- numeric()

    #Time point 1 (baseline variables take observed values)
      week[1] <- 1
      treatment[1] <- g_data$treatment
      eligibility[1] <- g_data$eligibility
      age[1] <- g_data$age
      BMI[1] <- g_data$BMI
      smoke[1] <- g_data$smoke
      #For t-varying variables, assign last value to be 0 in the first record
      nausea[1] <- g_data$nausea; nausea_last[1] <- 0
      bleed[1] <- g_data$bleed; bleed_last[1] <- 0
      compliance[1] <- g_data$compliance; compliance_last[1] <- 0
    
      pred_data <- data.frame(week=week[1], treatment=treatment[1], eligibility=eligibility[1], age=age[1], BMI=BMI[1],
                              smoke=smoke[1], nausea=nausea[1], nausea_last=nausea_last[1], bleed=bleed[1], bleed_last=bleed_last[1],
                              compliance=compliance[1], compliance_last=compliance_last[1])
      
    #Did they drop out?
      #Note: need to select correct model object from the list "mods", 
      #first based on assigned treatment and then on which variable we are predicting
      drop[1] <- as.numeric(predict(mods[[treatment[1]+1]][[4]], newdata=pred_data, type="response")>runif(1))
    
    #If censored, no outcome; else, see if there is an outcome
      if (drop[1]==1) {
        delta[1] <- 0
      } else {
        delta[1] <- as.numeric(predict(mods[[treatment[1]+1]][[5]], newdata=pred_data, type="response")>runif(1))
      }
      
      #Repeat for other time points until either drop out or event occurs
      for (t in 2:maxint) {
        if (delta[t-1]==0 & drop[t-1]==0) {
          week[t] <- t
          treatment[t] <- treatment[t-1]
          eligibility[t] <- eligibility[t-1]
          age[t] <- age[t-1]
          BMI[t] <- BMI[t-1]
          smoke[t] <- smoke[t-1]
          nausea_last[t] <- nausea[t-1]
          bleed_last[t] <- bleed[t-1]
          compliance_last[t] <- compliance[t-1]
          
          pred_data <- data.frame(week=week[t], treatment=treatment[t], eligibility=eligibility[t], age=age[t], BMI=BMI[t], 
                                  smoke=smoke[t], nausea_last=nausea_last[t], bleed_last=bleed_last[t], 
                                  compliance_last=compliance_last[t])
        
        #Predict nausea at t
          nausea[t] <- as.numeric(predict(mods[[treatment[t]+1]][[1]], newdata=pred_data, type="response")>runif(1))
          pred_data <- cbind(pred_data, nausea=nausea[t])
          
        #Predict bleeding at t
          bleed[t] <- as.numeric(predict(mods[[treatment[t]+1]][[2]], newdata=pred_data, type="response")>runif(1))
            pred_data <- cbind(pred_data, bleed=bleed[t])
          
        #Predict compliance at t
          compliance[t] <- as.numeric(predict(mods[[treatment[t]+1]][[3]], newdata=pred_data, type="response")>runif(1))
            pred_data <- cbind(pred_data, compliance=compliance[t])
            
        #Did they drop out?
          drop[t] <- as.numeric(predict(mods[[treatment[t]+1]][[4]], newdata=pred_data, type="response")>runif(1))
          
        #Did they have the event?
          if (drop[t]==1) {
            delta[t] <- 0
          } else {
            delta[t] <- as.numeric(predict(mods[[treatment[t]+1]][[5]], newdata=pred_data, type="response")>runif(1))
          }
        } else {
          break
        }
      }
      
      res <- data.frame(id=rep(g_data$id, max(week)), treatment=treatment, week=week, delta=delta, drop=drop)
      return(res)
  }
  
#Run the function for each individual in the MC resample
  #Resampled from all EAGeR women
  pred <- lapply(1:montecarlo, function(x){pgf(x, data=MC, maxint=K)})
  pred <- do.call(rbind, pred)
  pred$last <- as.numeric(!duplicated(pred$id, fromLast=T))
  surv_all <- survfit(Surv(week, delta) ~ 1, data=pred[pred$last==1, ])
    surv_gform_exp <- survfit(Surv(week, delta) ~ treatment, data=pred[pred$last==1, ])
  
  #Resampled from those assigned to aspirin
  pred <- lapply(1:montecarlo, function(x){pgf(x, data=MC_aspirin, maxint=K)})
  pred <- do.call(rbind, pred)
  pred$last <- as.numeric(!duplicated(pred$id, fromLast=T))
  surv_aspirin <- survfit(Surv(week, delta) ~ 1, data=pred[pred$last==1, ])
  
  #Resampled from those assigned to placebo
  pred <- lapply(1:montecarlo, function(x){pgf(x, data=MC_placebo, maxint=K)})
  pred <- do.call(rbind, pred)
  pred$last <- as.numeric(!duplicated(pred$id, fromLast=T))
  surv_placebo <- survfit(Surv(week, delta) ~ 1, data=pred[pred$last==1, ])

  
###################################################################################################################################  
# Visualize natural course
  
  thm <- theme_classic() +
    theme(
      legend.position = c(0.25,0.93),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.direction="vertical"
    )
  
  #Compare observed risk and g-comp estimated risk in all EAGeR women
  surv <- list(Observed=surv_obs, Gformula=surv_all)
  pdf("~/Documents/Pitt Projects/Natural Course/results/timevar_eager_nc.pdf", height=5, width=5)
  ggsurvplot(surv, combine=TRUE, fun="event", xlim=c(0,K), break.time.by=5, 
             ylab="Risk of hCG confirmed pregnancy", xlab="Time (weeks)",
             palette=c("dark gray", "black"),
             legend=c(0.2,0.93), font.legend=c(12, "plain", "black"),
             legend.title="", legend.labs=c("Observed", "G-computation")) 
  dev.off()
  
  #Compare among those assigned to aspirin
  surv_exposed <- data.frame(t=c(1:26),
                             observed=(1 - surv_exp$surv[27:52]), 
                             post_strat=(1 - surv_gform_exp$surv[27:52]), 
                             pre_strat=(1 - surv_aspirin$surv))
  
  cols<-c("Observed"="black", "Stratification before g-comp"="red", "Stratification after g-comp"="blue")
  pdf("~/Documents/Pitt Projects/Natural Course/results/timevar_eager_nc_aspirin.pdf", height=5, width=5)
  ggplot() + thm +
    geom_step(aes(x=surv_exposed$t, y=surv_exposed$observed, color="Observed")) +
    geom_step(aes(x=surv_exposed$t, y=surv_exposed$post_strat, color="Stratification after g-comp")) +
    geom_step(aes(x=surv_exposed$t, y=surv_exposed$pre_strat, color="Stratification before g-comp")) +
    scale_colour_manual(name="", values=cols) + 
    xlab("Time (weeks)") + scale_x_continuous(expand=c(0, 0)) +
    ylab("Risk of hCG confirmed pregnancy") + scale_y_continuous(expand=c(0,0))
  dev.off()
  
  #Compare among those assigned to placebo
  surv_unexposed <- data.frame(t=c(1:26),
                               observed=(1 - surv_exp$surv[1:26]), 
                               post_strat=(1 - surv_gform_exp$surv[1:26]), 
                               pre_strat=(1 - surv_aspirin$surv))
  
  pdf("~/Documents/Pitt Projects/Natural Course/results/timevar_eager_nc_placebo.pdf", height=5, width=5)
  ggplot() + thm +
    geom_step(aes(x=surv_unexposed$t, y=surv_unexposed$observed, color="Observed")) +
    geom_step(aes(x=surv_unexposed$t, y=surv_unexposed$post_strat, color="Stratification after g-comp")) +
    geom_step(aes(x=surv_unexposed$t, y=surv_unexposed$pre_strat, color="Stratification before g-comp")) +
    scale_colour_manual(name="", values=cols) + 
    xlab("Time (weeks)") + scale_x_continuous(expand=c(0, 0)) +
    ylab("Risk of hCG confirmed pregnancy") + scale_y_continuous(expand=c(0,0))
  dev.off()
  