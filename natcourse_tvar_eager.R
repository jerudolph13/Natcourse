
###################################################################################################################################  
#
# Project: g-formula natural course paper
#
# Purpose: Estimate natural course in EAGeR -- Is compliance associated with becoming preg within 6 months? (time-varying)
#
# Author: Jacqueline Rudolph
#
# Last Update: 16 Apr 2020
#
###################################################################################################################################  

setwd("~/Documents/Pitt Projects/Natural Course/data")

packages <- c("survival", "tidyverse", "flexsurv", "survminer", "splines")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

#What size Monte Carlo resample?
montecarlo <- 5000

#When do we admin censor?
K <- 26


###################################################################################################################################  
# Read in and manipulate the data

eager <- read.table(file="eager_tvar.txt", sep="\t", header=TRUE) %>% 
  mutate(
    age = scale(age),
    BMI = scale(BMI),
    compliance_last = if_else(week == 1, 0, as.double(lag(compliance))),
    bleed_last = if_else(week == 1, 0, as.double(lag(bleed))),
    last = as.numeric(!duplicated(id, fromLast=T))
  )

  #Create a version of data with only the last record
  eager_last <- filter(eager, last == 1)

  
###################################################################################################################################  
# What was the observed survival?
  
surv_obs <- survfit(Surv(week, delta) ~ 1, data=eager_last)


###################################################################################################################################  
# Run parametric g-formula
  
#Model exposure (A), bleeding (B), censoring (C), and outcome (T)
  #Pooled logistic models
  mod.a <- glm(compliance ~ compliance_last + bleed + bleed_last + treatment + bs(age, df=3) + bs(BMI, df=3) + smoke + eligibility + 
                 as.factor(week), family=binomial(link="logit"), data=eager)
  
  mod.b <- glm(bleed ~ compliance_last + bleed_last + treatment + bs(age, df=3) + bs(BMI, df=3) + smoke + eligibility + 
                 as.factor(week), family=binomial(link="logit"), data=eager)
  
  mod.c <- glm(drop ~ compliance + compliance_last + bleed + bleed_last + treatment + bs(age, df=3) + bs(BMI, df=3) + smoke + 
                 eligibility + as.factor(week), family=binomial(link="logit"), data=eager)
  
  mod.t <- glm(delta ~ compliance + compliance_last + bleed + bleed_last + treatment + bs(age, df=3) + bs(BMI, df=3) + smoke + 
                 eligibility + as.factor(week), family=binomial(link="logit"), data=eager)
  
#Complex implementation: MC resample and reconstruct follow-up
  #In the time-varying case, this is necessary due to complexity of the system
  #Take Monte Carlo resample
  MC0 <- eager %>% 
    filter(week == 1) %>%                     #Just want the first record
    select(-c(id, week, last, delta, drop))   #Keep only baseline variables
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  #Sample with replacement
  MC<-MC0[index,]
  MC$id<-1:montecarlo
  
  #Use parametric g-formula (pgf) function to reconstruct follow-up for each person
  pgf <- function(i, data, maxint){
    #Grab record for individual i
    g_data <- data[data$id==i, ]  
    #Initialize variables
    treatment<-eligibility<-age<-BMI<-smoke<-bleed<-bleed_last<-compliance<-compliance_last<-week<-drop<-delta<- numeric()

    #Time point 1 (baseline variables take observed values)
      week[1] <- 1
      treatment[1] <- g_data$treatment
      eligibility[1] <- g_data$eligibility
      age[1] <- g_data$age
      BMI[1] <- g_data$BMI
      smoke[1] <- g_data$smoke
      bleed[1] <- g_data$bleed; bleed_last[1] <- 0
      compliance[1] <- g_data$compliance; compliance_last[1] <- 0
    
      pred_data <- data.frame(week=week[1], treatment=treatment[1], eligibility=eligibility[1], age=age[1], BMI=BMI[1], 
                              smoke=smoke[1], bleed=bleed[1], bleed_last=bleed_last[1],
                              compliance=compliance[1], compliance_last=compliance_last[1])
      
    #Did they drop out?
      drop[1] <- as.numeric(predict(mod.c, newdata=pred_data, type="response")>runif(1))
    
    #If censored, no outcome; else, see if there is an outcome
      if (drop[1]==1) {
        delta[1] <- 0
      } else {
        delta[1] <- as.numeric(predict(mod.t, newdata=pred_data, type="response")>runif(1))
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
          bleed_last[t] <- bleed[t-1]
          compliance_last[t] <- compliance[t-1]
          
          pred_data <- data.frame(week=week[t], treatment=treatment[t], eligibility=eligibility[t], age=age[t], BMI=BMI[t], 
                                  smoke=smoke[t], bleed_last=bleed_last[t], compliance_last=compliance_last[t])
          
        #Predict bleeding at t
          bleed[t] <- as.numeric(predict(mod.b, newdata=pred_data, type="response")>runif(1))
            pred_data <- cbind(pred_data, bleed=bleed[t])
          
        #Predict compliance at t
          compliance[t] <- as.numeric(predict(mod.a, newdata=pred_data, type="response")>runif(1))
            pred_data <- cbind(pred_data, compliance=compliance[t])
            
        #Did they drop out?
          drop[t] <- as.numeric(predict(mod.c, newdata=pred_data, type="response")>runif(1))
          
        #Did they have the event?
          if (drop[t]==1) {
            delta[t] <- 0
          } else {
            delta[t] <- as.numeric(predict(mod.t, newdata=pred_data, type="response")>runif(1))
          }
        } else {
          break
        }
      }
      
      res <- data.frame(id=rep(g_data$id, max(week)), week=week, delta=delta, drop=drop)
      return(res)
  }
  
  
  #Run the function for each individual in the MC resample
  pred <- lapply(1:montecarlo, function(x){pgf(x, data=MC, maxint=K)})
  pred <- do.call(rbind, pred)
  pred$last <- as.numeric(!duplicated(pred$id, fromLast=T))
  surv_gform <- survfit(Surv(week, delta) ~ 1, data=pred[pred$last==1, ])


###################################################################################################################################  
# Compare natural course survival curves to observed
  
  surv <- list(Observed=surv_obs, Gformula=surv_gform)
  pdf("~/Documents/Pitt Projects/Natural Course/results/timevar_eager_plmc.pdf", height=5, width=5)
  ggsurvplot(surv, combine=TRUE, fun="event", xlim=c(0,K), break.time.by=5, 
             ylab="Natural Course Risk for Conception", xlab="Time (weeks)",
             palette=c("dark gray", "black"),
             legend=c(0.2,0.93), font.legend=c(12, "plain", "black"),
             legend.title="", legend.labs=c("Observed", "G-computation")) 
  dev.off()

  
  