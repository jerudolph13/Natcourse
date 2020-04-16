
###################################################################################################################################  
#
# Project: g-formula natural course paper
#
# Purpose: Estiamte natural course in EAGeR -- Is BMI associated with becoming preg within 6 months?
#
# Author: Jacqueline Rudolph
#
# Last Update: 07 Apr 2020
#
###################################################################################################################################  

  #Things to work on:
    # AFT models being very wrong

setwd("/Users/jacquelinerudolph/Documents/Pitt Projects/Natural Course/data")

packages <- c("survival", "tidyverse", "flexsurv", "survminer", "splines", "ggsci")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

#What size Monte Carlo resample?
montecarlo <- 10000

#When do we admin censor?
K <- 26

###################################################################################################################################  
# Read in and manipulate the data

eager <- read.table(file="eager_base.txt", sep="\t", header=TRUE)

  #Scale continuous covariates
  eager$age <- scale(eager$age)
  eager$BPS <- scale(eager$BPS)
  eager$BPD <- scale(eager$BPD)
  
  #Make BMI exposure categorical (1=underweight, 2=healthy weight, 3=overweight, 4=obese)
  eager$bmi_cat <- ifelse(eager$BMI<18.5, 1,       
                    ifelse(eager$BMI<25, 2,        
                      ifelse(eager$BMI<30, 3, 4)))
  #Make BMI binary (0=underweight and healthy weight; 1=overweight and obese)
  #eager$bmi_bin <- ifelse(eager$BMI<25, 0, 1)
  
  #Create a version of data with only the last record
  eager$last <- as.numeric(!duplicated(eager$id, fromLast=T))
  eager_last <- eager[eager$last==1, ]

  #Make vars for splines, to help with prediction step (?)
  # age_spline <- bs(eager_last$age, df=3)
  # age1 <- age_spline[ , 1]; age2 <- age_spline[ , 2]; age3 <- age_spline[ , 3]
  # eager_last <- cbind(eager_last, age1, age2, age3)
  # 
  # bps_spline <- bs(eager_last$BPS, df=3)
  # bps1 <- bps_spline[ , 1]; bps2 <- bps_spline[ , 2]; bps3 <- bps_spline[ , 3]
  # eager_last <- cbind(eager_last, bps1, bps2, bps3)
  # 
  # bpd_spline <- bs(eager_last$BPD, df=3)
  # bpd1 <- bpd_spline[ , 1]; bpd2 <- bpd_spline[ , 2]; bpd3 <- bpd_spline[ , 3]
  # eager_last <- cbind(eager_last, bpd1, bpd2, bpd3)
  
  
  
###################################################################################################################################  
# What was the observed survival?
  
surv_obs <- survfit(Surv(week, delta) ~ 1, data=eager_last)
  surv_obs_exp <- survfit(Surv(week, delta) ~ bmi_cat, data=eager_last)
  

###################################################################################################################################  
# Run parametric g-formula
  
#Model exposure and the outcome
  #Model for exposure
  mod.a1 <- glm(bmi_bin ~ bs(age, df=3) + bs(BPS, df=3) + bs(BPD, df=3) + white +
                  high_school + married + employed + as.factor(site) + as.factor(income) + as.factor(exercise) +
                  as.factor(alcohol) + smoke + treatment + as.factor(week), family=binomial(link="logit"), data=eager)

  #Specify T and C distributions parametrically
  mod.c1 <- flexsurvreg(Surv(week, drop) ~ treatment + as.factor(bmi_cat) + bs(age, df=3) + bs(BPS, df=3) + bs(BPD, df=3) + white +
                          high_school + married + employed + as.factor(site) + as.factor(income) + as.factor(exercise) +
                          as.factor(alcohol) + smoke + treatment,
                          dist="exp", data=eager_last)
  mod.t1 <- flexsurvreg(Surv(week, delta) ~ as.factor(bmi_cat) + bs(age, df=3) + bs(BPS, df=3) + bs(BPD, df=3) + white +
                          high_school + married + employed + as.factor(site) + as.factor(income) + as.factor(exercise) +
                          as.factor(alcohol) + smoke + treatment,
                          dist="exp", data=eager_last)
  
  #Pooled logistic for T and C (fit was slightly better without interaction)
  mod.c2 <- glm(drop ~ as.factor(bmi_cat) + bs(age, df=3) + bs(BPS, df=3) + bs(BPD, df=3) + white +
                  high_school + married + employed + as.factor(site) + as.factor(income) + as.factor(exercise) +
                  as.factor(alcohol) + smoke + treatment + as.factor(week),# + I(bmi_bin*bs(week, df=3)),
                  family=binomial(link="logit"), data=eager)
  
  mod.t2 <- glm(delta ~ as.factor(bmi_cat) + bs(age, df=3) + bs(BPS, df=3) + bs(BPD, df=3) + white +
                  high_school + married + employed + as.factor(site) + as.factor(income) + as.factor(exercise) +
                  as.factor(alcohol) + smoke + treatment + as.factor(week),# + I(bmi_bin*bs(week, df=3)),
                  family=binomial(link="logit"), data=eager)
  
#Simple implementation: use model predictions given observed exposure and covariates
  #Predict survival from AFT model
  c_pred1 <- predict(mod.c1, type="response")#*gamma(summary(mod.c1)$scale+1)
  t_pred1 <- predict(mod.t1, type="response")#*gamma(summary(mod.t1)$scale+1)
  pred1 <- data.frame(id=eager_last$id, t_pred1, c_pred1) %>% 
    mutate(tmin=if_else(c_pred1 < t_pred1, c_pred1, t_pred1), 
           delta=if_else(t_pred1==tmin & tmin<=26, 1, 0),
           t=if_else(tmin>26, 26, tmin))
  surv_gform1 <- survfit(Surv(t, delta) ~ 1, data=pred1)

  #Predict survival from pooled logistic model
    #Need to predict from a data set that includes all possible weeks for each woman
    eager_expand <- eager_last %>% 
      select(-week) %>% 
      slice(rep(1:n(), each=K))
    week <- rep(1:K, dim(eager_last)[1])
    eager_expand <- cbind(eager_expand, week)
    
  c_pred2 <- as.numeric(predict(mod.c2, newdata=eager_expand, type="response") > runif(dim(eager_expand)[1]))
  t_pred2 <- as.numeric(predict(mod.t2, newdata=eager_expand, type="response") > runif(dim(eager_expand)[1]))
  pred2 <- data.frame(id=eager_expand$id, week=eager_expand$week, t_pred2, c_pred2) %>% 
    group_by(id) %>% 
    mutate(cumt=cumsum(cumsum(t_pred2)), cumc=cumsum(cumsum(c_pred2)), totcum=cumt+cumc) %>% 
    filter(totcum<=1 | (totcum==2 & t_pred2==1 & c_pred2==1))
  pred2$delta <- ifelse(pred2$t_pred2==1, 1, 0)
  pred2$last <- as.numeric(!duplicated(pred2$id, fromLast=T))
  surv_gform2 <- survfit(Surv(week, delta) ~ 1, data=pred2[pred2$last==1, ])

  
#Complex implementation: MC resample and reconstruct follow-up
  #Take Monte Carlo resample
  MC0 <- select(eager_last, -id, -week, -last, -delta, -drop)
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)
  MC<-MC0[index,]
  MC$id<-1:montecarlo
  
  #Use pgf function to reconstruct follow-up for each person
  pgf <- function(i, data, cmodel, ymodel, aft, maxint){
    g_data <- data[data$id==i, ]
    week<-drop<-delta<- numeric()
    
    #If AFT model was used, just predict time to censoring and time to outcome
    if (aft==1) {
      c_time <- predict(cmodel, newdata=g_data, type="response")
      y_time <- predict(ymodel, newdata=g_data, type="response")
      #Observed time is the minimum of the above two times
      week <- min(c_time, y_time) 
      if (week==y_time & week<=26){
        delta<-1
        drop<- 0
      } else if (week==c_time & week<=26) {
        delta <- 0
        drop <- 1
      } else {
        delta <- 0
        drop <- 0
        time <- 26
      }
      
      res <- data.frame(id=g_data$id, week=week, delta=delta, drop=drop)
      return(res)
      
    } else {
    #If pooled logistic was used, loop through time
      #Time point 1
      g_data$week <- week[1] <- 1
      drop[1] <- as.numeric(predict(cmodel, newdata=g_data, type="response")>runif(1))
    
      #If censored, no outcome; else, see if there is an outcome
      if (drop[1]==1) {
        delta[1] <- 0
      } else {
        delta[1] <- as.numeric(predict(ymodel, newdata=g_data, type="response")>runif(1))
      }
      
      #Repeat for other time points until either drop out or event occurs
      for (t in 2:maxint) {
        if (delta[t-1]==0 & drop[t-1]==0) {
          g_data$week <- week[t] <- t
          drop[t] <- as.numeric(predict(cmodel, newdata=g_data, type="response")>runif(1))
          
          if (drop[t]==1) {
            delta[t] <- 0
          } else {
            delta[t] <- as.numeric(predict(ymodel, newdata=g_data, type="response")>runif(1))
          }
        } else {
          break
        }
      }
      
      res <- data.frame(id=rep(g_data$id, max(week)), week=week, delta=delta, drop=drop)
      return(res)
      
    }
  }
  
  #MC-implementation with AFT models
  pred3 <- lapply(1:montecarlo, function(x){pgf(x, data=MC, cmodel=mod.c1, ymodel=mod.t1, aft=1, maxint=K)})
  pred3 <- do.call(rbind, pred3)
  surv_gform3 <- survfit(Surv(week, delta) ~ 1, data=pred3)
  
  #MC-implementation with pooled logistic models
  pred4 <- lapply(1:montecarlo, function(x){pgf(x, data=MC, cmodel=mod.c2, ymodel=mod.t2, aft=0, maxint=K)})
  pred4 <- do.call(rbind, pred4)
  pred4$last <- as.numeric(!duplicated(pred4$id, fromLast=T))
  surv_gform4 <- survfit(Surv(week, delta) ~ 1, data=pred4[pred4$last==1, ])


###################################################################################################################################  
# Compare natural course survival curves to observed
  
  surv <- list(Observed=surv_obs, Gformula=surv_gform2)
  ggsurvplot(surv, combine=TRUE, xlim=c(0,K), ylim=c(0.25, 1.0), break.time.by=5, palette="JAMA",
             legend=c(0.70,0.93), font.legend=c(12, "plain", "black"),
             legend.title="", legend.labs=c("Observed natural course", "G-computation natural course")) 
  

###################################################################################################################################  
# Pooled logistic fit well, so let's use that approach to estimate counterfactual outcomes
  
  #Set exposure to 2 (everyone is healthy weight)
  eager2 <- eager_expand
  eager2$bmi_cat <- 2
    c2 <- as.numeric(predict(mod.c2, newdata=eager2, type="response") > runif(dim(eager2)[1]))
    t2 <- as.numeric(predict(mod.t2, newdata=eager2, type="response") > runif(dim(eager2)[1]))
    pred2 <- data.frame(id=eager2$id, week=eager2$week, t2, c2) %>% 
      group_by(id) %>% 
      mutate(cumt=cumsum(cumsum(t2)), cumc=cumsum(cumsum(c2)), totcum=cumt+cumc) %>% 
      filter(totcum<=1 | (totcum==2 & t2==1 & c2==1))
    pred2$delta <- ifelse(pred2$t2==1, 1, 0)
    pred2$last <- as.numeric(!duplicated(pred2$id, fromLast=T))
    surv_set2 <- survfit(Surv(week, delta) ~ 1, data=pred2[pred2$last==1, ])
    
    surv <- list(Observed=surv_obs, Gformula=surv_gform2, intervention=surv_set2)
    ggsurvplot(surv, combine=TRUE, xlim=c(0,K), ylim=c(0, 1.0), break.time.by=5, palette="JAMA",
               legend=c(0.70,0.93), font.legend=c(12, "plain", "black"),
               legend.title="", legend.labs=c("Observed natural course", "G-formula natural course", 
                                                                   "All healthy weight")) 
 

  
  
