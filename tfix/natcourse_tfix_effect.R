
###################################################################################################################################  
#
# Project: g-formula natural course paper
#
# Purpose: Obtain confidence intervals for EAGeR example
#
# Author: Jacqueline Rudolph
#
# Last Update: 09 Apr 2020
#
###################################################################################################################################  


setwd("~/Documents/Pitt Projects/Natural Course/data")

packages <- c("survival", "tidyverse", "flexsurv", "survminer", "splines", "ggsci")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

#How many bootstrap resamples?
nboot <- 200

#What size Monte Carlo resample?
montecarlo <- 5000

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
  
  
###################################################################################################################################  
# Start bootstrap loop
  
bootrep <- function(r) {
  
  set.seed(123 + r)
  
  firstobs <- eager[eager$week==1, ]
  samp <- table(firstobs[sample(1:nrow(firstobs),nrow(firstobs),replace=T), (names(eager) == "id")])
  
  boot <- NULL
  if (r==0) {
    boot <- eager %>% 
      rename(bid=id)
  } else {
    for(zzz in 1:max(samp)){ 
      cc <- eager[eager$id %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid<-paste0(cc$id,zzz)
      boot <- rbind(boot, cc)
    }
  }
  
  #Create a version of data with only the last record
  boot$last <- as.numeric(!duplicated(boot$bid, fromLast=T))
  boot_last <- boot[boot$last==1, ]
  
  
###################################################################################################################################  
# Estimate survival
  
#Observed survival
surv_obs <- survfit(Surv(week, delta) ~ 1, data=boot_last)

#Model censoring (C) and the outcome (T)
  #Pooled logistic for T and C
  mod.c2 <- glm(drop ~ as.factor(bmi_cat) + bs(age, df=3) + bs(BPS, df=3) + bs(BPD, df=3) + white +
                  high_school + married + employed + as.factor(site) + as.factor(income) + as.factor(exercise) +
                  as.factor(alcohol) + smoke + treatment + as.factor(week),
                  family=binomial(link="logit"), data=boot)
  
  mod.t2 <- glm(delta ~ as.factor(bmi_cat) + bs(age, df=3) + bs(BPS, df=3) + bs(BPD, df=3) + white +
                  high_school + married + employed + as.factor(site) + as.factor(income) + as.factor(exercise) +
                  as.factor(alcohol) + smoke + treatment + as.factor(week),
                  family=binomial(link="logit"), data=boot)

#MC-implementation
  #Take MC resample
  MC0 <- select(boot_last, -bid, -week, -last, -delta, -drop)  #Keep only exposure and baseline confounders
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)       #Sample with replacement
  MC<-MC0[index, ]
  MC$id<-1:montecarlo
  
  #Set exposure to bmi_cat==2
  MC$bmi_cat <- 2
  
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
  
  pgf.res <- lapply(1:montecarlo, function(x) {pgf(x, data=MC, cmodel=mod.c2, ymodel=mod.t2, aft=0, maxint=K)})
  pgf.res <- do.call(rbind, pgf.res)
  pgf.res$last <- as.numeric(!duplicated(pgf.res$id, fromLast=T))
  surv_set2 <- survfit(Surv(week, delta) ~ 1, data=pgf.res[pgf.res$last==1, ])
  
  if (r==0){
    surv <- list(Observed=surv_obs, Gformula=surv_set2)
    ggsurvplot(surv, combine=TRUE, fun="event", xlim=c(0,K), break.time.by=5, 
               xlab="Time (weeks)", ylab="Cumulative probability of conceiving",
               palette=c("dark gray", "black"), legend=c(0.3,0.93), font.legend=c(12, "plain", "black"),
               legend.title="", legend.labs=c("Observed natural course", "All healthy weight")) 
  }
  
  
###################################################################################################################################  
# Compare risks

r_obs <- 1 - summary(surv_obs)$surv
r_int <- 1 - summary(surv_set2)$surv
  rd <- r_int - r_obs
  week <- c(1:length(rd))
  boot <- rep(r, length(rd))

res <- data.frame(boot=boot, r_obs=r_obs, r_int=r_int, rd=rd, week=week)

return(res)

}

  
###################################################################################################################################  
# Summarize across bootstraps

all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})
  all.boot <- do.call(rbind, all.boot)

res <- all.boot[all.boot$boot==0, ]
summ <- all.boot %>% 
  group_by(week) %>% 
  summarise(sd=sd(rd))

res$lower <- res$rd - 1.96*summ$sd
res$upper <- res$rd + 1.96*summ$sd

#Visualize the risk difference (with 95% CI)
thm <- theme_classic() +
  theme(legend.position = "none")

pdf("~/Documents/Pitt Projects/Natural Course/results/timefix_eager_rd.pdf", height=5, width=5)
ggplot(res) + thm +
  geom_ribbon(aes(x=week, ymin=lower, ymax=upper, alpha=0.05)) + 
  geom_line(aes(x=week, y=rd))+
  geom_hline(yintercept=0, linetype="dashed") +
  xlab("Week") + scale_x_continuous(expand=c(0, 0), limits=c(1, 26)) +
  ylab("Risk Difference") + scale_y_continuous(expand=c(0,0), limits=c(-0.02, 0.06))
dev.off()

#Output results
write.table(res, file="~/Documents/Pitt Projects/Natural Course/results/eager_timefix_int.txt", sep="\t", row.names=FALSE)
#res <- read.table(file="~/Documents/Pitt Projects/Natural Course/results/eager_timefix_int.txt", sep="\t", header=TRUE)

