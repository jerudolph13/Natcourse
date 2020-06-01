
###################################################################################################################################  
#
# Project: g-formula natural course paper
#
# Purpose: Obtain point estimate and confidence intervals for time-varying EAGeR example
#
# Author: Jacqueline Rudolph
#
# Last Update: 19 May 2020
#
###################################################################################################################################  


setwd("~/Documents/Pitt Projects/Natural Course/results")

packages <- c("survival", "tidyverse", "survminer", "splines", "parallel", "doParallel")
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
# Read in data

eager <- read.table(file="~/Documents/Pitt Projects/Natural Course/data/eager_tvar.txt", sep="\t", header=TRUE) %>% 
  mutate(age = scale(age),
         BMI = scale(BMI))

#If doing the intervention only on those who were assigned to aspirin
eager <- eager %>% filter(treatment==1)
  
  
###################################################################################################################################  
# Start bootstrap loop
  
bootrep <- function(r) {
  
  set.seed(123 + r)
  
#Randomly sample with replacement women from EAGeR
  firstobs <- eager[eager$week==1, ] #Subset to the first week, so 1 record per woman
  samp <- table(firstobs[sample(1:nrow(firstobs),nrow(firstobs),replace=T), (names(eager) == "id")]) #Sample
  
#Grab the rest of their records past the first week
  #If r==0, grab the original sample, so we can compute point estimate
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
# What was the observed survival?
  
surv_obs <- survfit(Surv(week, delta) ~ 1, data=boot_last)


###################################################################################################################################  
# Run parametric g-formula
  
  mod.n <- glm(nausea ~ compliance_last + bleed_last + nausea_last + 
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=boot)
    
  mod.b <- glm(bleed ~ compliance_last + bleed_last + nausea + nausea_last + 
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=boot)
    
  mod.a <- glm(compliance ~ compliance_last + bleed + bleed_last + nausea + nausea_last +
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=boot)
    
  mod.c <- glm(drop ~ compliance + compliance_last + bleed + bleed_last + nausea + nausea_last +
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=boot)
    
  mod.t <- glm(delta ~ compliance + compliance_last + bleed + bleed_last + nausea + nausea_last +
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=boot)
    
#Store model objects in a list we will use in the pgf function below
mods <- list(mod.n, mod.b, mod.a, mod.c, mod.t)
  
#Take Monte Carlo resample
MC0 <- boot %>% 
  filter(week == 1) %>%  #Just want the first record
  select(-c(bid, week, last, delta, drop))   #Keep only baseline variables
index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  #Sample with replacement
MC<-MC0[index,]  #Grab the sampled records
MC$id<-1:montecarlo  #Assign new IDs
  
#Use parametric g-formula (pgf) function to reconstruct follow-up for each person
pgf <- function(i, data, maxint, exposure){
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
  
  if (is.null(exposure)) {
    compliance[1] <- g_data$compliance
  } else {
    compliance[1] <- exposure
  }
  compliance_last[1] <- 0
  
  pred_data <- data.frame(week=week[1], treatment=treatment[1], eligibility=eligibility[1], age=age[1], BMI=BMI[1],
                          smoke=smoke[1], nausea=nausea[1], nausea_last=nausea_last[1], bleed=bleed[1], bleed_last=bleed_last[1],
                          compliance=compliance[1], compliance_last=compliance_last[1])
  
  #Did they drop out?
  #Note: need to select correct model object from the list "mods", 
  #first based on assigned treatment and then on which variable we are predicting
  drop[1] <- as.numeric(predict(mods[[4]], newdata=pred_data, type="response")>runif(1))
  
  #If censored, no outcome; else, see if there is an outcome
  if (drop[1]==1) {
    delta[1] <- 0
  } else {
    delta[1] <- as.numeric(predict(mods[[5]], newdata=pred_data, type="response")>runif(1))
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
      nausea[t] <- as.numeric(predict(mods[[1]], newdata=pred_data, type="response")>runif(1))
      pred_data <- cbind(pred_data, nausea=nausea[t])
      
      #Predict bleeding at t
      bleed[t] <- as.numeric(predict(mods[[2]], newdata=pred_data, type="response")>runif(1))
      pred_data <- cbind(pred_data, bleed=bleed[t])
      
      #Predict or set compliance at t
      if (is.null(exposure)) {
        compliance[t] <- as.numeric(predict(mods[[3]], newdata=pred_data, type="response")>runif(1))
      } else {
        compliance[t] <- exposure
      }
      pred_data <- cbind(pred_data, compliance=compliance[t])
      
      #Did they drop out?
      drop[t] <- as.numeric(predict(mods[[4]], newdata=pred_data, type="response")>runif(1))
      
      #Did they have the event?
      if (drop[t]==1) {
        delta[t] <- 0
      } else {
        delta[t] <- as.numeric(predict(mods[[5]], newdata=pred_data, type="response")>runif(1))
      }
    } else {
      break
    }
  }
  
  res <- data.frame(id=rep(g_data$id, max(week)), treatment=treatment, week=week, delta=delta, drop=drop)
  return(res)
}
  
  #Estimate natural course using pgf
  pgf.nc <- lapply(1:montecarlo, function(x) {pgf(x, data=MC, maxint=K, exposure=NULL)})
    pgf.nc <- do.call(rbind, pgf.nc)
    pgf.nc$last <- as.numeric(!duplicated(pgf.nc$id, fromLast=T))
    #Estimate survival under the natural course
    surv_nc <- survfit(Surv(week, delta) ~ 1, data=pgf.nc[pgf.nc$last==1, ])
    
  #Run pgf, setting compliance to be 1
  pgf.int <- lapply(1:montecarlo, function(x) {pgf(x, data=MC, maxint=K, exposure=1)})
    pgf.int <- do.call(rbind, pgf.int)
    pgf.int$last <- as.numeric(!duplicated(pgf.int$id, fromLast=T))
    #Estimate survival under the intervention
    surv_set1 <- survfit(Surv(week, delta) ~ 1, data=pgf.int[pgf.int$last==1, ])
  
  #Run pgf, setting compliance to be 0
    pgf.int <- lapply(1:montecarlo, function(x) {pgf(x, data=MC, maxint=K, exposure=0)})
    pgf.int <- do.call(rbind, pgf.int)
    pgf.int$last <- as.numeric(!duplicated(pgf.int$id, fromLast=T))
    #Estimate survival under the intervention
    surv_set0 <- survfit(Surv(week, delta) ~ 1, data=pgf.int[pgf.int$last==1, ])
    
  # if (r==0){
  #   surv <- list(Observed=surv_obs, Gformula=surv_nc)
  #   pdf("~/Documents/Pitt Projects/Natural Course/results/timevar_eager_nc2.pdf", height=5, width=5)
  #   ggsurvplot(surv, combine=TRUE, fun="event", xlim=c(1,K), break.time.by=5, 
  #              xlab="Time (weeks)", ylab="Risk of hCG confirmed pregnancy",
  #              palette=c("dark gray", "black"), legend=c(0.2,0.93), font.legend=c(12, "plain", "black"),
  #              legend.title="", legend.labs=c("Observed", "G-computation")) 
  #   dev.off()
  #   
  #   surv <- list(Observed=surv_obs, Gformula=surv_set1)
  #   pdf("~/Documents/Pitt Projects/Natural Course/results/timevar_eager_set.pdf", height=5, width=5)
  #   ggsurvplot(surv, combine=TRUE, fun="event", xlim=c(1,K), break.time.by=5, 
  #              xlab="Time (weeks)", ylab="Risk of hCG confirmed pregnancy",
  #              palette=c("dark gray", "black"), legend=c(0.3,0.93), font.legend=c(12, "plain", "black"),
  #              legend.title="", legend.labs=c("Observed natural course", "Always compliant")) 
  #   dev.off()
  #   
  #   surv <- list(Observed=surv_set0, Gformula=surv_set1)
  #   pdf("~/Documents/Pitt Projects/Natural Course/results/timevar_ate_risk.pdf", height=5, width=5)
  #   ggsurvplot(surv, combine=TRUE, fun="event", xlim=c(1,K), break.time.by=5,
  #             xlab="Time (weeks)", ylab="Risk of hCG confirmed pregnancy",
  #             palette=c("dark gray", "black"), legend=c(0.3,0.93), font.legend=c(12, "plain", "black"),
  #             legend.title="", legend.labs=c("Never compliant", "Always compliant"))
  #   dev.off()
  # }
  
  
###################################################################################################################################  
# Compare risks

r_obs <- 1 - summary(surv_obs)$surv
  #In some bootstraps, these vectors were of different length, so carry forward last value if a shorter vector
  if (length(r_obs) < K) {
    r_obs[length(r_obs):K] <- r_obs[length(r_obs)]
  }
r_nc <- 1 - summary(surv_nc)$surv
  if (length(r_nc) < K) {
    r_nc[length(r_nc):K] <- r_nc[length(r_nc)]
  }
r_int1 <- 1 - summary(surv_set1)$surv
  if (length(r_int1) < K) {
    r_int1[length(r_int1):K] <- r_int1[length(r_int1)]
  }
r_int0 <- 1 - summary(surv_set0)$surv
  if (length(r_int0) < K) {
    r_int0[length(r_int0):K] <- r_int0[length(r_int0)]
  }

rd_nc <- r_int1 - r_obs
rd_ate <- r_int1 - r_int0
week <- c(1:K)
boot <- rep(r, K)

res <- data.frame(boot=boot, r_obs=r_obs, r_nc=r_nc, r_int1=r_int1, r_int0=r_int0, rd_nc=rd_nc, rd_ate=rd_ate, week=week)

return(res)

}

  
###################################################################################################################################  
# Summarize across bootstraps

#Run if only interested in the point estimate
#boot.0 <- bootrep(r=0)

#Run all bootstraps to get point estimate and 95% confidence intervals
cores <- detectCores()
all.boot <- mclapply(0:nboot, function(tt) {bootrep(tt)}, mc.cores=cores, mc.set.seed=FALSE)
#all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})
  all.boot <- do.call(rbind, all.boot)

res <- all.boot[all.boot$boot==0, ]
summ <- all.boot %>% 
  group_by(week) %>% 
  summarise(se_nc=sd(rd_nc), se_ate=sd(rd_ate))

res$lower_nc <- res$rd_nc - 1.96*summ$se_nc
res$upper_nc <- res$rd_nc + 1.96*summ$se_nc

res$lower_ate <- res$rd_ate - 1.96*summ$se_ate
res$upper_ate <- res$rd_ate + 1.96*summ$se_ate

#Output results
write.table(res, file="eager_timevar_results.txt", sep="\t", row.names=FALSE)
#res <- read.table(file="eager_timevar_results.txt", sep="\t", header=TRUE)


###################################################################################################################################  
# Visualize the risk difference (with 95% CI)

thm <- theme_classic() +
  theme(
    legend.position = "none",
    axis.text=element_text(size=12),
    axis.title=element_text(size=14)
  )

pdf("timevar_nc-effect_rd.pdf", height=5, width=5)
ggplot(res) + thm +
  geom_ribbon(aes(x=week, ymin=lower_nc, ymax=upper_nc, alpha=0.05)) +
  geom_line(aes(x=week, y=rd_nc)) +
  geom_hline(yintercept=0, linetype="dashed") +
  xlab("Time (weeks)") + scale_x_continuous(expand=c(0, 0), limits=c(1, 26)) +
  ylab("Risk difference") + scale_y_continuous(expand=c(0,0), limits=c(-0.03, 0.05), breaks=c(-0.02, 0.0, 0.02, 0.04))
dev.off()

pdf("timevar_ate_rd.pdf", height=5, width=5)
ggplot(res) + thm +
  geom_ribbon(aes(x=week, ymin=lower_ate, ymax=upper_ate, alpha=0.05)) +
  geom_line(aes(x=week, y=rd_ate)) +
  geom_hline(yintercept=0, linetype="dashed") +
  xlab("Time (weeks)") + scale_x_continuous(expand=c(0, 0), limits=c(1, 26)) +
  ylab("Risk difference") + scale_y_continuous(expand=c(0,0), limits=c(-0.07, 0.27), 
                                               breaks=c(-0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25))
dev.off()

