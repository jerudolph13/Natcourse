
###############################################################################################
#
# Project: Natural course paper
#
# Purpose: Estimate natural course in EAGeR
#
# Author: Jacqueline Rudolph
#
# Last Update: 24 May 2021
#
###############################################################################################


packages <- c("survival", "tidyverse", "survminer", "splines")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

source("./plot_thm.R")

# What size Monte Carlo resample?
montecarlo <- 5000

# When do we admin censor?
K <- 26

# Set seed for predictions
set.seed(123)


# Read in data ------------------------------------------------------------

eager <- read.table(file="../data/eager_tvar.txt", sep="\t", header=TRUE) %>% 
  mutate(last = as.numeric(!duplicated(id, fromLast=T)),
         age = scale(age),
         BMI = scale(BMI))

# Create a version of data with only the last record
  eager_last <- filter(eager, last == 1)


# Observed survival -------------------------------------------------------
  
surv_obs <- survfit(Surv(week, conception) ~ 1, data=eager_last)
  # Survival by assigned treatment
  surv_exp <- survfit(Surv(week, conception) ~ treatment, data=eager_last)
drop_obs <- survfit(Surv(week, drop) ~ 1, data=eager_last)
var_obs <- eager %>% 
  group_by(week) %>% 
  summarize(avg_comp = mean(compliance),
            avg_nausea = mean(nausea),
            avg_bleed = mean(bleed),
            avg_drop = mean(drop))

# Run parametric g-formula ------------------------------------------------
  
# Model nausea (N), bleeding (B), exposure (A), censoring (C), and outcome (T)
  # Run all models stratified by treatment arm
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
  
  mod.t <- glm(conception ~ compliance + compliance_last + bleed + bleed_last + nausea + nausea_last +
                 eligibility + bs(age, df=3) + bs(BMI, df=3) + smoke + as.factor(week), 
                 family=binomial(link="logit"), data=eager, subset=(treatment==k))
  
  return(list(mod.n, mod.b, mod.a, mod.c, mod.t))
}
  
# Store model objects in a list we will use in the pgf function below
mods <- lapply(0:1, function(x){pgf.mod(x)})
  
# Monte Carlo (MC) implementation
  # Resample from all women
  MC0 <- eager %>% 
    filter(week == 1) %>%  #Just want the first record
    select(-c(id, week, last, conception, drop))  # Keep only baseline variables
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  # Sample with replacement
  MC<-MC0[index,]  # Grab the sampled records
  MC$id<-1:montecarlo  # Assign new IDs
    
  # Resample from those assigned to aspirin
  aspirin <- filter(eager, treatment==1)
  MC0 <- aspirin %>% 
    filter(week == 1) %>%  
    select(-c(id, week, last, conception, drop))  
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  
  MC_aspirin <- MC0[index,]  
  MC_aspirin$id <- 1:montecarlo 
    
  # Resample from those assigned to placebo
  placebo <- filter(eager, treatment==0)
  MC0 <- placebo %>% 
    filter(week == 1) %>%  
    select(-c(id, week, last, conception, drop))  
  index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)  
  MC_placebo <- MC0[index,]  
  MC_placebo$id <- 1:montecarlo 
  
# Use parametric g-formula (pgf) function to reconstruct follow-up for each person
  pgf <- function(i, data, maxint){
    # Grab record for individual i
    g_data <- data[data$id==i, ] 
    
    # Initialize variables
    treatment<-eligibility<-age<-BMI<-smoke<-nausea<-nausea_last<-bleed<-bleed_last <- numeric()
    compliance<-compliance_last<-week<-drop<-conception <- numeric()

    # Time point 1 (baseline variables take observed values)
    week[1] <- 1
    treatment[1] <- g_data$treatment
    eligibility[1] <- g_data$eligibility
    age[1] <- g_data$age
    BMI[1] <- g_data$BMI
    smoke[1] <- g_data$smoke
    # For t-varying variables, assign last value to be 0 in the first record
    nausea[1] <- g_data$nausea; nausea_last[1] <- 0
    bleed[1] <- g_data$bleed; bleed_last[1] <- 0
    compliance[1] <- g_data$compliance; compliance_last[1] <- 0
    
    pred_data <- data.frame(week=week[1], treatment=treatment[1], eligibility=eligibility[1], age=age[1], BMI=BMI[1],
                            smoke=smoke[1], nausea=nausea[1], nausea_last=nausea_last[1], bleed=bleed[1], bleed_last=bleed_last[1],
                            compliance=compliance[1], compliance_last=compliance_last[1])
      
    # Did they drop out?
      # Note: need to select correct model object from the list "mods", 
      # First based on assigned treatment and then on which variable we are predicting
    drop[1] <- as.numeric(predict(mods[[treatment[1]+1]][[4]], newdata=pred_data, type="response")>runif(1))
    
    # If censored, no outcome; else, see if there is an outcome
    if (drop[1]==1) {
      conception[1] <- 0
    } else {
      conception[1] <- as.numeric(predict(mods[[treatment[1]+1]][[5]], newdata=pred_data, type="response")>runif(1))
    }
      
    # Repeat for other time points until either drop out or event occurs
    for (t in 2:maxint) {
        if (conception[t-1]==0 & drop[t-1]==0) {
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
        
        # Predict nausea at t
          nausea[t] <- as.numeric(predict(mods[[treatment[t]+1]][[1]], newdata=pred_data, type="response")>runif(1))
          pred_data <- cbind(pred_data, nausea=nausea[t])
          
        # Predict bleeding at t
          bleed[t] <- as.numeric(predict(mods[[treatment[t]+1]][[2]], newdata=pred_data, type="response")>runif(1))
            pred_data <- cbind(pred_data, bleed=bleed[t])
          
        # Predict compliance at t
          compliance[t] <- as.numeric(predict(mods[[treatment[t]+1]][[3]], newdata=pred_data, type="response")>runif(1))
            pred_data <- cbind(pred_data, compliance=compliance[t])
            
        # Did they drop out?
          drop[t] <- as.numeric(predict(mods[[treatment[t]+1]][[4]], newdata=pred_data, type="response")>runif(1))
          
        # Did they have the event?
          if (drop[t]==1) {
            conception[t] <- 0
          } else {
            conception[t] <- as.numeric(predict(mods[[treatment[t]+1]][[5]], newdata=pred_data, type="response")>runif(1))
          }
        } else {
          break
        }
      }
      
      res <- data.frame(id=rep(g_data$id, max(week)), treatment=treatment, week=week, conception=conception, drop=drop,
                        compliance=compliance, nausea=nausea, bleed=bleed)
      return(res)
  }
  
# Run the function for each individual in the MC resample
# Resampled from all EAGeR women
pred <- lapply(1:montecarlo, function(x){pgf(x, data=MC, maxint=K)})
pred <- do.call(rbind, pred)
pred$last <- as.numeric(!duplicated(pred$id, fromLast=T))
surv_all <- survfit(Surv(week, conception) ~ 1, data=pred[pred$last==1, ])
  surv_gform_exp <- survfit(Surv(week, conception) ~ treatment, data=pred[pred$last==1, ])
drop_all <- survfit(Surv(week, drop) ~ 1, data=pred[pred$last==1, ])
var_all <- pred %>% 
  group_by(week) %>% 
  summarize(avg_comp = mean(compliance),
            avg_nausea = mean(nausea),
            avg_bleed = mean(bleed),
            avg_drop = mean(drop))
  
# Resampled from those assigned to aspirin
pred <- lapply(1:montecarlo, function(x){pgf(x, data=MC_aspirin, maxint=K)})
pred <- do.call(rbind, pred)
pred$last <- as.numeric(!duplicated(pred$id, fromLast=T))
surv_aspirin <- survfit(Surv(week, conception) ~ 1, data=pred[pred$last==1, ])
drop_aspirin <- survfit(Surv(week, drop) ~ 1, data=pred[pred$last==1, ])
var_aspirin <- pred %>% 
  group_by(week) %>% 
  summarize(avg_comp = mean(compliance),
            avg_nausea = mean(nausea),
            avg_bleed = mean(bleed),
            avg_drop = mean(drop))

# Resampled from those assigned to placebo
pred <- lapply(1:montecarlo, function(x){pgf(x, data=MC_placebo, maxint=K)})
pred <- do.call(rbind, pred)
pred$last <- as.numeric(!duplicated(pred$id, fromLast=T))
surv_placebo <- survfit(Surv(week, conception) ~ 1, data=pred[pred$last==1, ])
drop_placebo <- survfit(Surv(week, drop) ~ 1, data=pred[pred$last==1, ])
var_placebo <- pred %>% 
  group_by(week) %>% 
  summarize(avg_comp = mean(compliance),
            avg_nausea = mean(nausea),
            avg_bleed = mean(bleed),
            avg_drop = mean(drop))  

# Visualize natural course ------------------------------------------------
  
cols<-c("Observed"="dark gray", "G-computation"="black")
  
# Compare observed risk and g-comp estimated risk in all EAGeR women
surv_fig <- data.frame(t=c(1:26),
                       observed = 1 - surv_obs$surv,
                       gform = 1 - surv_all$surv)

pdf("../figures/timevar_nc.pdf", height=5, width=5)
ggplot(surv_fig) + thm +
  labs(x="\nTime (weeks)", y="Cumulative incidence of pregnancy\n") +
  geom_step(aes(x=t, y=observed, color="Observed"), size=0.75) +
  geom_step(aes(x=t, y=gform, color="G-computation"), size=0.75) +
  scale_colour_manual(name="", values=cols) + 
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0.0, 0.8))
dev.off()
  
# Compare among those assigned to aspirin
surv_exposed <- data.frame(t=c(1:26),
                           observed=(1 - surv_exp$surv[27:52]), 
                           pre_strat=(1 - surv_aspirin$surv))
  
pdf("../figures/timevar_nc_aspirin.pdf", height=5, width=5)
ggplot(surv_exposed) + thm +
  labs(x="\nTime (weeks)", y="Cumulative incidence of pregnancy\n") +
  geom_step(aes(x=t, y=observed, color="Observed"), size=0.75) +
  geom_step(aes(x=t, y=pre_strat, color="G-computation"), size=0.75) +
  scale_colour_manual(name="", values=cols) + 
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0.0, 0.8))
dev.off()
  
# Compare among those assigned to placebo
surv_unexposed <- data.frame(t=c(1:26),
                             observed=(1 - surv_exp$surv[1:26]), 
                             pre_strat=(1 - surv_placebo$surv))
  
pdf("../figures/timevar_nc_placebo.pdf", height=5, width=5)
ggplot(surv_unexposed) + thm +
  labs(x="\nTime (weeks)", y="Cumulative incidence of pregnancy\n") +
  geom_step(aes(x=t, y=observed, color="Observed"), size=0.75) +
  geom_step(aes(x=t, y=pre_strat, color="G-computation"), size=0.75) +
  scale_colour_manual(name="", values=cols) + 
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0,0), limits=c(0.0, 0.8))
dev.off()

# Compare covariate distributions
dists <- bind_rows("Observed" = var_obs, "Predicted" = var_all, .id="Type")
write.table(dists, file="../results/dists_by_week.txt", sep="\t", row.names=FALSE)

pdf("../figures/comp_by_week.pdf", height=5, width=6)
ggplot(dists) + thm2 +
  labs(x="\nTime (weeks)", y="Proportion compliant\n") +
  geom_line(aes(x=week, y=avg_comp, color=Type), size=0.75) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0.3, 1.0))
dev.off()

pdf("../figures/drop_by_week.pdf", height=5, width=6)
ggplot(dists) + thm2 +
  labs(x="\nTime (weeks)", y="Proportion censored\n") +
  geom_line(aes(x=week, y=avg_drop, color=Type), size=0.75) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0.0, 0.025))
dev.off()

pdf("../figures/nausea_by_week.pdf", height=5, width=6)
ggplot(dists) + thm2 +
  labs(x="\nTime (weeks)", y="Proportion reporting nausea\n") +
  geom_line(aes(x=week, y=avg_nausea, color=Type), size=0.75) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0.2, 0.5))
dev.off()

pdf("../figures/bleed_by_week.pdf", height=5, width=6)
ggplot(dists) + thm2 +
  labs(x="\nTime (weeks)", y="Proportion reporting bleeding\n") +
  geom_line(aes(x=week, y=avg_bleed, color=Type), size=0.75) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0))#, limits=c(0.2, 0.5))
dev.off()
