
###################################################################################################
#
# Purpose: For g-formula natural course paper, time-fixed case 1: binary exposure,
#             3 binary confounders, exponential outcome, no LTFU
#
# Author: Jacqueline Rudolph
#
# Last Update: 16 apr 2020
#
##################################################################################################

# Things explored:
#     - Including more variables than are necessary, esp. interaction terms
#     - Leaving variables out

# Could be explored
#     - Changing distribution of AFT model (but how to properly predict T?)

packages <- c("survival", "tidyverse", "flexsurv", "survminer")

for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}


###################################################################################################################################  
# Simulate data

set.seed(123)
n <- 5000       #Sample size
K <- 1          #Maximum time

#Empty data frame to hold data
dat <- data.frame(X=rep(0, n),
                  Z1=rep(0, n),
                  Z2=rep(0, n),
                  Z3=rep(0, n),
                  Tv=rep(0, n),
                  Y=rep(0, n))

for (i in 1:n){
  #Binary confounders
  dat$Z1[i] <- rbinom(1, 1, 0.4)  #P(Z1)=0.4
  dat$Z2[i] <- rbinom(1, 1, 0.2)  #P(Z2)=0.2
  dat$Z3[i] <- rbinom(1, 1, 0.6)  #P(Z3)=0.6
  
  #Exposure, P(X)=0.5
  p_x <- 1/(1 + exp(-(-log(1/0.5 - 1) 
                      - log(2)*dat$Z1[i] + log(2)*0.4      #Effect for Z1: RR=2
                      + log(1.5)*dat$Z2[i] - log(1.5)*0.2  #Effect for Z2: RR=1.5  
                      - log(3)*dat$Z3[i] + log(3)*0.6)))   #Effect for Z3: RR=3
  dat$X[i] <- rbinom(1, 1, p_x)
  
  #Outcome, exponential with lambda=0.01
  p_t <- log(2)*dat$X[i] - log(2)*dat$Z1[i] + log(1.5)*dat$Z2[i] - log(3)*dat$Z3[i]
  dat$Tv[i] <- rexp(1, exp(0.01)*exp(p_t))
    dat$Y[i] <- ifelse(dat$Tv[i] > K, 0, 1)
    dat$Tv[i] <- ifelse(dat$Tv[i] > K, K, dat$Tv[i])
}


###################################################################################################################################  
# Observed survival

sum(dat$Y)/n #P(Y=1) = 0.4768
surv_obs <- survfit(Surv(Tv, Y) ~ 1, data=dat)


###################################################################################################################################  
# Run parametric g-formula

#Model exposure and the outcome
  #Exposure
  mod.x <- glm(X ~ Z1 + Z2 + Z3, family=binomial(link="logit"), data=dat)
  
  #Outcome: exponential AFT
    #Why AFT? Time is continuous, and we don't want to discretize.
    #Drawback: in real data, assuming time-to-event follows specified 
    #          parametric distribution may be too restrictive
  mod.t <- flexsurvreg(Surv(Tv, Y) ~ X + Z1 + Z2 + Z3, dist="exp", data=dat)

#Simple  g-formula implementation: just predict from outcome model
  #Predict from AFT model
  t_pred1 <- rep(0, n)
  for (i in 1:n) {
    desX <- dat[i , c("X", "Z1", "Z2", "Z3")]
    t_pred1[i]<-rexp(1, exp(coef(mod.t)[names(coef(mod.t))=="rate"])*
                    exp(coef(mod.t)[!names(coef(mod.t))=="rate"]%*%t(desX)))
  }
  t_pred1 <- ifelse(t_pred1>K, K, t_pred1)
  delta1 <- ifelse(t_pred1<K, 1, 0)

#Complex implementation: Monte Carlo (MC) resample and build follow-up
  #In this simple simulation, MC implementation is likely unnecessary because 
  #data and causal model are simple (time-fixed exposure, few confounders). 
  #However, let's show what it looks like here, so we can build to more complex scenarios.
  
  #Resample with replacement
  montecarlo <- 10000  #MC resample size
  MC0 <- dat[ , (names(dat) %in% c("X", "Z1", "Z2", "Z3"))]   #Pull out baseline variables
  MC <- MC0[sample(1:nrow(MC0), montecarlo, replace=TRUE), ]  #Sample with replacement
  MC$id<-1:montecarlo  #Assign new IDs
  
  #Use parametric g-formula (pgf) function to reconstruct follow-up for each person
  pgf <- function(i, data){
    #Pull out record for individual i
    g_data <- data[data$id==i, ]

    #Predict time using observed exposure and confounders
    desX <- g_data[ , c("X", "Z1", "Z2", "Z3")]
    time1<-rexp(1, exp(coef(mod.t)[names(coef(mod.t))=="rate"])*
                       exp(coef(mod.t)[!names(coef(mod.t))=="rate"]%*%t(desX)))
    
    #For time-fixed exposure, not necessary to predict exposure like this
    #But will be necessary in time-varying settings, so let's see what happens here
    g_data$X2 <- as.numeric(predict(mod.x, newdata=g_data, type="response")>runif(1))
    desX <- g_data[ , c("X2", "Z1", "Z2", "Z3")] %>% 
      rename(X=X2)
    time2<-rexp(1, exp(coef(mod.t)[names(coef(mod.t))=="rate"])*
                       exp(coef(mod.t)[!names(coef(mod.t))=="rate"]%*%t(desX)))
    
    time <- data.frame(time1, time2)
    return(time)
  }
  
  #Run pgf, looping over each observation in the MC resample
  T_pred <- lapply(1:montecarlo, function(x){pgf(x, MC)})
  T_pred <- do.call(rbind, T_pred)
    #Outcome using MC-sampled exposure
    delta2 <- ifelse(T_pred[ , 1] > K, 0, 1)
    t_pred2 <- ifelse(T_pred[ , 1] > K, K, T_pred[ , 1])
    #Outcome using predicted exposure
    delta3 <- ifelse(T_pred[ , 2] > K, 0, 1)
    t_pred3 <- ifelse(T_pred[ , 2] > K, K, T_pred[ , 2])

    
###################################################################################################################################    
# Estimate g-formula natural course survival

#Predict from model without MC resample 
sum(delta1)/n #P(Y=1) = 0.4738
surv_gform1 <- survfit(Surv(t_pred1, delta1) ~ 1)
    
#Using MC-sampled exposure
sum(delta2)/montecarlo #P(Y=1) = 0.4752
surv_gform2 <- survfit(Surv(t_pred2, delta2) ~ 1) 
  
#Using model-predicted exposure  
sum(delta3)/montecarlo #P(Y=1) = 0.471
surv_gform3 <- survfit(Surv(t_pred3, delta3) ~ 1) 

#Compare survival curves
surv <- list(Observed=surv_obs, Gformula1=surv_gform1, Gformula2=surv_gform2, Gformula3=surv_gform3)
ggsurvplot(surv, combine=TRUE, xlim=c(0,1), ylim=c(0.50, 1.00), break.time.by=0.25,
           legend=c(0.75,0.8), legend.title="", legend.labs=c("Observed natural course", 
                                                              "\n g-formula natural course \n (no MC-sampling)",
                                                              "\n g-formula natural course \n (MC, sampled exposure)",
                                                              "\n g-formula natural course \n (MC, model-predicted exposure)")) 

#As the above are all equivalent, in future choose faster no MC implementation
surv <- list(Observed=surv_obs, Gformula1=surv_gform1)
pdf("~/Documents/Pitt Projects/Natural Course/results/timefix_simple.pdf", height=5, width=5)
ggsurvplot(surv, combine=TRUE, fun="event", xlim=c(0,1), break.time.by=0.25, ylab="Natural Course Risk", 
           palette=c("dark gray", "black"), linetype="strata",
           legend=c(0.2,0.93), font.legend=c(12, "plain", "black"),
           legend.title="", legend.labs=c("Observed", "G-computation")) 
dev.off()

