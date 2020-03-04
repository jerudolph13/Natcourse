

###################################################################################################
#
# Purpose: For g-formula natural course paper, time-fixed case 1: binary exposure,
#             3 binary confounders, exponential outcome, no LTFU
#
# Author: Jacqueline Rudolph
#
# Last Update: 02 mar 2020
#
##################################################################################################

# Things I could explore:
#     - Including more variables than are necessary, esp. interaction terms
#     - Being flexible in outcome model

packages <- c("survival", "tidyverse", "flexsurv", "survminer")

for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}


### Simulate data

set.seed(123)
n <- 5000
montecarlo <- 50000
K <- 1

Z1<-Z2<-Z3<-X<-Tv<-Y <- rep(0, n)

for (i in 1:n){
  
  Z1[i] <- rbinom(1, 1, 0.4)
  Z2[i] <- rbinom(1, 1, 0.2)
  Z3[i] <- rbinom(1, 1, 0.6)
  
  p_x <- 1/(1 + exp(-(-log(1/0.5 - 1) - log(2)*Z1[i] + log(2)*0.4 + log(1.5)*Z2[i] - log(1.5)*0.2 - log(3)*Z3[i] + log(3)*0.6)))
  X[i] <- rbinom(1, 1, p_x)
  
  p_t <- log(2)*X[i] - log(2)*Z1[i] + log(1.5)*Z2[i] - log(3)*Z3[i]
  Tv[i] <- rexp(1, exp(0.01)*exp(p_t))
  Y[i] <- ifelse(Tv[i] > K, 0, 1)
  Tv[i] <- ifelse(Tv[i] > K, K, Tv[i])
}
dat <- data.frame(Z1, Z2, Z3, X, Tv, Y)


### Observed survival

surv_obs <- survfit(Surv(Tv, Y) ~ 1, data=dat)


### Implement parametric g-formula for natural course

mod.x <- glm(X ~ Z1 + Z2 + Z3, family=binomial(link="logit"), data=dat)
mod.t <- flexsurvreg(Surv(Tv, Y) ~ X + Z1 + Z2 + Z3, data=dat, dist="exp")

MC0 <- dat[ , (names(dat) %in% c("X", "Z1", "Z2", "Z3"))]
index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)
MC<-MC0[index,]
MC$id<-1:montecarlo


pgf <- function(i, data){
  g_data <- data[data$id==i, ]
  g_data$X2 <- as.numeric(predict(mod.x, newdata=g_data, type="response")>runif(1))

  desX <- g_data[ , c("X", "Z1", "Z2", "Z3")]
  time1 <- rexp(1, exp(coef(mod.t)[names(coef(mod.t))=="rate"])*
                    exp(coef(mod.t)[!names(coef(mod.t))=="rate"]%*%t(desX)))
  
  desX <- g_data[ , c("X2", "Z1", "Z2", "Z3")]
  time2 <- rexp(1, exp(coef(mod.t)[names(coef(mod.t))=="rate"])*
                  exp(coef(mod.t)[!names(coef(mod.t))=="rate"]%*%t(desX)))
  
  time <- data.frame(time1, time2)
  return(time)
}

T_pred <- lapply(1:montecarlo, function(x){pgf(x, MC)})
T_pred <- do.call(rbind, T_pred)
Y1 <- ifelse(T_pred[ , 1] > K, 0, 1)
Y2 <- ifelse(T_pred[ , 2] > K, 0, 1)
T1 <- ifelse(T_pred[ , 1] > K, K, T_pred[ , 1])
T2 <- ifelse(T_pred[ , 2] > K, K, T_pred[ , 2])

### Estimate g-formula NC survival

surv_gform1 <- survfit(Surv(T1, Y1) ~ 1) #Using MC-sampled X
surv_gform2 <- survfit(Surv(T2, Y2) ~ 1) #Using model-predicted X

### Compare survival curves

surv <- list(Observed=surv_obs, Gformula1=surv_gform1, Gformula2=surv_gform2)
ggsurvplot(surv, combine=TRUE, xlim=c(0,1), ylim=c(0.50, 1.00), break.time.by=0.25,
           legend=c(0.75,0.9), legend.title="", legend.labs=c("Observed natural course", 
                                                               "\n G-formula natural course \n (MC-sampled exposure)",
                                                               "\n G-formula natural course \n (model-predicted exposure)")) 


