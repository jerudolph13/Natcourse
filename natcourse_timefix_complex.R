

###################################################################################################
#
# Purpose: For g-formula natural course paper, time-fixed case 2: binary exposure,
#             2 binary confounders, 2 continuous confounders, interaction terms,
#             complex survival outcome (Weibull mixture model), LTFU
#
# Author: Jacqueline Rudolph
#
# Last Update: 05 mar 2020
#
##################################################################################################

# How concerned should I be about the following warning:
# glm.fit: fitted probabilities numerically 0 or 1 occurred 


packages <- c("survival", "tidyverse", "flexsurv", "survminer", "simsurv", "splines", "data.table")

for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}


### Simulate data

set.seed(123)
n <- 5000
montecarlo <- 10000
K <- 10

X<-Z1<-Z2<-Z3<-Z3_sq<-Z4<-Z4_cub<-C<-Tv<-Y <- rep(0, n)

for (i in 1:n){
  Z1[i] <- rbinom(1, 1, 0.4)
  Z2[i] <- rbinom(1, 1, 0.2)
  Z3[i] <- rnorm(1, 0, 1)
  Z3_sq[i] <- Z3[i]^2
  Z4[i] <- rgamma(1, shape=2, scale=1/2)
  Z4_cub[i] <- Z4[i]^3
  
  #Generate exposure
  p_x <- 1/(1 + exp(-(-log(1/0.5 - 1) - log(2)*Z1[i] + log(2)*0.4 + log(1.5)*Z2[i] - log(1.5)*0.2 
                        - log(2.5)*Z3[i] - log(2.5)*Z3_sq[i] + log(1.8)*Z4[i] + log(1.8)*Z4_cub[i]
                        + log(1.3)*Z1[i]*Z2[i] - log(1.3)*0.4*0.2)))
  X[i] <- rbinom(1, 1, p_x)
}

#Covariates and coefficients to simulate drop out
cov_c <- data.frame(Z1, Z2, Z3, Z3_sq, Z4, Z4_cub, int13=Z1*Z3, int13sq=Z1*Z3_sq)
betas_c <- c(Z1=log(2), Z2=log(1.5), Z3=log(2.5), Z3_sq=log(2.5), Z4=-log(1.5), Z4_cub=-log(1.5), int13=log(1.3), int13sq=log(1.3))

#Covariates and coefficients to simulate outcome
cov_y <- data.frame(X, Z1, Z2, Z3, Z3_sq, Z4, Z4_cub, intx2=X*Z2)
betas_y <- c(X=-log(3), Z1=log(2), Z2=log(1.5), Z3=log(2.5), Z3_sq=log(2.5), Z4=-log(1.5), Z4_cub=-log(1.5), intx2=-0.5)

out_sim <- simsurv(dist="weibull", lambdas=c(1.5, 0.1), gammas=c(3.0, 1.2),
              mixture=TRUE, pmix=0.05,
              betas = betas_y, 
              x = cov_y,                     
              maxt = K,                   
              interval=c(0,500))
out_sim <- out_sim[ , -3]
out_sim <- rename(out_sim, Tv=eventtime)
dat <- cbind(cov_y, out_sim)

out_sim <- simsurv(dist="weibull", lambdas=0.5, gammas=0.3,
                   betas = betas_c, 
                   x = cov_c,                     
                   maxt = K,                   
                   interval=c(0,500))
Cv <- out_sim[ , "eventtime"]
dat <- cbind(dat, Cv)

makelong <- function(i,data,maxint){
  
  newdat <- data[data$id==i, ]
  X<-Z1<-Z2<-Z3<-Z4<-C<-Int<-Y<- numeric()
  
  Int[1] <- 1
  X[1] <- newdat$X; Z1[1] <- newdat$Z1; Z2[1] <- newdat$Z2; Z3[1] <- newdat$Z3; Z4[1] <- newdat$Z4
  if (newdat$Cv <= 1 & newdat$Cv<=newdat$Tv) {
    C[1] <- 1
    Y[1] <- 0
    } else {
      C[1] <- 0
      if (newdat$Tv <= 1) {
        Y[1] <- 1
      } else {
        Y[1] <- 0
      }
        }
  
  for (t in 2:maxint) {
    if (Y[t-1] == 0 & C[t-1] == 0) {
      Int[t] <- t
      X[t] <- X[t-1]; Z1[t] <- Z1[t-1]; Z2[t] <- Z2[t-1]; Z3[t] <- Z3[t-1]; Z4[t] <- Z4[t-1]
      if (newdat$Cv > (t-1) & newdat$Cv <= t & newdat$Cv <= newdat$Tv & newdat$Cv != maxint) {
        C[t] <- 1
        Y[t] <- 0
      } else {
          C[t] <- 0
          if (newdat$Tv > (t-1) & newdat$Tv <= t & newdat$Tv!=maxint) {
            Y[t] <- 1
          } else {
              Y[t] <- 0
              }
        }
    } else {
      break
      }
  }
  
  ID <- rep(i, max(Int))
  long <- data.frame(ID, Int, X, Z1, Z2, Z3, Z4, C, Y)
  return(long)
}

long <- lapply(1:n, function (x) {makelong(x, data=dat, maxint=K)})
long <- do.call(rbind, long)


### Observed survival

surv_obs <- survfit(Surv(Int, Y) ~ 1, data=long)


### Implement parametric g-formula for natural course
    #Model 1: As flexible as possible
    #Model 2: No interaction terms
    #Model 3: No interaction terms and linear continuous variables

mod.c1 <- glm(C ~ Z1 + Z2 + bs(Z3, knots=2, degree=2) + bs(Z4, knots=2, degree=2) + I(Z1*bs(Z3, knots=2, degree=2)) + as.factor(Int), 
             family=binomial(link="logit"), data=long)
mod.c2 <- glm(C ~ Z1 + Z2 + bs(Z3, knots=2, degree=2) + bs(Z4, knots=2, degree=2) + as.factor(Int), 
              family=binomial(link="logit"), data=long)
mod.c3 <- glm(C ~ Z1 + Z2 + Z3 + Z4 + Int, 
              family=binomial(link="logit"), data=long)

mod.y1 <- glm(Y ~ X + Z1 + Z2 + bs(Z3, knots=2, degree=2) + bs(Z4, knots=2, degree=2) + I(X*Z2) + as.factor(Int), 
             family=binomial(link="logit"), data=long)
mod.y2 <- glm(Y ~ X + Z1 + Z2 + bs(Z3, knots=2, degree=2) + bs(Z4, knots=2, degree=2) + as.factor(Int), 
             family=binomial(link="logit"), data=long)
mod.y3 <- glm(Y ~ X + Z1 + Z2 + Z3 + Z4 + Int, 
             family=binomial(link="logit"), data=long)

MC0 <- long[long$Int == 1, c("X", "Z1", "Z2", "Z3", "Z4")]
index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)
MC <- MC0[index, ]
MC$id <- 1:montecarlo

pgf <- function(i, data, cmodel, ymodel, maxint){
  g_data <- data[data$id==i, ]
  X<-Z1<-Z2<-Z3<-Z4<-Int<-C<-Y<- numeric()
  
  Int[1] <- 1
  Z1[1] <- g_data$Z1; Z2[1] <- g_data$Z2; Z3[1] <- g_data$Z3; Z4[1] <- g_data$Z4
  X[1] <- g_data$X
  
  desX <- data.table(Z1=Z1[1], Z2=Z2[1], Z3=Z3[1], Z4=Z4[1], Int=Int[1])
  C[1] <- as.numeric(predict(cmodel, newdata=desX, type="response")>runif(1))
  
  if (C[1]==1) {
    Y[1] <- 0
  } else {
    desX <- data.table(X=X[1], Z1=Z1[1], Z2=Z2[1], Z3=Z3[1], Z4=Z4[1], Int=Int[1])
    Y[1] <- as.numeric(predict(ymodel, newdata=desX, type="response")>runif(1))
    }


  for (t in 2:maxint) {
    if (Y[t-1]==0 & C[t-1]==0) {
      Int[t] <- t
      X[t] <- X[t-1]; Z1[t] <- Z1[t-1]; Z2[t] <- Z2[t-1]; Z3[t] <- Z3[t-1]; Z4[t] <- Z4[t-1]
      
      desX <- data.table(Z1=Z1[t], Z2=Z2[t], Z3=Z3[t], Z4=Z4[t], Int=Int[t])
      C[t] <- as.numeric(predict(cmodel, newdata=desX, type="response")>runif(1))
      
      if (C[t]==1) {
        Y[t] <- 0
      } else {
        desX <- data.table(X=X[t], Z1=Z1[t], Z2=Z2[t], Z3=Z3[t], Z4=Z4[t], Int=Int[t])
        Y[t] <- as.numeric(predict(ymodel, newdata=desX, type="response")>runif(1))
      }
    } else {
      break
    }
  }
  
  ID <- rep(i, max(Int))
  Int0 <- Int - 1
  res <- data.frame(ID, Int, Int0, X, Y)
  return(res)
}

gform1 <- lapply(1:montecarlo, function(x){pgf(x, data=MC, cmodel=mod.c1, ymodel=mod.y1, maxint=K)})
gform1 <- do.call(rbind, gform1)

gform2 <- lapply(1:montecarlo, function(x){pgf(x, data=MC, cmodel=mod.c2, ymodel=mod.y2, maxint=K)})
gform2 <- do.call(rbind, gform2)

gform3 <- lapply(1:montecarlo, function(x){pgf(x, data=MC, cmodel=mod.c3, ymodel=mod.y3, maxint=K)})
gform3 <- do.call(rbind, gform3)

### Estimate g-formula NC survival

surv_gform1 <- survfit(Surv(Int, Y) ~ 1, data=gform1)
surv_gform2 <- survfit(Surv(Int, Y) ~ 1, data=gform2)
surv_gform3 <- survfit(Surv(Int, Y) ~ 1, data=gform3)

### Compare survival curves

surv <- list(Observed=surv_obs, Gformula=surv_gform1)
ggsurvplot(surv, combine=TRUE, xlim=c(0,K), ylim=c(0.85, 1.0), break.time.by=1,
           legend=c(0.75,0.95), legend.title="", legend.labs=c("Observed natural course", "G-formula natural course")) 

surv <- list(Observed=surv_obs, Gformula=surv_gform2)
ggsurvplot(surv, combine=TRUE, xlim=c(0,K), ylim=c(0.85, 1.0), break.time.by=1,
           legend=c(0.75,0.95), legend.title="", legend.labs=c("Observed natural course", "G-formula natural course"))

surv <- list(Observed=surv_obs, Gformula=surv_gform3)
ggsurvplot(surv, combine=TRUE, xlim=c(0,K), ylim=c(0.85, 1.0), break.time.by=1,
           legend=c(0.75,0.95), legend.title="", legend.labs=c("Observed natural course", "G-formula natural course"))


