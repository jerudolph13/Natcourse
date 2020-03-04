

###################################################################################################
#
# Purpose: For g-formula natural course paper, time-fixed case 2: binary exposure,
#             2 binary confounders, 2 continuous confounders, interaction terms,
#             complex survival outcome (Weibull mixture model), no LTFU
#
# Author: Jacqueline Rudolph
#
# Last Update: 02 mar 2020
#
##################################################################################################

#Things to explore:
#   - Incorrect but flexible model
#   - Moderately misspecified model
#   - Grossly misspecified model


packages <- c("survival", "tidyverse", "flexsurv", "survminer", "simsurv", "splines")

for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}


### Simulate data

set.seed(123)
n <- 5000
montecarlo <- 50000
K <- 1

Z1<-Z2<-Z3<-Z3_sq<-Z4<-X<-Tv<-Y <- rep(0, n)

for (i in 1:n){

  Z1[i] <- rbinom(1, 1, 0.4)
  Z2[i] <- rbinom(1, 1, 0.2)
  Z3[i] <- rnorm(1, 0, 1)
  Z3_sq[i] <- Z3[i]^2
  Z4[i] <- rgamma(1, shape=2, scale=1/2)
  
  p_x <- 1/(1 + exp(-(-log(1/0.5 - 1) - log(2)*Z1[i] + log(2)*0.4 + log(1.5)*Z2[i] - log(1.5)*0.2 
                        - log(3)*Z3[i] - log(3)*Z3_sq[i] + log(4)*Z4[i]
                        + log(1.8)*Z1[i]*Z2[i] - log(1.8)*0.4*0.2)))
  X[i] <- rbinom(1, 1, p_x)

}


cov <- data.frame(X, Z1, Z2, Z3, Z3_sq, Z4, int12=Z1*Z2, intx2=X*Z2, intx4=X*Z4)
betas <- c(X=log(2), Z1=log(2), Z2=log(1.5), Z3=log(3), Z3_sq=log(3), Z4=-log(1.5), 
           int12=log(1.8), intx2=1.5, intx4=-0.5)

out_sim <- simsurv(dist="weibull", lambdas=c(1.5, 0.1), gammas=c(3.0, 1.2),
              mixture=TRUE, pmix=0.2,
              betas = betas, 
              x = cov,                     
              maxt = 1,                   
              interval=c(0,500))
out_sim <- out_sim[ , -1]
out_sim <- rename(out_sim, Tv=eventtime, Y=status)

dat <- cbind(cov, out_sim)


### Observed survival

surv_obs <- survfit(Surv(Tv, Y) ~ 1, data=dat)


### Implement (flexible) parametric g-formula for natural course

mod.x <- glm(X ~ Z1 + Z2 + Z3 + Z3_sq + bs(Z4, knots=4, degree=2) + int12, family=binomial(link="logit"), data=dat)

mod.t <- flexsurvreg(Surv(Tv, Y) ~ X + Z1 + Z2 + Z3 + Z4 + int12 + intx2 + intx4, 
                     data=dat, dist="gengamma")

MC0 <- dat[ , -c("Tv", "Y")]
index <- sample(1:nrow(MC0), montecarlo, replace=TRUE)
MC<-MC0[index,]
MC$id<-1:montecarlo


pgf <- function(i, data){
  g_data <- data[data$id==i, ]
  g_data$X2 <- as.numeric(predict(mod.x, newdata=g_data, type="response")>runif(1))
  
  desX <- g_data[ , c("X", "Z1", "Z2", "Z3")]
  time1 <- 
  return(time)
}

T_pred <- lapply(1:montecarlo, function(x){pgf(x, MC)})
T_pred <- do.call(rbind, T_pred)
Y_pred <- ifelse(T_pred > K, 0, 1)
T_pred <- ifelse(T_pred > K, K, T_pred)

### Estimate g-formula NC survival

surv_gform <- survfit(Surv(T_pred, Y_pred) ~ 1)

### Compare survival curves

surv <- list(Observed=surv_crude, Gformula=surv_gform)
ggsurvplot(surv, combine=TRUE, xlim=c(0,1), ylim=c(0.50, 1.00), break.time.by=0.25,
           legend=c(0.75,0.95), legend.title="", legend.labs=c("Observed natural course", "G-formula natural course")) 


