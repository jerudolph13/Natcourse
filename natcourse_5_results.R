
###############################################################################################
#
# Project: Natural course paper
#
# Purpose: Summarize and visualize results
#
# Author: Jacqueline Rudolph
#
# Last Update: 01 Jun 2021
#
############################################################################################## 

library("tidyverse")

source("./plot_thm.R")


# Read in results ---------------------------------------------------------

res <- read.table(file="../results/eager_timevar_results.txt", sep="\t", header=TRUE)


# Summarize results -------------------------------------------------------

res2 <- res %>% filter(week==26)

paste0("Natural course effect NOT controlling for drop out = ", round(res2$rd_nc_drop, 3),
       " (95% CI: ", round(res2$lower_nc_drop, 3), ", ", round(res2$upper_nc_drop, 3), ")")

paste0("Natural course effect controlling for drop out = ", round(res2$rd_nc_nodrop, 3),
       " (95% CI: ", round(res2$lower_nc_nodrop, 3), ", ", round(res2$upper_nc_nodrop, 3), ")")

paste0("ACE NOT controlling for drop out = ", round(res2$rd_ate_drop, 3),
       " (95% CI: ", round(res2$lower_ate_drop, 3), ", ", round(res2$upper_ate_drop, 3), ")")

paste0("ACE controlling for drop out = ", round(res2$rd_ate_nodrop, 3),
       " (95% CI: ", round(res2$lower_ate_nodrop, 3), ", ", round(res2$upper_ate_nodrop, 3), ")")


# Visualize results -------------------------------------------------------

# Risk and RD functions for natural course effect
cols<-c("Natural Course"="dark gray", "Always Compliant"="black")
jpeg("../figures/timevar_nc-effect_risk.jpg", height=5, width=5, units="in", res=300)
ggplot(res) + thm1 +
  labs(x="\nTime (weeks)", y="Cumulative incidence of pregnancy\n", tag="A)") +
  geom_step(aes(x=week, y=r_obs, color="Natural Course"), size=0.75) +
  geom_step(aes(x=week, y=r_nodrop1, color="Always Compliant"), size=0.75) +
  scale_colour_manual(name="", values=cols) +
  scale_x_continuous(expand=c(0, 0), limits=c(1, 26)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 0.8))
dev.off()

jpeg("../figures/timevar_nc-effect_rd.jpg", height=5, width=5, units="in", res=300)
ggplot(res) + thm1 +
  labs(x="\nTime (weeks)", y="Risk difference\n", tag="B)") +
  geom_ribbon(aes(x=week, ymin=lower_nc_nodrop, ymax=upper_nc_nodrop), alpha=0.5) +
  geom_line(aes(x=week, y=rd_nc_nodrop)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(expand=c(0, 0), limits=c(1, 26)) +
  scale_y_continuous(expand=c(0, 0), limits=c(-0.03, 0.08), breaks=c(-0.02, 0.0, 0.02, 0.04, 0.06, 0.08))
dev.off()

# Risk and RD functions for the ATE
cols<-c("Never Compliant"="dark gray", "Always Compliant"="black")
jpeg("../figures/timevar_ate_risk.jpg", height=5, width=5, unit="in", res=300)
ggplot(res) + thm1 +
  labs(x="\nTime (weeks)", y="Cumulative incidence of pregnancy\n", tag="A)") +
  geom_step(aes(x=week, y=r_int0, color="Never Compliant"), size=0.75) +
  geom_step(aes(x=week, y=r_int1, color="Always Compliant"), size=0.75) +
  scale_colour_manual(name="", values=cols) +
  scale_x_continuous(expand=c(0, 0), limits=c(1, 26)) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 0.8))
dev.off()

jpeg("../figures/timevar_ate_rd.jpg", height=5, width=5, unit="in", res=300)
ggplot(res) + thm1 +
  labs(x="\nTime (weeks)", y="Risk difference\n", tag="B)") +
  geom_ribbon(aes(x=week, ymin=lower_ate_nodrop, ymax=upper_ate_nodrop), alpha=0.5) +
  geom_line(aes(x=week, y=rd_ate_nodrop)) +
  geom_hline(yintercept=0, linetype="dashed") +
  scale_x_continuous(expand=c(0, 0), limits=c(1, 26)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 0.4), breaks=seq(0, 0.4, 0.05))
dev.off()

