############################################################################################
############## Load packages
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(janitor)
library(tidyr)
library(zoo)
library(tibble)
library(lmerTest)
library(lme4)
library(lmerTest)
library(MCMCglmm)
library(coda)

###################################################################################
############## Load data
options(digits=14) 

FullPop_temp <- read.table("RTL_FEC_IgG_FullPop.csv", header=TRUE, sep=",")

# Rescaling variables 
FullPop_temp$Sex_R <- as.factor(FullPop_temp$Sex_R) 
FullPop_temp$AgeClass <- as.factor(FullPop_temp$AgeClass) 
FullPop_temp$RTL_z <- scale(FullPop_temp$RTL)
FullPop_temp$IgG_TC_sc <- scale(FullPop_temp$IgG_TC)
FullPop_temp$FEC_sc <- round(FullPop_temp$Strongyles/100)
FullPop_temp$Age_z <- scale(FullPop_temp$Age)
FullPop_temp <- as.data.frame(FullPop_temp)

###############################################################################################
###############################################################################################
# Prior
prior_TL_FEC_IgG_Surv <- list(G=list(G1=list(V=diag(4), nu=1, alpha.mu=c(0,0,0,0), alpha.V=diag(4)*1000), 
                                     G2=list(V=diag(4), nu=1, alpha.mu=c(0,0,0,0), alpha.V=diag(4)*1000),
                                     G3=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                                     G4=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000),
                                     G5=list(V=diag(1), nu=1, alpha.mu=0, alpha.V=diag(1)*1000)),
                              R=list(V=diag(4), nu=4.002, fix=4))

# Tetravariate model with survival
tetravar_TL_FEC_IgG_Survival <- MCMCglmm::MCMCglmm(cbind(RTL_z, FEC_sc, IgG_TC_sc, Survival) ~ trait - 1 +
                                                     at.level(trait, 1):Age_z +  at.level(trait, 1):AgeClass + 
                                                     at.level(trait, 2):AgeClass + at.level(trait, 2):poly(Age_z,2)*Sex_R + 
                                                     at.level(trait, 3):poly(Age_z,3) + at.level(trait, 3):AgeClass*Sex_R + 
                                                     at.level(trait, 4):poly(Age_z,3) + at.level(trait, 4):AgeClass*Sex_R,
                                                   data=FullPop_temp,
                                                   random = ~us(trait):SampleYear + us(trait):ID + us(at.level(trait, 1)):qPCRPlate + us(at.level(trait, 1)):qPCRRow + us(at.level(trait, 3)):IgPlateRef, 
                                                   rcov = ~us(trait):units,  
                                                   family = c("gaussian", "poisson", "gaussian", "threshold"),
                                                   prior = prior_TL_FEC_IgG_Surv, 
                                                   trunc=T,
                                                   nitt=795000,
                                                   thin=500,
                                                   burnin=95000)

###############################################################################################
###################################################################################################
