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
library(rstanarm)
library(MCMCglmm)
library(coda)

###################################################################################
############## Load data
options(digits=14) 

Adults3Yrs_temp <- read.table("RTL_FEC_IgG_3YrsAdults.csv", header=TRUE, sep=",")

########## Rescaling variables 
Adults3Yrs_temp$Sex_R <- as.factor(Adults3Yrs_temp$Sex_R) 
Adults3Yrs_temp$RTL_z <- scale(Adults3Yrs_temp$RTL)
Adults3Yrs_temp$IgG_TC_sc <- scale(Adults3Yrs_temp$IgG_TC)
Adults3Yrs_temp$FEC_sc <- round(Adults3Yrs_temp$Strongyles/100)
Adults3Yrs_temp$Age_z <- scale(Adults3Yrs_temp$Age)
Adults3Yrs_temp <- as.data.frame(Adults3Yrs_temp)

###############################################################################################
###################################################################################################
# Prior
prior_TL_FEC_IgG = list(G = list(G1 = list(V=diag(3), nu=3, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                                 G2 = list(V=diag(3), nu=3, alpha.mu=rep(0,3), alpha.V=diag(3)*1000),
                                 G3 = list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                                 G4 = list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000),
                                 G5 = list(V=diag(1), nu=1, alpha.mu=rep(0,1), alpha.V=diag(1)*1000)),
                                 R = list(V = diag(3), nu = 3.002))
                        


# Trivariate Model (3 Yr Adults and Above)
trivar_TL_FEC_IgG_num <- MCMCglmm::MCMCglmm(cbind(RTL_z, FEC_sc, IgG_TC_sc) ~ trait - 1 + 
                                              at.level(trait,1):Age_z + 
                                              at.level(trait,2):poly(Age_z,2) + at.level(trait,2):Sex_R + 
                                              at.level(trait,3):Age_z + at.level(trait,3):Sex_R,
                                            data=Adults3Yrs_temp,
                                            random = ~us(trait):ID + us(trait):SampleYear + us(at.level(trait,1)):qPCRPlate + us(at.level(trait,1)):qPCRRow + us(at.level(trait,3)):IgPlateRef, 
                                            rcov = ~us(trait):units,  
                                            family = c("gaussian", "poisson", "gaussian"),
                                            prior = prior_TL_FEC_IgG, 
                                            trunc=T,
                                            nitt=495000,
                                            thin=400,
                                            burnin=95000)

# Sex-specific costs of parasitism on RTL 
adults_fec_mod1  <- stan_lmer(RTL_z ~ FEC_sc * Sex_R + Age + (1|ID) + (1|SampleYear) + (1|qPCRPlate) + (1|qPCRRow),
                              cores=8, 
                              seed=12345,
                              iter=8000,
                              data= Adults3Yrs_temp)

head(print(summary(adults_fec_mod1, digits=5, probs=c(0.025, 0.975))),4)

tail(print(summary(adults_fec_mod1, digits=5, probs=c(0.025, 0.975))),7)

# Sex-specific costs of immune response on RTL 
adults_igg_mod1  <- stan_lmer(RTL_z ~ IgG_TC_sc * Sex_R + Age + (1|ID) + (1|SampleYear) + (1|qPCRPlate) + (1|qPCRRow),
                              cores=8, 
                              seed=12345,
                              iter=8000,
                              data= Adults3Yrs_temp)

head(print(summary(adults_igg_mod1, digits=5, probs=c(0.025, 0.975))),4)

tail(print(summary(adults_igg_mod1, digits=5, probs=c(0.025, 0.975))),7)

############################################################################################
############################################################################################
############################################################################################


# Models run using data of adults aged 1Yr and above
###################################################################################
############## Load data
Adults1Yrs_temp <- read.table("RTL_FEC_IgG_1YrsAdults.csv", header=TRUE, sep=",")

########## Rescaling variables 
Adults1Yrs_temp$Sex_R <- as.factor(Adults1Yrs_temp$Sex_R) 
Adults1Yrs_temp$RTL_z <- scale(Adults1Yrs_temp$RTL)
Adults1Yrs_temp$IgG_TC_sc <- scale(Adults1Yrs_temp$IgG_TC)
Adults1Yrs_temp$FEC_sc <- round(Adults1Yrs_temp$Strongyles/100)
Adults1Yrs_temp$Age_z <- scale(Adults1Yrs_temp$Age)
Adults1Yrs_temp <- as.data.frame(Adults1Yrs_temp)

###############################################################################################

# Trivariate Model (1 Yr Adults and Above)
trivar_TL_FEC_IgG_num <- MCMCglmm::MCMCglmm(cbind(RTL_z, FEC_sc, IgG_TC_sc) ~ trait - 1 + 
                                              at.level(trait,1):Age_z + 
                                              at.level(trait,2):poly(Age_z,2) + at.level(trait,2):Sex_R + 
                                              at.level(trait,3):Age_z + at.level(trait,3):Sex_R,
                                            data=Adults1Yrs_temp,
                                            random = ~us(trait):ID + us(trait):SampleYear + us(at.level(trait,1)):qPCRPlate + us(at.level(trait,1)):qPCRRow + us(at.level(trait,3)):IgPlateRef, 
                                            rcov = ~us(trait):units,  
                                            family = c("gaussian", "poisson", "gaussian"),
                                            prior = prior_TL_FEC_IgG, 
                                            trunc=T,
                                            nitt=495000,
                                            thin=400,
                                            burnin=95000)

###############################################################################################
###################################################################################################
