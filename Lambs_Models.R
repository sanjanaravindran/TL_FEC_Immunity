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
library(ggpubr)

###################################################################################
options(digits=14) 

############## Load datasets ######################################################
para_tel_lambs <- read.table("RTL_FEC_Lambs.csv", header=T, sep=",")
ig_tel_lambs <- read.table("RTL_IgG_Lambs.csv", header=T, sep=",")
ig_para_lambs <- read.table("RTL_FEC_IgG_Lambs.csv", header=T, sep=",")

####################################################################################################

##################
# FEC model ######
##################
lambs_fec_mod0 <- stan_lmer(RTL ~ Strongyles + (1|SampleYear) + (1|qPCRPlate) + (1|qPCRRow),
                            cores=8, 
                            seed=12345,
                            iter=8000,
                            data=para_tel_lambs)


# # Model Results 
head(print(summary(lambs_fec_mod0, digits=5, probs=c(0.025, 0.975))),4)

tail(print(summary(lambs_fec_mod0, digits=5, probs=c(0.025, 0.975))),7)

####################
# FEC and Sex model#
####################
lambs_fec_mod1  <- stan_lmer(RTL ~ Strongyles * Sex_R + (1|SampleYear) + (1|qPCRPlate) + (1|qPCRRow),
                             cores=8, 
                             seed=12345,
                             iter=8000,
                             data=para_tel_lambs)

# # Model Results 
head(print(summary(lambs_fec_mod1, digits=5, probs=c(0.025, 0.975))),4)

tail(print(summary(lambs_fec_mod1, digits=5, probs=c(0.025, 0.975))),7)

####################################################################################################

####################
# IgG model ########
####################
lambs_igg_mod0 <- stan_lmer(RTL ~ IgG_TC + (1|SampleYear) + (1|qPCRPlate) + (1|qPCRRow),
                            cores=8, 
                            seed=12345,
                            iter=8000,
                            data=ig_tel_lambs)


# Model Results 
head(print(summary(lambs_igg_mod0, digits=5, probs=c(0.025, 0.975))),4)

tail(print(summary(lambs_igg_mod0, digits=5, probs=c(0.025, 0.975))),6)

####################
# IgG and sex model#
####################
lambs_igg_mod1 <- stan_lmer(RTL ~ IgG_TC * Sex_R + (1|SampleYear) + (1|qPCRPlate) + (1|qPCRRow),
                            cores=8, 
                            seed=12345,
                            iter=8000,
                            data=ig_tel_lambs)

# Model Results 
head(print(summary(lambs_igg_mod1, digits=5, probs=c(0.025, 0.975))),4)

tail(print(summary(lambs_igg_mod1, digits=5, probs=c(0.025, 0.975))),6)

####################################################################################################

########################
# FEC, IgG and TL model#
########################
lambs_fec_igg_mod0 <- stan_lmer(RTL ~ Strongyles + IgG_TC + (1|SampleYear) + (1|qPCRPlate) + (1|qPCRRow),
                                cores=8, 
                                seed=12345,
                                iter=8000,
                                data=ig_para_lambs)


head(print(summary(lambs_fec_igg_mod0, digits=5, probs=c(0.025, 0.975))),4)
tail(print(summary(lambs_fec_igg_mod0, digits=5, probs=c(0.025, 0.975))),6)

######################################################################
