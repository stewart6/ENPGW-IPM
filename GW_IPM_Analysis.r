
# Analysis code for gray whale integrated population dynamics model
# Published in Stewart et al. 2023:
# Boom-bust cycles in gray whales associated with dynamic and changing Arctic conditions; Science Magazine
# Model developed by Joshua Stewart in Nimble
# Model improved and translated to Stan by M. Tim Tinker


# LOAD libraries
require(gtools)
require(ggplot2)
require(dplyr)
require(readxl)
require(bayesplot)
require(cmdstanr)
require(posterior)
require(parallel)
require(loo)
library(patchwork)

rstan::rstan_options(javascript=FALSE)

# Import data -----------------------------------
load("GW_IPM_Data_2022.RData")

# Import covariates (two options)
# Two covariates (ice & benthic biomass)
CovFile = "Env_Covariates_x2.xlsx"

# Three covariates (ice, benthic biomass & estimated zooplankton density)
# CovFile = "Env_Covariates_x3.xlsx"

# It is possible to include interactions among covariates in the model. No interactions are included in these formulations.
CovZ = read_excel(CovFile,sheet = 1)
Intxns = read_excel(CovFile,sheet = 2)
#
AddExtraYear = 1 ; # Add an extra year to time series to accommodate either abundance or other data inputs as final data year (as of publication date and with original dataset, use 1; including additional data years may require using 0)


# Set value for assumed proportional difference in 
#  detection rate of anthropogenic vs. natural deaths
#  (1 = identical detection rates, 2 = 2x more likely, etc.)
Anthmort_dtct_factor = 1


# Process data-----------------------------------
# Three stranding periods with separate detection probabilities estimated by the model
Strandings$ObsPeriodUpd <- 1
Strandings$ObsPeriodUpd[Strandings$Year%in%1990:2000] <- 2
Strandings$ObsPeriodUpd[Strandings$Year>2000] <- 3

# Format index years
Year1 = min(Abundance$Year1) ; Year0 = Year1-1
YearT =  max(Abundance$Year1) + AddExtraYear  
Years = seq(Year1,YearT)
Abund.year = Abundance$Year1 - Year0
# Reference vectors
Nyrs = length(Years)
yrm1 = seq(1,Nyrs-1)
Nyrm1 = length(yrm1)
# Data & observations
Abund = Abundance$Abundance
Abund.se = Abundance$SE
BC = filter(BodyCondition,Migration=="Southbound") %>% select(CorrectedYear,MWL)
Cond = BC$MWL
Cond.year = BC$CorrectedYear - Year0
CondNorth = filter(BodyCondition,Migration=="Northbound") %>% select(CorrectedYear,MWL)
CondN <- CondNorth$MWL
CondN.year = CondNorth$CorrectedYear - Year0
# Separate anthropogenic strandings from non-anthropogenic
Strand = Strandings$Count - Strandings$AnthOrigin
AnthStrand = Strandings$AnthOrigin
Strand.year = Strandings$Year - Year0
Strand.period = Strandings$ObsPeriodUpd
Calves = Production$Estimate
Calves.se = Production$SE
Calves.year = Production$Year - Year0
# Process covariates (if no interactions were included, there will be an error; ignore it and continue running code)
V_names = colnames(CovZ); V_names = V_names[-c(1,2)]
NV = length(V_names)
Nitxn = nrow(Intxns)
if(Nitxn>0){
  Intxns$V1 = numeric(length = Nitxn)
  Intxns$V2 = numeric(length = Nitxn)
  for(i in 1:Nitxn){
    Intxns$V1[i] = which(V_names==Intxns$Var1[i])
    Intxns$V2[i] = which(V_names==Intxns$Var2[i])
  }
  Intxmat = as.matrix(Intxns[,-c(1,2,3)])
}
V_z = matrix(0,nrow=Nyrm1,ncol = NV)
yrobs_V = matrix(0,nrow=Nyrm1,ncol = NV)
yrmis_V = matrix(0,nrow=Nyrm1,ncol = NV)
N_Vobs = numeric(length = NV)
N_Vmis = numeric(length = NV)
for(i in 1:NV){
  V_tmp = as.numeric(as.matrix(CovZ[,i+2]))
  Yr_adj = CovZ$WhaleYear-Year0
  for(t in 1:Nyrm1){
    ii = which(Yr_adj==t & !is.na(V_tmp))
    if(length(ii)==1){
      N_Vobs[i] = N_Vobs[i] + 1
      yrobs_V[N_Vobs[i],i] = t
      V_z[N_Vobs[i],i] = (V_tmp[ii] - mean(V_tmp,na.rm = T))/sd(V_tmp,na.rm = T)
    }else{
      N_Vmis[i] = N_Vmis[i] + 1
      yrmis_V[N_Vmis[i],i] = t
    }
  }
}
if(Nitxn>0){
  for(i in 1:Nitxn){
    V_names[NV+i] = paste0(V_names[Intxmat[i,1]],"-",V_names[Intxmat[i,2]])
  }
}

# Set K_hat (see materials and methods)
K_med_approx = 24500

# Data list for model
stan.data <- list(Nyrs = Nyrs,
                  Nyrm1=Nyrm1,
                  N_Abund = length(Abund),
                  N_Calves = length(Calves),
                  N_Strand = length(Strand),
                  N_Cond = length(Cond),
                  N_CondN = length(CondN),
                  NV=NV,
                  N_Vobs=N_Vobs,
                  N_Vmis=N_Vmis,
                  yrobs_V=t(yrobs_V),
                  yrmis_V=t(yrmis_V),
                  V_z=t(V_z),
                  Abund_year=Abund.year,
                  Abund = Abund,
                  Abund_se = Abund.se,
                  Calves = Calves,
                  Calves_se = Calves.se,
                  Calves_year = Calves.year,
                  Strand = Strand,
                  AnthStrand = AnthStrand,
                  Strand_year = Strand.year,
                  Strand_per = Strand.period,
                  Cond = Cond,
                  Cond_year=Cond.year,
                  CondN = CondN,
                  CondN_year = CondN.year,
                  N_init_pri = Abund[1],
                  r_max = 0.1,
                  K_med_approx = K_med_approx,
                  ADF=Anthmort_dtct_factor)
if(Nitxn>0){
  stan.data$Nitxn=Nitxn
  stan.data$Intxvars=Intxmat
}else{
  stan.data$Nitxn=0
  stan.data$Intxvars=numeric()
}

#
# Fit model-------------------------------------------------
#
fitmodel = "GW_modfit_v4.stan"
#
parms = c("sig_E","sig_A","sig_V","Cond_base_lgt","Brate_base_lgt","gammaN_base","gammaA_base",
          "Beta","phi","alpha","gamma","phi_star","alpha_star","gamma_star","K_mean","K_med",
          "NBAdjust","tau","tauN","r_max_exp","Ppn_dtct_N","Ppn_dtct_A","C_max","B_max","M_min","M",
          "B","N","BN","C","CN","K_inst","D_N","D_A","Env_eff","X","log_lik","eps")

nburnin = 500                    # number warm-up (burn-in) samples
nsamples = 5000                  # desired total number samples
cores = detectCores()
ncore = min(10,cores-3)
Niter = round(nsamples/ncore)
mod <- cmdstan_model(fitmodel)   # compiles model (if necessary)
suppressMessages(                # Suppress messages/warnings (if desired)
  suppressWarnings (
    fit <- mod$sample(
      data = stan.data,
      seed = 123,
      chains = ncore,
      parallel_chains = ncore,
      refresh = 100,
      iter_warmup = nburnin,
      iter_sampling = Niter
    )
  )
)
# Note: if error, run: fit$output(1)
# Create table of summary stats for params of interest:

source("cmdstan_sumstats.r")


# Summarize results ----------------------------------------
#
mcmc_trace(fit$draws("sig_E"))
mcmc_trace(fit$draws("sig_A"))
mcmc_trace(fit$draws("Beta"))
#
mcmc_areas(fit$draws(variables = c("Beta")),
           area_method="equal height",
           adjust = 2, prob = 0.8) + 
  ggtitle(paste0("Posterior distributions, Environmental effects")) +
  labs(x="Estimated Covariate Effect",y="Posterior density") +
  scale_y_discrete(labels=V_names) +
  theme_classic()
mean(fit$draws(variables="Beta[1]")>0) # p>0
mean(fit$draws(variables="Beta[2]")>0) # p>0
#mean(fit$draws(variables="Beta[3]")>0) #if using model with 3 covariates

# standardized D-D effects on vital rates:
#  - proportional reduction in Condition, birth rate and survival rate
#    given an increase from 0%K to 100% K
mcmc_areas(fit$draws(variables = c("phi_star","alpha_star","gamma_star")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle(paste0("Posterior distributions, Density-dependent effects")) +
  labs(x="Proportional change, 0% -> 100% K", y="Posterior density") +
  scale_y_discrete(labels=c("Body condition (phi*)",
                            "Calf production (alpha*)",
                            "Survival rate (gamma*)")) +
  theme_classic()
#

#Abundance trend
N_exp = sumstats[which(startsWith(vn,"N[")),1]
N_exp_lo = sumstats[which(startsWith(vn,"N[")),4]  
N_exp_hi = sumstats[which(startsWith(vn,"N[")),8]  
df_Trends = data.frame(Year = Years,N_exp=N_exp,
                       N_exp_lo=N_exp_lo,N_exp_hi=N_exp_hi)
df_Trends$Surveys = rep(NA,nrow(df_Trends))
df_Trends$Surveys[Abund.year] = Abund
df_Trends$Surveys_se = rep(NA,nrow(df_Trends))
df_Trends$Surveys_se[Abund.year] = Abund.se
K_mean = sumstats[which(startsWith(vn,"K_mean")),1]
K_mean_lo = sumstats[which(startsWith(vn,"K_mean")),5]
K_mean_hi = sumstats[which(startsWith(vn,"K_mean")),7]
#
plt_Trnd = ggplot(df_Trends,aes(x=Year,y=N_exp)) +
  geom_ribbon(aes(ymin=N_exp_lo,ymax=N_exp_hi),alpha=0.3) +
  geom_line() +
  geom_hline(yintercept = K_mean, linetype = 2) +
  geom_hline(yintercept = K_mean_lo, linetype = 3) +
  geom_hline(yintercept = K_mean_hi, linetype = 3) +
  geom_point(aes(x=Year,y=Surveys)) +
  geom_errorbar(aes(x=Year,ymin=Surveys-1.96*Surveys_se,
                    ymax=Surveys+1.96*Surveys_se)) +
  labs(x="Year",y="Population abundance") +
  ggtitle(paste0("Estimated trends"),
          subtitle = "Longterm mean K indicated by dashed line") +
  theme_classic()
print(plt_Trnd)
#

# Condition trend
C_exp = sumstats[which(startsWith(vn,"C[")),1]
C_exp_lo = sumstats[which(startsWith(vn,"C[")),4]  
C_exp_hi = sumstats[which(startsWith(vn,"C[")),8]  
df_BCond = data.frame(Year = Years[-length(Years)],C_exp=C_exp,
                       C_exp_lo=C_exp_lo,C_exp_hi=C_exp_hi)
df_BCond$Obs_C = rep(NA,nrow(df_BCond))
df_BCond$Obs_C_se = rep(NA,nrow(df_BCond))
Obs_Cond_sum = BC %>% group_by(CorrectedYear) %>% 
  summarise(BC_mn = mean(MWL),
            BC_sd = sd(MWL))
for(i in 1:nrow(Obs_Cond_sum)){
  df_BCond$Obs_C[df_BCond$Year==Obs_Cond_sum$CorrectedYear[i]] = Obs_Cond_sum$BC_mn[i]
  df_BCond$Obs_C_se[df_BCond$Year==Obs_Cond_sum$CorrectedYear[i]] = Obs_Cond_sum$BC_sd[i]
}  

NorthBC = BodyCondition %>% filter(Migration=="Northbound") %>%
  group_by(CorrectedYear) %>% summarize(BC_mn=mean(MWL),
                                        BC_sd=sd(MWL)) %>%
  rename(Year=CorrectedYear)
CN_exp = sumstats[which(startsWith(vn,"CN[")),1]
CN_exp_lo = sumstats[which(startsWith(vn,"CN[")),4]  
CN_exp_hi = sumstats[which(startsWith(vn,"CN[")),8]  
df_NBCond = data.frame(Year = Years[-length(Years)],CN_exp=CN_exp,
                      CN_exp_lo=CN_exp_lo,CN_exp_hi=CN_exp_hi)

plt_Cond = ggplot(df_BCond,aes(x=Year,y=C_exp)) +
  geom_ribbon(aes(ymin=C_exp_lo,ymax=C_exp_hi),alpha=0.3) +
  geom_line() +
  geom_point(aes(x=Year,y=Obs_C)) +
  geom_errorbar(aes(x=Year,ymin=Obs_C-Obs_C_se,ymax=Obs_C+Obs_C_se)) +
  labs(x="Year",y="Mean Condition") +
  ggtitle(paste0("Trends in Body Condition")) +
  theme_classic() +
  geom_ribbon(data=df_NBCond,aes(ymin=CN_exp_lo,ymax=CN_exp_hi),fill="red",alpha=0.3) +
  geom_line(data=df_NBCond,aes(x=Year,y=CN_exp),color="red",inherit.aes = F) +
  geom_point(data=NorthBC,aes(x=Year,y=BC_mn),color="red",inherit.aes=F) +
  geom_errorbar(data=NorthBC,aes(x=Year,ymin=BC_mn-BC_sd,ymax=BC_mn+BC_sd),color="red",inherit.aes=F)

print(plt_Cond)

# Condition plotted separately
SCond <- ggplot(df_BCond,aes(x=Year,y=C_exp)) +
  geom_ribbon(aes(ymin=C_exp_lo,ymax=C_exp_hi),alpha=0.3) +
  geom_line() +
  geom_point(aes(x=Year,y=Obs_C)) +
  geom_errorbar(aes(x=Year,ymin=Obs_C-Obs_C_se,ymax=Obs_C+Obs_C_se)) +
  labs(x="Year",y="Mean Southbound Condition") +
  ggtitle(paste0("Trends in Body Condition")) +
  theme_classic() 
NCond <- ggplot(df_NBCond,aes(x=Year,y=CN_exp)) +
  geom_ribbon(data=df_NBCond,aes(ymin=CN_exp_lo,ymax=CN_exp_hi),fill="red",alpha=0.3) +
  geom_line(data=df_NBCond,aes(x=Year,y=CN_exp),color="red",inherit.aes = F) +
  geom_point(data=NorthBC,aes(x=Year,y=BC_mn),color="red",inherit.aes=F) +
  geom_errorbar(data=NorthBC,aes(x=Year,ymin=BC_mn-BC_sd,ymax=BC_mn+BC_sd),color="red",inherit.aes=F) +
  labs(x="Year",y="Mean Northbound Condition") +
  ggtitle(paste0("Adjusted Northbound Body Condition")) +
  theme_classic() 

SCond / NCond

# Observation error around mean model fit for SBound and NBound condition
mcmc_areas(fit$draws(variables = c("tau","tauN")),
           area_method="equal height",
           prob = 0.8) + 
  ggtitle(paste0("Posterior distributions, Body Condition Observation Error")) +
  labs(x="Parameter value",y="Posterior density") +
  scale_y_discrete(labels=c("Southbound Error",
                            "Northbound Error")) +
  theme_classic()

# Mortality Trends
D_N_mn = sumstats[which(startsWith(vn,"D_N[")),1]
D_N_lo = sumstats[which(startsWith(vn,"D_N[")),4]  
D_N_hi = sumstats[which(startsWith(vn,"D_N[")),8]  
D_A_mn = sumstats[which(startsWith(vn,"D_A[")),1]
D_A_lo = sumstats[which(startsWith(vn,"D_A[")),4]  
D_A_hi = sumstats[which(startsWith(vn,"D_A[")),8]  
df_Mort = data.frame(Year = Years[-length(Years)],D_est_mn=D_N_mn,
                Source=rep("Natural",Nyrm1),
                D_est_lo=D_N_lo,D_est_hi=D_N_hi)
df_Mort = rbind(df_Mort,data.frame(Year = Years[-length(Years)],
                Source=rep("Anthropogenic",Nyrm1),D_est_mn=D_A_mn,
                D_est_lo=D_A_lo,D_est_hi=D_A_hi))
df_Mort$Source <- factor(df_Mort$Source,levels=c("Natural","Anthropogenic"))
plt_Mort = ggplot(df_Mort,aes(x=Year,y=D_est_mn,group=Source)) +
  geom_ribbon(aes(ymin=D_est_lo,ymax=D_est_hi,fill=Source),alpha=0.3) +
  geom_line(aes(color=Source)) +
  labs(x="Year",y="Estimated Mortality") +
  ggtitle(paste0("Trends in Mortality")) +
  theme_classic()
print(plt_Mort)

# Annual K trend (K_t)
Env_eff = sumstats[which(startsWith(vn,"Env_eff[")),1]
K_inst = sumstats[which(startsWith(vn,"K_inst[")),1]
K_inst_lo = sumstats[which(startsWith(vn,"K_inst[")),4]  
K_inst_hi = sumstats[which(startsWith(vn,"K_inst[")),8]  
df_Kinst = data.frame(Year = Years[1:Nyrm1],K_inst=K_inst,
                      K_inst_lo=K_inst_lo,K_inst_hi=K_inst_hi)
plt_Kinst = ggplot(df_Kinst,aes(x=Year,y=K_inst)) +
  geom_ribbon(aes(ymin=K_inst_lo,ymax=K_inst_hi),alpha=0.3) +
  geom_line() +
  ylim(c(0,150000))+
  geom_hline(yintercept=25000,lty=2) +
  geom_hline(yintercept=10000,lty=3) +
  labs(x="Year",y="Instantaneous K") +
  ggtitle(paste0("Trends in Instantaneous K")) +
  theme_classic()
print(plt_Kinst)


# Covariate 1 state space model fit
X1_exp = sumstats[which(startsWith(vn,"X[") & endsWith(vn, ",1]")),1]
X1_exp_lo = sumstats[which(startsWith(vn,"X[") & endsWith(vn, ",1]")),4]  
X1_exp_hi = sumstats[which(startsWith(vn,"X[") & endsWith(vn, ",1]")),8]  
df_X1 = data.frame(Year = Years[1:Nyrm1],X1_exp=X1_exp,
                   X1_exp_lo=X1_exp_lo,X1_exp_hi=X1_exp_hi)

ggplot(df_X1,aes(x=Year,y=X1_exp)) +
  geom_ribbon(aes(ymin=X1_exp_lo,ymax=X1_exp_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Covariate Z-Score") +
  ggtitle(paste0("Crustacean Biomass")) +
  #geom_vline(xintercept=2020)+
  theme_classic()

# Covariate 2 state space model fit
X2_exp = sumstats[which(startsWith(vn,"X[") & endsWith(vn, ",2]")),1]
X2_exp_lo = sumstats[which(startsWith(vn,"X[") & endsWith(vn, ",2]")),4]  
X2_exp_hi = sumstats[which(startsWith(vn,"X[") & endsWith(vn, ",2]")),8]  
df_X2 = data.frame(Year = Years[1:Nyrm1],X2_exp=X2_exp,
                   X2_exp_lo=X2_exp_lo,X2_exp_hi=X2_exp_hi)

ggplot(df_X2,aes(x=Year,y=X2_exp)) +
  geom_ribbon(aes(ymin=X2_exp_lo,ymax=X2_exp_hi),alpha=0.3) +
  geom_line() +
  labs(x="Year",y="Covariate Z-Score") +
  ggtitle(paste0("Ice Access")) +
  theme_classic()

# Estimated proportion of mortalities detected as strandings, by time period
mcmc_areas(fit$draws(variables = c("Ppn_dtct_N")),
           area_method="equal height",
           adjust = 2, prob = 0.8) + 
  ggtitle(paste0("Stranding Detection Proportion")) +
  labs(x="Proportion of Deaths Detected as Strandings",y="Posterior density") +
  scale_y_discrete(labels=c("Pre-1990","1990-2000","Post-2000")) +
  theme_classic()

# Birth / calf production trend
Bexp = sumstats[which(startsWith(vn,"BN[")),1]
Bexp_lo = sumstats[which(startsWith(vn,"BN[")),4]
Bexp_hi = sumstats[which(startsWith(vn,"BN[")),8]
df_Brate = data.frame(Year = Years[1:Nyrm1],Bexp=Bexp,
                      Bexp_lo=Bexp_lo,Bexp_hi=Bexp_hi)
df_Brate$Surveys = rep(NA,nrow(df_Brate))
df_Brate$Surveys[Calves.year] = Calves
df_Brate$Surveys_se = rep(NA,nrow(df_Brate))
df_Brate$Surveys_se[Calves.year] = Calves.se
colnames(df_Brate)[2:4] <- c("Bexp","Bexp_lo","Bexp_hi")

plt_B = ggplot(df_Brate,aes(x=Year,y=Bexp)) +
  geom_ribbon(aes(ymin=Bexp_lo,ymax=Bexp_hi),alpha=0.3) +
  geom_line() +
  geom_point(aes(x=Year,y=Surveys)) +
  geom_errorbar(aes(x=Year,ymin=Surveys-1.96*Surveys_se,
                    ymax=Surveys+1.96*Surveys_se)) +
  labs(x="Year",y="Number of Births") +
  ggtitle(paste0("Estimated Births")) +
  theme_classic()
print(plt_B)


# Birth rate
Brexp = sumstats[which(startsWith(vn,"B[")),1]
Brexp_lo = sumstats[which(startsWith(vn,"B[")),4]
Brexp_hi = sumstats[which(startsWith(vn,"B[")),8]
df_Brater = data.frame(Year = Years[1:Nyrm1],Bexp=Brexp,
                      Bexp_lo=Brexp_lo,Bexp_hi=Brexp_hi)
colnames(df_Brater)[2:4] <- c("Bexp","Bexp_lo","Bexp_hi")

plt_Br = ggplot(df_Brater,aes(x=Year,y=Bexp)) +
  geom_ribbon(aes(ymin=Bexp_lo,ymax=Bexp_hi),fill="#5A8C47",alpha=0.5) +
  geom_line() +
  labs(x="Year",y="Birth Rate") +
  ggtitle(paste0("Estimated Birth Rate")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))
print(plt_Br)

# Maximum birth rate
mcmc_areas(fit$draws(variables = c("B_max")),
           area_method="equal height",
           adjust = 2, prob = 0.8) + 
  ggtitle(paste0("Posterior distribution, Maximum Birth Rate")) +
  labs(x="Parameter value",y="Posterior density") +
  scale_y_discrete(labels=V_names) +
  theme_classic()
median(fit$draws(variables="B_max"))



# Realized mortality rate (combined hazards)
Mexp = sumstats[which(startsWith(vn,"M[")),1]
Mexp_lo = sumstats[which(startsWith(vn,"M[")),4]
Mexp_hi = sumstats[which(startsWith(vn,"M[")),8]
df_Mrate = data.frame(Year = Years[1:Nyrm1],Mexp=Mexp,
                      Mexp_lo=Mexp_lo,Mexp_hi=Mexp_hi)
colnames(df_Mrate)[2:4] <- c("Mexp","Mexp_lo","Mexp_hi")

plt_Mr = ggplot(df_Mrate,aes(x=Year,y=Mexp)) +
  geom_ribbon(aes(ymin=Mexp_lo,ymax=Mexp_hi),fill="#E6433B",alpha=0.5) +
  geom_line() +
  labs(x="Year",y="Mortality Rate") +
  ggtitle(paste0("Estimated Annual Total Mortality Rate (Combined Hazards)")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))
print(plt_Mr)

# Minimum mortality rate
mcmc_areas(fit$draws(variables = c("M_min")),
           area_method="equal height",
           adjust = 2, prob = 0.8) + 
  ggtitle(paste0("Posterior distribution, Minimum Mortality Rate")) +
  labs(x="Parameter value",y="Posterior density") +
  scale_y_discrete(labels=V_names) +
  theme_classic()
median(fit$draws(variables="M_min"))




# Manuscript Figures

library(ggpattern)
library(ggplot2)
library(grid)
library(png)
library(patchwork)
library(jpeg)
library(lubridate)
library(tidyverse)

whalepic <- readPNG(paste0(getwd(),"/GrayWhale-lightened.png"))
g <- rasterGrob(whalepic, interpolate=TRUE)


Abund_Trnd = ggplot(df_Trends,aes(x=Year,y=N_exp)) +
  geom_ribbon(aes(ymin=N_exp_lo,ymax=N_exp_hi),fill="gray50",alpha=0.5) +
  geom_line() +
  geom_point(aes(x=Year,y=Surveys)) +
  geom_errorbar(aes(x=Year,ymin=Surveys-1.96*Surveys_se,
                    ymax=Surveys+1.96*Surveys_se)) +
  labs(x="Year",y="Estimated Abundance") +
  xlim(c(1967,2023))+
    ggtitle(paste0("A) Abundance")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))+
  annotation_custom(g, xmin=-Inf, xmax=1990, ymin=-Inf, ymax=60000)


Birth_Trnd = ggplot(df_Brate,aes(x=Year,y=Bexp)) +
  geom_ribbon(aes(ymin=Bexp_lo,ymax=Bexp_hi),fill="#5A8C47",alpha=0.5) +
  geom_line() +
  geom_point(aes(x=Year,y=Surveys)) +
  geom_errorbar(aes(x=Year,ymin=Surveys-1.96*Surveys_se,
                    ymax=Surveys+1.96*Surveys_se)) +
  labs(x="Year",y="Estimated Births") +
  xlim(c(1968,2022))+
  ylim(c(0,3500))+
  ggtitle(paste0("B) Births")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))


Mort_Trnd = ggplot(df_Mort,aes(x=Year,y=D_est_mn,group=Source)) +
  geom_ribbon(aes(ymin=D_est_lo,ymax=D_est_hi,fill=Source,alpha=Source)) +
  scale_fill_manual(values=c("#E6433B","#E6843B"))+
  scale_alpha_manual(values=c(0.4,0.8))+
  geom_line() +
  labs(x="Year",y="Estimated Mortality") +
  xlim(c(1968,2022))+
  ylim(c(0,3500))+
  ggtitle(paste0("C) Deaths")) +
  theme_classic() +
  theme(legend.position=c(0.8,0.9),legend.key.size=unit(0.55,'cm'))+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))

SCond <- ggplot(df_BCond,aes(x=Year,y=C_exp)) +
  geom_hline(yintercept=0.175,lty=2,color="gray50") +
  geom_ribbon(aes(ymin=C_exp_lo,ymax=C_exp_hi),fill="#3F83CF",alpha=0.6) +
  geom_line() +
  geom_point(aes(x=Year,y=Obs_C)) +
  geom_errorbar(aes(x=Year,ymin=Obs_C-Obs_C_se,ymax=Obs_C+Obs_C_se)) +
  labs(x="Year",y="Estimated Condition") +
  ylim(c(0.13,0.21))+
  xlim(c(1968,2022))+
  ggtitle(paste0("D) Body Condition (Southbound)")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))

NCond <- ggplot(df_NBCond,aes(x=Year,y=CN_exp)) +
  geom_hline(yintercept=0.175,lty=2,color="gray50") +
  geom_ribbon(data=df_NBCond,aes(ymin=CN_exp_lo,ymax=CN_exp_hi),fill="#3F83CF",alpha=0.25) +
  geom_line(data=df_NBCond,aes(x=Year,y=CN_exp),color="black",inherit.aes = F) +
  geom_point(data=NorthBC,aes(x=Year,y=BC_mn),color="black",inherit.aes=F) +
  geom_errorbar(data=NorthBC,aes(x=Year,ymin=BC_mn-BC_sd,ymax=BC_mn+BC_sd),color="black",inherit.aes=F) +
  xlim(c(1968,2022))+
  ylim(c(0.13,0.21))+
  labs(x="Year",y="Estimated Condition") +
  ggtitle(paste0("E) Body Condition (Northbound)")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))



# Figure 1
Abund_Trnd / (Birth_Trnd + Mort_Trnd) / ( SCond + NCond) + plot_layout(heights=c(3,1,1))

# Figure 1 alternative orientation (stacked)
Abund_Trnd / Birth_Trnd / Mort_Trnd / SCond / NCond + plot_layout(heights=c(2,1,1,1,1))



# 
Annual_K = ggplot(df_Kinst,aes(x=Year,y=K_inst)) +
  geom_hline(yintercept=25000,lty=2,color="gray50") +
  geom_hline(yintercept=10000,lty=3,color="gray50") +
  geom_ribbon(aes(ymin=K_inst_lo,ymax=K_inst_hi),fill="gray50",alpha=0.5) +
  geom_line() +
  ylim(c(7000,160000))+
  labs(x="Year",y="Number of Whales") +
  xlim(c(1968,2022))+
  ggtitle(expression(paste("A) Annual Carrying Capacity (", italic(K[t]),")"))) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))

Crust_Cov = ggplot(df_X1,aes(x=Year,y=X1_exp)) +
  geom_hline(yintercept=0,lty=2,color="gray50") +
  geom_ribbon(aes(ymin=X1_exp_lo,ymax=X1_exp_hi),fill="#fe987b",alpha=0.4) +
  geom_line() +
  labs(x="Year",y="Scaled Anomaly") +
  xlim(c(1968,2022))+
  ylim(c(-3,3.5))+
  ggtitle(paste0("C) Crustacean Biomass")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))

Ice_Cov = ggplot(df_X2,aes(x=Year,y=X2_exp)) +
  geom_hline(yintercept=0,lty=2,color="gray50") +
  geom_ribbon(aes(ymin=X2_exp_lo,ymax=X2_exp_hi),fill="#a2d2df",alpha=0.5) +
  geom_line() +
  labs(x="Year",y="Scaled Anomaly") +
  xlim(c(1968,2022))+
  ylim(c(-3,3.5))+
  ggtitle(paste0("B) Ice Access")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))


Benthos <- read.csv("AllStation_BenthicInfauna_CLASS_ALL3_082522_Grebmeier.csv",na.strings=c(-1000,-999.99))

Chirikov <- Benthos %>% filter(Latitude > 63.3 & Latitude < 68.4 & Longitude > -171.4 & Longitude < -165.5) %>%
  mutate(DBO_Station_Name = "Chirikov")

BenthDecline <- ggplot(Chirikov,aes(x=abd_Crustacea,y=gc_Crustacea,color=as.factor(DataYear))) + 
  geom_smooth(method="lm", formula=y ~ x+0, se=F) + 
  geom_point() + 
  scale_color_viridis_d(option="magma",guide = guide_coloursteps(),name="Year",labels=c("1970",rep("",7),"1988",rep("",7),"2003",rep("",6),"2010",rep("",7),"2019")) +
  ylab("Crustacean Biomass (gC)") + xlab("Crustacean Abundance") +
  ggtitle("D) Benthic Prey Per-Capita Biomass") +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))


#Figure 2:

#(Annual_K + BenthDecline) / (Ice_Cov + Crust_Cov )
(Annual_K + Ice_Cov) / (Crust_Cov + BenthDecline)


#Supplemental Figures

#Beta posteriors
mcmc_areas(fit$draws(variables = c("Beta")),
           area_method="equal height",
           adjust = 2, prob = 0.9) + 
  ggtitle(paste0("Covariate Effects on K")) +
  labs(x="Estimated Covariate Effect",y="Posterior density") +
  #scale_y_discrete(labels=c("Crustacean Biomass","Ice Access","Zooplankton Density")) +
  scale_y_discrete(labels=c("Crustacean Biomass","Ice Access")) +
  theme_classic() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))

#Density dependent effects on vital rates
mcmc_areas(fit$draws(variables = c("phi_star","alpha_star","gamma_star")),
           area_method="equal height",
           prob = 0.9) + 
  ggtitle(paste0("Density-dependent Effects on Vital Rates")) +
  labs(x="Proportional Decrease, 0% -> 100% K", y="Posterior density") +
  scale_y_discrete(labels=c("Body condition",
                            "Calf production",
                            "Survival rate")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))

# Stranding detection proportion
mcmc_areas(fit$draws(variables = c("Ppn_dtct_N")),
           area_method="equal height",
           adjust = 2, prob = 0.9) + 
  ggtitle(paste0("Stranding Detection Proportion")) +
  labs(x="Proportion of Deaths Detected as Strandings",y="Posterior density") +
  scale_y_discrete(labels=c("Pre-1990","1990-2000","Post-2000")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))

# Birth rate
Brexp = sumstats[which(startsWith(vn,"B[")),1]
Brexp_lo = sumstats[which(startsWith(vn,"B[")),4]
Brexp_hi = sumstats[which(startsWith(vn,"B[")),8]
df_Brater = data.frame(Year = Years[1:Nyrm1],Bexp=Brexp,
                       Bexp_lo=Brexp_lo,Bexp_hi=Brexp_hi)
colnames(df_Brater)[2:4] <- c("Bexp","Bexp_lo","Bexp_hi")

plt_Br = ggplot(df_Brater,aes(x=Year,y=Bexp)) +
  geom_ribbon(aes(ymin=Bexp_lo,ymax=Bexp_hi),fill="#5A8C47",alpha=0.5) +
  geom_line() +
  labs(x="Year",y="Birth Rate") +
  ggtitle(paste0("Estimated Birth Rate")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))
print(plt_Br)

# Realized mortality rate (combined hazards)
Mexp = sumstats[which(startsWith(vn,"M[")),1]
Mexp_lo = sumstats[which(startsWith(vn,"M[")),4]
Mexp_hi = sumstats[which(startsWith(vn,"M[")),8]
df_Mrate = data.frame(Year = Years[1:Nyrm1],Mexp=Mexp,
                      Mexp_lo=Mexp_lo,Mexp_hi=Mexp_hi)
colnames(df_Mrate)[2:4] <- c("Mexp","Mexp_lo","Mexp_hi")

plt_Mr = ggplot(df_Mrate,aes(x=Year,y=Mexp)) +
  geom_ribbon(aes(ymin=Mexp_lo,ymax=Mexp_hi),fill="#E6433B",alpha=0.5) +
  geom_line() +
  labs(x="Year",y="Mortality Rate") +
  ggtitle(paste0("Estimated Annual Total Mortality Rate (Combined Hazards)")) +
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))
print(plt_Mr)


# Density dependent effects on growth rate
df_Trends <- df_Trends %>% mutate(lambda=N_exp/lag(N_exp,1))

DensDep <- data.frame(N_exp=df_Trends$N_exp[1:Nyrm1], LogLambda = log(df_Trends$lambda[2:Nyrs]))

ggplot(DensDep,aes(x=N_exp,y=LogLambda)) +
  geom_point() +
  geom_smooth(method='lm',color='black')+
  xlab("Median Estimated Annual Abundance") +
  ylab("Log(Median Realized Annual Growth Rate)")+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title=element_text(size=16))
