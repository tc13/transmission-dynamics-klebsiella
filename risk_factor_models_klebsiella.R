#Risk Factor Models for Acquisition of ESBL K. pneumoniae among Cambodian Neonates
#Generalised linear models, fitted with HMCMC in Stan
#Author: Thomas Crellen, Postdoctoral Researcher, University of Oxford, thomas.crellen@ndm.ox.ac.uk

#Depending on your computing resources, allow approximately 6 hours for this script to run
#The final RDATA file will be approximately 85Mb in size

#Required packages
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#install.packages(c("loo", "lubridate"), dependencies = TRUE)

#clear global environment
rm(list=ls())

#load libraries
library(rstan) 
#note: ensure rstan is correctly configured for your system see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(loo)
library(lubridate)

#Packages for plotting
library(ggplot2)
library(bayesplot)

#rstan setting
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

#Set working directory as this folder - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#logistic function (inverse logit)
logistic <- function(x){
  odds <- exp(x)
  prob <- odds/(odds+1)
  return(prob)
}

#Read in dataset
#Each row is a patient day in the neonatal unit. Records have been removed for infants colonised at entry and swab intervals with no outcome recorded
#Records are shown for patients up to the first positive rectal swab for ESBL K. pneumoniae
kp <- read.csv("klebsiella_acquisition_reproducible.csv", header=T)

#seperate df with line per swabbing interval (nrow=402) for interval outcome
kp.interval <- kp[kp$interval_day==kp$interval_length,]

#make binary variables for regression
kp.interval$outcome <- with(kp.interval, ifelse(k_pneumoniae.colonised==1,1,0))
kp$sex <- with(kp, ifelse(sex=="Male",1,0))
kp$probio <- with(kp, ifelse(probio=="Yes",1,0))
kp$breast.milk <- with(kp, ifelse(breast.milk=="Yes",1,0))
kp$premature <- with(kp, ifelse(gestation.cat=="Preterm",1,0))

#index each patient interval
counter <- 0
interval_index <- c()
id <- ""
interval <- 0
for(i in 1:nrow(kp)){
  id_new <- kp$id_num[i]
  interval_new <-kp$swab_interval[i] 
  if(id_new!=id|interval_new!=interval){
    counter<-counter+1
    id <- id_new
    interval <- interval_new}
  interval_index[i] <-counter}

#Variable if colonised with E. coli in previous swab interval
ecoli_colonised_interval <- c()
id_marker <- ""
interval_outcome <- 0
interval_counter <- 0

for(i in 1:nrow(kp)){
  #update variables for new interval
  if(kp$swab_interval[i]!=interval_counter){
    interval_counter <- kp$swab_interval[i]
    interval_outcome <- kp$e_coli.colonised[i]
  }
  #update variables for new patient
  if(kp$id[i]!=id_marker){
    id_marker <- kp$id[i]
    interval_outcome <- kp$e_coli.colonised[i]
    interval_counter <- kp$swab_interval[i]
    #if colonised "on entry" then presumed colonised from admission
    if(kp$e_coli.colonised_entry[i]==1){
      interval_outcome <- 1}
  }
  #update vector with e_coli colonised values 
  ecoli_colonised_interval[i] <- interval_outcome
}

kp$ecoli_colonised_interval <- ecoli_colonised_interval

#Variables for month and year 
month_year_levels <- unique(as.factor(kp$month_year))
month_num  <- sapply(kp$month_year, function(x) match(x, month_year_levels), USE.NAMES = F)

#list of variables for stan input, 
#note the outcome variable is a different length (402) to the total number of observations (871)
stan_input <- list(
  N=nrow(kp),                          #all observations (patient days, n=871)
  K=max(interval_index),               #number of intervals (n=402)     
  index=interval_index,                #indexes each interval
  outcome= kp.interval$outcome,        #outcome from each interval (0/1),n=402
  L= kp.interval$interval_length,      #gives num of days in each interval
  #time invarying variables
  sex=kp$sex,                          #sex: male=1, female=0
  age_entry=kp$age.entry,              #age in days at first ward entry
  probio=kp$probio,                    #probiotic taken =1, not taken =0
  breast_milk=kp$breast.milk,          #breast milk fed =1, not fed =0
  severe=kp$severe,                    #any severe symptoms (see paper) =1, no symptoms =0
  premature=kp$premature,              #born prematurely (see paper) =1, born term =0
  #time varying covariates
  e_coli = ecoli_colonised_interval,   #E.coli colonised in previous swab interval =1, not colonised =0
  nurses = kp$nurses_NU,               #Number of nurses on the ward
  #antibiotics taken within 48hours
  ampicillin_gentamicin_48 = kp$ampicillin_gentamicin.48hrs, 
  ampicillin_mono_48 = kp$ampicillin_mono.48hrs, 
  cloxacillin_oral_48 = kp$cloxacillin_oral.48hrs,
  ceftriaxone_mono_48 = kp$ceftriaxone_mono.48hrs, 
  cloxacillin_gentamicin_48 = kp$cloxacillin_gentamicin.48hrs,
  imipenem_mono_48 = kp$imipenem_mono.48hrs,
  #antibiotics taken within 96 hours
  ampicillin_gentamicin_96 = kp$ampicillin_gentamicin.96hrs, 
  ampicillin_mono_96 = kp$ampicillin_mono.96hrs, 
  cloxacillin_oral_96 = kp$cloxacillin_oral.96hrs,
  ceftriaxone_mono_96 = kp$ceftriaxone_mono.96hrs, 
  cloxacillin_gentamicin_96 = kp$cloxacillin_gentamicin.96hrs,
  imipenem_mono_96 = kp$imipenem_mono.96hrs,
  #colonisation pressure from other patients
  foi = kp$k_pneumoniae.n_infected,
  #month
  month=month_num,
  N_month = length(month_year_levels))

#########################
## Running stan models ##
#########################

#compile models by reading in .stan files
rf.mod.A <- stan_model(file = "risk_factor_mod_A.stan",
                  model_name="rfA")

rf.mod.B <- stan_model(file= "risk_factor_mod_B.stan",
                      model_name="rfB")

rf.mod.C <- stan_model(file= "risk_factor_mod_C.stan",
                       model_name="rfC")

rf.mod.D <- stan_model(file = "risk_factor_mod_D.stan",
                       model_name="rfD")

#Sampling. Assumes a four core processor, change `cores` argument if fewer cores are available 

rf.fit.A <- sampling(rf.mod.A, data=stan_input, iter=5000, warmup=2000,
                seed=5, chains=4, thin=5, cores=4,
                pars=c("alpha","beta", "interval_p", "log_lik"),
                control=list(adapt_delta = 0.95))

#Run Model A
rf.fit.B <- sampling(rf.mod.B, data=stan_input, seed = 6,
                 iter=5000,warmup=2000,chains=4,thin=5, cores=4,
                 pars=c("alpha","beta", "log_lik"),
                 control=list(adapt_delta = 0.95))

#Run Model C
rf.fit.C <- sampling(rf.mod.C, data=stan_input, seed = 7,
                 iter=5000,warmup=2000,chains=4,thin=5, cores=4,
                 pars=c("alpha","beta", "log_lik"),
                 control=list(adapt_delta = 0.95))

#Run Model D
rf.fit.D <- sampling(rf.mod.D, data=stan_input, seed = 10, 
                 iter=22000,warmup=10000,chains=4,thin=15, cores=4,
                 pars=c("alpha","beta", "mu","sigma", "log_lik"),
                 control=list(adapt_delta = 0.99))

############################
## Convergence assessment ##
############################

#Examine convergence metrics for Tail ESS, Bulk ESS and R-hat
mon.A <- monitor(rf.fit.A)
summary(mon.A[1:15, 'Tail_ESS'])
summary(mon.A[1:15, 'Bulk_ESS'])
summary(mon.A[1:15, 'Rhat'])

#extract posterior as an array (retains chain-specific information)
post.A.chains <- as.array(rf.fit.A)

#theme set for bayesplot
theme_set(bayesplot::theme_default(base_family = "sans"))
color_scheme_set("viridis")

#Plot chains
mcmc_trace(post.A.chains, 
           pars = c("alpha", "beta[1]", "beta[2]", "beta[3]", 
                    "beta[4]", "beta[5]", "beta[6]", "beta[7]", 
                    "beta[8]", "beta[9]", "beta[10]", "beta[11]", 
                    "beta[12]", "beta[13]", "beta[14]"), 
           facet_args = list(ncol=3, strip.position="left"))

#####################################
## Extract posterior distributions ##
#####################################

post.A <- extract(rf.fit.A)
post.B <- extract(rf.fit.B)
post.C <- extract(rf.fit.C)
post.D <- extract(rf.fit.D)

#############################
## WAIC from fitted models ##
#############################

waic.A <- waic(post.A$log_lik)
waic.B <- waic(post.B$log_lik)
waic.C <- waic(post.C$log_lik)
waic.D <- waic(post.D$log_lik)
compare(waic.A, waic.B, waic.C, waic.D)

##############################
## Parameter Interpretation ##
##############################

#Models are fitted on the log-odds scale
#To obtain odds ratios from Model A:
kleb.odds <- apply(post.A$beta, 2, exp)

#vector of variable names
kleb.post.labels <- c("Male", "Age at Admission (Unit)", "Probiotic Taken", "Breast Milk Fed", "Severe", "Born Premature", "E. coli Colonised", "Nurses (Unit)", "Ampicillin", "Ampicillin + Gentamicin", "Ceftriaxone", "Oral Cloxacillin", "Cloxacillin + Gentamicin", "Imipenem")

#Create data.frame - continuous variables show unit increase
#odds ratios reported in results 
kleb.post.df <- data.frame(
  male = kleb.odds[,1],
  age = kleb.odds[,2],
  probiotic = kleb.odds[,3],
  breast_milk = kleb.odds[,4],
  severe = kleb.odds[,5],
  premature = kleb.odds[,6],
  e_coli = kleb.odds[,7],
  nurse = kleb.odds[,8],
  amp = kleb.odds[,9],
  amp_gent = kleb.odds[,10],
  ceft = kleb.odds[,11],
  clox = kleb.odds[,12],
  clox_gent = kleb.odds[,13],
  imip = kleb.odds[,14])

#Quantiles from posterior distribution
kleb.quant <- as.data.frame(t(apply(kleb.post.df, 2, function(x) quantile(x, probs=c(0.025, 0.1, 0.5, 0.9, 0.975)))))
colnames(kleb.quant) <- c("lower.95", "lower.80", "median", "upper.80", "upper.95")
kleb.quant$variable <- kleb.post.labels

####################################################
# Figure 2A: Posterior Distribution of covariates 
# Odds Ratio for ESBL K. pneumo acquisition
####################################################

#ggplot theme set
theme_set(theme_bw(base_size = 28))

ggplot(kleb.quant, aes(x=median, y=variable, group=variable))+
  geom_vline(xintercept = 1, linetype="dashed")+
  geom_errorbarh(aes(xmin=lower.95, xmax=upper.95), colour="black", size=1, height=0)+
  geom_errorbarh(aes(xmin=lower.80, xmax=upper.80), colour="#4575b4", size=3, height=0)+
  geom_point(size=4)+
  scale_y_discrete(limits=rev(kleb.post.labels), name="")+
  xlab("Odds Ratio for Acquisition")+
  coord_cartesian(xlim=c(0,6))+
  scale_x_continuous(expand=c(0,0))

##########################################################
## Cumulative Risk of Colonisation for Specific Infants ##
##########################################################

#5 day old, non-severe, non-premature girl. 3 nurses per day
#No probiotic, breast milk fed, not on antibiotics. No ESBL E. coli

low_risk <- logistic(post.A$alpha + post.A$beta[,1]*0 + post.A$beta[,2]*4 + post.A$beta[,3]*0+ post.A$beta[,4]*1+ post.A$beta[,5]*0+ post.A$beta[,6]*0+ post.A$beta[,7]*0+ post.A$beta[,8]*3 + post.A$beta[,9]*0+ post.A$beta[,10]*0 + post.A$beta[,11]*0 + post.A$beta[,12]*0+ post.A$beta[,13]*0 + post.A$beta[,14]*0)
lr_list <- list()

#obtain daily risk over 8 days
for(i in 1:8){
  lr_list[[i]] <- quantile((low_risk*(1-low_risk)^(i-1)), probs = c(0.1, 0.5, 0.9))}
lr_quant <- t(bind_cols(lr_list))

#obtain cumulative risk over 8 days
lr_cum_list <- list()
for(i in 1:8){
  lr_cum_list[[i]] <- quantile((1-(1-low_risk)^i), probs = c(0.1, 0.5, 0.9))}
lr_cum <- as.data.frame(t(bind_cols(lr_cum_list)))
colnames(lr_cum) <- c("lower","median","upper")
lr_cum$day <- seq(1:nrow(lr_cum))

#add zeros
lr_cum <- rbind(c(0,0,0,0), lr_cum)
lr_cum$scenario <- "1"

#same as low risk, but takes Amp-Gent
lr_amp_gent <- logistic(post.A$alpha + post.A$beta[,1]*0 + post.A$beta[,2]*4 + post.A$beta[,3]*0+ post.A$beta[,4]*1+ post.A$beta[,5]*0+ post.A$beta[,6]*0+ post.A$beta[,7]*0+ post.A$beta[,8]*3 + post.A$beta[,9]*0+ post.A$beta[,10]*1 + post.A$beta[,11]*0 + post.A$beta[,12]*0+ post.A$beta[,13]*0 + post.A$beta[,14]*0)

#if infant switches onto amp-gent on day 3
#matrix of daily probabilities
lr_amp_gent_mat <- rbind(
  low_risk, low_risk, lr_amp_gent, lr_amp_gent, lr_amp_gent, lr_amp_gent, lr_amp_gent, lr_amp_gent)

#cumulative probability
m <- t(1-apply(1-lr_amp_gent_mat, 2, FUN=cumprod))

#quantiles
m.q <- as.data.frame(t(apply(m, 2, function(x) quantile(x, probs = c(0.1, 0.5, 0.9)))))
colnames(m.q) <- c("lower", "median", "upper")
m.q$day <- seq(1:8)
m.q <- rbind(c(0,0,0,0), m.q)
m.q$scenario <- "2"

#bind 2 different scenarios
cum_risk <- rbind(lr_cum, m.q)

########################################################
## Figure 2C                                          ##
## Cumulative Risk of Acquisition under Two Scenarios ##
########################################################

ggplot(cum_risk, aes(y=median, x=day, group=scenario, 
                     colour=scenario, fill=scenario))+
  geom_point(size=5)+geom_line(size=1.5)+
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  scale_y_continuous(name="Cumulative Probability of Acquisition",
                     limits = c(0,1), expand = c(0,0))+
  scale_x_continuous(name="Days in Neonatal Unit", limits = c(0,8.2),
                     expand=c(0,0), breaks = c(0,1,2,3,4,5,6,7,8))+
  theme(legend.position = c(0.35,0.95), 
        legend.justification = c(0.9,0.9),
        legend.key.size = unit(0.8, "cm"))+
  scale_fill_viridis_d(name="Scenario", option="A", begin = 0.25, end=0.65, labels=c("No Antibiotics", "Amp + Gent"))+
  scale_color_viridis_d(name="Scenario", option="A",begin = 0.25, end=0.65, labels=c("No Antibiotics", "Amp + Gent"))

