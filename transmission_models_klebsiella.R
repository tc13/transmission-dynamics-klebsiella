#Transmission Models for Acquisition of ESBL K. pneumoniae among Cambodian Neonates
#Linear models, fitted with HMCMC in Stan
#Author: Thomas Crellen, Postdoctoral Researcher, University of Oxford, thomas.crellen@ndm.ox.ac.uk

#Depending on your computing resources, allow approximately 24 hours for this script to run
#The final RDATA file will be approximately 1Gb in size

#Required packages
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
#note: ensure rstan is correctly configured for your system see https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
#install.packages(c("loo", "lubridate", "reshape2", "ggplot2", "bayesplot"), dependencies = TRUE)

#clear global environment
rm(list=ls())

#load libraries
require(rstan) 
require(loo)
require(lubridate)
require(reshape2)
require(ggplot2)
require(bayesplot)

#rstan setting
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

#Set working directory as this folder - only works in RStudio
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#set seed
set.seed(1234)

#Read in dataset for transmission models
#each line refers to a patient day where an infant is at risk of acquiring a specific sequence type (ST) of ESBL K. pneumoniae
#patients colonised at entry and intervals where no outcome is recorded have been removed
#patients are present in the data until they become colonised, at which point they are no longer at risk
#ST information is obtained from the assemblies of the K. pnuemoniae genomes (see methods in paper)

tm <- read.csv("klebsiella_ST_transmission_reproducible_covariates.csv", header=T)

#index each patient interval
counter <- 0
interval_index <- c()
id <- ""
interval <- 0
for(i in 1:length(tm$id)){
  id_new <- tm$id[i]
  interval_new <- tm$swab_interval[i] 
  if(id_new!=id | interval_new!=interval){
    counter <- (counter+1)
    id <- id_new
    interval <- interval_new
  }
  interval_index[i] <- counter
}

#Outcome data.frame
tm.swabs <- tm[tm$rectal.swab==1,]

#get numeric ID for each ST
ST.names.all <- unique(tm$ST)
tm$ST_numeric <- match(tm$ST, ST.names.all)

#modify coavariate
tm$sex <- ifelse(tm$sex=="Male",1,0)

#list of variables for stan input
stan_tm_input <- list(
  N = nrow(tm),
  K = max(interval_index),
  L = tm.swabs$interval_length,
  outcome = tm.swabs$first_colonised,
  patients_colonised = tm$patients_colonised,
  patients_colonised_lag = tm$patients_colonised_lag,
  N_ST = length(unique(tm$ST)),
  ST = tm$ST_numeric,
  #additional covariates
  sex= tm$sex,
  probio= tm$probiotic,
  breast_milk= tm$breast_milk,
  premature= tm$premature,
  e_coli= tm$e_coli,
  nurses= tm$nurses,
  ampicillin_gentamicin_96= tm$ampicillin_gentamicin_96hrs,
  ceftriaxone_mono_96= tm$ceftriaxone_mono_96)

#Additional variables for transmission model 3 to account for the number of patients colonised on previous days
tm$date <- ymd(tm$date)
dates_windows <- seq.Date(from=min(tm$date), to=max(tm$date), by = "1 day")
stan_tm_input$N_days <- length(dates_windows)
stan_tm_input$day <- match(tm$date, dates_windows) 

#matrix where rows are STs, columns are dates, and cells are the number colonized
colonized_matrix <- matrix(nrow=length(ST.names.all), ncol=length(dates_windows), data=0)
rownames(colonized_matrix) <- ST.names.all
colnames <- dates_windows

#fill up colonized matrix
for(i in 1:length(ST.names.all)){
  st <- ST.names.all[i]
  df <- tm[tm$ST==st&tm$colonised==0,]
  index <- match(dates_windows, df$date)
  colonized <- df$patients_colonised[index]
  colonized_zero <- ifelse(is.na(colonized), 0,colonized)
  colonized_matrix[i,] <- colonized_zero
}

#add values to stan input list
stan_tm_input$colonised_matrix <- colonized_matrix

##############################
## Initialisation functions ##
##############################

#Functions to randomly assign initial values within a valid range
n_chains <- 4

init_fn_1 <- function(chain_id) {
  list(alpha=runif(1, 0.001,0.1))
}

init_fn_2 <- function(chain_id) {
  list(alpha=runif(1, 0.001, 0.1), beta=runif(1, 0.01, 0.05))
}

#initalise function for mod 2 with covariates
init_fn_2_cov <- function(chain_id) {
  list(alpha=runif(1, 0.001, 0.1), beta=runif(4, 0.001, 0.025))
}

init_fn_3 <- function(chain_id) {
  list(alpha=runif(1, 0.001,0.1), beta=runif(1, 0.001,0.05), gamma=runif(1, 0.001,0.05), 
       lambda=runif(1, 0.1, 8))
}

init_fn_4 <- function(chain_id) {
  list(alpha=runif(1, 0.001,0.1), beta=runif(stan_tm_input$N_ST, 0.0001,0.05), 
       scale=runif(1, 0.1, 2), location=runif(1, 0.1, 2))
}

init_fn_5 <- function(chain_id) {
  list(beta=runif(1, 0.001,0.1), alpha=runif(stan_tm_input$N_ST, 0.0001,0.05), 
       scale=runif(1, 0.1, 2), location=runif(1, 0.1, 2))
}

#then assign initial values to lists
init_ll_1 <- lapply(1:n_chains, function(id) init_fn_1(chain_id = id))
init_ll_2 <- lapply(1:n_chains, function(id) init_fn_2(chain_id = id))
init_ll_2_cov <- lapply(1:n_chains, function(id) init_fn_2_cov(chain_id = id))
init_ll_3 <- lapply(1:n_chains, function(id) init_fn_3(chain_id = id))
init_ll_4 <- lapply(1:n_chains, function(id) init_fn_4(chain_id = id))
init_ll_5 <- lapply(1:n_chains, function(id) init_fn_5(chain_id = id))

#########################
## Running stan models ##
#########################

#compile models from .stan files
tm.mod.1 <- stan_model(file="transmission_mod_1.stan", model_name = "tm1")

tm.mod.2 <- stan_model(file="transmission_mod_2.stan", model_name = "tm2")

#model 2 with additional covariates
tm.mod.2.cov <- stan_model(file="transmission_mod_2_covariates.stan", model_name="tm2cov")

tm.mod.3 <- stan_model(file="transmission_mod_3.stan", model_name = "tm3")

#model 3 with alternate priors on lambda
tm.mod.3p <- stan_model(file="transmission_mod_3_priors.stan", model_name = "tm3p")

tm.mod.4 <- stan_model(file= "transmission_mod_4.stan", model_name = "tm4")

tm.mod.5 <- stan_model(file="transmission_mod_5.stan", model_name = "tm5")

#Sampling. Assumes a four core processor, change `cores` argument if fewer cores are available 
tm.fit.1 <- sampling(tm.mod.1, data=stan_tm_input, seed=1,
               iter=4000,warmup=2000,chains=4, thin=6, cores=4,
               pars=c("alpha","log_lik"), init = init_ll_1,
               control=list(adapt_delta = 0.95))

tm.fit.2  <- sampling(tm.mod.2, data=stan_tm_input, seed=2,
               iter=4000,warmup=2000,chains=4, thin=6, cores=4,
               pars=c("alpha", "beta", "log_lik"), init= init_ll_2,
               control=list(adapt_delta = 0.95))

#model 2 with additional covariates
tm.fit.2.cov  <- sampling(tm.mod.2.cov, data=stan_tm_input, seed=9,
                      iter=4000,warmup=2000,chains=4, thin=6, cores=4,
                      pars=c("alpha", "beta", "log_lik"), init= init_ll_2_cov,
                      control=list(adapt_delta = 0.95))

tm.fit.3 <- sampling(tm.mod.3, data=stan_tm_input, seed=3,
               iter=4000,warmup=2000,chains=4,thin=6, cores=4,
               pars=c("alpha", "beta", "gamma", "lambda", "log_lik"), 
               init = init_ll_3, control=list(adapt_delta = 0.95))

tm.fit.3p <- sampling(tm.mod.3p, data=stan_tm_input, seed=3,
                     iter=4000,warmup=2000,chains=4,thin=6, cores=4,
                     pars=c("alpha", "beta", "gamma", "lambda", "log_lik"), 
                     init = init_ll_3, control=list(adapt_delta = 0.90))

tm.fit.4 <- sampling(tm.mod.4, data=stan_tm_input, seed=4,
                 iter=16000,warmup=8000,chains=4,thin=12, cores=4,
                 pars=c("alpha", "beta", "scale", "location", "log_lik"),
                 init = init_ll_4, control=list(adapt_delta = 0.85))

tm.fit.5 <- sampling(tm.mod.5, data=stan_tm_input, seed=4,
                     iter=16000,warmup=8000,chains=4,thin=12, cores=4,
                     pars=c("alpha", "beta", "scale", "location", "log_lik"),
                     init = init_ll_5, control=list(adapt_delta = 0.85))

#############################
## Convergence Diagnostics ##
#############################

#Examine convergence metrics for Tail ESS, Bulk ESS and R-hat
mon.1 <- monitor(tm.fit.1)
summary(mon.1[1, 'Tail_ESS'])
summary(mon.1[1, 'Bulk_ESS'])
summary(mon.1[1, 'Rhat'])

mon.2 <- monitor(tm.fit.2)
summary(mon.2[1:2, 'Tail_ESS'])
summary(mon.2[1:2, 'Bulk_ESS'])
summary(mon.2[1:2, 'Rhat'])

mon.2.cov <- monitor(tm.fit.2.cov)
summary(mon.2.cov[1:5, 'Tail_ESS'])
summary(mon.2.cov[1:5, 'Bulk_ESS'])
summary(mon.2.cov[1:5, 'Rhat'])

mon.3 <- monitor(tm.fit.3)
summary(mon.3[1:4, 'Tail_ESS'])
summary(mon.3[1:4, 'Bulk_ESS'])
summary(mon.3[1:4, 'Rhat'])

mon.3p <- monitor(tm.fit.3p)
summary(mon.3p[1:4, 'Tail_ESS'])
summary(mon.3p[1:4, 'Bulk_ESS'])
summary(mon.3p[1:4, 'Rhat'])

mon.4 <- monitor(tm.fit.4)
summary(mon.4[1:65, 'Tail_ESS'])
summary(mon.4[1:65, 'Bulk_ESS'])
summary(mon.4[1:65, 'Rhat'])

mon.5 <- monitor(tm.fit.5)
summary(mon.5[1:65, 'Tail_ESS'])
summary(mon.5[1:65, 'Bulk_ESS'])
summary(mon.5[1:65, 'Rhat'])

######################
## Posterior chains ##
######################

#extract posterior as an array (retains chain-specific information)
post.1.chains <- as.array(tm.fit.1, pars=c("alpha"))
post.2.chains <- as.array(tm.fit.2, pars=c("alpha", "beta"))
post.2.cov.chains <- as.array(tm.fit.2.cov, pars=c("alpha", "beta"))
post.3.chains <- as.array(tm.fit.3, pars=c("alpha", "beta", "gamma", "lambda"))
post.4.chains <- as.array(tm.fit.4, pars=c("alpha", "beta", "scale", "location"))
post.5.chains <- as.array(tm.fit.5, pars=c("alpha", "beta", "scale", "location"))

#theme set for bayesplot
theme_set(bayesplot::theme_default(base_family = "sans"))
color_scheme_set("viridis")

#Plot chains
mcmc_trace(post.1.chains, 
           pars = c("alpha"), 
           facet_args = list(ncol=1, strip.position="left"))

mcmc_trace(post.2.chains, 
           pars = c("alpha", "beta"), 
           facet_args = list(ncol=2, strip.position="left"))

mcmc_trace(post.2.cov.chains,
           pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]"), 
           facet_args = list(ncol=2, nrow=2, strip.position="left"))



mcmc_trace(post.3.chains, 
           pars = c("alpha", "beta", "gamma", "lambda"), 
           facet_args = list(ncol=2, nrow=2, strip.position="left"))

mcmc_trace(post.4.chains, 
           pars = c("alpha", "scale", "location"), 
           facet_args = list(ncol=3, strip.position="left"))

mcmc_trace(post.5.chains, 
           pars = c("beta", "scale", "location"), 
           facet_args = list(ncol=3, strip.position="left"))

################################
## extract posteriors as list ##
################################

post.1 <- extract(tm.fit.1)
post.2 <- extract(tm.fit.2)
post.2.cov <- extract(tm.fit.2.cov, pars=c("alpha", "beta"))
post.3 <- extract(tm.fit.3, pars=c("alpha", "beta", "gamma", "lambda"))
post.3p <- extract(tm.fit.3p ,pars=c("alpha", "beta", "gamma", "lambda"))
post.4 <- extract(tm.fit.4)
post.5 <- extract(tm.fit.5)

#############################
## WAIC from fitted models ##
#############################

waic.1 <- waic(post.1$log_lik)
waic.2 <- waic(post.2$log_lik)
waic.3 <- waic(post.3$log_lik)
waic.3p <- waic(post.3$log_lik)
waic.4 <- waic(post.4$log_lik)
waic.5 <- waic(post.5$log_lik)

compare(waic.1, waic.2, waic.3, waic.3p, waic.4, waic.5)

##############################
## Parameter Interpretation ##
##############################

#function to obtain median and 95% credible interval from a vector
CrI <- function(x) signif(quantile(x, probs = c(0.025, 0.5, 0.975)),2)

#values for table 3
CrI(post.1$alpha)

CrI(post.2$alpha)
CrI(post.2$beta)

CrI(post.3$alpha)
CrI(post.3$beta)
CrI(post.3$gamma)
CrI(post.3$lambda)

CrI(post.4$alpha)
CrI(post.4$scale)
CrI(post.4$location)

CrI(post.5$beta)
CrI(post.5$scale)
CrI(post.5$location)

#####################################################################
## Plotting priors and posterior together for transmission model 3 ##
#####################################################################

x.prob <- seq(0, 1, by=0.001)
x.pos <-  seq(0, 20, by=0.01)

y1 <- dbeta(x=x.prob, shape1=2, shape2=8)
y2 <- dnorm(x=x.pos, mean = 1, sd = 2)
y3 <- dnorm(x=x.pos, mean=0, sd=5)

###################################
#lambda estimates original prior ##
###################################

png("/Users/Tom/Google Drive/MORU/manuscript/Supp-Fig-4-prior1-lambda.png", width=7, height=5, units="in", res=280)

plot(x.pos, y2, type="l", lwd=2, ylim=c(0,0.5), xlim=c(0,10),
     xaxs="i", yaxs="i", ylab="Probability Density", 
     xlab="Parameter Estimate", col="red")

lines(density(post.3$lambda, adjust=2), lwd=2, col="blue")

abline(v=median(post.3$lambda), col="blue",lwd=2, lty=2)
#title and legend
title(main = expression(paste("Transmission Model 3, environmental decay (", lambda, ") prior normal(", mu, "=1, ", sigma, "=2)")))

legend(7, 0.45, legend=c("Prior", "Posterior", "Posterior Median"),
       col=c("red", "blue", "blue"), lty=c(1,1,2), cex=0.8)

dev.off()

##############################
#lambda estimates new prior ##
##############################

png("/Users/Tom/Google Drive/MORU/manuscript/Supp-Fig-4-prior2-lambda.png", width=7, height=5, units="in", res=280)

plot(x.pos, y3, type="l", lwd=2, ylim=c(0,0.5), xlim=c(0,10),
     xaxs="i", yaxs="i", ylab="Probability Density", 
     xlab="Parameter Estimate", col="red")

lines(density(post.3p$lambda, adjust=2), lwd=2, col="blue")

abline(v=median(post.3p$lambda), col="blue",lwd=2, lty=2)
#title and legend
title(main = expression(paste("Transmission Model 3, environmental decay (", lambda, ") prior normal(", mu, "=0, ", sigma, "=5)")))

legend(7, 0.45, legend=c("Prior", "Posterior", "Posterior Median"),
       col=c("red", "blue", "blue"), lty=c(1,1,2), cex=0.8)
 
dev.off()

###################################################
#comparing estimates for environmental half-life ##
###################################################

png("/Users/Tom/Google Drive/MORU/manuscript/Supp-Fig-4-half-life.png", width=7, height=5, units="in", res=280)

plot(density(log(2)/post.3$lambda*24), type="l", lwd=2, ylim=c(0,0.5), xlim=c(0,15),
     xaxs="i", yaxs="i", ylab="Probability Density", 
     xlab="Environmental Half Life (Hours)", col="orange", main="")

lines(density(log(2)/post.3p$lambda*24), lwd=2, col="purple")

abline(v=median(log(2)/post.3$lambda*24), lwd=2, lty=2, col="orange")
abline(v=median(log(2)/post.3p$lambda*24), lwd=2, lty=2, col="purple")

abline(v=12, lwd=2, lty=4, col="darkgreen")
title(main = expression(paste("Transmission Model 3, environmental half-life estimates")))

legend(6, 0.45, legend=c("Posterior: Prior Normal(1, 2)", "Posterior Median", "Posterior: Prior Normal(0, 5)", "Posterior Median", "Published Estimate"),
       col=c("orange", "orange", "purple", "purple", "darkgreen"), lty=c(1,2,1,2,4), cex=0.8)

dev.off()

####################
# genomic metadata #
####################

meta <- read.table("klebsiella_genomic_metadata.txt", 
                   header=T, sep = "\t")

carriage <- meta[meta$Isolate_from=="carriage",]
ST_table <- table(carriage$ST)
ST_levels <- names(ST_table)[order(ST_table)]

#extract values from model 4
alpha <- mean(post.4$alpha) 
beta_ST <- apply(post.4$beta, 2, mean)

#predictions from model - force of infection / colonisation pressure
s <- seq(0, 6, by=1)
pred <- matrix(ncol = length(s), nrow = length(beta_ST), data=0) 

for(i in 1:length(s)){
  pred[,i]<- alpha+beta_ST*s[i]
}

#format data.frame for plotting
pred.df <- as.data.frame(pred)
colnames(pred.df) <- s
pred.df$ST <- ST.names.all

#melt to create long format data.frame
pred_long <- melt(pred.df, id.vars=c("ST"), variable.name = "x", value.name = "p")
pred_long$x <- as.numeric(as.character(pred_long$x))
pred_long$ST <- as.character(pred_long$ST)
pred_long$ST_select <- ifelse(pred_long$ST %in% rev(ST_levels)[1:4], pred_long$ST, "other")

#aethetics
col_vec <- c("#7a0177", "#225ea8", "#ec7014", "#bd0026", "gray60")
alpha_ST <- c(0.5,1,1,1,1,1)
size_ST <- c(2.25,3.5,3.5,3.5,3.5)

########################################
## Figure 3C Force of infection by ST
########################################

#ggplot theme set
theme_set(theme_bw(base_size = 28))

ggplot(data=pred_long, aes(x=x, y=p, group=ST, alpha=ST_select, size=ST_select, colour=ST_select)) + 
  geom_jitter(height=0, width=0.04) +
  geom_line(size=1)+
  scale_y_continuous(name = "Force of Infection",
                     limits=c(0,0.31), expand=c(0,0))+
  scale_x_continuous(name="Number of Infants Colonised with ST", 
                     labels = c(0,1,2,3,4,5,6),
                     breaks= c(0,1,2,3,4,5,6))+
  scale_colour_manual(name=expression(atop(italic("K. pneumoniae s.l.")," Sequence Types")), 
                      breaks=c("ST334", "ST101", "ST1074-1LV", "ST45", "other"), 
                      values=rev(col_vec), labels=c("ST 334","ST 101", "ST 1074", "ST 45", "Other"))+
  scale_alpha_manual(values=alpha_ST, guide=FALSE)+
  scale_size_manual(values=size_ST, guide=FALSE)+
  #format legend
  theme(legend.position=c(0.02,0.98), legend.justification=c(0, 1),
        legend.key.size = unit(0.8, 'cm'))
