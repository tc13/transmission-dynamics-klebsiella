//Klebsiella risk factor model C
//implemented in rstan, author: Thomas Crellen

data{
  int<lower=1> N;                                     //number of patient days (observations)
  int<lower=1> K;                                     //number of intervals
  int<lower=1> L[K];                                  //interval lengths
  int<lower=0,upper=1> outcome[K];                    //outcome for interval (n=402)
  int<lower=0,upper=1> sex[N];                        //sex (male=1)
  int<lower=0> age_entry[N];                          //age in days at entry
  int<lower=0,upper=1> probio[N];                     //probio (taken=1)
  int<lower=0,upper=1> breast_milk[N];                //breast_milk (taken=1)
  int<lower=0,upper=1> severe[N];                     //severe (yes=1, defined in Turner et al. 2016 PIDJ)
  int<lower=0,upper=1> premature[N];                  //premature (yes=1, defined in Turner et al. 2016 PIDJ)
  int<lower=0,upper=1> e_coli[N];                     //colonised with e_coli
  int<lower=0> nurses[N];                             //nurses in NU per day
  int<lower=0,upper=1> ampicillin_gentamicin_96[N];      //taken ampicillin + gentamicin in past 48hrs/96hrs
  int<lower=0,upper=1> ampicillin_mono_96[N];            //taken ampicillin in past 48hrs/96hrs
  int<lower=0,upper=1> cloxacillin_oral_96[N];           //taken oral cloxacillin in past 48hrs/96hrs
  int<lower=0,upper=1> ceftriaxone_mono_96[N];           //taken ceftriaxone in past 48hrs/96hrs
  int<lower=0,upper=1> cloxacillin_gentamicin_96[N];     //taken cloxacillin + gentamicin in past 48hrs/96hrs
  int<lower=0,upper=1> imipenem_mono_96[N];              //taken imipenem in past 48hrs/96hrs
  int<lower=0> foi[N];                                   //number of other patients colonised with ESBL K. pneumoniae
}

parameters{
  real alpha;                                         //intercept
  real beta[15];                                      //slopes for covariates 
}

model{
  vector[K] p;                                        //vector for interval probability (n=402)
  vector[N] s;                                        //vector for day probability (n=817)
  int counter;                                        //set up counter
  //priors (weakly informative)
  alpha ~ normal(0, 10);
  beta ~ normal(0, 5);
  //build likelihood function
  for (i in 1:N) {
    //store values for each day
    real day_odds;
    real day_prob;
    real day_prob_inv;
    
    //get odds per day from linear predictors
    day_odds = exp(alpha + beta[1]*sex[i] + beta[2]*age_entry[i] + beta[3]*probio[i] + beta[4]*breast_milk[i] + beta[5]*severe[i] + beta[6]*premature[i] + beta[7]*e_coli[i] + beta[8]*nurses[i] + beta[9]*ampicillin_mono_96[i] + beta[10]*ampicillin_gentamicin_96[i] + beta[11]*ceftriaxone_mono_96[i] + beta[12]*cloxacillin_oral_96[i] + beta[13]*cloxacillin_gentamicin_96[i] + beta[14]*imipenem_mono_96[i] + beta[15]*foi[i]);
    
    day_prob = day_odds / (day_odds + 1);               //convert odds to probability
    day_prob_inv = 1-day_prob;                          //1-probability 
    s[i] = day_prob_inv;                                //store daily 1-probablity in s vector
  }
  
  counter = 0;                                          //start counter at zero
  for(k in 1:K){
    vector[(L[k])] interval_temp;                       //vector with length of interval k
    for(j in 1:(L[k])){                                 //for day in interval
      interval_temp[j] = s[(counter+j)];                //daily probability is the jth element of interval_prob 
    }
    //multiply the inverse probability by taking the sum of the logs
    p[k] = 1-exp(sum(log(interval_temp)));
    counter=counter+(L[k]);                             //increase counter
  }
  outcome ~ bernoulli(p);                               //regression
}

generated quantities{                                   //stores probabilites from the model                          
  vector[K] interval_p;                                 //probability of colonisation per interval
  vector[N] day_p;                                      //probability of colonisation per day
  vector[K] log_lik;                                    //stores log-likelihood of the model    
  
  int counter;                                        
  for (i in 1:N) {
    real day_odds;
    real day_prob;
    real day_prob_inv;
    
    day_odds = exp(alpha + beta[1]*sex[i] + beta[2]*age_entry[i] + beta[3]*probio[i] + beta[4]*breast_milk[i] + beta[5]*severe[i] + beta[6]*premature[i] + beta[7]*e_coli[i] + beta[8]*nurses[i] + beta[9]*ampicillin_mono_96[i] + beta[10]*ampicillin_gentamicin_96[i] + beta[11]*ceftriaxone_mono_96[i] + beta[12]*cloxacillin_oral_96[i] + beta[13]*cloxacillin_gentamicin_96[i] + beta[14]*imipenem_mono_96[i] + beta[15]*foi[i]);    
    
    day_prob = day_odds / (day_odds + 1);         
    day_prob_inv = 1-day_prob;                    
    day_p[i] = day_prob_inv;                      
  }
  counter = 0;                                      
  for(k in 1:K){
    vector[(L[k])] interval_temp;                 
    for(j in 1:(L[k])){                           
      interval_temp[j] = day_p[(counter+j)];    
    }
    interval_p[k] = 1-exp(sum(log(interval_temp)));       
    counter=counter+(L[k]);                             
  }
  for(i in 1:K){
    log_lik[i] = bernoulli_lpmf(outcome[i] | interval_p[i]);
  }
}