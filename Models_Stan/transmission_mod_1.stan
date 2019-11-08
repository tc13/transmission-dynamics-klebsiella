//Klebsiella transmission model 1
//implemented in rstan. Author Thomas Crellen

data{
  int<lower=1> N;                         //Number of observations
  int<lower=1> K;                         //Number of intervals
  int<lower=1> L[K];                      //length of each interval
  int<lower=0,upper=1> outcome[K];        //outcome for each cluster
  int<lower=0> patients_colonised[N];     //number of patients colonised with ST
}
parameters{
  real<lower=0,upper=1> alpha;
}
model{
  vector[K] p;                              //vector for interval probability (n=402)
  vector[N] s;                              //vector for day probability (n=817)
  int counter;                              //set up counter
  
  //priors (weakly informative)
  alpha ~ beta(2, 6);
  
  //build likelihood function
  for (i in 1:N) {
    //store values for each day
    real day_prob;
    real day_prob_inv;
    
    //get odds per day from linear predictors
    day_prob= alpha;
    
    day_prob_inv = 1-day_prob;                      //1-probability 
    s[i] = day_prob_inv;                            //store daily 1-probablity in s vector
  }
  
  counter = 0;                                        //start counter at zero
  for(k in 1:K){
    vector[(L[k])] interval_temp;                   //vector with length of interval k
    for(j in 1:(L[k])){                             //for day in interval
      interval_temp[j] = s[(counter+j)];          //daily probability is the jth element of interval_prob 
    }
    //multiply the inverse probability by taking the sum of the logs
    p[k] = 1-exp(sum(log(interval_temp)));
    counter=counter+(L[k]);                         //increase counter
  }
  outcome ~ bernoulli(p);                             //regression
}
generated quantities{                                   //stores probabilites from the model                          
  vector[K] interval_p;                                 //probability of colonisation per interval
  vector[N] day_p;                                      //probability of colonisation per day
  vector[K] log_lik;                                    //stores log-likelihood of the model    
  
  int counter;                                        
  for (i in 1:N) {
    real day_prob;
    real day_prob_inv;
    day_prob = alpha;
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