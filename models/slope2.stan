// This is a Stan implementation of the bbsBayes slope model

data {
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;

  int<lower=0> count[ncounts];              // count observations
  int<lower=1> strat[ncounts];               // strata indicators
  int<lower=1> year[ncounts]; // year index
  int<lower=1> fixedyear; // centering value for years
  
  int<lower=0> firstyr[ncounts]; // first year index
  
  int<lower=1> obser[ncounts];              // observer indicators
  int<lower=1> nobservers;
  
}

parameters {
  vector[ncounts] noise_raw;             //non centered observation-level random effect to model over-dispersion
  real lambda[ncounts];             // Poisson means
  
  vector[nstrata] beta_p; //non centered random effect slopes
  real BETA; //Slope hyperparameter

  vector[nstrata] strata_p; //non centered random effect intercepts
  real STRATA; //intercept hypterparameter

  real eta; //first-year observer intercept (observer startup effect)
  
  matrix[nstrata,nyears] yeareffect_raw; //random year-effects modeling departures from slope-based trajectory

  vector[nobservers] obs_raw;    // sd of year effects

  real<lower=0> sdnoise;    // sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta;    // sd of slopes 
  real<lower=0> sdstrata;    // sd of intercepts
  real<lower=0> sdyear[nstrata];    // sd of year effects

  
}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  vector[ncounts] noise;           // extra-Poisson log-normal variance
  vector[nstrata] beta;
  vector[nstrata] strata;
  vector[nobservers] obs; //observer effects
  matrix[nstrata,nyears] yeareffect;

     obs = sdobs*obs_raw;
     noise = sdnoise*noise_raw;
      

// covariate effect on intercepts and slopes
  beta = (sdbeta*beta_p) + BETA;
  strata = (sdstrata*strata_p) + STRATA;
  
   for(s in 1:nstrata){
     yeareffect[s,] = sdyear[s]*yeareffect_raw[s,];
   }
  

  for(i in 1:ncounts){
    E[i] =  beta[strat[i]] * (year[i]-fixedyear) + strata[strat[i]] + yeareffect[strat[i],year[i]] + obs[obser[i]] + eta*firstyr[i] + noise[i];
  }
  
  }
  
model {

  sdnoise ~ normal(0,1); //prior on scale of extra Poisson log-normal variance
  noise_raw ~ normal(0,1); //non centered prior normal tailed extra Poisson log-normal variance
  
  sdobs ~ normal(0,1); //prior on sd of gam hyperparameters
  sdyear ~ normal(0,1); // prior on sd of yeareffects - stratum specific
  obs_raw ~ normal(0,1); //non centered prior on observer effects
 
  
  //nu ~ gamma(2,0.1); // alternate prior on df for t-distribution of heavy tailed 
 for(s in 1:nstrata){

  yeareffect_raw[s,] ~ normal(0,1);
  sum(yeareffect_raw[s,]) ~ normal(0,0.001*nyears);// soft sum to zero constraint on year-effects within strata
  
 }
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope - hyperparameter
  STRATA ~ normal(0,1);// prior on fixed effect mean intercept - hyperparameter
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  sdstrata ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.1); //prior on sd of slope variation

  beta_p ~ normal(0,1); //non centered prior on stratum-level slopes
  strata_p ~ normal(0,1); //non centered prior on stratum-level slopes

  //sum to zero constraints
  sum(strata_p) ~ normal(0,0.001*nstrata);
 
}
