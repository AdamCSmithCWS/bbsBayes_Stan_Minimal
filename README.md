# bbsBayes_Stan_Minimal

This is an minimum working example of an early attempt to use Stan to fit a relatively simple hierarchical model applied to data from the North American Breeding Bird Survey. I'm trying to build a Stan version of the JAGS models currently used to analyse the North American Breeding Bird Survey (BBS) monitoring data [@sauer2011][@smith2020]. The model is effectively a hierarchical, over-dispersed Poisson regression, with a number of intercept parameters to account for observer-effects, route-level variation in abundance, and random annual fluctuations around the slope-line.

This repo includes the code, model, and data for a minimal working example of one of the simplest models, using a species with pretty representative count structure and relatively few data (few strata and routes). It's currently set up to use only the most recent 20-years of data, but the continental scale monitoring program includes \> 50 years of data. The data are observations (summed counts) of individual birds observed during annual surveys of BBS "routes".

The model works. It generates reasonable estimates, no divergent transitions, effective sample sizes are generally \~80% of the total samples. If the max_treedepth argument is increased beyond 10.

## The problem

The problem is that the sampler exceeds the maximum tree depth, even if the limit is set to \~15-17 (haven't tried higher yet because life is short). The high max_treedepth limit means that the sampler takes a really long time, even for this species which has very few data compared to many in the dataset.

For example, at tree_depth 15, the n_leapfrog steps is \> 30,000. So even a limited initial run of 500 (400 warmup) requires \~10 hours. If max_treedepth at 17, it's closer to 24 hours for the same small run (500 transitions).

## Help!

Is there an alternative parameterization that might make this posterior easier to sample?

I'm really new to Stan (this is my 3rd bespoke model), and I'm no mathematician. So it's entirely possible that I'm doing something dumb.

## What I've tried so far

The priors are pretty reasonable. They're not strongly informative, but they are tuned to restrict the prior values to reasonable ranges.

The random effects all use an uncentered parameterization. I've fit versions with centered parameterizations, and they are notably worse (much lower n_eff values).

Three of the random effects also include soft sum-to-zero constraints `{stan, eval = FALSE} sum(strata_p) ~ normal(0,0.001*nstrata)`. There are some inherent collinearities in the BBS data, and these constraints should (and seem to) help estimate parameters that have some correlation, such as the various stratum-level parameters that model and account for temporal changes in the observations:

-   slopes (rates of change in abundance over time),

-   intercepts (mean abundance in a region),

-   year-effects (annual, random fluctuations in abundance, in addition to the changes due to the slope)

I've tried versions without these constraints and the constraints help, in particular they improve the sampling of the hyperparameters on slope and intercept.

## The data

```{r dataload, message = FALSE}

library(rstan)
rstan_options(auto_write = TRUE)
library(shinystan)
load("data/prepared_BBS_data_Pacific_Wren.RData")
str(stan_data)

```

I've prepared the data using the package bbsBayes [www.github.com/BrandonEdwards/bbsBayes](www.github.com/BrandonEdwards/bbsBayes), which downloads the raw BBS data and generates species-specific datafiles designed for a customized JAGS model. I've then modified the data object to suit the Stan model. The process used to download, create and modify, the BBS data is outlined at the end of this document.

The data object is a list with the following components:

-   count - the counts of individual birds on each route and year (the example data here has \~4000 counts)

-   strat - the integer indicator for each of the 18 geographic strata included in this species' range (some species in the BBS database cover \~ 150 strata)

-   year - the integer indicating which year of the time-series

-   firstyr - 1 if this route-year combination represents the observers first year on this particular route

-   obser - integer indicating which observer conducted the survey

-   ncounts - number of observations

-   nstrata - number of geographic strata

-   nyears - number of years included

-   nobservers - number of unique observer-route combinations in the data

-   fixedyear - midpoint of the time-series, used to center the year-values

## The model

```{stan, eval = FALSE}
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
  
  
 for(s in 1:nstrata){

  yeareffect_raw[s,] ~ normal(0,1);
  sum(yeareffect_raw[s,]) ~ normal(0,0.001*nyears);// soft sum to zero constraint
  
 }
  count ~ poisson_log(E); //vectorized count likelihood with log-transformation
  
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope - hyperparameter
  STRATA ~ normal(0,1);// prior on fixed effect mean intercept - hyperparameter
  eta ~ normal(0,1);// prior on first-year observer effect
  
  
  sdstrata ~ normal(0,1); //prior on sd of intercept variation
  sdbeta ~ normal(0,0.1); //prior on sd of slope variation

  beta_p ~ normal(0,1); //non centered prior on stratum-level slopes
  strata_p ~ normal(0,1); //non centered prior on stratum-level intercepts

  //sum to zero constraints
  sum(strata_p) ~ normal(0,0.001*nstrata);
  sum(beta_p) ~ normal(0,0.001*nstrata);
  
}

```

Here's an example of the settings for a model run that required \~10 hours to complete

```{r model_run, eval = FALSE}

mod.file = "models/slope.stan"

parms = c("sdnoise",
          "sdyear",
          "sdobs",
          "beta_p",
          "sdbeta",
          "strata_p",
          "sdstrata",
          "BETA",
          "STRATA",
          "eta")

## compile model
slope_model = stan_model(file=mod.file)

# run sampler on model, data
slope_stanfit <- sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=500,
                               warmup=400,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15))


# save(list = c("slope_stanfit","species","mod.file","model","strat"),
#      file = "output/Pacific Wren_slope_stan_saved_output.RData")
```

```{r summary,eval = TRUE}
load("output/Pacific Wren_slope_stan_saved_output.RData")

get_elapsed_time(slope_stanfit)/3600 ## in hours
          warmup   sample
chain:1 7.798528 2.105014
chain:2 8.093528 2.224678
chain:3 7.706306 2.151717

#print(slope_stanfit)

conv_diag = get_sampler_params(slope_stanfit, inc_warmup = FALSE)

signif(sapply(conv_diag,simplify = TRUE, function(x) colMeans(x)),2)
                  [,1]     [,2]     [,3]
accept_stat__  8.8e-01  8.7e-01  8.2e-01
stepsize__     3.8e-03  4.0e-03  4.6e-03
treedepth__    1.5e+01  1.5e+01  1.5e+01
n_leapfrog__   3.3e+04  3.3e+04  3.3e+04
divergent__    0.0e+00  0.0e+00  0.0e+00
energy__      -6.6e+04 -6.6e+04 -6.6e+04

#launch_shinystan(slope_stanfit) 


```





## Original data prep for JAGS model

```{r earlier_data_prep, eval = FALSE}
# This is setup that is not required to run this example ------------------
#install.packages("bbsBayes")
#fetch_bbs_data()
#yes

library(bbsBayes)



# load and stratify PAWR data ---------------------------------------------
species = "Pacific Wren"
strat = "bbs_usgs"
model = "slope"

#stratify the BBS data based on the standard USGS stratification,
# intersection of states_provinces and Bird Conservation Regions
strat_data = stratify(by = strat)

#prepare BBS data for Pacific Wren, retaining only counts since 1999
jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             min_year = 1999)

# exporting the original JAGS version of the model
# some version of this has been used since 2011 to produce status and trend estimates by the Canadian and US federal agencies
model_to_file(model = model, filename = "Original_JAGS_model.R")
# the priors in this model are different, but generally only in the sense that these
# JAGS priors are mostly "uninformative" and suggest far flatter prior distributions than is realistic
# the priors in the Stan model called below are not strongly informative, but they should be more reasonable




#select out portions of the
stan_data = jags_data[c("ncounts", #number of counts (observations)
                         "nstrata", #number of groups (geographic strata) for the random-effect slopes
                         "count", #counts of individuals observed at a site (BBS survey route) in a year
                         "strat", #group indicators (strata)
                         "year", #years scaled as 1 through nyears.
                         "firstyr", #indicator variable for the first year an observer surveyed a particular route - nuisance parameter (start-up effect)
                         "fixedyear")] #mid-point year in the time-series, used to center the year variables
stan_data[["nyears"]] <- max(jags_data$year) #number of years in the time-series
stan_data[["nobservers"]] <- sum(jags_data$nobservers)#number of observer-route combinations in each stratum - random-effect nuisance parameter
stan_data[["obser"]] <- as.integer(factor(paste(jags_data$strat,jags_data$obser,sep = "_"))) #observer-route indicators


save(list = c("stan_data",
              "species",
              "strat",
              "jags_data",
              "model"),
    file = "data/prepared_BBS_data_Pacific_Wren.RData")




```

## References
