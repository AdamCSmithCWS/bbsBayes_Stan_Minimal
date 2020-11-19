## Trying to build a Stan version of the JAGS models currently used to 
## analyse the North American Breeding Bird Survey monitoring data
## This repo includes a minimal example of one of the simplest models
## for a species that provides a reasonably representative example
## set up to use only the most recent 20-years of data as a simplification.




# This is setup that is not required to run this example ------------------
#install.packages("bbsBayes")
#fetch_bbs_data()
#yes

library(bbsBayes)



# # load and stratify PAWR data ---------------------------------------------
# species = "Pacific Wren"
# strat = "bbs_usgs"
# model = "slope"
# 
# #stratify the BBS data based on the standard USGS stratification,
# # intersection of states_provinces and Bird Conservation Regions
# strat_data = stratify(by = strat)
# 
# #prepare BBS data for Pacific Wren, retaining only counts since 1999
# jags_data = prepare_jags_data(strat_data = strat_data,
#                              species_to_run = species,
#                              model = model,
#                              min_year = 1999)
# 
# # exporting the original JAGS version of the model
# # some version of this has been used since 2011 to produce status and trend estimates by the Canadian and US federal agencies
# model_to_file(model = model, filename = "Original_JAGS_model.R")
# # the priors in this model are different, but generally only in the sense that these
# # JAGS priors are mostly "uninformative" and suggest far flatter prior distributions than is realistic
# # the priors in the Stan model called below are not strongly informative, but they should be more reasonable
# 
# 
# 
# 
# #select out portions of the
# stan_data = jags_data[c("ncounts", #number of counts (observations)
#                          "nstrata", #number of groups (geographic strata) for the random-effect slopes
#                          "count", #counts of individuals observed at a site (BBS survey route) in a year
#                          "strat", #group indicators (strata)
#                          "year", #years scaled as 1 through nyears.
#                          "firstyr", #indicator variable for the first year an observer surveyed a particular route - nuisance parameter (start-up effect)
#                          "fixedyear")] #mid-point year in the time-series, used to center the year variables
# stan_data[["nyears"]] <- max(jags_data$year) #number of years in the time-series
# stan_data[["nobservers"]] <- sum(jags_data$nobservers)#number of observer-route combinations in each stratum - random-effect nuisance parameter
# stan_data[["obser"]] <- as.integer(factor(paste(jags_data$strat,jags_data$obser,sep = "_"))) #observer-route indicators
# 
# 
# save(list = c("stan_data",
#               "species",
#               "strat",
#               "jags_data",
#               "model"),
#     file = "data/prepared_BBS_data_Pacific_Wren.RData")



load("data/prepared_BBS_data_Pacific_Wren.RData")
library(rstan)
rstan_options(auto_write = TRUE)
library(shinystan)



# replace the real counts with simulated counts ---------------------------
nstrata = stan_data$nstrata
ncounts = stan_data$ncounts
strat = stan_data$strat
year = stan_data$year
fixedyear = stan_data$fixedyear
nyears = stan_data$nyears
nobservers = stan_data$nobservers
obser = stan_data$obser
firstyr = stan_data$firstyr

sdnoise = 0.25
sdobs = 1.1
sdyear = runif(nstrata,0.05,1)
sdstrata = 1
sdbeta = 0.01


STRATA = -1
B = -0.02
eta = -0.1
beta = rnorm(nstrata,B,sdbeta)
strata = rnorm(nstrata,STRATA,sdstrata)
yeareffect = matrix(NA,nrow = nstrata,ncol = nyears)
for(s in 1:nstrata){
  yeareffect[s,1:nyears] = rnorm(nyears,0,sdyear[s])
}
obs = rnorm(nobservers,0,sdobs)
noise = rnorm(ncounts,0,sdnoise)
E = vector(mode = "numeric",length = ncounts)
simcount =vector(mode = "integer",length = ncounts)
for(i in 1:ncounts){
  E[i] =  beta[strat[i]] * (year[i]-fixedyear) + strata[strat[i]] + 
    yeareffect[strat[i],year[i]] + obs[obser[i]] + eta*firstyr[i] + noise[i]
  simcount[i] = rpois(1,exp(E[i]))
}



# loo capable version -----------------------------------------------------


parms = c("sdnoise",
          "sdyear",
          "sdobs",
          "beta_p",
          "sdbeta",
          "strata_p",
          "sdstrata",
          "BETA",
          "STRATA",
          "eta",
          "log_lik")
mod.file = "models/slope_alt.stan"

stan_data_sim = stan_data
stan_data_sim$count = simcount

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
slope_stanfit <-sampling(slope_model,
                         data=stan_data_sim,
                         verbose=TRUE, refresh=100,
                         chains=3, iter=500,
                         warmup=400,
                         cores = 3,
                         pars = parms,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 15))

print(slope_stanfit)
get_elapsed_time(slope_stanfit)/3600 ## in hours

library(loo)
library(tidyverse)
log_lik_1 <- extract_log_lik(slope_stanfit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 10)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 10)
print(loo_1)

plot(loo_1$pointwise[,"influence_pareto_k"],stan_data$count)

loo2 = data.frame(loo_1$pointwise)

 loo2$flag = cut(loo2$influence_pareto_k,breaks = c(0,0.5,0.7,1,Inf))
 dts = data.frame(count = stan_data$count,
                  obser = stan_data$obser,
                  strat = stan_data$strat,
                  year = stan_data$year)
  loo2 = cbind(loo2,dts)

  plot(log(loo2$count+1),loo2$influence_pareto_k)
  
  obserk = loo2 %>% group_by(obser) %>% 
    summarise(n = n(),
              mean_k = mean(influence_pareto_k),
              max_k = max(influence_pareto_k),
              sd_k = sd(influence_pareto_k),
              strat = mean(strat),
              sd = sd(strat))
  plot(obserk$n,obserk$max_k)
  plot(obserk$n,obserk$sd_k)
  
  
  yeark = loo2 %>% group_by(year) %>% 
    summarise(n = n(),
              mean_k = mean(influence_pareto_k),
              q90 = quantile(influence_pareto_k,0.9),
              max_k = max(influence_pareto_k),
              sd_k = sd(influence_pareto_k),
              strat = mean(strat),
              sd = sd(strat))
  plot(yeark$year,yeark$max_k)
  plot(yeark$year,yeark$mean_k)
  plot(yeark$year,yeark$sd_k)
  plot(yeark$year,yeark$q90)
  
  stratk = loo2 %>% group_by(strat) %>% 
    summarise(n = n(),
              mean_k = mean(influence_pareto_k),
              q90_k = quantile(influence_pareto_k,0.9),
              max_k = max(influence_pareto_k),
              sd_k = sd(influence_pareto_k),
              strat = mean(strat),
              sd = sd(strat))
  plot(stratk$strat,stratk$max_k)
  plot(stratk$n,stratk$mean_k)
  
  plot(stratk$strat,stratk$mean_k)
  plot(stratk$strat,stratk$sd_k)
  plot(stratk$strat,stratk$q90_k)
  
  
save(list = c("slope_stanfit","species","mod.file","model","strat"),
     file = "output/Pacific Wren_slope_stan_saved_output_reparam3.RData")


# heavy-tailed version with loo -------------------------------------------


parms = c("sdnoise",
          "sdyear",
          "sdobs",
          "beta_p",
          "sdbeta",
          "strata_p",
          "sdstrata",
          "BETA",
          "STRATA",
          "eta",
          #"nu",
          "log_lik")
mod.file = "models/slope_heavy.stan"


## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
slope_stanfit <-sampling(slope_model,
                         data=stan_data,
                         verbose=TRUE, refresh=100,
                         chains=3, iter=500,
                         warmup=400,
                         cores = 3,
                         pars = parms,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 15))

print(slope_stanfit)
get_elapsed_time(slope_stanfit)/3600 ## in hours

library(loo)
library(tidyverse)
log_lik_1 <- extract_log_lik(slope_stanfit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 10)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 10)
print(loo_1)







load("output/Pacific Wren_slope_stan_saved_output.RData")
print(slope_stanfit)
get_elapsed_time(slope_stanfit)/3600 ## in hours

launch_shinystan(slope_stanfit) 

