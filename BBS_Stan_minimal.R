## Trying to build a Stan version of the JAGS models currently used to 
## analyse the North American Breeding Bird Survey monitoring data
## This repo includes a minimal example of one of the simplest models
## for a species that provides a reasonably representative example
## set up to use only the most recent 20-years of data as a simplification.




# # This is setup that is not required to run this example ------------------
# #install.packages("bbsBayes")
# #fetch_bbs_data()
# #yes
# 
# library(bbsBayes)
# 
# 
# 
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
          # "n",
          # "nsmooth",
          "eta")

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
stime = system.time(slope_stanfit <-
                      sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=500,
                               warmup=400,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15)))

stime[[3]]/3600

save(list = c("stime","slope_stanfit","species","mod.file","model","strat"),
     file = paste(species,mod.file,model,"saved_output.RData",sep = "_"))

launch_shinystan(slope_stanfit) 



