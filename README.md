---
bibliography: references.bib
---

# bbsBayes_Stan_Minimal

This is an minimum working example of an early attempt to use Stan to fit a relatively simple hierarchical model applied to data from the North American Breeding Bird Survey. I'm trying to build a Stan version of the JAGS models currently used to analyse the North American Breeding Bird Survey (BBS) monitoring data [@sauer2011][@smith2015][@smith2020]. This repo includes the code, model, and data for a minimal working example of one of the simplest models, using a species with pretty representative data. It's currently set up to use only the most recent 20-years of data, but the continental scale monitoring program includes \> 50 years of data. The data are observations (summed counts) of individual birds observed during annual surveys of BBS "routes".

The model works. It generates reasonable estimates, no divergent transitions, effective sample sizes are generally \~80% of the total samples. If the max_treedepth argument is increased.

## The problem

The problem is that the sampler exceeds the maximum tree depth, even if the limit is set to \~15-17 (haven't tried higher because life is short). The high max_treedepth limit means that the sampler takes a really long time, even for this species which has very few data compared to many in the dataset.

## What I've tried so far

The priors are pretty reasonable. They're not strongly informative, but they are tuned to restrict the prior values to reasonable ranges.

The random effects all use an uncentered parameterization. I've fit versions with centered parameterizations, and they are notably worse (much lower n_eff values).

Three of the random effects also include soft sum-to-zero constraints \` sum(strata_p) \~ normal(0,0.001\*nstrata) \`. There are some inherent collinearities in the BBS data, and these constraints help to estimate parameter-combinations such as the various stratum-level parameters that model and account for temporal changes in the observations:

-   slopes (rates of change in abundance over time),

-   intercepts (mean abundance in a region),

-   year-effects (annual, random fluctuations in abundance, in addition to the changes due to the slope)

I've tried versions without these constraints and the constraints help, in particular they improve the estimates of the hyperparameters on slope and intercept.
