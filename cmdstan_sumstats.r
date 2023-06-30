# Script to generate sumstats and mcmc matrix similar to rstan,
#  after fitting a model using cmdstan 
#  NOTE: following lines show sample of cmdstan fitting code
# ------------------------------------------------------
# require(parallel)
# require(cmdstanr)
# require(posterior)
# nburnin = 500                    # number warm-up (burn-in) samples 
# nsamples = 10000                 # desired total number samples
# fitmodel = c("filename.stan")    # name of file with stan code
# parms = c("Par1","Par2")         # vector of parameters to save
# stan.data <- list(N=N,X=X,Y=Y)   # list of input data variables
# cores = detectCores()
# ncore = min(40,cores-4)
# Niter = round(nsamples/ncore)
# mod <- cmdstan_model(fitmodel)   # compiles model (if necessary) 
# suppressMessages(                # Suppress messages/warnings (if desired)
#   suppressWarnings ( 
#     fit <- mod$sample(
#       data = stan.data,
#       seed = 1234,
#       chains = ncore,
#       parallel_chains = ncore,
#       refresh = 100,
#       iter_warmup = nburnin,
#       iter_sampling = Niter,
#       max_treedepth = 12,
#       adapt_delta = 0.8
#     )
#   )
# )
# source("cmdstan_sumstats.r")
# ------------------------------------------------------
#
sumstats = fit$draws(parms) %>%
  summarise_draws(mean, mcse = mcse_mean, sd, 
                  ~quantile(.x, probs = c(0.025, 0.05, .5, .95, .975),na.rm=T),
                  N_eff = ess_bulk, rhat)
sumstats = as.data.frame(sumstats)
row.names(sumstats) = sumstats$variable; sumstats = sumstats[,-1] 
#
mcmc = as_draws_matrix(fit$draws(variables = parms))
Nsims = nrow(mcmc)
mcmc_array = as_draws_array(fit$draws(variables = parms))
vn = colnames(mcmc)
# rm(mod,fit)       # can also remove mod & fit from memory if no longer needed 
