#!/usr/bin/env Rscript

###################################
## filename: infer-titers.R
## author: Dylan H. Morris (dylanhmorris.com)
## description: infer titers from raw
## well data
####################################

#'
#' infer_titers
#'
#' Bayesian titer inference with Stan
#'
#' @param model_src_path path to the stan model from
#' which to draw samples.
#' 
#' @param well_status binary array of whether inoculated
#' wells were positive (1 / TRUE) or negative (0 / FALSE)
#' 
#' @param dilution numeric array of log10 well dilution factors
#' (0 = no dilution, -1 = 10-fold dilution, etc).
#' 
#' @param sample_id which sample the well corresponds to
#' 
#' @param titer_prior_mean Normal prior mean for log10 titers
#' 
#' @param titer_prior_sd Normal prior sd for log10 titers
#' 
#' @param prior_check whether to perform a prior predictive check
#' (fit no data, generating predictive checks based solely on the
#' prior distributions given) Default FALSE.
#' 
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' 
#' @return An object of class `stanfit` returned by `rstan::sampling`
infer_titers <- function(model_src_path,
                         well_status,
                         dilution,
                         sample_id,
                         titer_prior_mean,
                         titer_prior_sd,
                         prior_check = FALSE,
                         debug = FALSE,
                         ...) {
    n_samples <- max(sample_id)
    n_obs <- length(well_status)
    
    ## perform data quality checks
    check_data_quality(n_obs,
                       dilution,
                       list("sample" = sample_id))

    print(unique(sample_id))
    
    standata <- list(
        well_status = well_status,
        dilution = dilution,
        sample_id = sample_id,
        n_total_datapoints = n_obs,
        n_used_datapoints = n_obs,
        n_samples = n_samples,
        n_used_samples = n_samples,
        titer_prior_mean = titer_prior_mean,
        titer_prior_sd = titer_prior_sd,
        debug = debug)
    
    if(prior_check){
        standata[["n_used_datapoints"]] = 0
    }

    fit <- rstan::stan(
        model_src_path,
        data = standata,
        ...)
    
    return(fit)
}

suppressPackageStartupMessages(library(varstab))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(parallel))

## read command line args
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
hyperparam_path <- args[2]
mcmc_model_path <- args[3]
mcmc_output_path <- args[4]
prior_check <- grepl("prior-check", mcmc_output_path)


cat('reading in data from file ', data_path, ' ...\n')
dat <- read_tsv(data_path,
                col_types = cols())
cat('data loaded successfully!\n')

cat('reading in hyperparameters from file ', hyperparam_path, ' ...\n')
source(hyperparam_path)
cat('hyperparameters loaded successfully!\n')

###############################
## Compile, fit, and save model
###############################
## pass stan data and hyperparams
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)
rstan_options(auto_write = TRUE)
fixed_seed <- 234


if("virus_detect" %in% names(dat)){
    well_status = dat$virus_detect
} else if ("well_status" %in% names(dat)) {
    well_status = dat$well_status
} else {
    stop("Data must contain well statuses")
}

if(prior_check){
    cat("Performing prior check; not fitting to observations...\n")
}

print(unique(dat$sample_id))



fit_titers <- infer_titers(
    mcmc_model_path,
    well_status = well_status,
    dilution = dat$dilution,
    sample_id = dat$sample_id,
    titer_prior_mean = titer_hypers[["titer_prior_mean"]],
    titer_prior_sd = titer_hypers[["titer_prior_sd"]],
    seed = fixed_seed,
    prior_check = prior_check)


## check that sampled correctly and only save
## if so
sampler_params <- get_sampler_params(fit_titers, inc_warmup = FALSE)

if(is.null(sampler_params)){
    stop("Stan model failed to sample")
}

divs <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))

if(strict){
    if(divs > 0){
        stop("Divergent transitions when sampling. Check priors and parametrization")
    }
}

cat("Saving results to ", mcmc_output_path, "...\n")

saveRDS(fit_titers, mcmc_output_path)

warnings()
