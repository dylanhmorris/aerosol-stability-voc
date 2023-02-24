#!/usr/bin/env Rscript

###################################
## filename: infer-titers.R
## author: Dylan H. Morris (dylanhmorris.com)
## description: infer titers from raw
## well data
####################################

#' infer_timeseries_halflives()
#' 
#' function to fit half-life inference
#' model using Stan when the titer observations
#' form true timeseries
#' 
#' @param model_src_path path to the stan model from
#' which to draw samples.

#' @param well_status binary array of whether inoculated
#' wells were positive (1 / TRUE) or negative (0 / FALSE)
#'
#' @param dilution numeric array of log10 well dilution factors
#' (0 = no dilution, -1 = 10-fold dilution, etc).
#'
#' @param time time since initial virus deposition
#' that the sample in the well was taken
#' 
#' @param sample_id which sample the well corresponds to
#'
#' @param experiment_id which experiment the well corresponds to
#' 
#' @param timeseries_id which timeseries the well corresponds to
#'
#' @param log_change_virus_genomes Log change in virus
#' genetic material quantity
#' relative to the intercept
#' (e.g. calculated from Ct values
#' or estimated copy numbers)
#' 
#' @param log_hl_prior_mean Normal prior mean for the
#' log (base e) of the half-life
#'
#' @param log_hl_prior_sd Normal prior sd for the
#' log (base e) of the half-life
#'
#' @param intercept_prior_mean Normal prior mean
#' for mean initial log10 titer (mean value at t = 0)
#'
#' @param intercept_prior_sd Normal prior sd
#' for mean initial log10 titers (mean value at t = 0)
#'
#' @param sd_intercept_prior_mode Half-Normal prior mode
#' for the estimated sd of the individual initial titers
#' (i.e. the titer values at t = 0) about their
#' shared mean for a given experiment
#' 
#' @param sd_intercept_prior_sd Half-Normal prior sd
#' for the estimated sd of the individual initial titers
#' (i.e. the titer values at t = 0) about their
#' shared mean for a given experiment

#' @param prior_check whether to perform a prior predictive check
#' (fit no data, generating predictive checks based solely on the
#' prior distributions given) Default FALSE.
#' 
#' @param debug whether to use debugging settings (default FALSE)

#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' 
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' 
#' @export
infer_timeseries_halflives <- function(model_src_path,
                                       well_status,
                                       dilution,
                                       time,
                                       sample_id,
                                       experiment_id,
                                       timeseries_id,
                                       log_change_virus_genomes,
                                       log_hl_prior_mean,
                                       log_hl_prior_sd,
                                       intercept_prior_mean,
                                       intercept_prior_sd,
                                       sd_intercept_prior_mode,
                                       sd_intercept_prior_sd,
                                       prior_check = FALSE,
                                       debug = FALSE,
                                       ...) {
    n_samples <- max(sample_id)
    n_obs <- length(well_status)
    n_experiments <- max(experiment_id)
    n_timeseries <- max(timeseries_id)
    
    ## perform data quality checks
    check_data_quality(n_obs,
                       dilution,
                       list(
                           "sample" = sample_id,
                           "experiment" = experiment_id,
                           "timeseries" = timeseries_id))
    
    sample_dat <- tibble::tibble(
                              sample_id = sample_id,
                              experiment_id = experiment_id,
                              timeseries_id = timeseries_id,
                              time = time,
                              log_change_virus_genomes =
                                  log_change_virus_genomes)
    
    sample_dat <- dplyr::distinct(
                             sample_dat,
                             sample_id,
                             .keep_all = TRUE) |>
        arrange(sample_id)
    
    timeseries_dat <- dplyr::distinct(
                                 sample_dat,
                                 timeseries_id,
                                 .keep_all = TRUE) |>
        arrange(timeseries_id)
    
    
    standata <- list(
        well_status = well_status,
        dilution = dilution,
        sample_id = sample_id,
        experiment_id = experiment_id,
        n_total_datapoints = n_obs,
        n_used_datapoints = n_obs,
        n_titers = n_samples,
        n_used_samples = n_samples,
        n_experiments = n_experiments,
        n_runs = n_timeseries,
        sample_times = sample_dat$time,
        sample_experiment_id = sample_dat$experiment_id,
        sample_run_id = sample_dat$timeseries_id,
        log_change_virus_genomes =
            sample_dat$log_change_virus_genomes,
        run_experiment_id = timeseries_dat$experiment_id,
        log_hl_prior_mean = log_hl_prior_mean,
        log_hl_prior_sd = log_hl_prior_sd,
        intercept_prior_mean = intercept_prior_mean,
        intercept_prior_sd = intercept_prior_sd,
        sd_intercept_prior_mode =
            sd_intercept_prior_mode,
        sd_intercept_prior_sd =
            sd_intercept_prior_sd,
        debug = debug)
    
    if(prior_check){
        standata[["n_used_datapoints"]] = 0
    }
    
    fit <- rstan::stan(model_src_path,
                       data = standata,
                       ...)
    
    return(fit)
}


suppressPackageStartupMessages(library(varstab))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))

## read command line args
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
hyperparam_path <- args[2]
mcmc_model_path <- args[3]
mcmc_output_path <- args[4]

prior_check <- grepl("prior-check", mcmc_output_path)
omit_pcr <- grepl("no-pcr", mcmc_output_path)

## read and process data
cat('reading in data from file ', data_path, ' ...\n')
dat <- read_tsv(data_path,
                col_types = cols())
cat('data loaded successfully!\n')

cat('reading in hyperparameters from file ', hyperparam_path, ' ...\n')
source(hyperparam_path)
cat('hyperparameters loaded successfully!\n')

## optionally do not PCR correct for physical loss
if(omit_pcr){
    dat <- dat |>
        mutate(log_change_virus_genomes = 0)
} else {
    dat <- dat |>
        mutate(log_change_virus_genomes = log10_ct_rel_t0)
}

## remove unmodeled pre-spray t = -1
dat <- dat |>
    filter(!pre_spray)

print(min(dat$time))
print(min(dat$sample_id))
print(max(dat$sample_id))
print(length(unique(dat$sample_id)))

###############################
## Compile, fit, and save model
###############################
## pass stan data and hyperparams
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)

fit_halflives <- infer_timeseries_halflives(
    mcmc_model_path,
    well_status = dat$virus_detect,
    dilution = dat$dilution,
    time = dat$time,
    sample_id = dat$sample_id,
    experiment_id = dat$experiment_id,
    timeseries_id = dat$run_id,
    log_change_virus_genomes = dat$log_change_virus_genomes,
    log_hl_prior_mean = aerosol_hypers[["aerosol_log_hl_prior_mean"]],
    log_hl_prior_sd = aerosol_hypers[["aerosol_log_hl_prior_sd"]],
    intercept_prior_mean =
        aerosol_hypers[["aerosol_intercept_prior_mean"]],
    intercept_prior_sd =
        aerosol_hypers[["aerosol_intercept_prior_sd"]],
    sd_intercept_prior_mode =
        aerosol_hypers[["aerosol_sd_intercept_prior_mode"]],
    sd_intercept_prior_sd =
        aerosol_hypers[["aerosol_sd_intercept_prior_sd"]],
    seed = 23421,
    prior_check = prior_check,
    control = list(max_treedepth = 13,
                   adapt_delta = 0.95))


## check that sampled correctly and only save
## if so
sampler_params <- get_sampler_params(fit_halflives, inc_warmup = FALSE)

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

saveRDS(fit_halflives, mcmc_output_path)

warnings()
