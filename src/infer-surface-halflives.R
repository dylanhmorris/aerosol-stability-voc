#!/usr/bin/env Rscript

###################################
## filename: infer-surface-halflives.R
## author: Dylan H. Morris (dylanhmorris.com)
## description: infer halflives from raw
## well data
####################################

#' infer_halflives
#'
#' Infer halflives when individual datapoints
#' do not form a true timeseries (destructively
#' sampled, as in surface experiments).
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
#' @param time time since initial virus deposition
#' that the sample in the well was taken
#' 
#' @param sample_id which sample the well corresponds to
#' 
#' @param experiment_id which experiment the well corresponds to
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
#'
#' @param prior_check whether to perform a prior predictive check
#' (fit no data, generating predictive checks based solely on the
#' prior distributions given) Default FALSE.
#' 
#' @param debug whether to use debugging settings (default FALSE)
#' 
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' 
#' @return An object of class `stanfit` returned by `rstan::sampling`
#' 
infer_halflives <- function(model_src_path,
                            well_status,
                            dilution,
                            time,
                            sample_id,
                            experiment_id,
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

    

    if(length(dilution) != n_obs)
        stop(paste0("Must have same number of dilution factors ",
                    "as well status observations"))
    if(length(sample_id) != n_obs)
        stop(paste0("Must have same number of sample_ids ",
                    "as well status observations"))
    if(length(experiment_id) != n_obs)
        stop(paste0("Must have same number of experiment_ids ",
                    "as well status observations"))
    if(length(time) != n_obs)
        stop(paste0("Must have same number of times ",
                    "as well status observations"))
    
    sample_dat_long <- tibble::tibble(
        sample_id = sample_id,
        experiment_id = experiment_id,
        time = time)

    test_summary <- sample_dat_long |>
        dplyr::group_by(sample_id) |>
        dplyr::reframe(
                   time_test = time == first(time),
                   experiment_test = experiment_id == first(experiment_id))
    
    if(!all(test_summary$time_test)){
        stop(paste0("All observations with the same sample_id",
                    "must have the same time"))
    }
    if(!all(test_summary$experiment_test)){
        stop(paste0("All observations with the same sample_id",
                    "must have the same experiment_id"))
    }
    
    sample_dat <- sample_dat_long |>
        dplyr::distinct(sample_id,
                        .keep_all = TRUE) |>
        dplyr::arrange(sample_id)
    ## this is critical in case sample_ids are
    ## supplied out of order in the main dataset
    
    standata <- list(
        well_status = well_status,
        dilution = dilution,
        sample_id = sample_id,
        experiment_id = experiment_id,
        n_total_datapoints = n_obs,
        n_used_datapoints = n_obs,
        n_samples = n_samples,
        n_used_samples = n_samples,
        n_experiments = n_experiments,
        sample_times = sample_dat$time,
        sample_experiment_id = sample_dat$experiment_id,
        log_hl_prior_mean = log_hl_prior_mean,
        log_hl_prior_sd = log_hl_prior_sd,
        intercept_prior_mean = intercept_prior_mean,
        intercept_prior_sd = intercept_prior_sd,
        sd_intercept_prior_mode = sd_intercept_prior_mode,
        sd_intercept_prior_sd = sd_intercept_prior_sd,
        debug = debug)

    if(prior_check){
        standata[["n_used_datapoints"]] = 0
    }


    fit <- rstan::stan(
        model_src_path,
        data = standata,
        ...)

    return (fit)
}



suppressPackageStartupMessages(library(varstab))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(dplyr))

## read command line args
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
hyperparam_path <- args[2]
mcmc_model_path <- args[3]
mcmc_output_path <- args[4]
prior_check <- grepl("prior-check", mcmc_output_path)

cat('reading in titration data from file ',
    data_path,
    '...\n')

dat <- read_tsv(data_path,
                col_types = cols())

cat('data loaded successfully!\n')


cat('reading in hyperparameters from file ',
    hyperparam_path, ' ...\n')
source(hyperparam_path)
cat('hyperparameters loaded successfully!\n')

strict <- TRUE
fixed_seed <- 532
###############################
## Compile, fit, and save model
###############################
## pass stan data and hyperparams
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)


fit_halflives <- infer_halflives(
    mcmc_model_path,
    well_status = dat$well_status,
    dilution = dat$dilution,
    time = 1.0 * dat$time,
    sample_id = dat$sample_id,
    experiment_id = dat$experiment_id,
    log_hl_prior_mean =
        surface_hypers[["surface_log_hl_prior_mean"]],
    log_hl_prior_sd =
        surface_hypers[["surface_log_hl_prior_sd"]],
    intercept_prior_mean =
        surface_hypers[["surface_intercept_prior_mean"]],
    intercept_prior_sd =
        surface_hypers[["surface_intercept_prior_sd"]],
    sd_intercept_prior_mode =
        surface_hypers[["surface_sd_intercept_prior_mode"]],
    sd_intercept_prior_sd =
        surface_hypers[["surface_sd_intercept_prior_sd"]],
    prior_check = prior_check,
    seed = fixed_seed,
    control = list(max_treedepth = 12,
                   adapt_delta = 0.8))


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
