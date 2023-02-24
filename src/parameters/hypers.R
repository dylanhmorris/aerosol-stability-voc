#!/usr/bin/env Rscript

###################################
## filename: hypers.R
## author: Dylan H. Morris (dylanhmorris.com)
## description: contains hyperparameters for
## Bayesian models of decay of aerosolized
## virus
####################################

titer_hypers <- list(
    titer_prior_mean = 3,
    titer_prior_sd = 3)

aerosol_hypers <- list(
    aerosol_log_hl_prior_mean = log(5),
    aerosol_log_hl_prior_sd = log(4),
    aerosol_intercept_prior_mean = 2,
    aerosol_intercept_prior_sd = 2,
    aerosol_sd_intercept_prior_mode = 0.4,
    aerosol_sd_intercept_prior_sd = 0.2
)

surface_hypers <- list(
    surface_log_hl_prior_mean = log(5),
    surface_log_hl_prior_sd = log(4),
    surface_intercept_prior_mean = 3,
    surface_intercept_prior_sd = 2,
    surface_sd_intercept_prior_mode = 0.4,
    surface_sd_intercept_prior_sd = 0.2
)

flags <- list(
    debug = FALSE)

hyperparam_list <- c(
    titer_hypers,
    aerosol_hypers,
    surface_hypers,
    flags)


strict <- TRUE
