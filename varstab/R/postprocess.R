#!/usr/bin/env Rscript


########################################
## filename: postprocess.R
## author: Dylan H. Morris (dylanhmorris.com)
## description: Postprocessing of MCMC output
## for variant stability
#######################################

#' add_titer_metadata
#'
#' add metadata to a set of MCMC chains
#' ordered by titer_id
#'
#' @param mcmc_draws mcmc chains in spread_draws tibble
#' format from tidybayes
#' @param fitting_data data used to fit the mcmc model
#' as a tibble with a titer id column
#' @param id_column name of the id column, default "titer_id"
#' @param multiple = {"all", NULL, "error"}. Flag passed
#' to dplyr::inner_join(). Whether to expect multiple
#' rows of the output table per row of the index table.
#' As we almost always do, this defaults to
#' "all" to avoid a dplyr warning.
#' @return tibble of draws with metadata added
#' @export
add_titer_metadata <- function(mcmc_draws,
                               fitting_data,
                               id_column = "titer_id",
                               verbose_errors = FALSE,
                               multiple = "all"){
    
    filt_dat <- dplyr::distinct(fitting_data,
                                fitting_data[id_column],
                                .keep_all = TRUE)
    
    result <- dplyr::inner_join(filt_dat,
                                mcmc_draws,
                                by = id_column,
                                multiple = multiple)

    if(!dim(result)[1] == dim(mcmc_draws)[1]) {
        if ( verbose_errors ) {
            cat("Draws:\n")
            print(mcmc_draws)
            cat("\n\nIndex created from data:\n")
            print(filt_dat)
            cat("\n\nResult:\n")
            print(result)
        }

        stop(sprintf(paste0(
            "Error: joining data to draws on column %s ",
            "failed to produce set of draws ",
            "with the same number of rows, ",
            "suggesting draws were lost; ",
            "check your data to make sure the ",
            "join column has all the needed values"),
            id_column))
    }
    return(result)
}


#' format_est_ci
#'
#' format a central estimate (e.g. a mode)
#' and its bounding interval (CI, e.g. a
#' credible interval or a confidence 
#' interval) for in-text display
#' as:
#'
#' central [lower CI, upper CI]
#'
#' @param estimate the central estimate
#' @param lower_ci the lower bound of the CI
#' @param upper_ci the upper bound of the CI
#' @param significant_digits how many digits to treat as
#' significant in rounding / truncating
#' 
#' @export
format_est_ci <- function(estimate,
                          lower_ci,
                          upper_ci,
                          significant_digits = 3){
    sprintf(
        "%g [%g, %g]",
        estimate |>
        signif(digits = significant_digits),
        lower_ci |>
        signif(digits = significant_digits),
        upper_ci |>
        signif(digits = significant_digits))
}



#' to_decay_rate
#'
#' convert a an nth-life
#' (i.e. how long until the fraction in
#' question is all that remains)
#' into a decay rate in logs per hour
#' (default log base 10)
#'
#' @param fractional_life fractional life to convert
#' in a given time unit
#' @param fraction used (default 1 / 2: half-life)
#' @param decay_rate_log_base base of the logarithm
#' used to calculate the decay rate (default 10)
#' @return decay_ rate in logs of the given
#' base per time unit
#' @export
to_decay_rate <- function(fractional_life,
                          fraction = 1 /2,
                          decay_rate_log_base = 10){

    log_frac <- -log(fraction) / log(decay_rate_log_base)
    return(log_frac / fractional_life)
}


#' to_fractional_life
#'
#' convert a decay rate to an nth-life
#' (i.e. how long until the fraction in
#' question is all that remains)
#'
#' @param decay_rate exponential decay rate
#' @param fraction fraction to calculate for
#' @param decay_rate_log_base base of the logarithm
#' used to calculate the decay rate (default 10)
#' @return fractional life for the given fraction
#' @export
to_fractional_life <- function(decay_rate,
                               fraction,
                               decay_rate_log_base = 10){

    log_frac <- -log(fraction) / log(decay_rate_log_base)
    return(log_frac / decay_rate)
}


#' to_half_life
#'
#' convert a decay rate to a half_life
#'
#' @param decay_rate exponential decay rate
#' @param decay_rate_log_base base of the logarithm
#' used to calculate the decay rate (default 10)
#' @return fractional life for the given fraction
#' @export
to_half_life <- function(decay_rate,
                         decay_rate_log_base = 10){
    fraction <- 1 / 2
    return (to_fractional_life(decay_rate,
                               fraction,
                               decay_rate_log_base))
}



#' get_random_draws
#'
#' Given a tidy dataframe of
#' MCMC results, get a random
#' subsample of draws, returned
#' as another tidy dataframe
#' with the same scheme but fewer
#' rows
#'
#' This is important to have since the same draw
#' will occur on multiple lines of of the data frame
#' 
#' @param mcmc_draws mcmc chains in spread_draws tidy
#' tibble format from tidybayes
#' @param n_samples how many draws to sample
#' @param draw_id_column name of the draw id column,
#' default ".draws" (as is default in tidybayes)
#' @param replace whether to sample with replacement.
#' Default FALSE.
#' @return tibble of subsampled draws
#' @export
get_random_draws <- function(mcmc_draws,
                             n_samples,
                             draw_id_column = ".draw",
                             replace = FALSE) {

    chosen_draws <- sample(unique(mcmc_draws[[draw_id_column]]),
                           n_samples,
                           replace = replace)

    return (
        mcmc_draws |>
        dplyr::filter(get(draw_id_column) %in%
                      chosen_draws)
    )

}
