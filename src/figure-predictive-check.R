#!/usr/bin/env Rscript

########################################
## filename: figure-predictive-check.R
## author: Dylan Morris <dylanhmorris.com>
## plot predictive checks for titer inference
#######################################

suppressPackageStartupMessages(library(varstab))   # plotting style
suppressPackageStartupMessages(library(ggplot2))   # plotting
suppressPackageStartupMessages(library(dplyr))     # SQL joins
suppressPackageStartupMessages(library(tidyr))     # unnest()
suppressPackageStartupMessages(library(cowplot))   # save_plot()
suppressPackageStartupMessages(library(readr))     # I/O
suppressPackageStartupMessages(library(tidybayes)) # spread_draws()


## read command line args
args <- commandArgs(trailingOnly = TRUE)

n_args <- length(args)

outpath <- args[n_args]
check_results_path <- args[n_args - 1]
titers_path <- args[n_args - 2]
data_path <- args[1]

#################################
## overall plot styling
#################################
set.seed(34634) # reproducible! (since we use random draws)
LOD_log10_per_ml <- 0.50
n_draws_to_sample <- 50 
n_intercepts_to_sample <- 6

#################################
# read in needed data / style
#################################

cat("reading data (this may take a while)...\n")
check_chains <- readRDS(check_results_path)
ests_chains <- readRDS(titers_path)
dat <- read_tsv(data_path,
                col_types = cols())

## get parameters / manipulate data
surface <- (check_chains@model_name == "infer_halflives")
aerosol <- (check_chains@model_name ==
            "infer_timeseries_halflives")
if( !(surface | aerosol) ){
    stop(sprintf("Unknown model name %s",
                 check_chains@model_name))
}
is_prior <- grepl("prior-check", outpath)


if(! "virus_detect" %in% names(dat))
    dat <- dat |> rename(virus_detect = well_status)
if(! "cell" %in% names(dat))
    dat <- dat |> mutate(cell = "Vero_E6")

if(! "lineage" %in% names(dat))
    dat <- dat |> mutate(lineage = lineage_name)

if(! "pre_spray" %in% names(dat))
    dat <- dat |> mutate(pre_spray = FALSE)
shape_scale <- scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))



## define plot-type specific variables
if(surface){
    intercept_id_col <- "sample_id"
    log10_ct_rel_t0 <- 0  ## override any data value
    titer_ylab <- expression("Virus titer (TCID"[50] * "/mL media)")
    convert_units <- 1
    scale_x_checks <- scale_x_continuous(
        breaks = c(0, 24, 48, 72, 96),
        expand = c(0.075, 0.075))
    xmax <- 110
    xmax_disp <- 100
} else if (aerosol) {
    intercept_id_col <- "run_id"
    log10_ct_rel_t0 <- NA # only get a value if present in data
    titer_ylab <- expression("Virus titer (TCID"[50] * "/L air)")
    convert_units <- 10/3
    scale_x_checks <- scale_x_continuous(
        expand = c(0.075, 0.075))
    xmax <- 9
    xmax_disp <- 8
}

if(aerosol){
} else if (surface) {
}

if (is_prior) {
    ylim <- c(1e-2, 1e9)
} else if (surface) {
    ylim <- c(1e-1, 1e6) * convert_units
} else {
    ylim <- c(1e-1, 1e4) * convert_units
}

xlim <- c(0, xmax_disp)
autoscale <- sqrt(xmax / 2)
line_alpha = 0.15


if ( grepl("tmprss2-reg", outpath) ) {
    cell_line <- "Vero_E6-TMPRSS2-T2A-ACE2"
    ylim[2] <- convert_units * 1e6  ## this one spills over a bit
} else if( grepl("tmprss2-rml", outpath) ) {
    cell_line <- "Vero-TMPRSSII-RML"
} else {
    cell_line <- "Vero_E6"
}


dat <- dat |>
    mutate(log_copy_adjustment = log10_ct_rel_t0)

cat("data read succesfully!\n")



##################################################
## calculate posterior draws for regression lines
##################################################

cat("extracting draws for decay rates / intercepts (this may also take a while)...\n")


tidy_draws <- check_chains |>
    spread_draws(intercept_pred[!!intercept_id_col])


tidy_draws <- tidy_draws |>
    add_titer_metadata(dat, intercept_id_col) |>
    ungroup() |>
    mutate(
        lineage_name = factor(lineage_name,
                              ordered = TRUE,
                              levels = variant_order()),
        line_id = interaction(.draw, get(intercept_id_col)))


hl_draws <- check_chains |>
    spread_draws(decay_rate[experiment_id]) |>
    select(
        decay_rate,
        .draw,
        experiment_id)


tidy_draws <- tidy_draws |>
    inner_join(
        hl_draws,
        by = c(".draw", "experiment_id")) |>
    filter(cell == cell_line)




cat("extracting positive wells...\n")
pos_wells <- dat |>
    group_by(sample_id) |>
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

cat('extracting titer estimates...\n')

titer_ests_draws <- ests_chains |>
    spread_draws(sampled_titer[sample_id])

## get human readable names and detectability
titer_ests_draws <- titer_ests_draws |>
    add_titer_metadata(dat, "sample_id") |>
    inner_join(pos_wells,
               by = "sample_id") |>
    mutate(detectable = n_pos > 0)

## sort by time
titer_ests_draws <- titer_ests_draws |>
    arrange(time)

## calculate LOD titers
titer_ests_draws <- titer_ests_draws |>
    mutate(
        lineage_name = factor(lineage_name,
                              ordered = TRUE,
                              levels = variant_order()),
        log10_titer_per_ml = ifelse(
            detectable,
            sampled_titer + 1 - log_copy_adjustment,
            LOD_log10_per_ml)) |>
    filter(cell == cell_line) |>
    filter(!pre_spray)




cat("Randomly downsampling predictive checks...\n")

downsampled_draws <- tidy_draws |>
    get_random_draws(n_draws_to_sample) |> 
    group_by(.draw, experiment_id) |> # now get random intercept
    dplyr::slice_sample(n = n_intercepts_to_sample) |>
    ungroup()


panel <- downsampled_draws |>
    ggplot(aes(
        x = time,
        fill = lineage_name)) +
    geom_hline(
        aes(yintercept = convert_units * 10^LOD_log10_per_ml),
        linewidth = 2,
        linetype = "dotted") +
    stat_exp_curve(        
        aes(intercept = convert_units * 10^(1 + intercept_pred),
            rate = -decay_rate,
            color = lineage_name,
            group = line_id),
        base = 10,
        xmin = 0,
        xmax = xmax,
        alpha = line_alpha) +
    titer_pointinterval(
        mapping = aes(x = time,
                      y = convert_units * 10^log10_titer_per_ml,
                      shape = detectable,
                      point_fill = lineage_name,
                      group = sample_id),
        data = titer_ests_draws,
        size = 17) +
    scale_fill_variant() +
    scale_fill_variant(aesthetics = "point_fill") +
    scale_color_variant() +
    shape_scale +
    facet_wrap(~lineage_name) +
    scale_y_log10_mathformat() +
    scale_x_checks +
    coord_cartesian(xlim = xlim,
                    ylim = ylim) +
    theme_project(base_size = 20) +
    theme(legend.position = "none") +
    labs(x = "Time (h)",
         y = titer_ylab)


####################################
## compose full figure from panels
####################################


margin <- theme(
    plot.margin = margin(b = 0.5, t = 0.5, l = 1, r = 1, unit = "cm"))
    
cat('making full figure...\n')

full_fig <- panel + margin

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = 7,
          base_asp = 1.6)
warnings()
