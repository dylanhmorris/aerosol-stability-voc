#!/usr/bin/env Rscript

########################################
## filename: figure-main.R
## author: Dylan H. Morris (dylanhmorris.com)
## plot main figure
#######################################


script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggdist',     # stat_halfeye.
    'ggplot2',    # for plotting
    'tidyr',      # for crossing()
    'cowplot',    # publication ready ggplot
    'varstab'    # plotting style and functions
)


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


## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
hl_path <- args[2]
titer_path <- args[3]
outpath <- args[4]

figure <- grepl("figure", outpath)
table <- grepl("table", outpath)

hl_chains <- readRDS(hl_path)
titer_chains <- readRDS(titer_path)
dat <- read_tsv(data_path,
                col_types = cols())

#################################
## overall plot styling
#################################

set.seed(989327) # reproducible! (since we use random draws)
n_draws_to_sample <- 50 
n_intercepts_to_sample <- 6
LOD_log10_per_ml <- 0.5
LOD <- 10^LOD_log10_per_ml

line_alpha <- 0.05
regression_linewidth <- 1
density_alpha <- plotting_params("density_alpha")
density_interval_size <- 15
density_point_size <- 12
figure_theme_base <- (
    theme_project(base_size = 45) + 
    theme(legend.position = "none",
          panel.spacing.x = unit(2.5, "lines"),
          panel.spacing.y = unit(2, "lines"),
          plot.margin = margin(l = 0.05, r = 0.1,
                               unit = "in")))

# make sure point of halfeye plot pointintervals are visible
density_y_scale <- scale_y_continuous(expand = c(0.125, 0))

surface <- (hl_chains@model_name == "infer_halflives")
aerosol <- (hl_chains@model_name ==
            "infer_timeseries_halflives")
omit_pcr <- grepl("no-pcr", outpath)
surface_linear <- grepl("linear", outpath)

if( !(surface | aerosol) ) {
    stop(sprintf("Unknown model name %s",
                 hl_chains@model_name))
}

cat("Processing fitting data...\n")
## draw n_lines random regression lines
hl_x_scale <- scale_x_continuous(
    expand = c(0, 0),
    breaks = c(1, 2, 4, 8, 16, 32),
    trans = "log2")
hl_xlim <- c(1/2, 64)

ratio_xlim <- c(1/12, 12)

if(surface_linear){
    hl_x_scale <- scale_x_continuous(
        expand = c(0, 0),
        breaks = c(0, 2, 4, 6, 8))
    hl_xlim <- c(0, 8)
}


if (aerosol) {
    intercept_id_col <- "run_id"
    titer_ylab <- expression("Virus titer (TCID"[50] * "/L air)")
    convert_units <- 10/3
    log10_ct_rel_t0 <- NA # only get a value if present in data
    xmax <- 8
    ylim <- c(
        LOD * convert_units / 10,
        convert_units * 1e4)
    scale_x_regression <- scale_x_continuous()
} else if (surface) {
    intercept_id_col <- "sample_id"
    titer_ylab <- expression("Virus titer (TCID"[50] * "/mL media)")
    convert_units <- 1
    log10_ct_rel_t0 <- 0  ## override any data value
    xmax <- 96
    scale_x_regression <- scale_x_continuous(
        breaks = c(0, 24, 48, 72, 96))
    ylim <- c(
        LOD * convert_units / 10,
        convert_units * 1e6)
}


if ( grepl("tmprss2-reg", outpath) ) {
    cell_line <- "Vero_E6-TMPRSS2-T2A-ACE2"
    ratio_xlim <- c(1/20, 20) ## this one spills over a bit
    ylim[2] <- convert_units * 1e6  ## this one spills over a bit
} else if( grepl("tmprss2-rml", outpath) ) {
    cell_line <- "Vero-TMPRSSII-RML"
} else {
    cell_line <- "Vero_E6"
}

ref_lineage <- "WA1"


if( !("cell" %in% names(dat)) ){
    dat <- dat |>
        mutate(
            cell = "Vero_E6")
}

if( !("well_status" %in% names(dat)) ){
    dat <- dat |>
        mutate(
            well_status = virus_detect)
}

if (!("lineage" %in% names(dat))){
    dat <- dat |>
        mutate(
            lineage = lineage_name)
}
if (!("pre_spray" %in% names(dat))){
    dat <- dat |>
        mutate(
            pre_spray = FALSE)
}

## handle qPCR adjustment
dat <- dat |>
    mutate(log_copy_adjustment = log10_ct_rel_t0)

if(omit_pcr) {
    dat <- dat |>
        mutate(log_copy_adjustment = 0)
}

dat <- dat |>
    mutate(
        lineage_name = ordered(lineage_name,
                               levels = variant_order()),
        variant_display_name = as.character(lineage_name)) |>
    arrange(lineage_name) |>
    mutate(variant_display_name = ordered(
               variant_display_name,
               levels = unique(variant_display_name)))


if(! (dat$sample_id |> unique() |> length() == max(dat$sample_id))) {
    print(dat |> select(sample_id) |> arrange() | unique(),
          n = 1000)
    stop("Missing sample ids!")
}

##################################################
## calculate posterior draws for half lives
##################################################
cat("Processing MCMC draws...\n")

hl_draws <- hl_chains |>
    spread_draws(log_half_life[experiment_id]) |>
    add_titer_metadata(dat, "experiment_id",
                       verbose_errors = TRUE)




## get needed draws and add human readable names

if(aerosol){
    draws <- hl_chains |>
        spread_draws(run_intercept[run_id]) |>
        rename(intercept = run_intercept) |>
        add_titer_metadata(dat, id_column = "run_id") |>
        mutate(intercept_id = run_id)
} else if (surface) {
    draws <- hl_chains |>
        spread_draws(intercept[sample_id]) |>
        add_titer_metadata(dat, id_column = "sample_id") |>
        mutate(intercept_id = sample_id)
}


draws <- hl_chains |>
    spread_draws(log_half_life[experiment_id]) |>
    inner_join(draws,
               by = c("experiment_id", ".draw"),
               multiple = "all") |>
    mutate(decay_rate = to_decay_rate(exp(log_half_life)))

pos_wells <- dat |>
    group_by(sample_id) |>
    summarise(
        n_wells = n(),
        n_pos = sum(well_status))    


titer_draws <- titer_chains |>
    spread_draws(sampled_titer[sample_id]) |>
    add_titer_metadata(dat, id_column = "sample_id") |>
    inner_join(pos_wells,
               by = "sample_id",
               multiple = "all") |>
    mutate(detectable = n_pos >= 1) |>
    filter(cell == cell_line) |>
    filter(!pre_spray)


## convert from TCID50/(0.1mL) to TCID50/mL
## and visualize 0 positive well titers at
## the traditional LOD, and
## qPCR adjust.
titer_draws <- titer_draws |>
    mutate(log10_titer_per_mL = ifelse(
               detectable,
               sampled_titer + 1 - log_copy_adjustment,
               LOD_log10_per_ml)) |>
    ungroup() |>
    arrange(desc(time))

largest_titers <- titer_draws |>
    group_by(sample_id) |>
    summarise(
        q95 = quantile(log10_titer_per_mL, 0.95, na.rm=TRUE))





###################################
## plot panel showing half-life
## estimates
###################################
hl_plot <- hl_draws |>
    filter(cell == cell_line) |>
    ggplot(aes(x= exp(log_half_life),
               fill = lineage_name)) +
    stat_default_eye(
        side="top") +
    coord_cartesian(xlim = hl_xlim) +
    hl_x_scale +
    density_y_scale +
    scale_fill_variant() +
    ylab("Posterior density") +
    xlab("Half-life (h)\n") +
    facet_grid(
        variant_display_name~.) + 
    figure_theme_base + 
    theme(strip.text = element_blank())

###################################
## plot panel showing raw
## data and regression lines
###################################
cat('plotting fit of regression to raw data...\n')
shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))


draws <- draws |>
    get_random_draws(n_draws_to_sample) |> 
    group_by(.draw, experiment_id) |> # now get random intercept
    dplyr::slice_sample(n = n_intercepts_to_sample) |>
    ungroup()

to_plot <- draws |>
    mutate(unlogged_intercept = (
        10^(1 + intercept) *
        convert_units)) |>
    filter(cell == cell_line)


fit_panel <- to_plot |>
    ggplot(aes(x = time,
               group = interaction(.draw, intercept_id))) +
    geom_hline(aes(yintercept = LOD * convert_units),
               linewidth = 2,
               linetype = "dotted") +
    stat_exp_curve(        
        aes(intercept = unlogged_intercept,
            rate = -decay_rate,
            color = lineage_name),
        base = 10,
        xmin = 0,
        xmax = xmax,
        alpha = line_alpha,
        linewidth = regression_linewidth) +
    titer_pointinterval(
        mapping = aes(
            x = time,
            y = (10^log10_titer_per_mL) * convert_units,
            shape = detectable,
            point_fill = lineage_name,
            group = sample_id),
        data = titer_draws) +
    scale_fill_variant() +
    scale_fill_variant(aesthetics = "point_fill") +
    scale_color_variant() +
    shape_scale + 
    scale_y_log10_mathformat() +
    scale_x_regression +
    coord_cartesian(
        ylim = ylim,
        xlim = c(0, xmax)) +
     facet_grid(
        variant_display_name~.)

# styling: no facet labels because is internal
fit_panel <- fit_panel +
    figure_theme_base + 
    theme(legend.position = "none") +
    xlab("Time (h)\n") +
    ylab(titer_ylab) +
    theme(strip.text = element_blank())


hl_rat <- hl_draws |>
    inner_join(
        hl_draws |> filter(lineage == ref_lineage) |>
        select(.draw, cell, ref_log_hl = log_half_life),
        by = c(".draw", "cell")) |>
    ungroup() |>
    mutate(hl_diff = exp(log_half_life) - exp(ref_log_hl),
           hl_fold_change = exp(log_half_life - ref_log_hl))


ratio_plot <- hl_rat |>
    filter(cell == cell_line) |>
    ggplot(aes(x = hl_fold_change,
               fill = lineage_name)) +
    stat_default_eye(
        side = "top") +
    geom_vline(
        linewidth = 2,
        linetype = "dashed",
        xintercept = 1) +
    density_y_scale +
    scale_x_continuous(expand = c(0, 0),
                       breaks = c(1/8, 1/4, 1/2, 1, 2, 4, 8),
                       labels = c(
                           expression(scriptstyle(frac(1, 8))),
                           expression(scriptstyle(frac(1, 4))),
                           expression(scriptstyle(frac(1, 2))),
                           1, 2, 4, 8),
                       trans = "log2") +
    coord_cartesian(xlim = ratio_xlim) +
    scale_fill_variant() +
    figure_theme_base + 
    ylab("") +
    xlab("Half-life fold change\nrelative to WA1") +
    facet_grid(
        variant_display_name~.)


cat('making full figure...\n')
full_fig <- plot_grid(
    fit_panel +
    theme(plot.margin = margin(
              l = 0.1,
              r = 0.15,
              b = 0.1,
              t = 0.35,
              unit = "in")),
    hl_plot +
    theme(plot.margin = margin(l = 0.15,
                               unit = "in")),
    ratio_plot +
    theme(plot.margin = margin(l = 0,
                               unit = "in")),
    label_size = 50,
    align = "h",
    axis = "bt",
    ncol = 3,
    labels = c('A', 'B', 'C'))

cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = 20,
          base_asp = 1.3)
warnings()

