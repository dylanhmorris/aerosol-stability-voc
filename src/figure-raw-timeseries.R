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


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly = TRUE)

n_args <- length(args)

outpath <- args[n_args]
titers_path <- args[n_args - 1]
data_path <- args[1]
LOD_log10_per_ml <- 0.50
LOD <- 10^LOD_log10_per_ml

#################################
## overall plot styling
#################################
set.seed(234) # reproducible! (since we use random draws)

intercept_id_col <- "run_id"
log10_ct_rel_t0 <- NA # only get a value if present in data
titer_ylab <- expression("Virus titer (TCID"[50] * "/L air)")
convert_units <- 10/3
scale_x <- scale_x_continuous(
    expand = c(0.075, 0.075),
    breaks = c(0, 2, 4, 6, 8))
xlim <- c(0, 8)

## read data / style files

cat("reading data (this may take a while)...\n")
ests_chains <- readRDS(titers_path)

dat <- read_tsv(data_path,
                col_types = cols())

dat <- dat |>
    mutate(log_copy_adjustment = log10_copies_rel_t0)

cat("extracting positive wells...\n")
pos_wells <- dat |>
    group_by(sample_id) |>
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

## conditional styling
pre_spray <- grepl("pre-spray", outpath)
if ( grepl("tmprss2-reg", outpath) ) {
    cell_line <- "Vero_E6-TMPRSS2-T2A-ACE2"
} else if( grepl("tmprss2-rml", outpath) ) {
    cell_line <- "Vero-TMPRSSII-RML"
} else {
    cell_line <- "Vero_E6"
}

if(pre_spray){
    xlim[1] <- -2
}

cat("data read succesfully!\n")



#######################
## get posterior draws
#######################

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
## adding one to convert to per mL from per 0.1 mL
plot_draws <- titer_ests_draws |>
    mutate(
        lineage_name = factor(lineage_name,
                              ordered = TRUE,
                              levels = variant_order()),
        log10_titer_per_ml = ifelse(
            detectable,
            sampled_titer + 1 - log_copy_adjustment,
            LOD_log10_per_ml)) |>
    filter(cell == cell_line)

if(!pre_spray){
    plot_draws <- plot_draws |>
        filter(!pre_spray)
}

largest_titers <- plot_draws |>
    group_by(sample_id) |>
    summarise(
        q95 = quantile(log10_titer_per_ml, 0.95, na.rm=TRUE))

largest_titer <- max(largest_titers$q95)

## y limits are dynamic based on largest titer
ylim <- c(
    LOD * convert_units / 10,
    convert_units * 10^largest_titer)


shape_scale <- scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))

panel <- plot_draws |>
    ggplot(aes(
        x = time,
        y = convert_units * 10^log10_titer_per_ml,
        fill = lineage_name)) +
    geom_hline(
        aes(yintercept = convert_units * 10^LOD_log10_per_ml),
        linewidth = 2,
        linetype = "dotted") +
    titer_pointinterval(
        mapping = aes(
            group = sample_id,
            shape = detectable,
            point_fill = lineage_name),
        size = 17) +
    stat_summary(
        mapping = aes(group = run_id),
        fun = median,
        geom = "line",
        color = "grey") +
    stat_summary(
        mapping = aes(group = run_id),
        fun = median,
        geom = "point",
        color = "grey") +
    scale_fill_variant() +
    scale_fill_variant(aesthetics = "point_fill") +
    scale_color_variant() +
    shape_scale +
    facet_wrap(~lineage_name) +
    scale_y_log10_mathformat() +
    scale_x +
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
