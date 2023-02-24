#!/usr/bin/env Rscript

########################################
## filename: figure-ratios.R
## author: Dylan H. Morris (dylanhmorris.com)
## plot figure showing ratios of half-lives
## compared to wild-type
#######################################


script_packages <- c(
    'rstan',      # stan interface
    'readr',      # tsv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # for plotting
    'tidyr',      # for crossing()
    'cowplot',    # publication ready ggplot
    'varstab'     # project functions
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only = TRUE))
}

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
omit_pcr <- grepl("no-pcr", outpath)


hl_chains <- readRDS(hl_path)
dat <- read_tsv(data_path,
                col_types = cols())
if (!("lineage" %in% names(dat))){
    dat <- dat |>
        mutate(
            lineage = lineage_name)
}

if( !("cell" %in% names(dat)) ){
    dat <- dat |>
        mutate(
            cell = "Vero_E6")
}


dat <- dat |>
    mutate(
        cell_display_name = ordered(
            case_when(
                cell == "Vero_E6" ~ "Vero E6",
                cell == "Vero_E6-TMPRSS2-T2A-ACE2" ~ "+TMPRSS2",
                cell == "Vero-TMPRSSII-RML" ~ "+TMPRSS2 (RML)",
                TRUE ~ NA),
            levels = rev(c(
                "Vero E6",
                "+TMPRSS2",
                "+TMPRSS2 (RML)"))))

#################################
## overall plot styling
#################################
ref_lineage <- "WA1"

##################################################
## calculate posterior draws for regression lines
##################################################

hl_draws <- hl_chains |>
    spread_draws(log_half_life[experiment_id]) |>
    add_titer_metadata(dat, "experiment_id") |>
    mutate(
        lineage_name = factor(lineage_name,
                              ordered = TRUE,
                              levels = variant_order()),
        variant_display_name =
            ifelse(
                lineage == lineage_name,
                lineage,
                paste0(lineage_name,
                       " (",
                       lineage,
                       ")"))) |>
    arrange(desc(lineage_name)) |>
    mutate(variant_display_name = factor(
               variant_display_name,
               ordered = TRUE,
               levels = unique(variant_display_name)))


hl_draws <- hl_draws |>
    inner_join(
        hl_draws |> filter(lineage == ref_lineage) |>
        select(.draw, cell, ref_log_hl = log_half_life),
        by = c(".draw", "cell")) |>
    ungroup() |>
    mutate(hl_diff = exp(log_half_life) - exp(ref_log_hl),
           hl_fold_change = exp(log_half_life - ref_log_hl))
    

ratio_plot <- hl_draws |>
    filter(lineage != ref_lineage) |>
    ggplot(aes(y = cell_display_name,
               x = hl_fold_change,
               fill = lineage_name)) +
    stat_default_eye(size = 8,
                     height = 1.5) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       trans = "log2") +
    geom_vline(xintercept = 1,
               linewidth = 2.5,
               linetype = "dashed",
               color = "black") +
    coord_cartesian(xlim = c(1/10, 10)) +
    scale_fill_variant() +
    theme_project() +
    theme(legend.position = "none") +
    ylab("Cell line") +
    xlab("Half-life fold change relative to WA1") +
    facet_grid(lineage_name~.)


ratio_tab <- hl_draws |>
    group_by(variant_display_name, cell) |>
    summarise(med_diff = median(hl_diff),
              q025_diff = quantile(hl_diff, 0.025),
              q975_diff = quantile(hl_diff, 0.975),
              med_fold_change = median(hl_fold_change),
              q025_fold_change = quantile(hl_fold_change, 0.025),
              q975_fold_change = quantile(hl_fold_change, 0.975),
              formatted_fold_change = format_est_ci(
                  med_fold_change,
                  q025_fold_change,
                  q975_fold_change))

              

####################################
## compose full figure from panels
####################################

if(figure){
    ## compose full figure from panels
    cat('making full figure...\n')

    ## save the plot to outpath
    cat('saving figure to ', outpath, '...\n')
    save_plot(outpath,
              ratio_plot,
              base_width = 10,
              base_height = 15)
} else if (table) {
    write_tsv(ratio_tab, outpath)
} else {
    stop("Ambiguous output format. Should be either figure or table")
}
