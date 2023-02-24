#!/usr/bin/env Rscript

########################################
## filename: figure-halflives.R
## author: Dylan H. Morris (dylanhmorris.com)
## plot figure showing half-lives
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

##################################################
## calculate posterior draws for half lives
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

hl_plot <- hl_draws |>
    ggplot(aes(y = cell_display_name,
               x = exp(log_half_life),
               fill = lineage_name)) +
    stat_default_eye(size=8,
                     height=1.5) +
    coord_cartesian(xlim = c(1/2, 40)) +
    scale_x_continuous(expand = c(0, 0.25),
                       breaks = c(1, 2, 4, 8, 16, 32),
                       trans = "log2") +
    scale_fill_variant() +
    theme_project() +
    theme(legend.position = "none") +
    ylab("Cell line") +
    xlab("Half-life (h)") +
    facet_grid(lineage_name~.)

hl_tab <- hl_draws |>
    group_by(variant_display_name, cell) |>
    summarise(med = median(exp(log_half_life)),
              q025 = quantile(exp(log_half_life), 0.025),
              q975 = quantile(exp(log_half_life), 0.975),
              formatted_hl = format_est_ci(
                  med,
                  q025,
                  q975))

              

if(figure){
    ## compose full figure from panels
    cat('making full figure...\n')

    ## save the plot to outpath
    cat('saving figure to ', outpath, '...\n')
    save_plot(outpath,
              hl_plot,
              base_width = 10,
              base_height = 15)
} else if (table) {
    write_tsv(hl_tab, outpath)
} else {
    stop("Ambiguous output format. Should be either figure or table")
}
warnings()
