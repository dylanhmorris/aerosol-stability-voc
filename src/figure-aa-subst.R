#!/usr/bin/env Rscript

########################################
## filename: figure-aa-subst.R
## author: Dylan H. Morris (dylanhmorris.com)
## plot figure showing ratios of half-lives
## compared to wild-type versus
## amino acid substitutions rel
## to wild-type
#######################################


script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # for plotting
    'tidyr',      # for crossing()
    'cowplot',    # publication ready ggplot
    'varstab'     # plotting style
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
mutation_path <- args[3]
outpath <- args[4]

cat("Reading in data...\n")
hl_chains <- readRDS(hl_path)
dat <- read_tsv(data_path,
                col_types = cols())

mutation_dat <- read_tsv(mutation_path,
                         col_types = cols())

mutation_counts <- mutation_dat |>
    group_by(pango_lineage) |>
    summarise(substitution_count = n()) |>
    rename(lineage = pango_lineage) |>
    bind_rows(
        tibble(lineage="WA1", substitution_count=0))


#################################
## overall plot styling
#################################
ref_lineage = "WA1"

##################################################
## calculate posterior draws for regression lines
##################################################

hl_draws <- hl_chains |>
    spread_draws(log_half_life[experiment_id]) |>
    add_titer_metadata(dat, "experiment_id") |>
    inner_join(mutation_counts,
               by = "lineage") |>
    arrange(substitution_count, lineage) |>
    mutate(
        lineage_name = ordered(
            lineage_name,
            levels = unique(lineage_name))
    )
print(hl_draws |> select(lineage_name) |> unique())

hl_draws <- hl_draws |>
    inner_join(
        hl_draws |> filter(lineage == ref_lineage) |>
        select(.draw, cell, ref_log_hl = log_half_life),
        by = c(".draw", "cell")) |>
    ungroup() |>
    mutate(hl_diff = exp(log_half_life) - exp(ref_log_hl),
           hl_fold_change = exp(log_half_life - ref_log_hl))



hl_plot <- hl_draws |>
    filter(cell == "Vero_E6") |>
    filter(lineage != ref_lineage) |>
    ggplot(aes(x = substitution_count,
               y = hl_fold_change,
               fill = lineage_name)) +
    stat_default_eye() +
    geom_hline(yintercept = 1,
               linewidth = 2,
               linetype = "dashed") +
    scale_x_continuous() +
    scale_y_continuous(expand = c(0, 0),
                       breaks = c(0.25, 0.5, 1, 2, 4),
                       trans = "log2") +
    coord_cartesian(xlim = c(0, 50),
                    ylim = c(1/10, 10)) +
    scale_fill_variant(name="") +
    theme_project(base_size = 35) +
    theme(legend.position = "bottom",
          legend.key.width = unit(1, "mm"),
          legend.text = element_text(size=25),
          legend.direction = "horizontal") +
    guides(fill = guide_legend(
               override.aes = list(size = 25))) + 
    xlab("Substitutions relative to WA1") +
    ylab("Half-life fold change relative to WA1")


####################################
## compose full figure from panels
####################################

## compose full figure from panels
cat('making full figure...\n')

    ## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          hl_plot,
          base_width = 14,
          base_height = 12)
