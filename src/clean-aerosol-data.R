#!/usr/bin/env Rscript

#####################################
## name: clean-data.R
## author: Dylan Morris <dylanhmorris.com>
##
## process raw data and save cleaned
## data
####################################


suppressPackageStartupMessages(library(readr))   # delim I/O
suppressPackageStartupMessages(library(dplyr))   # SQL joins etc
suppressPackageStartupMessages(library(stringr)) # str_replace methods
suppressPackageStartupMessages(library(varstab)) # pango_to_who

## functions

clean_data <- function(dat_path,
                       pcr_path,
                       outpath,
                       without_pcr = FALSE,
                       quantify_initial_loss = FALSE) {

    cat("Reading data...\n")
    dat <- read_delim(dat_path,
                      delim = ";")
    pcr <- read_delim(pcr_path,
                      delim = ";")


    ## for now use only one PCR run per sample
    pcr <- pcr |>
        group_by(sample_id) |>
        filter(replicate_PCR == max(replicate_PCR)) |>
        ungroup()

    ## strip whitespace
    dat <- dat |>
        mutate(sample_id = str_replace_all(sample_id, " ", ""))
    pcr <- pcr |>
        mutate(sample_id = str_replace_all(sample_id, " ", ""))
    
    if(!without_pcr){
        dat <- dat |>
            inner_join(pcr |> select(sample_id,
                                     CT,
                                     quantity_per_ml),
                       by = "sample_id")
    } else {
        dat <- dat |> mutate(CT = 0)
    }
    
    cat("Processing data...\n")
    dat <- dat |>
        filter(!is.na(replicate))

    ## clean variant names
    dat <- dat |>
        mutate(lineage_name = pango_to_who(
                   lineage,
                   pango_if_fail = TRUE))

    ## fix typos in cell line names
    dat <- dat |>
        mutate(cell = case_when(
                   cell == "VERO E6" ~ "Vero_E6",
                   cell == "VERO TMPRSII" ~ "Vero_E6-TMPRSS2-T2A-ACE2",
                   cell == "VERO TMPRSII-RML" ~ "Vero-TMPRSSII-RML",
                   TRUE ~ NA))

    dat <- dat |>
        group_by(virus,
                 lineage,
                 xp,
                 temperature,
                 cell,
                 humidity,
                 material) |>
        mutate(experiment_id = cur_group_id()) |>
        ungroup()
    
    dat <- dat |>
        arrange(material) |>
        group_by(material) |>
        mutate(material_id = cur_group_id()) |>
        ungroup()
    
    dat <- dat |>
        arrange(cell) |>
        group_by(cell) |>
        mutate(cell_id = cur_group_id()) |>
        ungroup()

    dat <- dat |>
        mutate(
            pre_spray = time < 0)
    
    dat <- dat |>
        arrange(desc(pre_spray),
                lineage,
                cell_id,
                replicate,
                time) |>
        group_by(pre_spray,
                 experiment_id,
                 replicate,
                 time) |>
        mutate(sample_id = cur_group_id(),
               titration_run_id =
                   as.numeric(factor(rank(replicate_titration,
                                          ties.method = "min"),
                                     ordered=TRUE))) |>
        ungroup()
    

    dat <- dat |>
        arrange(lineage, cell_id, replicate, time) |>
        group_by(experiment_id,
                 replicate) |>
        mutate(run_id = cur_group_id()) |>
        ungroup()

    dat <- dat |>
        group_by(run_id) |>
        arrange(pre_spray, time, .by_group = TRUE) |>
        mutate(log10_ct_rel_t0 = log10(2) * (first(CT) - CT),
               log10_copies_rel_t0 = log10(quantity_per_ml) -
                   log10(first(quantity_per_ml))) |>
        ungroup()
    
    dat <- dat |> arrange(sample_id, titration_run_id)
    
    print(dat)
    
    cat("Saving cleaned data to ", outpath)
    write_tsv(dat,
              outpath)

}

## load and clean data
args <- commandArgs(trailingOnly = TRUE)
raw_data_path <- args[1]
pcr_data_path <- args[2]
out_path <- args[3]

without_pcr = grepl("no_pcr", out_path)
quantify_initial_loss = grepl("initial_loss", out_path)

clean_data(raw_data_path,
           pcr_data_path,
           out_path,
           without_pcr = without_pcr,
           quantify_initial_loss = quantify_initial_loss)

warnings()
