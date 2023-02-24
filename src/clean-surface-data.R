library(readr)   # delim I/O
library(dplyr)   # SQL joins etc
library(stringr) # str_ functions
library(tidyr)   # pivot_ functions

process_dilution_row <- function(row) {
    if(is.null(row)) {
        result <- rep(FALSE, 4)
    } else if (grepl("Empty row flag", row)) {
        result <- rep(NA, 4)
    } else {
        result <- sapply(
            1:4,
            function(x){grepl(x, row)})
    }
    
    return (result)
}


process_plate <- function(plate){
    
    output <- tibble(row_name = rep(NA, 8),
                     col1 = NA,
                     col2=NA,
                     col3=NA,
                     col4=NA,
                     sample_id=NA)

    dilution_rows <- plate %>%
        select(starts_with(
            "Input 32 well plate section"))
    sample_id <- plate["input sample id"]
    for (i_row in 1:length(dilution_rows)){
        row <- dilution_rows[i_row]
        row_name <- str_match(
            names(row),
            "\\[([A-Z])\\]")[2]
        row_wells <- process_dilution_row(row)
        
        output[i_row, ] <- c(row_name, row_wells, sample_id)
    }
    
    return (output)
}


args <- commandArgs(trailingOnly = TRUE)
data_path <- args[1]
out_path <- args[2]


if (!str_ends(data_path, ".tsv")) {
    stop("Must provide .tsv input path")
}
if (!str_ends(out_path, ".tsv")) {
    stop("Must provide .tsv output path")
}


data <- read_tsv(data_path)

data_as_plates <- list()

for(i_plate in 1:dim(data)[1]) {
    plate <- data[i_plate, ]
    sample_id <- as.character(plate["input sample id"])
    data_as_plates[[sample_id]] <- process_plate(plate)
}


## fix samples with issues
invalid_wells <- list(
    c("WA1_T0_4", 2, 2),
    c("WA1_T0_4", 2, 3),
    c("WA1_T0_4", 2, 4),
    c("WA1_T48_4", 1, 3),
    c("Beta_T72_4", 1, 2),
    c("Beta_T72_4", 1, 3))

for (well in invalid_wells){
    cat("Removing invalid well", well, "\n\n")
    name <- well[1]
    row <- as.numeric(well[2])
    col <- as.numeric(well[3]) + 1

    data_as_plates[[name]][row, col] <- NA
}


## melt plate form data into long-form
results <- tibble()
for(plate in data_as_plates){
    results <- results %>%
        bind_rows(
            plate %>%
            pivot_longer(starts_with("col")))
}


results <- results |>
    rowwise() |>
    mutate(
        dilution = -(which(
            LETTERS == row_name)))

results <- results |>
    mutate(
        sample_name = sample_id,
        split_name = str_split(
            sample_name,
            "_",
            simplify = TRUE),
        variant = split_name[, 1],
        variant = ifelse(variant == "RML7",
                         "B.1",
                         variant),
        time = as.numeric(
            str_sub(split_name[, 2], 2)),
        replicate = as.numeric(
            split_name[, 3])) %>%
    rename(well_status = value)


filtered <- results |> 
    filter(!is.na(well_status))


dat <- filtered |>
    arrange(variant, time, replicate) |>
    group_by(sample_name) |>
    mutate(sample_id = cur_group_id()) |>
    ungroup() |>
    arrange(variant) |>
    group_by(variant) |>
    mutate(experiment_id = cur_group_id()) |>
    ungroup()


dat <- dat |>
    group_by(sample_id, dilution) |>
    mutate(CPE_full = all(well_status))


summary_dat <- dat |>
    group_by(variant, replicate, time, sample_id) |>
    summarise(
        last_full_CPE = min(max(dilution) + 1,
                            min(dilution[CPE_full])),
        wells_past = sum(well_status[dilution <= last_full_CPE]))

print(summary_dat,
      n = 1000)

write_tsv(
    summary_dat,
    paste0(out_path, "-summary.tsv"))

dat <- dat |>
    select(
        sample_name,
        sample_id,
        experiment_id,
        lineage_name = variant,
        time,
        replicate,
        dilution,
        well_status,
        row_name,
        col_name = name,
        CPE_full)

write_tsv(
    dat,
    out_path)
