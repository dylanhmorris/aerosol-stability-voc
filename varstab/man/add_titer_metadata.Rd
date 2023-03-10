% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/postprocess.R
\name{add_titer_metadata}
\alias{add_titer_metadata}
\title{add_titer_metadata}
\usage{
add_titer_metadata(
  mcmc_draws,
  fitting_data,
  id_column = "titer_id",
  verbose_errors = FALSE,
  multiple = "all"
)
}
\arguments{
\item{mcmc_draws}{mcmc chains in spread_draws tibble
format from tidybayes}

\item{fitting_data}{data used to fit the mcmc model
as a tibble with a titer id column}

\item{id_column}{name of the id column, default "titer_id"}

\item{multiple}{= {"all", NULL, "error"}. Flag passed
to dplyr::inner_join(). Whether to expect multiple
rows of the output table per row of the index table.
As we almost always do, this defaults to
"all" to avoid a dplyr warning.}
}
\value{
tibble of draws with metadata added
}
\description{
add metadata to a set of MCMC chains
ordered by titer_id
}
