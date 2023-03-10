% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R, R/plotting_style.R
\name{titer_pointinterval}
\alias{titer_pointinterval}
\title{titer_pointinterval}
\usage{
titer_pointinterval(
  size = 28,
  point_alpha = 0.75,
  interval_alpha = 0.75,
  stroke_mult = 0.1,
  .width = 0.95,
  ...
)

titer_pointinterval(
  size = 28,
  point_alpha = 0.75,
  interval_alpha = 0.75,
  stroke_mult = 0.1,
  .width = 0.95,
  ...
)
}
\arguments{
\item{size}{overall size of the geom (default 30)}

\item{point_alpha}{transparency of the point estimate
point (default 0.75)}

\item{interval_alpha}{transparency of the interval lines
(default 0.75)}

\item{stroke_mult}{linewidth for the point estimate point
border relative to the size (default 0.1, so a width of
2.8 with the default size of 28).}

\item{.width}{width of the interval default 0.95}

\item{...}{other keyword arguments passed to stat_pointinterval()}

\item{interval_size_ratio}{how much greater should the
linewidth of the smaller (default 68\%) credible interval
be compared to that of the larger (default 95\%) credible
interval? Default 2.5}
}
\value{
a stat_pointinterval object that
can be added to a ggplot object

a stat_pointinterval object that
can be added to a ggplot object
}
\description{
default stat_pointinterval for estimated
virus titers

default stat_pointinterval for estimated
virus titers
}
