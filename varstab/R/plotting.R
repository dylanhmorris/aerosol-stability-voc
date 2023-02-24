
##################################################
## filename: plotting.R
## author: Dylan H. Morris <dylan@dylanhmorris.com>
##
## description: functions for visualizing
## variant aerosol stability model output
##################################################

#' StatExpCurve
#' 
#' Stat class for exponential growth/decay
#' curves in ggplot
StatExpCurve <- ggplot2::ggproto(
  "StatExpCurve",
  ggplot2::Stat, 
  required_aes =
      c("intercept", "rate"),
  compute_group =
      function(data,
               scales,
               params,
               xmin = 0,
               xmax = 1,
               n = 50,
               base = exp(1)) {
          grid <- data.frame(
              x = seq(xmin, xmax, length = n))
          
          raw_yval <- data$intercept * exp(log(base) * (data$rate * grid$x))
          ytrans <- scales$y$transform
          if (is.null(ytrans)) {
              ytrans <- function(x){x}
          }
          grid$y <- ytrans(raw_yval)
          
          return (grid)
      }
  )


#' Plot an exponential decay curve using
#' ggplot
#'
#' @param n number of distinct points to plot (default 50)
#' @param xmin first x-value to plot (default 0)
#' @param xmax last x-value to plot (default 1)
#' @param base base of the exponential (default e)
#'
#' @export
stat_exp_curve <- function(mapping = NULL,
                           data = NULL,
                           geom = "line",
                           position = "identity",
                           na.rm = FALSE,
                           show.legend = NA, 
                           inherit.aes = TRUE,
                           n = 50,
                           xmin = 0,
                           xmax = 1,
                           base = exp(1),
                           ...) {
  layer(
      stat = StatExpCurve,
      data = data,
      mapping = mapping,
      geom = geom, 
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(n = n,
                    xmin = xmin,
                    xmax = xmax,
                    base = base,
                    na.rm = na.rm,
                    ...)
  )
}



#' stat_default_eye
#'
#' default eyeplot for figures in this project
#'
#' @param size size argument for the eyeplot, default 20
#' @param interval_size_ratio how much greater should the
#' linewidth of the smaller (default 68%) credible interval
#' be compared to that of the larger (default 95%) credible
#' interval? Default 2.5
#' @param ... other parameters passed to stat_eye()
#' 
#' @return list of stats that can be added to a ggplot object
#' @export
stat_default_eye <- function(size = 15,
                             interval_size_ratio = 2.5,
                             ...){

    return (
        list(
            ggdist::stat_eye(
                        .width = plotting_params("interval_width")[1],
                        interval_alpha = 1,
                        size = interval_size_ratio * size,
                        point_alpha = 1,
                        slab_alpha = plotting_params("density_alpha"),
                        ...),
            ggdist::stat_eye(
                        .width = plotting_params("interval_width")[2],
                        slab_alpha = 0,
                        point_alpha = 0,
                        interval_alpha = 1,
                        size = size,
                        ...)
        )
    )
}


#' titer_pointinterval
#'
#' default stat_pointinterval for estimated
#' virus titers
#'
#' @param size overall size of the geom (default 30)
#' @param interval_size_ratio how much greater should the
#' linewidth of the smaller (default 68%) credible interval
#' be compared to that of the larger (default 95%) credible
#' interval? Default 2.5
#' @param point_alpha transparency of the point estimate
#' point (default 0.75)
#' @param interval_alpha transparency of the interval lines
#' (default 0.75)
#' @param stroke_mult linewidth for the point estimate point
#' border relative to the size (default 0.1, so a width of
#' 2.8 with the default size of 28).
#' @param .width width of the interval default 0.95
#' @param ... other keyword arguments passed to stat_pointinterval()
#' @return a stat_pointinterval object that
#' can be added to a ggplot object
#' @export
titer_pointinterval <- function(size = 28,
                                point_alpha = 0.75,
                                interval_alpha = 0.75,
                                stroke_mult = 0.1,
                                .width = 0.95,
                                ...){
    return (
        ggdist::stat_pointinterval(
                    .width = .width,
                    interval_alpha = interval_alpha,
                    size = size,
                    stroke = stroke_mult * size,
                    point_alpha = point_alpha,
                    slab_alpha = 0,   # workaround for
                    show_slab = TRUE, # error with
                    ...)              # show_slab = FALSE
    )
}


