#!/usr/bin/env Rscript

###################################
## shared styling for all plots
##
####################################

## styling for text

power10_breaks <- scales::trans_breaks("log10",
                                       function(x) 10^x)
power10_format <- scales::trans_format("log10",
                                       scales::math_format(10^.x))


#' variant_order
#' 
#' order in which to label virus variants
#'
#' @export
variant_order <- function() {

    return (
        c(
            "WA1",
            "B.1",
            "Alpha",
            "Beta",
            "Gamma",
            "Delta",
            "Omicron")
    )
}

#' plotting_params
#'
#' @param param_to_get which parameter is desired
#'
#' @return parameter value
#' @export
plotting_params <- function(param_to_get){

    param_list <- list(
        density_alpha = 0.75,
        interval_width = c(0.68, 0.95),
        interval_size_domain = c(2, 12)
    )

    if(! param_to_get %in% names(param_list) ){
        stop(sprintf("Unknown parameter %s",
                     param_to_get))
    }
    
    return (param_list[[param_to_get]])
    
}


#' scale_x_log10_mathformat
#' 
#' convenience scale function
#' for putting x breaks on
#' powers of 10 in ggplot
#'
#' @param ... arguments passed to scale_x_continuous()
#'
#' @return a scale function
#' 
#'@export
scale_x_log10_mathformat <- function(...){
    return (ggplot2::scale_x_continuous(trans = "log10",
                                        breaks = power10_breaks,
                                        label = power10_format,
                                        ...))
}

#'
#'
#' convenience scale function
#' for putting y breaks on
#' powers of 10 in ggplot
#'
#' @param ... arguments passed to scale_y_continuous()
#'
#' @return a scale function
#' 
#'@export
scale_y_log10_mathformat <- function(...){
    return (ggplot2::scale_y_continuous(trans = "log10",
                                        breaks = power10_breaks,
                                        label = power10_format,
                                        ...))
}


#' theme_project
#'
#' variant of theme_classic() ggplot theme for
#' this project
#' 
#' @param base_size base font size for the theme, passed to
#' theme_classic (default 30)
#' @param ... other parameters passed to theme_classic()
#' @return theme
#' @export
theme_project <- function(base_size = 30,
                          ...){
    x_axis_margin <- ggplot2::margin(t = base_size / 2)
    y_axis_margin <- ggplot2::margin(r = base_size / 2)

    ggplot2::theme_classic(
                 base_size = base_size,
                 ...) +
    cowplot::background_grid(major = "xy",
                             minor = "none",
                             size.major = 0.5) +
    ggplot2::theme(axis.title.x = element_text(
                       margin = x_axis_margin),
                   axis.title.y = element_text(
                       margin = y_axis_margin),
                   strip.background = element_blank(),
                   panel.border = element_blank()
                   )
}


#' pango_to_who
#'
#' converts a vector of PANGO lineage designations
#' to a vector of WHO Greek letter variant names
#'
#' @param pango_lineage PANGO lineage(s) of the variant(s)
#' @param retain_pango flag for whether to retain the
#' PANGO lineage in parentheses after to WHO name. Default FALSE.
#' @param pango_if_fail Flag for what to do if lookup fails. If TRUE
#' (default) return the PANGO lineage inputed. If FALSE,
#' return NA.
#' 
#' @return vector of WHO variant names, possibly with
#' PANGO lineage in parentheses 
#' @export
pango_to_who <- function(pango_lineage,
                         pango_if_fail = FALSE){
    
    Greek_letters <- dplyr::case_when(
                                grepl("B.1.1.7", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Alpha",
                            

                                grepl("B.1.351", pango_lineage,
                                      ignore.case = TRUE) ~
                            
                                "Beta",


                                grepl("P.1", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Gamma",
                            

                                grepl("B.1.617.2", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Delta",
                            

                                grepl("B.1.42[7,9]", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Epsilon",
                            

                                grepl("P.2", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Zeta",
                            
                            
                                grepl("B.1.525", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Eta",

                            
                                grepl("P.3", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Theta",
                            
                                
                                grepl("B.1.526", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Iota",
                            

                                grepl("B.1.617.1", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Kappa",
                            

                                grepl("C.37", pango_lineage,
                                      ignore.case = TRUE) ~
                                "Lambda",
                            

                                grepl("B.1.621", pango_lineage,
                                      ignore.case = TRUE) ~
                                    "Mu",

                                grepl("B.1.1.529", pango_lineage,
                                      ignore.case = TRUE) ~
                                    "Omicron")
    
    if(pango_if_fail) {
        Greek_letters <- ifelse(is.na(Greek_letters),
                                pango_lineage,
                                Greek_letters)
    }

    return (Greek_letters)
}



#' variant_colors
#'
#' color scheme for variants
#' of concern
#'
#' Palette generated with Colorgorical
#' http://vrl.cs.brown.edu/color
#'
variant_colors <- list(
    "B.1" = "#2f7842",
    "Alpha" = "#b51936",
    "Gamma" = "pink",
    "Beta" = "#0362a0",
    "Delta" = "#f76015",
    "Omicron" = "#b543a6",
    "WA1" = "#969696")


#' scale_fill_variant
#'
#' default fill scale for variants
#' for this project
#'
#' @param ... other parameters passed to scale_fill_manual()
#' 
#' @return fill scale
#' @export
scale_fill_variant <- function(...) {
    
    return( ggplot2::scale_fill_manual(
                         values = unlist(variant_colors),
                         ...) )
}


#' scale_color_variant
#'
#' default color scale for variants
#' for this project
#'
#' @param ... other parameters passed to scale_color_manual()
#' 
#' @return color scale
#' @export
scale_color_variant <- function(...) {
    
    return( ggplot2::scale_color_manual(
                         values = unlist(variant_colors),
                         ...) )
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
                    ...)             # show_slab = FALSE
    )
}


