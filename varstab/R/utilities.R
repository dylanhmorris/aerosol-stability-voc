##################################################
## filename: utilities.R
## author: Dylan H. Morris <dylanhmorris.com>
##
## description: miscellaneous utilities for
## variant stability analysis 
##################################################

#' check_data_quality
#'
#' Verify that data does not have obvious errors
#' before passing to Stan fitting functions
#'
#' @param n_obs number of observations
#' 
#' @param dilution numeric array of log10 well dilution factors
#' (0 = no dilution, -1 = 10-fold dilution, etc).
#' 
#' @param id_columns list of columns in the data containing
#' ids for each observation (all of which must be integer vectors
#' of length equal to the number of observations)
#' 
#' @return TRUE 
#' @export
check_data_quality <- function(n_obs,
                               dilution,
                               id_columns){
    
    if(length(dilution) != n_obs)
        stop(paste0("Must have same number of dilution factors ",
                    "as well status observations"))

    for(id_column_name in names(id_columns)) {
        id_col <- id_columns[[id_column_name]]
        if(length(id_col) != n_obs){
            stop(paste0("Must have same number of ",
                        id_column_name,
                        " ids ",
                        "as well status observations"))
        }

        if(any(id_col != floor(id_col))) {
            stop(paste0(id_column_name,
                        " ids must be ",
                        "positive integers"))
        }
    }

    return (TRUE)

}
                               
