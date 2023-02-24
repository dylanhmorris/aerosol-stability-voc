#!/usr/bin/env Rscript

#####################################
## name: install.R
## author: Dylan H. Morris <dylanhmorris.com>
##
## installs project package and its
## dependencies
####################################

install_if_absent <- function(package_name){
    if (!(package_name %in% installed.packages())){
        cat(sprintf("Attempting to install package %s\n",
                    package_name))
        install.packages(pkgs = package_name,
                         repos = "http://cloud.r-project.org")
    } else {
        cat(sprintf("Package %s already installed\n", package_name))
    }
}

install_local_force <- function(package_name){
    cat(sprintf("Installing local package %s...\n",
                package_name))
    remotes::install_local(package_name, force=TRUE)

}


args <- commandArgs(trailingOnly=TRUE)

## install CRAN packages
install_if_absent("remotes")

cat("Attempting to install project R package and any missing dependencies...\n")
install_local_force("varstab")
