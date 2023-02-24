#!/usr/bin/env Rscript

#####################################
## name: uninstall.R
## author: Dylan H. Morris <dylanhmorris.com>
##
## uninstalls project package and
## destroys virtual environment,
## if they are installed / exist
####################################

uninstall_if_present <- function(package_name){
    if (package_name %in% installed.packages())
        remove.packages(pkgs = package_name)
    else
        cat(sprintf("Package %s not installed\n", package_name))
}

to_uninstall <- c(
    "varstab")

cat("Uninstalling project package(s)...\n")
for (package in to_uninstall)
    uninstall_if_present(package)
