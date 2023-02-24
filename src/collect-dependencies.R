packages <- renv::dependencies()

packages <- unique(packages$Package)

excludes <- c(
    "varstab",
    "usethis",
    "renv")

usethis::proj_set("varstab")

for(package in packages){
    if (!(package %in% excludes))
        usethis::use_package(package, min_version = TRUE)
}
