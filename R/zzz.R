.onLoad <- function(libname, pkgname) {
    repos <- getOption("repos")
    repos["r-multiverse"] <- "https://community.r-multiverse.org"
    options(repos = repos)
    invisible(NULL)
}
