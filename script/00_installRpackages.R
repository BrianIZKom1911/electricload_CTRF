# Required packages
pkg_needed <- c(
    "here",         # import here() to locate the main directory where .Rproj is
    "tidyverse",    # dplyr, ggplot2, ... common operations
    "AER",          # sandwich, lmtest, used in regression summary
    "data.table",   # import frollmean() to get moving average
    "suncalc",      # used for sunlight calculation
    "broom",        # used to print clean output
    "pbapply", "R6" # used to build the R6 class; wrapped in my package "CTRF"
)
# Check if your computer has installed these packages
install_missing <- function(pkgnames) {
    installed <- pkgnames %in% rownames(installed.packages())
    if (any(!installed)) {
        cat("Installing packages:\n")
        cat(pkgnames[!installed], sep = ", ")
        cat("\n")
        install.packages(pkgnames[!installed])
    } else {
        cat("All required packages already exist.\n")
    }
}

install_missing(pkd_needed)
