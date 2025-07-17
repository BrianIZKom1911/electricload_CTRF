# Required packages
pkg_needed <- c("tidyverse", "AER", "suncalc", "here", "data.table", "broom", "pbapply", "R6")

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