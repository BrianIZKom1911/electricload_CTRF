# Install remotes and devtools if needed
pkg_need <- c("remotes", "devtools")
for (pkg in pkg_need) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        install.packages(pkg)
    }
}

# Install my package using remotes
remotes::install_github("BrianIZKom1911/CTRF@v0.1.0") 
# ... or devtools if you really want
#devtools::install_github("BrianIZKom1911/CTRF@v0.1.0")

# Then simply call 
library(CTRF)
