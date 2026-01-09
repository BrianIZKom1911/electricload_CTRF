# This script generates scatterplot and TRF each hour (24) for four
# regions, which can be compressed into (4) GIFs.  
rm(list = ls())

library(here)
library(dplyr)
library(ggplot2)
library(R6)
library(pbapply)
library(CTRF) # Mine

regions <- c("Coast", "NC", "SC", "South")
# Build an table to retrieve region-specific information in the Loop
df_info <- as.data.frame(matrix(NA, nrow=4, ncol=3))
colnames(df_info) <- c("Region", "pq", "xlim")
df_info$Region <- regions
df_info$pq <- list(c(3, 2), c(3, 2), c(2, 3), c(2, 3))
df_info$xlim <- list(c(-10, 42), c(-10, 42), c(-10, 42), c(-8, 40))

# Define my color scheme
df_colors <- data.frame(
    Hour=seq(0, 23),
    Color=c(
        "royalblue1", "dodgerblue", "deepskyblue", "skyblue1", 
        "lightblue1", "lightcyan", "mediumturquoise", "mediumaquamarine", 
        "aquamarine", "chartreuse", "greenyellow", "yellow", 
        "gold", "goldenrod1", "orange", "sienna1", 
        "tomato", "lightcoral", "palevioletred1", "orchid2", 
        "purple", "slateblue1", "blue", "mediumblue"
    )
)

# Part I ------------------------------------------------------
# Draw separate TRF and scatterplots with fixed axes and an time tag,
# used to generate GIF

# Looped for all four regions
for (region in regions) {
    cat("Working on ", region, "Region:\n")
    # Save to each region's directory
    savedir <- here("extraoutput", region)
    dir.create(savedir, showWarnings=FALSE)
    
    # Import data ready to run
    dt_run <- readRDS(here("script", "S1_models", paste0(region, "_run.RDS")))
    dt_run$logy <- log(dt_run$load)
    # Remove Ike Hurricane influenced period for Coast only
    if (region=="Coast") {
        dt_run <- dt_run[dt_run$Ike==0, ]
    }
    
    # Prepare for the Loop
    # Get the needed terms
    pq <- df_info$pq[[ match(region, df_info$Region) ]]
    ## predetermined degrees (from preliminary results)
    p <- pq[1]; q <- pq[2]
    nom1 <- as.vector(outer(c("ts"), 1:p, paste0)); nom1[1] <- "ts"
    nom2 <- as.vector(outer(c("sin", "cos"), 1:q, paste0))
    x_terms <- c(nom1, nom2)
    # Set the four limits for the plot
    ya <- min(dt_run$logy); yb <- max(dt_run$logy)
    xlim <- df_info$xlim[[ match(region, df_info$Region) ]] # x-axis limits
    ylim <- c((round(ya, 1)-0.1), round(yb, 1)) # y-axis limits acc.t log load
    ## round ya and yb to 1 digit after decimal, then ya.rounded-0.1
    # Trim temp extremes for the other regions
    iftrim <- ifelse(region=="NC", FALSE, TRUE)
    
    hrs <- 0:23
    for (hr in hrs) {
        dt_hr <- data_finalize(dt_run, hr, iftrim)
        # Find the temp bounds used in simulation
        a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
        x_range <- c(ceiling(a), floor(b)) # make the range an subset of (a, b)
        # Keep only needed variables
        dt_b <- dt_hr |>
            dplyr::select(tidyselect::all_of(c("logy", x_terms))) |>
            na.omit()
        # Run TRF bootstrap as the R6 class in my package
        trf <- TRF_Model$new(dt_b, "logy", x_terms)
        trf$create_sim_matrix(x_range, c(a, b), pq)
        cat("Hour ", hr, "Progress:\n")
        system.time({
            trf$run_bootstrap(1000)
        })
        # Draw plot
        trf$aggregate_results()
        color_hr <- df_colors$Color[df_colors$Hour==hr]
        plot_hr <- draw_tempplot(dt_run, hr, "logy", "Log load (MW)", color_hr, xlim, ylim)
        trf_plot <- trf$plot(base_plot=plot_hr)
        trf_plot +
            annotate(
                "label", x=xlim[2]-7, y=ylim[1]+0.2, 
                label=paste0(hr, ":00"), fill="white", color="black"
            )
        ggsave(
            file=file.path(savedir, paste0("trf_plot_bt_h", hr, ".png")), 
            bg="white", width=8, height=0.618*8
        )
    }
}
#End of script.
