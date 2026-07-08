# Draw for four regions, three seasons, all hours (24), including PRCP.
# This script generates more CTRF graphs. They can be compressed into (48) GIFs 
# or just kept for record.
rm(list = ls())

library(here)
library(ggplot2)
library(dplyr)
library(broom)
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

# Run Hourly CTRF on Subsamples -----------------------------------------------
# Workday and seasons

DO_THIS_CTRF <- function(saveto, region, season, workdayonly=FALSE) {
    #'
    #' This is an fixed procedure to do all CTRF plots for given region and season.
    #' Dataset and parameters for the Bootstrap class are fixed, e.g. df_into, group_vars _orders.
    #' If you want to use different parameters, must change the procedure hereafter.
    #' 
    dir.create(saveto, showWarnings = FALSE)
    data <- readRDS(here("script", "S1_models", paste0(region, "_run.RDS")))
    data$logy <- log(data$load)
    # Remove Ike Hurricane influenced period for Coast only
    if (region=="Coast") {
        data <- data[data$Ike==0, ]
    }
    # Get the needed terms
    pq <- df_info$pq[[ match(region, df_info$Region) ]]
    ## predetermined degrees (from preliminary results)
    p <- pq[1]; q <- pq[2]
    nom1 <- as.vector(outer(c("ts"), 1:p, paste0)); nom1[1] <- "ts"
    nom2 <- as.vector(outer(c("sin", "cos"), 1:q, paste0))
    xb_terms <- c(nom1, nom2)
    
    #### These parameters are fixed
    xc_terms <- c(
        "ntd", "s_t", "sin1_t", "cos1_t",
        "prcp", "s_p", "sin1_p", "cos1_p",
        "rhum", "s_h", "sin1_h", "cos1_h",
        "wsp", "s_w", "sin1_w", "cos1_w",
        "skc", "s_c", "sin1_c", "cos1_c"
    )
    group_vars <- list(
        vars0 = xb_terms,
        vars1 = xc_terms[1:4],
        vars2 = xc_terms[5:8],
        vars3 = xc_terms[9:12],
        vars4 = xc_terms[13:16],
        vars5 = xc_terms[17:20]
    )
    group_orders <- list(pq, c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1)) 
    ##### END
    
    # Select months
    if (season=="summer") { # Summer: June, July, August
        smons <- c(6, 7, 8)
    } else if (season=="winter") { # Winter: January, February, December
        smons <- c(1, 2, 12)
    } else { # Spring and Fall
        smons <- c(3, 4, 5, 9, 10, 11)
    }
    if (workdayonly) { # whether exclude holidays&weekends
        dt_sub <- data[(data$workday==1) & (data$Month %in% smons), ]
    } else {
        dt_sub <- data[data$Month %in% smons, ]
    }
    # Trim temp extremes for the other regions
    iftrim <- region!="NC"
    
    hrs <- 0:23
    for (hr in hrs) {
        dt_hr <- data_finalize(dt_sub, hr, iftrim)
        # Find the temp bounds used in simulation
        a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
        x_range <- c(ceiling(a), floor(b)) # make the range an subset of (a, b)
        # Generate and keep needed variables
        dt_b <- dt_hr |>
            dplyr::mutate(
                s_t = ts*ntd, sin1_t = sin1*ntd, cos1_t = cos1*ntd,
                s_p = ts*prcp, sin1_p = sin1*prcp, cos1_p = cos1*prcp,
                s_h = ts*rhum, sin1_h = sin1*rhum, cos1_h = cos1*rhum,
                s_w = ts*wsp, sin1_w = sin1*wsp, cos1_w = cos1*wsp,
                s_c = ts*skc, sin1_c = sin1*skc, cos1_c = cos1*skc
            ) |>
            dplyr::select(tidyselect::all_of(c("logy", xb_terms, xc_terms))) |>
            na.omit()
        # Run CTRF Bootstrap
        ctrf <- CTRF_Model$new(
            data = dt_b, 
            y_var = "logy", 
            x_vars = c(xb_terms, xc_terms),
            var_groups = group_vars,
            orders = group_orders
        )
        ctrf$create_sim_matrices(x_range, c(a, b))
        cat("Hour ", hr, "Progress:\n")
        system.time({
            ctrf$run_bootstrap(1000)
        })
        # Draw plot with different y_lim for different variables
        ctrf$aggregate_results()
        df_ctrf_2 <- ctrf$aggregated_data[[3]] # prcp cTRF coefs
        df_ctrf_3 <- ctrf$aggregated_data[[4]] # rhum cTRF coefs
        df_ctrf_4 <- ctrf$aggregated_data[[5]] # wsp cTRF coefs
        df_ctrf_5 <- ctrf$aggregated_data[[6]] # skc cTRF coefs
        draw_save_CTRF(
            df_ctrf_2, "PRCP", paste0(hr, ":00"), 
            suffix=paste0("2_h", hr, "_", season), save_dir=saveto, 
            x_breaks=seq(-10, 40, 5), y_lim=c(-0.1, 0.11)
        )
        draw_save_CTRF(
            df_ctrf_3, "RH", paste0(hr, ":00"), 
            suffix=paste0("3_h", hr, "_", season), save_dir=saveto, 
            x_breaks=seq(-10, 40, 5), y_lim=c(-0.1, 0.13)
        )
        draw_save_CTRF(
            df_ctrf_4, "WSP", paste0(hr, ":00"), 
            suffix=paste0("4_h", hr, "_", season), save_dir=saveto, 
            x_breaks=seq(-10, 40, 5), y_lim=c(-0.1, 0.11)
        )
        draw_save_CTRF(
            df_ctrf_5, "SKC", paste0(hr, ":00"), 
            suffix=paste0("5_h", hr, "_", season), save_dir=saveto, 
            x_breaks=seq(-10, 40, 5), y_lim=c(-0.1, 0.11)
        )
    }
    invisible()
}

# Test: NC workdays, every season (3)
region <- regions[2]
for (sn in c("summer", "winter", "SnF")) {
    savedir <- here(paste0("output_", region), sn) # NC workdays in main output
    DO_THIS_CTRF(savedir, region, sn, workdayonly = TRUE)
}

# Loop: every region (4), every season (3)
for (region in regions) {
    cat("Working on ", region, "Region:\n")
    for (sn in c("summer", "winter", "SnF")) {
        savedir <- here("extraoutput", region, sn) # more in extra output
        DO_THIS_CTRF(savedir, region, sn, workdayonly = FALSE)
    }
}
# End of script.
