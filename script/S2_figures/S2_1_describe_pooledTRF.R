rm(list = ls())

library(here)
library(ggplot2)
library(dplyr)
library(broom)
library(pbapply)
library(R6)
library(CTRF) # import my package

# This is simplified from the one in CTRF; whose name ends in 0 to differentiate
data_finalize0 <- function(df) { # create final variables
    df <- df |>
    CTRF::fff33(colname = "temperature") |> # FFF terms
    dplyr::mutate( # normalize/ standardize variables
        prcp = precipitation, # exception
        rhum = CTRF::stdize(relative_humidity), 
        wsp = CTRF::stdize(wind_speed, median=TRUE),
        skc = CTRF::stdize(skycover, median=TRUE),
        ntd = CTRF::nmlize(n_day)
    )
    df[df$prcp > 1, "prcp"] <- 1 # winsorize precipitation
    df$prcp <- df$prcp*10
    return(df)
}

regions <- c("Coast", "NC", "SC", "South")
df_info <- as.data.frame(matrix(NA, nrow=4, ncol=3))
colnames(df_info) <- c("Region", "pq")
df_info$Region <- regions
df_info$pq <- list(c(3, 2), c(3, 2), c(2, 3), c(3, 2))

# Description ---------------------------------------------------

# Loop: plot elec and temp series
for (region in regions) {
    savedir <- here(paste0("output_", region))
    dir.create(savedir, showWarnings = FALSE)
    dt_run <- readRDS(here("script", "S1_models", paste0(region, "_run.RDS")))
    if (region=="Coast") { # Remove Ike influenced period for Coast only
        dt_run <- dt_run[dt_run$Ike==0, ]
    }
    cat("Drawing pictures for ", region, " Region:\n")
    ggplot(data=dt_run, aes(x=Date))+
        geom_line(aes(y=load), color="orangered", linewidth=0.25, show.legend=FALSE)+
        geom_line(aes(y=load_ma), color="grey", linewidth=0.5, show.legend=FALSE)+
        scale_y_continuous(name="Megawatt")+
        labs(x="")
    ggsave(file=file.path(savedir , "plot_series_load.png"), width=8.5, height=0.618*6.5)

    if (region=="NC") { # draw temp series for NC only
        ggplot(dt_run, aes(x=Date))+
            geom_line(aes(y=temperature), color="turquoise3", linewidth=0.25, show.legend=FALSE)+
            scale_y_continuous(name="degree Celsius")+
            labs(x="")
        ggsave(file=file.path(savedir , "plot_series_temp.png"), width=8.5, height=0.618*6.5)    
    }
}

# Distribution. See Appendix scripts.

# All regions ------------------------------------------------
DO_THIS_POOL <- function(region, y_trim=FALSE) {
    #'
    #' This is an fixed procedure for the specific use in this script.
    #' Dataset and parameters are fixed. 
    #' Models run inside are for visual illustration only.
    #' If you want to use different parameters, must change the procedure hereafter.
    #'
    if (region=="Coast") { # Remove Ike influenced period for Coast only
        dt_run <- dt_run[dt_run$Ike==0, ]
        # Add moving average and demeaned load after removing Ike
        n <- sum(dt_run$Year==2002) + 1
        dt_run$load_ma <- data.table::frollmean(dt_run$load, n, align="right", na.rm=TRUE)
        avg2002 <- mean(dt_run$load[1:(n-1)])
        dt_run$load_ma[1:(n-1)] <- avg2002
        dt_run$load_dt <- dt_run$load - dt_run$load_ma
    }
    # I) Scatter plot: elec and temp
    plot_xy_hr <- ggplot()+ # Plot all points, even if there are frowning outliers
        geom_point(
            data = dt_run, 
            mapping = aes(x=temperature, y=load_dt, color=factor(Hour)), 
            alpha=0.3
        )+
        labs(x="Temperature C", y="Load (MW)", color="hour")+
        theme_classic()
    # Remove extremes, four-sided or two-sided
    tl <- floor(quantile(dt_run$temperature, prob=0.001, na.rm=TRUE))
    th <- ceiling(quantile(dt_run$temperature, prob=0.999, na.rm=TRUE))
    if (y_trim) {
        yl <- floor(quantile(dt_run$load_dt, prob=0.001, na.rm=TRUE))
        yh <- ceiling(quantile(dt_run$temperature, prob=0.999, na.rm=TRUE))
    } else {
        yl <- min(dt_run$load_dt); yh <- max(dt_run$load_dt)
    }
    dt_pool <- dt_run |>
        dplyr::filter(temperature>=tl & temperature<=th & load_dt>=yl & load_dt<=yh) |>
        data_finalize0() |>
        dplyr::mutate(logy = log(load))
    
    # II) Let Bootstrap TRF class run
    # Get the needed terms
    pq <- df_info$pq[[ match(region, df_info$Region) ]]
    ## predetermined degrees (from preliminary results)
    p <- pq[1]; q <- pq[2]
    nom1 <- as.vector(outer(c("ts"), 1:p, paste0)); nom1[1] <- "ts"
    nom2 <- as.vector(outer(c("sin", "cos"), 1:q, paste0))
    x_terms <- c(nom1, nom2, "ntd", "prcp", "rhum", "wsp", "skc")
    # (1) --- Y is load_dt
    dt_b <- dt_pool |> 
        dplyr::select(tidyselect::all_of(c("load_dt", "logy", x_terms))) |>
        na.omit()
    trf <- TRF_Model$new(dt_b, "load_dt", x_terms)
    trf$create_sim_matrix(c(tl, th), c(tl, th), c(3, 2))
    ## First argument is the range of plotting curve and must be a subset of
    ## the second argument which must equal or contain the range of computing 
    ## Fourier transform
    system.time({
        trf$run_bootstrap(1000)
    })
    trf$aggregate_results()
    trf_plot <- trf$plot(base_plot=plot_xy_hr)
    ggsave(file=file.path(savedir , "trf_plot_loaddt.bt.png"), bg="white", width=10, height=0.618*10)
    cat(region, " -- Scatterplot and TRF0 are done\n")
    # (2) --- Y is logy. Regress out weather covariates for regional comparison
    trf <- TRF_Model$new(dt_b, "logy", x_terms)
    trf$create_sim_matrix(c(tl, th), c(tl, th), pq)
    system.time({
        trf$run_bootstrap(1000)
    })
    trf$aggregate_results()
    df_trf_plot <- trf$aggregated_data
    df_trf_plot <- df_trf_plot[, c("temperature", "CI_lb", "fv_median", "CI_ub")]
    write.csv(df_trf_plot, file=file.path(savedir , paste0(region, "_trf_bt2.csv")), row.names=FALSE)
    cat(region, " -- TRF+ BT values have been saved.\n")
    invisible()
}

region <- regions[1] # Coast
DO_THIS_POOL(region, TRUE)
region <- regions[2] # NC
DO_THIS_POOL(region, FALSE)
region <- regions[3] # SC
DO_THIS_POOL(region, TRUE)
region <- regions[4] # South
DO_THIS_POOL(region, FALSE)
# End of script.
