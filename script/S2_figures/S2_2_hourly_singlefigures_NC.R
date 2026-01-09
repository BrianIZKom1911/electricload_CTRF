# North Central, as shown in the main text
# Draw eleven hours in the TRF (1) and base TRF (1) figures:
# 0, 3, 6, 8, 10, 12, 14, 16, 18, 20, 22
rm(list = ls())

library(here)
library(ggplot2)
library(dplyr)
library(broom)
library(R6)
library(pbapply)
library(CTRF) # import my package

regions <- c("Coast", "NC", "SC", "South")
region <- regions[2] # take NC as the major instance
savedir <- here(paste0("output_", region), "hourly")
dir.create(savedir, showWarnings = FALSE)
dt_run <- readRDS(here("script", "S1_models", paste0(region, "_run.RDS")))
dt_run$logy <- log(dt_run$load)

# Run TRF ---------------------------------------------------
ya <- min(dt_run$logy); yb <- max(dt_run$logy)
xlim <- c(-10, 42) # x-axis limits
ylim <- c((round(ya, 1)-0.1), round(yb, 1)) # y-axis limits acc.t log load
## round ya and yb to 1 digit after decimal, then ya.rounded-0.1
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
#############################
# One Figure for all curves #
# To make separate figures and further make a GIF, see Appendix scripts
y_var <- "logy"
x1 <- c("ts", "ts2", "ts3")
x2 <- c("sin1", "cos1", "sin2", "cos2", "sin3", "cos3")
hrs <- c(0, 3, 6, 8, 10, 12, 14, 16, 18, 20, 22)
dfs <- vector(mode="list", length=length(hrs))
# Use this loop to store all the dfs:
# Allow different choices of (p, q)
# 6 | (1, 2)
# 3 | (2, 2)
# 16, 18, 20 | (2, 3)
# 0, 8, 10, 22 | (3, 2)
# 12, 14 | (3, 3)
for (i in 1:length(hrs)) {
    hr <- hrs[i]
    dt_hr <- data_finalize(dt_run, hour=hr)
    a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
    x_range <- c(
        round(quantile(dt_hr$temperature, probs=0.005)),
        floor(b)
    ) # set the range of ploting the curve
    # Select variables
    if (hr==6) {
        p <- 1; q <- 2
    } else if (hr==3) {
        p <- 2; q <- 2
    } else if (hr %in% c(16, 18, 20)) {
        p <- 2; q <- 3
    } else if (hr %in% c(0, 8, 10, 22)) {
        p <- 3; q <- 2
    } else if (hr==12 | hr==14) {
        p <- 3; q <- 3
    }
    order <- c(p, q)
    x_terms <- c(x1[1:p], x2[1:(2*q)], "ntd", "workday")
    dt_b <- dt_hr |>
        dplyr::select(tidyselect::all_of(c("logy", x_terms))) |>
        na.omit()
    # Run TRF Bootstrap
    trf <- TRF_Model$new(dt_b, "logy", x_terms)
    trf$create_sim_matrix(x_range, c(a, b), order)
    cat("Hour ", hr, "Progress:\n")
    system.time({
        trf$run_bootstrap(1000)
    })
    # Draw plot
    trf$aggregate_results()
    df <- trf$aggregated_data
    dfs[[i]] <- df
}
# rbind all dfs
df_all <- do.call(rbind, lapply(seq_along(dfs), function(i) {
    df <- dfs[[i]]; hr <- hrs[i]
    df$Hour <- hr # add identifier column
    df
})
)
df_all <- merge(df_all, df_colors, by="Hour")
# Plot all in one 
ggplot(df_all, aes(x = temperature, y = fv_median))+
    geom_line(aes(color = factor(Hour)))+ # plot the median line
    geom_ribbon( # add confidence bands
        aes(ymin = CI_lb, ymax = CI_ub, fill = factor(Hour)), 
        alpha = 0.2, color = NA
    )+
    labs(x = "Temperature (C)", y = "Load (log MW)")+
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    scale_color_manual(values = unique(df_all$Color), name="Hour")+
    # use existing color # legend title
    scale_fill_manual(values = unique(df_all$Color), name="Hour")+ 
    theme_classic()
ggsave(file.path(savedir, "trf_plot_bt_all11.png"), bg="white", width=8, height=0.618*8)

# Run Hourly CTRF ---------------------------------------------------
# Get single Figure for all curves at selected hours.
# Similar as above, it adds cTRFs for comparison.
xb_terms <- c("ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2")
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
group_orders <- list(c(3, 2), c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1))

hrs <- c(0, 3, 6, 8, 10, 12, 14, 16, 18, 20, 22)
dfs_0 <- vector(mode="list", length=length(hrs))
dfs_01 <- vector(mode="list", length=length(hrs))
dt_sub <- dt_run[(dt_run$workday==1), ] # only use workday data!
# Use this loop to store the base dfs and summed dfs: 
for (i in 1:length(hrs)) {
    hr <- hrs[i]
    dt_hr <- data_finalize(dt_sub, hour=hr)
    a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
    x_range <- c( # set the range of plotting the curve
        round(quantile(dt_hr$temperature, probs=0.005)),
        floor(b)
    )
    dt_b <- dt_hr |> # Generate and keep needed variables
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
        x_vars = c(xb_terms, xc_terms), # month dummies not controlled.
        ## Controlling month dummies is unnecessary. See the difference in subsamples
        var_groups = group_vars,
        orders = group_orders
    )
    ctrf$create_sim_matrices(x_range, c(a, b))
    cat("Hour ", hr, "Progress:\n")
    system.time({
        ctrf$run_bootstrap(1000)
    })
    ctrf$aggregate_results()
    ctrf$plot(savedir, hrsuffix=paste0("_h", hr)) # plot CTRF graphs
    
    ## Base TRF plus time
    df_ctrf_0 <- ctrf$aggregated_data[[1]] # base TRF coefs
    df_ctrf_1 <- ctrf$aggregated_data[[2]] # time cTRF coefs
    df_ctrf_01 <- df_ctrf_0
    df_ctrf_01[, c(2:1001)] <- df_ctrf_0[, c(2:1001)] + df_ctrf_1[, c(2:1001)]
    df_ctrf_01$CI_lb <- apply(df_ctrf_01[, -1], 1, quantile, probs=0.025, na.rm=TRUE)
    df_ctrf_01$fv_median <- apply(df_ctrf_01[, -1], 1, median, na.rm=TRUE)
    df_ctrf_01$CI_ub <- apply(df_ctrf_01[, -1], 1, quantile, probs=0.975, na.rm=TRUE)
    
    dfs_0[[i]] <- df_ctrf_0 # base TRF BT values
    dfs_01[[i]] <- df_ctrf_01 # base+time cTRF BT values
}
# Plot all in one twice: First for base TRF, second for base+time cTRF
# 1) rbind all dfs_0
df_all <- do.call(rbind, lapply(seq_along(dfs_0), function(i) {
    df <- dfs_0[[i]]; hr <- hrs[i]
    df$Hour <- hr # add identifier column
    df
})
)
write.csv(df_all[, c(1002:1005)], file=file.path(savedir, "ctrf_00_bt.csv"), row.names=FALSE)
## save medians of df_bt
df_all <- merge(df_all, df_colors, by="Hour")
# Draw picture #1
ggplot(df_all, aes(x = temperature, y = fv_median))+
    geom_line(aes(color = factor(Hour)))+
    geom_ribbon(
        aes(ymin = CI_lb, ymax = CI_ub, fill = factor(Hour)), 
        alpha = 0.2, color = NA
    )+
    annotate("text", x=15, y=10.1, label="Base TRF (2002)", color="black")+
    labs(x = "Temperature (C)", y = "Load (log MW)")+
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    scale_color_manual(values = unique(df_all$Color), name="Hour")+
    scale_fill_manual(values = unique(df_all$Color), name="Hour")+ 
    theme_classic()
ggsave(file.path(savedir, "ctrf_00_all11.png"), bg="white", width=7, height=0.618*7)

# 2) rbind all dfs_01
df_all <- do.call(rbind, lapply(seq_along(dfs_01), function(i) {
    df <- dfs_01[[i]]; hr <- hrs[i]
    df$Hour <- hr # add identifier column
    df
})
)
write.csv(df_all[, c(1002:1005)], file=file.path(savedir, "ctrf_01_bt.csv"), row.names=FALSE)
## save medians of df_bt
df_all <- merge(df_all, df_colors, by="Hour")
# Draw picture #2
ggplot(df_all, aes(x = temperature, y = fv_median))+
    geom_line(aes(color = factor(Hour)))+
    geom_ribbon(
        aes(ymin = CI_lb, ymax = CI_ub, fill = factor(Hour)), 
        alpha = 0.2, color = NA
    )+
    annotate("text", x=15, y=10.2, label="Base TRF + Time (2024)", color="black")+
    labs(x = "Temperature (C)", y = "Load (log MW)")+
    scale_x_continuous(limits = xlim)+
    scale_y_continuous(limits = ylim)+
    scale_color_manual(values = unique(df_all$Color), name="Hour")+
    scale_fill_manual(values = unique(df_all$Color), name="Hour")+ 
    theme_classic()
ggsave(file.path(savedir, "ctrf_01_all11.png"), bg="white", width=7, height=0.618*7)
# End of script.
