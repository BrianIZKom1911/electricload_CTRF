# Selected hours: 9:00--23:00
# Draw eight hours: 9, 11, 13, 15, 17, 19, 21, 23
rm(list = ls())

library(here)
library(ggplot2)
library(dplyr)
library(broom)

fff33 <- function(data, colname, a=NULL, b=NULL){
  a <- ifelse(is.null(a), min(data[[colname]]), a)
  b <- ifelse(is.null(b), max(data[[colname]]), b)
  data %>%
    mutate(ts = (.data[[colname]]-a)/(b-a), 
           ts2 = ts^2, ts3 = ts^3, 
           sin1 = sin(2*pi*ts), cos1 = cos(2*pi*ts), 
           sin2 = sin(4*pi*ts), cos2 = cos(4*pi*ts), 
           sin3 = sin(6*pi*ts), cos3 = cos(6*pi*ts)) -> data
  return(data)
}
stdize <- function(x, mean=NULL, median=FALSE){
  mean <- ifelse(is.null(mean), mean(x, na.rm=TRUE), mean) 
  ## default mean function or provided number
  m <- ifelse(median, median(x, na.rm=TRUE), mean)
  ## output of previous step (default) or median
  stdev <- sd(x, na.rm=TRUE)
  xstd <- (x - m)/stdev
  return(xstd)
}
nmlize <- function(x, a=NULL, b=NULL){
  a <- ifelse(is.null(a), min(x), a)
  b <- ifelse(is.null(b), max(x), b)
  xnml <- (x - a)/(b - a)
  return(xnml)
}

data_finalize <- function(df, hour) { # create final variables
  df %>%
    filter(Hour==hour) %>%
    fff33(., "temperature") %>% # FFF terms
    mutate(prcp = precipitation,
           rhum = stdize(relative_humidity), # normalize/ standardize variables
           wsp = stdize(wind_speed, median=TRUE),
           skc = stdize(skycover, median=TRUE),
           ntm = nmlize(n_month),
           ntd = nmlize(n_day)) -> df
  df[df$prcp > 1, "prcp"] <- 1 # winsorize precipitation
  df$prcp <- df$prcp*10
  return(df)
}

md <- "D:\\OneDrive - University of Missouri\\transfer_desktop\\MU\\2025S_Mod_chapter2"
regions <- c("NC", "SC", "Coast", "South")
region <- regions[1]
savedir <- file.path(md, paste0("output_", region), "hourly")
dt_run <- readRDS(file.path(md, "S1_modelanalysis", paste0(region, "_run.RDS")))
dt_run$logy <- log(dt_run$load)

# Source the OOP script
source(file.path(md, "S2_makefigures", "1_boots_R6class.R"))
# Run TRF ---------------------------------------------------
ya <- min(dt_run$logy); yb <- max(dt_run$logy)
xlim <- c(-10, 42)
ylim <- c((round(ya, 1)-0.1), round(yb, 1)) # check the range of x and y
df_colors <- data.frame(
  Hour=seq(0, 23),
  Color=c("royalblue1", "dodgerblue", "deepskyblue", "skyblue1", 
          "lightblue1", "lightcyan", "mediumturquoise", "mediumaquamarine", 
          "aquamarine", "chartreuse", "greenyellow", "yellow", 
          "gold", "goldenrod1", "orange", "sienna1", 
          "tomato", "lightcoral", "palevioletred1", "orchid2", 
          "purple", "slateblue1", "blue", "mediumblue")
)
###############################
# 1 Separate figure each hour #
x_terms <- c("ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2")
hrs <- c(8, 9, 10, 11, 14, 17, 20, 23)
for (hr in hrs) {
  dt_hr <- data_finalize(dt_run, hour=hr)
  a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
  x_range <- c(ceiling(a), floor(b)) # find the limits
  # Select variables
  dt_hr %>%
    select(all_of(c("logy", x_terms))) %>%  # keep only needed variables
    na.omit() -> dt_b
  # Run TRF boots
  trf <- TRF_Model$new(dt_b, "logy", x_terms)
  trf$create_sim_matrix(x_range, c(a, b), c(3, 2))
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
    annotate("label", x=15, y=10.4, label=paste0(hr, ":00"), fill="white", color="black")
  ggsave(file=file.path(savedir, paste0("trf_plot_bt_h", hr, ".png")), bg="white", width=8, height=0.618*8)
}

########################
# 2 One Figure for all #
y_var <- "logy"
x1 <- c("ts", "ts2", "ts3")
x2 <- c("sin1", "cos1", "sin2", "cos2", "sin3", "cos3")
hrs <- c(0, 3, 6, 8, 10, 12, 14, 16, 18, 20, 22)
dfs <- vector(mode="list", length=length(hrs))
# Use this loop to store all the dfs:
# Allow different choices of (p, q)
#6:00 | (1, 2)
#3:00 | (2, 2)
#16, 18, 20 | (2, 3)
#0, 8, 10, 22 | (3, 2)
#12, 14 | (3, 3)

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
  dt_hr %>%
    select(all_of(c("logy", x_terms))) %>% # keep only needed variables
    na.omit() -> dt_b
  # Run TRF boots
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
  geom_ribbon(aes(ymin = CI_lb, ymax = CI_ub, fill = factor(Hour)), 
              alpha = 0.2, color = NA)+ # add confidence bands
  labs(x = "Temperature (C)", y = "Load (log MW)")+
  scale_x_continuous(limits = xlim)+
  scale_y_continuous(limits = ylim)+
  scale_color_manual(values = unique(df_all$Color), name="Hour")+
  # use existing color # legend title
  scale_fill_manual(values = unique(df_all$Color), name="Hour")+ 
  theme_classic()
ggsave(file.path(savedir, "trf_plot_bt_all11.png"), bg="white", width=8, height=0.618*8)

# Run CTRF ---------------------------------------------------
# Controlling month dummies is unnecessary. See the difference in subsamples.
#months <- c("m01", "m02", "m03", "m04", "m05", "m06", "m07", "m08", "m09", "m10", "m11")

# (04/20): don't include workday. Split the sample into two.
xb_terms <- c("ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2")
xc_terms <- c(
  "ntd", "s_t", "sin1_t", "cos1_t",
  "prcp", "s_p", "sin1_p", "cos1_p",
  "rhum", "s_h", "sin1_h", "cos1_h",
  "skc", "s_c", "sin1_c", "cos1_c",
  "wsp", "s_w", "sin1_w", "cos1_w"
)
group_vars <- list(
  vars0 = xb_terms,
  vars1 = xc_terms[1:4],
  vars2 = xc_terms[5:8],
  vars3 = xc_terms[9:12],
  vars4 = xc_terms[13:16],
  vars5 = xc_terms[17:20]
)
group_orders <- list(
  c(3, 2), c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1)
)
hrs <- c(0, 3, 6, 8, 10, 12, 14, 16, 18, 20, 22)
dfs_0 <- vector(mode="list", length=length(hrs))
dfs_01 <- vector(mode="list", length=length(hrs))
dt_sub <- dt_run[(dt_run$workday==1), ] # only use workday data!
# Use this loop to store the base dfs and summed dfs: 
for (i in 1:length(hrs)) {
  hr <- hrs[i]
  dt_hr <- data_finalize(dt_sub, hour=hr)
  a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
  x_range <- c(
    round(quantile(dt_hr$temperature, probs=0.005)),
    floor(b)
  ) # set the range of ploting the curve
  dt_hr %>%
    mutate(s_t = ts*ntd, sin1_t = sin1*ntd, cos1_t = cos1*ntd,
           s_p = ts*prcp, sin1_p = sin1*prcp, cos1_p = cos1*prcp,
           s_h = ts*rhum, sin1_h = sin1*rhum, cos1_h = cos1*rhum,
           s_w = ts*wsp, sin1_w = sin1*wsp, cos1_w = cos1*wsp,
           s_c = ts*skc, sin1_c = sin1*skc, cos1_c = cos1*skc) %>%
    select(all_of(c("logy", xb_terms, xc_terms))) %>%  # keep only needed variables
    na.omit() -> dt_b
  # Run CTRF boots
  ctrf <- CTRF_Model$new(
    data = dt_b, 
    y_var = "logy", 
    x_vars = c(xb_terms, xc_terms), # month dummies not controlled
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
  
  dfs_0[[i]] <- df_ctrf_0
  dfs_01[[i]] <- df_ctrf_01
}
# Plot all in one twice 
# rbind all dfs
df_all <- do.call(rbind, lapply(seq_along(dfs_0), function(i) {
  df <- dfs_0[[i]]; hr <- hrs[i]
  df$Hour <- hr # add identifier column
  df
})
)
write.csv(df_all[, c(1002:1005)], file=file.path(savedir, "ctrf_00_bt.csv"), row.names=FALSE)
## save medians of df_bt
df_all <- merge(df_all, df_colors, by="Hour")
ggplot(df_all, aes(x = temperature, y = fv_median))+
  geom_line(aes(color = factor(Hour)))+
  geom_ribbon(aes(ymin = CI_lb, ymax = CI_ub, fill = factor(Hour)), 
              alpha = 0.2, color = NA)+
  annotate("text", x=15, y=10.1, label="Base TRF (2002)", color="black")+
  labs(x = "Temperature (C)", y = "Load (log MW)")+
  scale_x_continuous(limits = xlim)+
  scale_y_continuous(limits = ylim)+
  scale_color_manual(values = unique(df_all$Color), name="Hour")+
  scale_fill_manual(values = unique(df_all$Color), name="Hour")+ 
  theme_classic()
ggsave(file.path(savedir, "ctrf_00_all11.png"), bg="white", width=7, height=0.618*7)

df_all <- do.call(rbind, lapply(seq_along(dfs_01), function(i) {
  df <- dfs_01[[i]]; hr <- hrs[i]
  df$Hour <- hr # add identifier column
  df
})
)
write.csv(df_all[, c(1002:1005)], file=file.path(savedir, "ctrf_01_bt.csv"), row.names=FALSE)
## save medians of df_bt
df_all <- merge(df_all, df_colors, by="Hour")
ggplot(df_all, aes(x = temperature, y = fv_median))+
  geom_line(aes(color = factor(Hour)))+
  geom_ribbon(aes(ymin = CI_lb, ymax = CI_ub, fill = factor(Hour)), 
              alpha = 0.2, color = NA)+
  annotate("text", x=15, y=10.2, label="Base TRF + Time (2024)", color="black")+
  labs(x = "Temperature (C)", y = "Load (log MW)")+
  scale_x_continuous(limits = xlim)+
  scale_y_continuous(limits = ylim)+
  scale_color_manual(values = unique(df_all$Color), name="Hour")+
  scale_fill_manual(values = unique(df_all$Color), name="Hour")+ 
  theme_classic()
ggsave(file.path(savedir, "ctrf_01_all11.png"), bg="white", width=7, height=0.618*7)

# Subsamples ---------------------------------------------------
# Workday and seasons
savedir <- file.path(md, paste0("output_", region), "subsample")
if (dir.exists(savedir)) {
  savedir <- file.path(md, paste0("output_", region), "subsample")
} else {
  dir.create(savedir)
}

## Summer: June, July, August
dt_sub <- dt_run[(dt_run$workday==1) & (dt_run$Month %in% c(6,7,8)), ]
xb_terms <- c("ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2")
xc_terms <- c(
  "ntd", "s_t", "sin1_t", "cos1_t",
  "prcp", "s_p", "sin1_p", "cos1_p",
  "rhum", "s_h", "sin1_h", "cos1_h",
  "skc", "s_c", "sin1_c", "cos1_c",
  "wsp", "s_w", "sin1_w", "cos1_w"
)
group_vars <- list(
  vars0 = xb_terms,
  vars1 = xc_terms[1:4],
  vars2 = xc_terms[5:8],
  vars3 = xc_terms[9:12],
  vars4 = xc_terms[13:16],
  vars5 = xc_terms[17:20]
)
group_orders <- list(
  c(3, 2), c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1)
)
hrs <- c(0, 3, 6, 8, 9, 10, 11, 12, 14, 16, 17, 18, 20, 22, 23)
hrs <- c(0, 3, 6, 12, 16, 18, 22)
for (hr in hrs) {
  dt_hr <- data_finalize(dt_sub, hour=hr)
  a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
  x_range <- c(ceiling(a), floor(b)) # find the limits
  dt_hr %>%
    mutate(s_t = ts*ntd, sin1_t = sin1*ntd, cos1_t = cos1*ntd,
           s_p = ts*prcp, sin1_p = sin1*prcp, cos1_p = cos1*prcp,
           s_h = ts*rhum, sin1_h = sin1*rhum, cos1_h = cos1*rhum,
           s_w = ts*wsp, sin1_w = sin1*wsp, cos1_w = cos1*wsp,
           s_c = ts*skc, sin1_c = sin1*skc, cos1_c = cos1*skc) %>%
    select(all_of(c("logy", xb_terms, xc_terms))) %>%  # keep only needed variables
    na.omit() -> dt_b
  # Run CTRF boots
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
  df_ctrf_3 <- ctrf$aggregated_data[[4]] # rhum cTRF coefs
  df_ctrf_4 <- ctrf$aggregated_data[[5]] # wsp TRF coefs
  df_ctrf_5 <- ctrf$aggregated_data[[6]] # skc cTRF coefs
  drawsave_CTRF(
    df_ctrf_3, "RH", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("4_h", hr, "_summer"), 
    save_dir=savedir, y_lim=c(-0.08, 0.1)
  )
  drawsave_CTRF(
    df_ctrf_4, "WSP", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("5_h", hr, "_summer"), 
    save_dir=savedir, y_lim=c(-0.05, 0.05)
  )
  drawsave_CTRF(
    df_ctrf_5, "SKC", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("6_h", hr, "_summer"), 
    save_dir=savedir, y_lim=c(-0.03, 0.07)
  )
}

## Winter: December, January, February 
dt_sub <- dt_run[(dt_run$workday==1) & (dt_run$Month %in% c(1,2,12)), ]
for (hr in hrs) {
  dt_hr <- data_finalize(dt_sub, hour=hr)
  a <- min(dt_hr$temperature); b <- max(dt_hr$temperature)
  x_range <- c(ceiling(a), floor(b)) # find the limits
  dt_hr %>%
    mutate(s_t = ts*ntd, sin1_t = sin1*ntd, cos1_t = cos1*ntd,
           s_p = ts*prcp, sin1_p = sin1*prcp, cos1_p = cos1*prcp,
           s_h = ts*rhum, sin1_h = sin1*rhum, cos1_h = cos1*rhum,
           s_w = ts*wsp, sin1_w = sin1*wsp, cos1_w = cos1*wsp,
           s_c = ts*skc, sin1_c = sin1*skc, cos1_c = cos1*skc) %>%
    select(all_of(c("logy", xb_terms, xc_terms))) %>%  # keep only needed variables
    na.omit() -> dt_b
  # Run CTRF boots
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
  df_ctrf_3 <- ctrf$aggregated_data[[4]] # rhum cTRF coefs
  df_ctrf_4 <- ctrf$aggregated_data[[5]] # wsp TRF coefs
  df_ctrf_5 <- ctrf$aggregated_data[[6]] # skc cTRF coefs
  drawsave_CTRF(
    df_ctrf_3, "RH", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("4_h", hr, "_winter"), 
    save_dir=savedir, y_lim=c(-0.1, 0.1)
  )
  drawsave_CTRF(
    df_ctrf_4, "WSP", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("5_h", hr, "_winter"), 
    save_dir=savedir, y_lim=c(-0.075, 0.075)
  )
  drawsave_CTRF(
    df_ctrf_5, "SKC", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("6_h", hr, "_winter"), 
    save_dir=savedir, y_lim=c(-0.08, 0.12)
  )
}
# End of script.