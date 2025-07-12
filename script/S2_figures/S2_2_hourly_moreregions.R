# Draw for other regions
# Only CTRF graphs from 2 to 5
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
           #ntm = nmlize(n_month),
           ntd = nmlize(n_day)) -> df
  df[df$prcp > 1, "prcp"] <- 1 # winsorize precipitation
  df$prcp <- df$prcp*10
  return(df)
}

md <- "D:\\OneDrive - University of Missouri\\transfer_desktop\\MU\\2025S_Mod_chapter2"
regions <- c("NC", "SC", "Coast", "South")
region <- regions[3]
#region <- regions[4]
savedir <- file.path(md, paste0("output_", region), "hourly")
dt_run <- readRDS(file.path(md, "S1_modelanalysis", paste0(region, "_run.RDS")))
dt_run$logy = log(dt_run$load)
if (region=="Coast"){
  dt_run <- dt_run[dt_run$Ike==0, ]
}

# Source the OOP script
source(file.path(md, "S2_makefigures", "1_boots_R6class.R"))

# Run CTRF in Subsamples -----------------------------------------------
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
#xb_terms <- c("ts", "ts2", "sin1", "cos1", "sin2", "cos2", "sin3", "cos3")
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
hrs <- c(0, 3, 6, 8, 9, 10, 11, 12, 14, 16, 17, 18, 20, 22, 23) # may be changed as you wish
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
    save_dir=savedir, y_lim=c(-0.1, 0.1)
  )
  drawsave_CTRF(
    df_ctrf_4, "WSP", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("5_h", hr, "_summer"), 
    save_dir=savedir, y_lim=c(-0.1, 0.1)
  )
  drawsave_CTRF(
    df_ctrf_5, "SKC", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("6_h", hr, "_summer"), 
    save_dir=savedir, y_lim=c(-0.1, 0.1)
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
    save_dir=savedir, y_lim=c(-0.1, 0.1)
  )
  drawsave_CTRF(
    df_ctrf_5, "SKC", v_breaks=seq(-10, 40, 5), 
    suffix=paste0("6_h", hr, "_winter"), 
    save_dir=savedir, y_lim=c(-0.12, 0.12)
  )
}
# End of script.