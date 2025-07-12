rm(list = ls())

library(here)
library(tidyverse)
library(AER)
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
# Create variables
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

print_output <- function(..., filepath, text=NULL){
  models <- list(...)
  sink(filepath)
  for (i in seq_along(models)) {
    cat("Model ", i, " summary:\n")
    print(lmtest::coeftest(models[[i]], vcov = sandwich::vcovHC(models[[i]], type="HC3")))
    print(broom::glance(models[[i]]))
    cat("\n")
  }
  cat(text)
  sink(file = NULL)
  invisible()
}

regions <- c("NC", "SC", "Coast", "South")
region <- regions[1]
savedir <- here(paste0("output_", region), "hourly")
if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}
dt_run <- readRDS(here("script", "S1_models", paste0(region, "_run.RDS")))

# Sum stats -----------------------------------------------------
# Try using workday as a covariate w_7.
# Sunlight is embedded in different hours. 
# Only consider 15 hours: Working hours 9:00--18:00 and off-work active hours 19:00--23:00
hrs <- seq(8, 23)

# Summary statistics for all hours selected
w_names <- c("load", "temperature", "dew_point_temperature", "wet_bulb_temperature", 
            "station_level_pressure", "precipitation", "relative_humidity", "wind_speed", "skycover")
col_names <- c("Variable", "Min", "Q1", "Median", "Mean", "Q3", "Max", "SD", "nobs", "Hour")
all_tab_stats <- list()
for (j in 1:length(hrs)) {
  hr <- hrs[j]
  tab_stats <- as.data.frame(matrix(NA, nrow=length(w_names), ncol=length(col_names)))
  colnames(tab_stats) <- col_names
  tab_stats[, 10] <- hr
  for (i in 1:length(w_names)) {
    w <- w_names[i]
    tab_stats[i, 2:7] <- summary(dt_run[dt_run$Hour==hr, w])
    tab_stats[i, 8] <- sd(dt_run[dt_run$Hour==hr, w], na.rm=TRUE)
    tab_stats[i, 9] <- sum(is.na(dt_run[dt_run$Hour==hr, w])==FALSE)
  }
  all_tab_stats[[j]] <- tab_stats
}
stacked_tab <- do.call(rbind, all_tab_stats)
stacked_tab$Variable <- w_names
write.csv(stacked_tab, file.path(savedir, "tab_stats_byhour.csv"), row.names=FALSE)

# Correlation table (03/27)
## Correlation coefficients between one and another variables, varying by hour of day
dfs_corr <- list()
for (j in 1:length(hrs)) {
  hr <- hrs[j]
  dt_w <- dt_run[dt_run$Hour==hr, w_names]
  tab_corr <- cor(dt_w, use = "pairwise.complete.obs")
  df_cor <- as.data.frame(tab_corr)
  df_cor$Variable <- rownames(df_cor)
  rownames(df_cor) <- NULL
  df_cor$Hour <- hr
  dfs_corr[[j]] <- df_cor
}
stacked_df_corr <- do.call(rbind, dfs_corr)
write.csv(stacked_df_corr, file.path(savedir, "tab_corr_byhour.csv"), row.names=FALSE)

# CTRF models --------------------------------------------------
for (hr in hrs) {
  savepath <- file.path(savedir, paste0("ctrf_h", hr, ".txt"))
  dt_hr <- data_finalize(dt_run, hour=hr)
  # workday's additional response as a constant
  ctrf_0 <- lm(
    log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+workday
    +ntd+ntd:ts+ntd:sin1+ntd:cos1
    +prcp+prcp:ts+prcp:sin1+prcp:cos1
    +rhum+rhum:ts+rhum:sin1+rhum:cos1
    +wsp+wsp:ts+wsp:sin1+wsp:cos1
    +skc+skc:ts+skc:sin1+skc:cos1, data=dt_hr
  )
  # workday's additional response as a function of temp
  ctrf_d <- lm(
    log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
    +workday+workday:ts+workday:sin2+workday:cos1
    +ntd+ntd:ts+ntd:sin1+ntd:cos1
    +prcp+prcp:ts+prcp:sin1+prcp:cos1
    +rhum+rhum:ts+rhum:sin1+rhum:cos1
    +wsp+wsp:ts+wsp:sin1+wsp:cos1
    +skc+skc:ts+skc:sin1+skc:cos1, data=dt_hr
  )
  print_output(ctrf_0, ctrf_d, filepath=savepath)
}

# Complement --------------------------------
# 1) allow  different (p, q)
# 2) compute fitted values at temp in seq(0, 30, 5)
# 3) compute minimum temp

# Define functions to compute fitted values and find minimum from TRF
create_sim_mx = function(x_seq, tmax, tmin, order) {
  # order=c(p, q) the numbers of pairs and poly; ncol = 2q+p+1
  p <- order[1]; q <- order[2]
  r <- x_seq
  s <- (r - tmin)/(tmax - tmin)
  mx1 <- matrix(sapply(1:p, function(i) s^i), ncol=p, byrow=FALSE)
  mx2 <- matrix(sapply(1:q, function(i) c(sin(2*i*pi*s), cos(2*i*pi*s))), ncol=2*q, byrow=FALSE)
  ones <- matrix(1, nrow=length(s), ncol=1)
  cbind(ones, mx1, mx2)
}
compute_gfit <- function(TRF, order, x_seq, t_range){
  # TRF model must match the order
  X_mx <- create_sim_mx(x_seq, t_range[2], t_range[1], order)
  p <- order[1]; q <- order[2]
  vars_1 <- c("ts", "ts2", "ts3", "ts4")
  vars_2 <- c("sin1", "cos1", "sin2", "cos2", "sin3", "cos3", "sin4", "cos4")
  vars <- c("(Intercept)", vars_1[1:p], vars_2[1:(2*q)])
  coef <- TRF$coefficients[vars]
  g_fit <- X_mx %*% coef
  df <- data.frame(temp=x_seq, y_est=g_fit)
  return(df)
}

get_BPT <- function(TRF, order, t_range){
  a <- t_range[1]; b <- t_range[2]
  p <- order[1]; q <- order[2]
  vars_1 <- c("ts", "ts2", "ts3", "ts4")
  vars_2 <- c("sin1", "cos1", "sin2", "cos2", "sin3", "cos3", "sin4", "cos4")
  vars <- c("(Intercept)", vars_1[1:p], vars_2[1:(2*q)])
  coef <- TRF$coefficients[vars]
  
  trf_s <- function(s){
    x1 <- c(s, I(s^2), I(s^3), I(s^4))
    x2 <- c(sin(2*pi*s), cos(2*pi*s), sin(4*pi*s), cos(4*pi*s), 
            sin(6*pi*s), cos(6*pi*s), sin(8*pi*s), cos(8*pi*s))
    x <- c(1, x1[1:p], x2[1:(2*q)])
    g <- coef %*% x
    return(g[1])
  }
  ans <- optimize(trf_s, interval=c(0, 1))
  s_o <- ans$minimum
  r_o <- s_o*(b-a)+a
  y_min <- ans$objective
  list_min <- list(s_opt = s_o, 
                   BPT = r_o,
                   BPE = y_min)
  return(list_min)
}

hrs <- c(0, 3, 6, 8, 10, 12, 14, 16, 18, 20, 22)
df_lims <- data.frame(Hour=hrs, a=NA, b=NA)
for (i in c(1:length(hrs))){
  hr <- hrs[i]
  tl <- min(dt_run[dt_run$Hour==hr, "temperature"])
  th <- max(dt_run[dt_run$Hour==hr, "temperature"])
  df_lims[i, 2:3] <- c(tl, th)
}
print(df_lims)
df_lims$BPT <- NA; df_lims$BPE <- NA # create new columns to store optima
temps <- seq(-5, 30, 5) # common t range and gap
dfs_fit <- vector(mode="list", length(hrs)) # list of df_fit

# Divide hours by choice of (p, q)
#0, 8, 10, 22 | (3, 2)
#3:00 | (2, 2)
#6:00 | (1, 2)
#12, 14 | (3, 3) 
#16, 18, 20 | (2, 3)

for (i in c(1, 2, 3, 4, 5, 11)){ #0, 3, 6, 8, 10, 22
  hr <- hrs[i]
  savepath <- file.path(savedir, paste0("trf_h", hr, ".txt"))
  dt_hr <- data_finalize(dt_run, hour=hr)
  if (i==2) { #3
    trf <- lm(log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+ntd+workday, data=dt_hr)
    order <- c(2, 2)
  } else if (i==3) { #6
    trf <- lm(log(load) ~ ts+sin1+cos1+sin2+cos2+ntd+workday, data=dt_hr)
    order <- c(1, 2)
  } else { #0, 8, 10, 22
    trf <- lm(log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd+workday, data=dt_hr)
    order <- c(3, 2)
  }
  print_output(trf, filepath=savepath, text=paste0("BIC preferred at h", hr))
  t_lims <- unlist(df_lims[df_lims$Hour==hr, c("a", "b")])
  df_fit_hr <- compute_gfit(trf, order, temps, t_lims)
  colnames(df_fit_hr)[colnames(df_fit_hr)=="y_est"] <- paste0("y_est_", hr)
  dfs_fit[[i]] <- df_fit_hr
  df_lims$BPT[i] <- get_BPT(trf, order, t_lims)$BPT
  df_lims$BPE[i] <- get_BPT(trf, order, t_lims)$BPE
}

for (i in 6:10){ #12, 14, 16, 18, 20
  hr <- hrs[i]
  savepath <- file.path(savedir, paste0("trf_h", hr, ".txt"))
  dt_hr <- data_finalize(dt_run, hour=hr)
  if (i <= 7){ #12, 14
    trf <- lm(log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+sin3+cos3+ntd+workday, data=dt_hr)
    order <- c(3, 3)
  } else { #16, 18, 20
    trf <- lm(log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+sin3+cos3+ntd+workday, data=dt_hr)
    order <- c(2, 3)
  }
  print_output(trf, filepath=savepath, text="BIC preferred")
  
  t_lims <- unlist(df_lims[df_lims$Hour==hr, c("a", "b")])
  df_fit_hr <- compute_gfit(trf, order, temps, t_lims)
  colnames(df_fit_hr)[colnames(df_fit_hr)=="y_est"] <- paste0("y_est_", hr)
  dfs_fit[[i]] <- df_fit_hr
  df_lims$BPT[i] <- get_BPT(trf, order, t_lims)$BPT
  df_lims$BPE[i] <- get_BPT(trf, order, t_lims)$BPE
}

df_all <- Reduce(function(x, y) merge(x, y, by="temp", all=TRUE), dfs_fit)
write.csv(df_all, file.path(savedir, "tab_fitvalues.csv"), row.names=FALSE)
write.csv(df_lims, file.path(savedir, "tab_optima.csv"), row.names=FALSE)

## Workday subsamples ----------
savedir <- here(paste0("output_", region), "subsample")
dt_sub <- dt_run[dt_run$workday==1, ]
hrs <- c(6, 8, 10, 14, 18, 22)

# Divide hours by choice of (p, q)
hr <- hrs[1] #6:00
savepath <- file.path(savedir, paste0("ctrf_h", hr, ".txt"))
dt_hr <- data_finalize(dt_sub, hour=hr)
ctrf_12 <- lm(log(load) ~ ts+sin1+cos1+sin2+cos2
             +ntd+ntd:ts+ntd:sin1+ntd:cos1
             +prcp+prcp:ts+prcp:sin1+prcp:cos1
             +rhum+rhum:ts+rhum:sin1+rhum:cos1
             +wsp+wsp:ts+wsp:sin1+wsp:cos1
             +skc+skc:ts+skc:sin1+skc:cos1
             +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
             data=dt_hr)
print_output(ctrf_12, filepath=savepath, text="BIC preferred")

for (i in c(2:6)) { #8, 10, 14, 18, 22
  hr <- hrs[i]
  savepath <- file.path(savedir, paste0("ctrf_h", hr, ".txt"))
  dt_hr <- data_finalize(dt_sub, hour=hr)
  if (i <= 3){ #10:00
    ctrf <- lm(
      log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
      +ntd+ntd:ts+ntd:sin1+ntd:cos1
      +prcp+prcp:ts+prcp:sin1+prcp:cos1
      +rhum+rhum:ts+rhum:sin1+rhum:cos1
      +wsp+wsp:ts+wsp:sin1+wsp:cos1
      +skc+skc:ts+skc:sin1+skc:cos1
      +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
      data=dt_hr
    )
  } else {
    ctrf <- lm(
      log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+sin3+cos3
      +ntd+ntd:ts+ntd:sin1+ntd:cos1
      +prcp+prcp:ts+prcp:sin1+prcp:cos1
      +rhum+rhum:ts+rhum:sin1+rhum:cos1
      +wsp+wsp:ts+wsp:sin1+wsp:cos1
      +skc+skc:ts+skc:sin1+skc:cos1
      +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
      data=dt_hr
    )
  }
  print_output(ctrf, filepath=savepath, text="BIC preferred")
}

## Sum stats -----------
hrs <- seq(0, 23)
w_names <- c("load", "temperature", "dew_point_temperature", "wet_bulb_temperature", 
             "station_level_pressure", "precipitation", "relative_humidity", "wind_speed", "skycover")
col_names <- c("Variable", "Min", "Q1", "Median", "Mean", "Q3", "Max", "SD", "nobs", "Hour", "Season")

all_tab_stats1 <- list()
all_tab_stats2 <- list()
dt_sub1 <- dt_sub[dt_sub$Month %in% c(6, 7, 8), ]
dt_sub2 <- dt_sub[dt_sub$Month %in% c(1, 2, 12), ]
for (j in 1:length(hrs)) {
  hr <- hrs[j]
  tab_stats1 <- as.data.frame(matrix(NA, nrow=length(w_names), ncol=length(col_names)))
  colnames(tab_stats1) <- col_names
  tab_stats2 <- tab_stats1
  
  tab_stats1[, 11] <- "Summer"
  tab_stats1[, 10] <- hr
  for (i in 1:length(w_names)) {
    w <- w_names[i]
    tab_stats1[i, 2:7] <- summary(dt_sub1[dt_sub1$Hour==hr, w])
    tab_stats1[i, 8] <- sd(dt_sub1[dt_sub1$Hour==hr, w], na.rm=TRUE)
    tab_stats1[i, 9] <- sum(is.na(dt_sub1[dt_sub1$Hour==hr, w])==FALSE)
  }
  all_tab_stats1[[j]] <- tab_stats1
  
  tab_stats2[, 11] <- "Winter"
  tab_stats2[, 10] <- hr
  for (i in 1:length(w_names)) {
    w <- w_names[i]
    tab_stats2[i, 2:7] <- summary(dt_sub2[dt_sub2$Hour==hr, w])
    tab_stats2[i, 8] <- sd(dt_sub2[dt_sub2$Hour==hr, w], na.rm=TRUE)
    tab_stats2[i, 9] <- sum(is.na(dt_sub2[dt_sub2$Hour==hr, w])==FALSE)
  }
  all_tab_stats2[[j]] <- tab_stats2
}
stacked_tab <- do.call(rbind, c(all_tab_stats1, all_tab_stats2))
stacked_tab$Variable <- w_names
write.csv(stacked_tab, file.path(savedir, "tab_stats_byhour.sub.csv"), row.names=FALSE)

## Corr. -----------
hrs <- seq(0, 23)
w_names <- c("load", "temperature", "dew_point_temperature", "wet_bulb_temperature", 
             "station_level_pressure", "precipitation", "relative_humidity", "wind_speed", "skycover")
dfs_corr1 <- list()
dfs_corr2 <- list()
dt_sub1 <- dt_sub[dt_sub$Month %in% c(6, 7, 8), ]
dt_sub2 <- dt_sub[dt_sub$Month %in% c(1, 2, 12), ]
for (j in 1:length(hrs)) {
  hr <- hrs[j]
  dt_w <- dt_sub1[dt_sub1$Hour==hr, w_names]
  tab_corr <- cor(dt_w, use = "pairwise.complete.obs")
  df_cor <- as.data.frame(tab_corr)
  df_cor$Variable <- rownames(df_cor)
  rownames(df_cor) <- NULL
  df_cor$Hour <- hr
  df_cor$Season <- "Summer"
  dfs_corr1[[j]] <- df_cor
  
  dt_w <- dt_sub2[dt_sub2$Hour==hr, w_names]
  tab_corr <- cor(dt_w, use = "pairwise.complete.obs")
  df_cor <- as.data.frame(tab_corr)
  df_cor$Variable <- rownames(df_cor)
  rownames(df_cor) <- NULL
  df_cor$Hour <- hr
  df_cor$Season <- "Winter"
  dfs_corr2[[j]] <- df_cor
}
stacked_df_corr <- do.call(rbind, c(dfs_corr1, dfs_corr2))
write.csv(stacked_df_corr, file.path(savedir, "tab_corr_byhour.sub.csv"), row.names=FALSE)
# End.