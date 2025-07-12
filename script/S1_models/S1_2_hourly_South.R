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
region <- regions[4]
savedir <- here(paste0("output_", region), "hourly")
if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}
dt_run <- readRDS(here("script", "S1_models", paste0(region, "_run.RDS")))

# TRF models --------------------------------------------------
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
