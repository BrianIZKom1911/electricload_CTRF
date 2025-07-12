# Continue with '_allregions.R' and regional running data
# Model (C)TRF, but I need to choose models seperately for each region
# Moreover, only Coast's load shows impact of Ike Hurricane 
rm(list=ls())

library(here)
library(dplyr)
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

# NC -----------------------------------------------------
region <- regions[1]
savedir <- here(paste0("output_", region))
if (!dir.exists(savedir)) {
  dir.create(savedir)
}
dt_run <- readRDS(paste0(region, "_run.RDS"))
# Get FFF terms
# Winsorize precipitation
dt_pool$prcp <- dt_pool$precipitation
dt_pool[dt_pool$prcp > 1, "prcp"] <- 1
# Get normalized/standardized variables
dt_pool$prcp <- dt_pool$prcp*10
## prcp = (rainfall - 0)/0.1. prcp=1 <=> rainfall=0.1; P(rainfall>0.1)=7.3% (NC)
dt_pool$rhum <- stdize(dt_pool$relative_humidity)
dt_pool$wsp <- stdize(dt_pool$wind_speed, median=TRUE)
dt_pool$skc <- stdize(dt_pool$skycover, median=TRUE)
dt_pool$ntm <- nmlize(dt_pool$n_month)
dt_pool$ntd <- nmlize(dt_pool$n_day)

## TRF -----------
## (0) detrended, no controls
trf_32.dt <- lm(load_dt ~ ts+ts2+ts3+sin1+cos1+sin2+cos2, data=dt_pool)
## (0.1) log(y), timetrend
trf_32 <- lm(log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd, data=dt_pool)
## (1) Month effects
trf_32.m <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, # Dec omitted
  data=dt_pool
)
print_output(
  trf_32.dt, trf_32, trf_32.m, 
  filepath=file.path(savedir, "trf_32.txt"), 
  text="\n(1) Demeaned load\n(2) Linear timetrend controlled\n(3) Month effects"
)

## (2) diurnal factors
trf_32_dl <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd+workday+daylight
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, # Dec omitted
  data=dt_pool
)
trf_32_dli <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd
  +workday*daylight # interaction
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, # Dec omitted
  data=dt_pool
)
trf_32_h <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd+workday
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11 # Dec omitted
  +h00+h01+h02+h03+h04+h05+h06+h07+h08+h10+h11 # 9:00 omitted
  +h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23,
  data=dt_pool
)
print_output(
  trf_32_dl, trf_32_dli, trf_32_h, 
  filepath=file.path(savedir, "trf_32_dayhour.txt"), 
  text="\n(1) Workday and daytime\n(2) Interaction effects\n(3) Workday and hour effects"
)

## (3) Weather covariates
trf_32_dwli <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
  +ntd+workday+prcp+rhum+skc+wsp+daylight+daylight:skc
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, # Dec omitted
  data=dt_pool
)
trf_32_dwh <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
  +ntd+workday+prcp+rhum+skc+wsp
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11 # Dec omitted
  +h00+h01+h02+h03+h04+h05+h06+h07+h08+h10+h11 # 9:00 omitted
  +h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23,
  data=dt_pool
)
print_output(
  trf_32_dwli, trf_32_dwh, 
  filepath=file.path(savedir, "trf_32_dwl.txt"), 
  text="Diurnal factors and weather covariates"
)
# It shows that skycover only affects daytime use

## CTRF -----------
ctrf_m <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+workday+daylight
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, # Dec omitted
  data=dt_pool
)
ctrf_md <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+daylight
  +workday+workday:ts+workday:sin2+workday:cos1
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, # Dec omitted 
  data=dt_pool
)
sink(file.path(savedir, "ctrf_d.txt"))
cat("(1) Workday effect as a constant\n")
coeftest(ctrf_m, vcov = sandwich::vcovHC(ctrf_m, type="HC3"))
glance(ctrf_m)
cat("\n(2) Workday effect as a function of temp\n")
coeftest(ctrf_md, vcov = sandwich::vcovHC(ctrf_md, type="HC3"))
glance(ctrf_md)
sink(file = NULL)

## Subsample -------------
# Workdays: non-holioday weekdays.
# Active daytime, 8:00--18:00; active nighttime 19:00--23:00
# Daylight hours are identified by location and day
savedir <- here(paste0("output_", region), "subsample")
if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}
dt_sub <- dt_pool[(dt_pool$workday==1) & (dt_pool$Hour %in% c(8:23)), ]

trf_32 <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd+daylight
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_sub
)
trf_32_h <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11
  +h08+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23, # 9:00 omitted
  data=dt_sub
)
print_output(
  trf_32, trf_32_h, 
  filepath = file.path(savedir, "trf_32.txt"), 
  text = "BIC preferred.\nSubsample: Workday active hours"
)

ctrf_32 <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+daylight
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_sub
)
ctrf_32_h <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11
  +h08+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23, # 9:00 omitted
  data=dt_sub
)
print_output(
  ctrf_32, ctrf_32_h, 
  filepath = file.path(savedir, "ctrf_32.txt"), 
  text = "With/without hour effects\nSubsample: Workday active hours"
)

# South --------------------------------------------------
region <- regions[4]
savedir <- here(paste0("output_", region))
if (!dir.exists(savedir)) {
  dir.create(savedir)
}
dt_run <- readRDS(paste0(region, "_run.RDS"))
# Get FFF terms
dt_pool <- fff33(dt_run, "temperature")
# Winsorize precipitation
dt_pool$prcp <- dt_pool$precipitation
dt_pool[dt_pool$prcp > 1, "prcp"] <- 1
# Get normalized/standardized variables
dt_pool$prcp <- dt_pool$prcp*10
## prcp = (rainfall - 0)/0.1. prcp=1 <=> rainfall=0.1; P(rainfall>0.1)=6% (South)
dt_pool$rhum <- stdize(dt_pool$relative_humidity)
dt_pool$wsp <- stdize(dt_pool$wind_speed, median=TRUE)
dt_pool$skc <- stdize(dt_pool$skycover, median=TRUE)
dt_pool$ntd <- nmlize(dt_pool$n_day)

## TRF -----------
trf_32 <- lm(log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd, data=dt_pool)

## (1) Month effects
trf_32.m <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, # Dec omitted
  data=dt_pool
)
print_output(trf_32, trf_32.m, filepath=file.path(savedir, "trf_32.txt"), text="Month effects")

## (2) Diurnal factors
trf_32_dl <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
  +ntd+workday+daylight
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11, 
  data=dt_pool
)
trf_32_h <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd+workday
  +h00+h01+h02+h03+h04+h05+h06+h07+h08+h10+h11 # 9:00 omitted
  +h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_pool
)
trf_32_dli <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd
  +workday*daylight # interaction
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_pool
)
print_output(
  trf_32_dl, trf_32_h, trf_32_dli, 
  filepath=file.path(savedir, "trf_32_dayhour.txt"), 
  text="Diurnal factors"
)

## (3) Weather covariates
trf_32_dwl <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
  +ntd+workday+daylight+prcp+rhum+skc+wsp
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_pool
)
trf_32_dwh <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
  +ntd+workday+prcp+rhum+skc+wsp
  +h00+h01+h02+h03+h04+h05+h06+h07+h08+h10+h11 # 9:00 omitted
  +h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_pool
)
trf_32_dwli <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2
  +ntd+workday+prcp+rhum+skc+wsp+daylight+daylight:skc
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_pool
)
print_output(
  trf_32_dwl, trf_32_dwh, trf_32_dwli, 
  filepath=file.path(savedir, "trf_32_dwl.txt"), 
  text="Diurnal factors and weather covariates"
)

## CTRF -----------
ctrf_m <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+workday+daylight
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_pool
)
ctrf_md <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+daylight
  +workday+workday:ts+workday:sin2+workday:cos1
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_pool
)
sink(file.path(savedir, "ctrf_d.txt"))
cat("(1) Workday effect as a constant\n")
coeftest(ctrf_m, vcov = sandwich::vcovHC(ctrf_m, type="HC3"))
glance(ctrf_m)
cat("\n(2) Workday effect as a function of temp\n")
coeftest(ctrf_md, vcov = sandwich::vcovHC(ctrf_md, type="HC3"))
glance(ctrf_md)
sink(file = NULL)

## Subsample ------------
savedir <- here(paste0("output_", region), "subsample")
if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}
dt_sub <- dt_pool[(dt_pool$workday==1) & (dt_pool$Hour %in% c(8:23)), ]
# TRF
trf_23 <- lm(
  log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+sin3+cos3+ntd+daylight
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_sub
)
trf_23_h <- lm(
  log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+sin3+cos3+ntd
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11
  +h08+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23, # 9:00 omitted
  data=dt_sub
)
print_output(
  trf_23, trf_23_h, 
  filepath = file.path(savedir, "trf_23.txt"), 
  text = "BIC preferred.\nSubsample: Workday active hours"
)

# CTRF
ctrf_23 <- lm(
  log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+sin3+cos3+daylight
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_sub
)
ctrf_23_h <- lm(
  log(load) ~ ts+ts2+sin1+cos1+sin2+cos2++sin3+cos3+
  +ntd+ntd:ts+ntd:sin1+ntd:cos1
  +prcp+prcp:ts+prcp:sin1+prcp:cos1
  +rhum+rhum:ts+rhum:sin1+rhum:cos1
  +wsp+wsp:ts+wsp:sin1+wsp:cos1
  +skc+skc:ts+skc:sin1+skc:cos1
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11
  +h08+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23, # 9:00 omitted
  data=dt_sub
)
print_output(
  ctrf_23, ctrf_23_h, 
  filepath = file.path(savedir, "ctrf_23.txt"), 
  text = "With/without hour effects\nSubsample: Workday active hours"
)

# SC -----------------------------------------------------------
region <- regions[2]
dt_run <- readRDS(paste0(region, "_run.RDS"))
# Get FFF terms
dt_pool <- fff33(dt_run, "temperature")
# Winsorize precipitation
dt_pool$prcp <- dt_pool$precipitation
dt_pool[dt_pool$prcp > 1, "prcp"] <- 1
# Get normalized/standardized variables
dt_pool$prcp <- dt_pool$prcp*10
## prcp = (rainfall - 0)/0.1. prcp=1 <=> rainfall=0.1; P(rainfall>0.1)=5.6% (SC)
dt_pool$rhum <- stdize(dt_pool$relative_humidity)
dt_pool$wsp <- stdize(dt_pool$wind_speed, median=TRUE)
dt_pool$skc <- stdize(dt_pool$skycover, median=TRUE)
dt_pool$ntd <- nmlize(dt_pool$n_day)

###########################################
# Only consider the subsample of workdays #
savedir <- here(paste0("output_", region), "subsample")
if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}
dt_sub <- dt_pool[(dt_pool$workday==1) & (dt_pool$Hour %in% c(8:23)), ]
# TRF
trf_23 <- lm(
  log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+sin3+cos3+ntd+daylight
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_sub
)
trf_23_h <- lm(
  log(load) ~ ts+ts2+sin1+cos1+sin2+cos2+sin3+cos3+ntd
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11
  +h08+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23, # 9:00 omitted
  data=dt_sub
)
print_output(
  trf_23, trf_23_h, 
  filepath = file.path(savedir, "trf_23.txt"), 
  text = "BIC preferred.\nSubsample: Workday active hours"
)

# Coast --------------------------------------------------------
region <- regions[3]
dt_run <- readRDS(paste0(region, "_run.RDS"))
# Drop Ike Hurricane
dt_run <- dt_run[dt_run$Ike==0, ] # Nobs=200874
# Get FFF terms
dt_pool <- fff33(dt_run, "temperature")
# Winsorize precipitation
dt_pool$prcp <- dt_pool$precipitation
dt_pool[dt_pool$prcp > 1, "prcp"] <- 1
# Get normalized/standardized variables
dt_pool$prcp <- dt_pool$prcp*10
## prcp = (rainfall - 0)/0.1. prcp=1 <=> rainfall=0.1; P(rainfall>0.1)=7.8% (Coast)
dt_pool$rhum <- stdize(dt_pool$relative_humidity)
dt_pool$wsp <- stdize(dt_pool$wind_speed, median=TRUE)
dt_pool$skc <- stdize(dt_pool$skycover, median=TRUE)
dt_pool$ntd <- nmlize(dt_pool$n_day)

savedir <- here(paste0("output_", region), "subsample")
if (!dir.exists(savedir)) {
  dir.create(savedir, recursive = TRUE)
}
dt_sub <- dt_pool[(dt_pool$workday==1) & (dt_pool$Hour %in% c(8:23)), ]

# TRF
trf_32 <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd+daylight
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11,
  data=dt_sub
)
trf_32_h <- lm(
  log(load) ~ ts+ts2+ts3+sin1+cos1+sin2+cos2+ntd
  +m01+m02+m03+m04+m05+m06+m07+m08+m09+m10+m11
  +h08+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+h21+h22+h23, # 9:00 omitted
  data=dt_sub
)
print_output(
  trf_32, trf_32_h, 
  filepath = file.path(savedir, "trf_32.txt"), 
  text = "BIC preferred.\nSubsample: Workday active hours"
)
# End of script.