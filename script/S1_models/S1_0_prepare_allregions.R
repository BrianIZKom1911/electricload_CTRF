# This script generate final data for analysis and give descriptive summary from various aspects
rm(list=ls())

library(here)
library(tidyverse)
library(data.table) 

geten_dummies <- function(data, colname){
  values <- unique(data[[colname]])
  for (val in values) {
    data[[paste0(colname, "_", val)]] <- as.numeric(data[[colname]] == val)
  }
  return(data)
}

regions <- c("NC", "SC", "Coast", "South")

# Import main data and process --------------------------------------
# Import Federal and Texas holidays
df_hds <- readRDS("holiday_dates.RDS")

for (region in regions) {
  cat("Finalizing ", region, " data...\n")
  dt_run <- read_csv(here("data_clean", paste0(region, "_main.csv")))
  # Process datetime
  dt_run$Date <- as.Date(dt_run$Date) # date in Central Time
  dt_run$datetime_UTC <- as.POSIXct(dt_run$datetime_UTC, tz="UTC")
  dt_run$datetime_CPT <- as.character(dt_run$datetime_CPT)
  dt_run$datetime_CPT <- gsub(" UTC", "", dt_run$datetime_CPT) # remove extra characters
  
  # Add moving average and demeaned load
  n <- sum(dt_run$Year==2002) + 1
  dt_run$load_ma <- data.table::frollmean(dt_run$load, n, align="right", na.rm=TRUE)
  avg2002 <- mean(dt_run$load[1:(n-1)])
  dt_run$load_ma[1:(n-1)] <- avg2002
  dt_run$load_dt <- dt_run$load - dt_run$load_ma
  # Add variables related to days and hours
  dt_run$n_day <- as.numeric(dt_run$Date - min(dt_run$Date)) + 1
  dt_run$day_week <- weekdays(dt_run$Date) # day of the week
  dt_run <- geten_dummies(dt_run, "Month")
  colnames(dt_run)[grepl("Month_", colnames(dt_run))] <- sprintf("m%02d", 1:12)
  dt_run <- geten_dummies(dt_run, "Hour")
  colnames(dt_run)[grepl("Hour_", colnames(dt_run))] <- c(sprintf("h%02d", 1:23), "h00")
  
  # Add workdays and daytime hours
  dt_run <- merge(dt_run, df_hds, by.x=c("Year", "Date"), by.y=c("year", "date"), all.x=TRUE)
  # Import sunrise and sunset times
  df_sun <- readRDS(paste0(region, "_sunlighttimes.RDS"))
  dt_run <- merge(dt_run, df_sun, by="Date")
  dt_run %>%
    mutate(offworkday = as.numeric(!is.na(holiday) | day_week %in% c("Saturday", "Sunday")),
           workday = 1 - offworkday,
           daylight = as.numeric(datetime_UTC >= sunrisetime & datetime_UTC <= sunsettime)) -> dt_run
  # Save as RDS
  dt_run <- dt_run[order(dt_run$datetime_UTC), ]
  row.names(dt_run) <- NULL
  saveRDS(dt_run, file=file.path(md, "S1_modelanalysis", paste0(region, "_run.RDS")))
  cat("Done.\n")
}

# Mark Ike Hurricane for Coast Region
dt_run <- readRDS(file.path(md, "S1_modelanalysis", "Coast_run.RDS"))
dt_run$Ike <- as.numeric(dt_run$Year==2008 & dt_run$m09==1)
saveRDS(dt_run, file=file.path(md, "S1_modelanalysis", "Coast_run.RDS"))

# Summarize data ------------------------------------------------
w_names <- c("load", "temperature", "dew_point_temperature", "wet_bulb_temperature", 
             "station_level_pressure", "precipitation", "relative_humidity", "wind_speed", "skycover")
x <- c(0, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10) # precipitation level
q <- c(0.005, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.995) # n-quantile for temp

for (region in regions) {
  savedir <- file.path(md, paste0("output_", region))
  if (dir.exists(savedir)) {
    cat(region, " Summary:\nSave directory exists.\n")
  } else {
    dir.create(savedir)
    cat(region, " Summary:\nSave directory created.\n")
  }
  dt_run <- readRDS(paste0(region, "_run.RDS"))
  
  if (region=="Coast"){
    dt_run <- dt_run[dt_run$Ike==0, ] # remove Ike for Coast
  }
  # Summary statistics: load, temp, prcp, rhum, wsp, skc  
  tab_stats <- as.data.frame(matrix(NA, nrow=length(w_names), ncol=9))
  colnames(tab_stats) <- c("Variables", "Min", "Q1", "Median", "Mean", "Q3", "Max", "SD", "nobs")
  for (i in 1:length(w_names)){
    w <- w_names[i]
    tab_stats[i, 2:7] <- summary(dt_run[[w]])
    tab_stats[i, 8] <- sd(dt_run[[w]], na.rm=TRUE)
    tab_stats[i, 9] <- sum(is.na(dt_run[[w]])==FALSE)
  }
  tab_stats$Variables <- w_names
  write.csv(tab_stats, file.path(savedir, "tab_statsum.csv"), row.names=FALSE)
  cat("Summary statistics saved.\n")
  
  # Correlation table
  dt_w <- dt_run[, w_names]
  tab_corr <- cor(dt_w, use = "pairwise.complete.obs")
  df_corr <- as.data.frame(tab_corr)
  df_corr$Variable <- rownames(df_corr)
  rownames(df_corr) <- NULL
  write.csv(df_corr, file.path(savedir, "tab_corr.csv"), row.names=FALSE)
  cat("Correlation table saved.\n")
  
  ## Precipitation frequency
  ## needed to determine the threshold of winsorizing precipitation 
  freq_prcp <- data.frame(threshold=x)
  freq_prcp$count <- sapply(x, function(a) sum(dt_run$precipitation>a))
  n <- sum(is.na(dt_run$precipitation)==FALSE)
  freq_prcp$frequncy <- (freq_prcp$count/n)*100
  write.csv(freq_prcp, file.path(savedir, "freq_prcp.csv"), row.names=FALSE)
  cat("Precipitation frequency saved.\n")
  
  ## Temperature quantiles (04/05)
  dt_run %>%
    group_by(Month, Hour) %>%
    summarise(
      tmin = min(temperature), 
      tmax = max(temperature),
      quantiles = list(map(q, ~quantile(temperature, prob=.x))),
      .groups = "drop"
    ) %>%
    unnest_wider(quantiles, names_sep = "_") -> df_qtemp
  colnames(df_qtemp)[grepl("quantiles_", colnames(df_qtemp))] <- sprintf("q_%03d", q*1000)
  write.csv(df_qtemp, file.path(savedir, "qn_temp.m_hr.csv"), row.names=FALSE)
  cat("Temeprature quantiles saved.\n")
}
# End of script.