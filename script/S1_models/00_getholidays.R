# Help to define workdays in the dataset
# Workdays are consistent over the state of Texas
rm(list=ls())
library(lubridate)
library(dplyr)

# Function to find the n-th weekday of a month
nth_weekday <- function(year, month, weekday, n) {
  first_day <- as.Date(paste(year, month, "01", sep = "-"))
  first_weekday <- first_day + (weekday - lubridate::wday(first_day)) %% 7
  first_weekday + (n - 1) * 7
}
# Function to find the last weekday of a month (since we don't know it's 4th or 5th)
last_weekday <- function(year, month, weekday) {
  date_01 <- as.Date(paste(year, month, "01", sep = "-"))
  numberofdays <- lubridate::days_in_month(date_01)
  last_day <- as.Date(paste(year, month, numberofdays, sep = "-"))
  date_last_weekday <- last_day - (lubridate::wday(last_day) - weekday) %% 7
  return(date_last_weekday)
}

# Function to find all holidays in a year
get_holidays <- function(year) {
  df <- data.frame(
    holiday = c("New Year's Day", "Martin Luther King Jr. Day", "Presidents' Day",
                 "Memorial Day", "Independence Day", "Labor Day", "Veterans Day", 
                 "Thanksgiving Day", "Christmas Day", "Texas Independence Day", 
                 "Juneteenth", "San Jacinto Day"),
    date = c(
      as.Date(paste(year, "01", "01", sep = "-")), # Fixed: New Year's Day
      nth_weekday(year, 1, 2, 3),                  # MLK: 3rd Monday of January
      nth_weekday(year, 2, 2, 3),                  # Presidents' Day: 3rd Monday of February
      last_weekday(year, 5, 2),                    # Memorial Day: Last Monday of May
      as.Date(paste(year, "07", "04", sep = "-")), # Fixed: Independence Day
      nth_weekday(year, 9, 2, 1),                  # Labor Day: 1st Monday of September
      as.Date(paste(year, "11", "11", sep = "-")), # Fixed: Veterans Day
      nth_weekday(year, 11, 5, 4),                 # Thanksgiving: 4th Thursday of November
      as.Date(paste(year, "12", "25", sep = "-")), # Fixed: Christmas Day
      as.Date(paste(year, "03", "02", sep = "-")), # Fixed: Texas Independence Day
      as.Date(paste(year, "06", "19", sep = "-")), # Fixed: Juneteenth
      as.Date(paste(year, "04", "21", sep = "-"))  # Fixed: San Jacinto Day
    )
  )
  return(df)
}

# Dataframe of all these holidays and dates (2000-2025)
years <- c(2000:2025)
df_hds <- data.frame()
for (y in years) {
  df <- get_holidays(y)
  df_hd <- cbind(data.frame(year=rep(y, 12)), df)
  df_hds <- rbind(df_hds, df_hd)
  print(y)
}

df_hds <- df_hds[order(df_hds$date), ]
rownames(df_hds) <- NULL
saveRDS(df_hds, file="holiday_dates.RDS")
