#install.packages("suncalc")
rm(list = ls())
library(suncalc)
#library(dplyr)

# NC: Dallas, 32.7792, -96.8089
# SC: San Antonio, 29.4250, -98.4939
# Coast: Houston, 29.7628, -95.3831
# South: McAllen, 26.2164	-98.2364

regions <- c("NC", "SC", "Coast", "South")
cities <- c("Dallas", "San Antonio", "Houston", "McAllen")
city_lats <- c(32.7792, 29.4250, 29.7628, 26.2164)
city_lons <- c(-96.8089, -98.4939, -95.3831, -98.2364)
df_info <- data.frame(regions, cities, city_lats, city_lons)

date_range <- seq(from=as.Date("2002-01-01"), to=as.Date("2025-01-01"), by="day")
for (region in regions){
  city <- df_info[df_info$regions==region, "cities"]
  city_lat <- df_info[df_info$regions==region, "city_lats"]
  city_lon <- df_info[df_info$regions==region, "city_lons"]
  df_sun <- data.frame(
    Date = date_range,
    sunrisetime = suncalc::getSunlightTimes(
      date=date_range, lat=city_lat, lon=city_lon, tz="UTC", keep="sunrise"
    )$sunrise, 
    sunsettime = suncalc::getSunlightTimes(
      date=date_range, lat=city_lat, lon=city_lon, tz="UTC", keep="sunset"
    )$sunset
  )
  saveRDS(df_sun, file=paste0(region, "_sunlighttimes.RDS"))
  print(paste0(city, " sunlighttime is saved."))
}
