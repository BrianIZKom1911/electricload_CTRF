# Region: North Central
# Draw Ridgeline plots for met.variables--PRCP*, RH, WSP, and SKC by seasons
# PRCP*: conditional on PRCP>0
library(here)
library(dplyr)
library(ggplot2)
library(ggridges)
library(viridis)
library(hrbrthemes)

# Locate
savepath <- here("extraoutput", "appfigures")
dir.create(savepath, showWarnings=FALSE)
# Import Region: North Central
data <- readRDS(file=here("script", "S1_models", "NC_run.RDS"))
# Add Month- and Hourfactors based on numeric values
data$month_nom <- month.name[data$Month]
data$month_nom <- factor(
    data$month_nom, 
    levels=c(
        "January", "February", "March", "April", "May", "June",  
        "July", "August", "September", "October", "November", "December"
    ), 
    ordered=TRUE
)
data$Hour <- sprintf("%02d", data$Hour)
data$Hour <- factor(data$Hour, level=sprintf("%02d", 0:23), ordered=TRUE)

#--------------------------------------------------
# Outcome: Load detrended, by Month
ggplot(data, aes(x=load_dt, y=month_nom, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="viridis")+
    scale_y_discrete(limits=rev, name=NULL)+ # flip the order
    scale_x_continuous(name="Load (relative to MA)")+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines"),
        strip.text.x=element_text(size=8)
    )
ggsave(file=file.path(savepath, "rdg_loaddt.mon.png"), bg="white", height=6.18, width=6.18)
# By Hour
ggplot(data, aes(x=load_dt, y=Hour, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="viridis")+
    scale_y_discrete(limits=rev, name="Time")+ # flip the order
    scale_x_continuous(name="Load (relative to MA)")+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines"),
        strip.text.x=element_text(size=8)
    )
ggsave(file=file.path(savepath, "rdg_loaddt.hr.png"), bg="white", height=8, width=4.6)

## Temperature, by Month ------------
ggplot(data, aes(x=temperature, y=month_nom, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="plasma")+
    #scale_fill_gradientn(colours = heat.colors(100))+
    scale_y_discrete(limits=rev, name=NULL)+
    scale_x_continuous(breaks=seq(-20, 40, 10), name="Temperature (C)")+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines"),
        strip.text.x=element_text(size=8)
    )
ggsave(file=file.path(savepath, "rdg_temp.mon.png"), bg="white", height=6.18, width=6.18)
# By Hour
ggplot(data, aes(x=temperature, y=Hour, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="plasma")+
    scale_y_discrete(limits=rev, name="Time")+
    scale_x_continuous(breaks=seq(-20, 40, 10), name="Temperature (C)")+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines"),
        strip.text.x=element_text(size=8)
    )
ggsave(file=file.path(savepath, "rdg_temp.hr.png"), bg="white", height=8, width=4.6)

## Relative Humidity, by Month
ggplot(data, aes(x=relative_humidity, y=month_nom, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="magma", direction=-1)+ # reverse the colors' order
    scale_y_discrete(limits=rev, name=NULL)+
    scale_x_continuous(breaks=seq(0, 100, 20), name="Relative Humidity (%)")+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines")
    )
ggsave(file=file.path(savepath, "rdg_rhum.mon.png"), bg="white", height=6.18, width=6.18)
# By Hour
ggplot(data, aes(x=relative_humidity, y=Hour, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="magma", direction=-1)+
    scale_y_discrete(limits=rev, name="Time")+
    scale_x_continuous(breaks=seq(0, 100, 20), name="Relative Humidity (%)")+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines")
    )
ggsave(file=file.path(savepath, "rdg_rhum.hr.png"), bg="white", height=8, width=4.6)

## Windspeed, by Month
ggplot(data, aes(x=wind_speed, y=month_nom, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="magma", direction=-1)+ # reverse the colors' order
    scale_y_discrete(limits=rev, name=NULL)+
    scale_x_continuous(breaks=seq(0, 12, 2), name="Windspeed (m/s)")+
    coord_cartesian(xlim=c(0, 12))+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines")
    )
ggsave(file=file.path(savepath, "rdg_wsp.mon.png"), bg="white", height=6.18, width=6.18)
# By Hour
ggplot(data, aes(x=wind_speed, y=Hour, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="magma", direction=-1)+
    scale_y_discrete(limits=rev, name="Time")+
    scale_x_continuous(breaks=seq(0, 14, 2), name="Windspeed (m/s)")+
    coord_cartesian(xlim=c(0, 14))+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines")
    )
ggsave(file=file.path(savepath, "rdg_wsp.hr.png"), bg="white", height=8, width=4.6)

## Skycover, by Month -------------
ggplot(data, aes(x=skycover, y=month_nom, fill=..x..))+
    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
    scale_fill_viridis(option="magma", direction=-1)+ # reverse the colors' order
    scale_y_discrete(limits=rev, name=NULL)+
    scale_x_continuous(name="Skycover")+
    coord_cartesian(xlim=c(0.0, 1.0))+
    theme_ipsum()+
    theme(
        legend.position="none",
        panel.spacing=unit(0.1, "lines")
    )
ggsave(file=file.path(savepath, "rdg_skc.mon.png"), bg="white", height=6.18, width=6.18)
# By Hour
#ggplot(data, aes(x=skycover, y=Hour, fill=..x..))+
#    geom_density_ridges_gradient(scale=3, rel_min_height=0.01)+
#    scale_fill_viridis(option="magma", direction=-1)+
#    scale_y_discrete(limits=rev, name="Time")+
#    scale_x_continuous(name="Skycover")+
#    coord_cartesian(xlim=c(0.0, 1.0))+
#    theme_ipsum()+
#    theme(
#        legend.position="none",
#        panel.spacing=unit(0.1, "lines")
#    )

#--------------------------------------------------
# More Figures: Hourly by seasons
# Add an season column, select only winter and summer
data <- data |>
    dplyr::mutate(season = ifelse(Month %in% c(6, 7, 8), "Summer", NA)) |>
    dplyr::mutate(season = ifelse(Month %in% c(1, 2, 12), "Winter", season))
    
ggplot(data=subset(data, !is.na(season)), aes(x=skycover, y=Hour, fill=season))+
    geom_density_ridges(scale=3, rel_min_height=0.01, alpha=0.6)+
    scale_fill_manual(values = c("Summer"="orange", "Winter"="steelblue1"), name=NULL)+
    scale_y_discrete(limits=rev, name="Time")+
    scale_x_continuous(name="Skycover")+
    coord_cartesian(xlim=c(0.0, 1.0))+
    theme_ipsum()+
    theme(legend.position="bottom", panel.spacing=unit(0.1, "lines"))
ggsave(file=file.path(savepath, "rdg_skc.2season.png"), bg="white", height=8.5, width=4.5)
