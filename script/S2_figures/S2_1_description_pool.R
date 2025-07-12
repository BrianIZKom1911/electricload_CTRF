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
data_finalize <- function(df) { # create final variables
  df %>%
    fff33(., "temperature") %>% # FFF terms
    mutate(prcp = precipitation,
           rhum = stdize(relative_humidity), # normalize/ standardize variables
           wsp = stdize(wind_speed, median=TRUE),
           skc = stdize(skycover, median=TRUE),
           ntd = nmlize(n_day)) -> df
  df[df$prcp > 1, "prcp"] <- 1 # winsorize precipitation
  df$prcp <- df$prcp*10
  return(df)
}

regions <- c("NC", "SC", "Coast", "South")
source(file.path("script", "S2_figures", "1_boots_R6class.R")) # Source the OOP script
# NC ---------------------------------------------------
region <- regions[1]
dt_run <- readRDS(here("script", "S1_modelanalysis", paste0(region, "_run.RDS")))
savedir <- here(paste0("output_", region))

## Description --------------
## Series: temp and elec
ggplot(data=dt_run, aes(x=Date))+
  geom_line(aes(y=load), color="orangered", linewidth=0.25, show.legend=FALSE)+
  geom_line(aes(y=load_ma), color="grey", linewidth=0.5, show.legend=FALSE)+
  scale_y_continuous(name="Megawatt")+
  labs(x="")
ggsave(file=file.path(savedir , "plot_series_load.png"), width=8.5, height=0.618*6.5)

ggplot(dt_run, aes(x=Date))+
  geom_line(aes(y=temperature), color="turquoise3", linewidth=0.25, show.legend=FALSE)+
  scale_y_continuous(name="degree Celsius")+
  labs(x="")
ggsave(file=file.path(savedir , "plot_series_temp.png"), width=8.5, height=0.618*6.5)

## Distribution ------------
## Omit this part when doing the other regions. ##
library(viridis) #new
library(hrbrthemes)

# Load, by month
ggplot(data=dt_run, aes(x=factor(Month), y=load_dt, fill=factor(Month)))+
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  scale_y_continuous(name="Megawatt")+
  theme_ipsum() +
  theme(legend.position="none") +
  xlab("month")
ggsave(file=file.path(savedir , "box_load_dt.mon.png"), bg="white", width=8, height=0.618*8)
# Load, by hour
ggplot(data=dt_run, aes(x=factor(Hour), y=load_dt, fill=factor(Hour)))+
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  scale_y_continuous(name="Megawatt")+
  theme_ipsum() +
  theme(legend.position="none") +
  xlab("hour")
ggsave(file.path(savedir , "box_load_dt.hr.png"), bg="white", width=8, height=0.618*6.5)

# temperature, by month
ggplot(data=dt_run, aes(x=factor(Month), y=temperature, fill=factor(Month)))+
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  scale_y_continuous(name="degree C")+
  theme_ipsum() +
  theme(legend.position="none") +
  xlab("month")
ggsave(file=file.path(savedir , "box_temp.mon.png"), bg="white", width=8, height=0.618*8)
# temperature, by hour
ggplot(data=dt_run, aes(x=factor(Hour), y=temperature, fill=factor(Hour)))+
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  scale_y_continuous(name="degree C")+
  theme_ipsum() +
  theme(legend.position="none") +
  xlab("hour")
ggsave(file.path(savedir , "box_temp.hr.png"), bg="white", width=8, height=0.618*6.5)

# relative_humidity, by month
ggplot(data=dt_run, aes(x=factor(Month), y=relative_humidity, fill=factor(Month)))+
  geom_boxplot()+
  scale_fill_viridis(discrete=TRUE, alpha=0.6, option="A")+
  #geom_jitter(color="black", size=0.4, alpha=0.3)+
  scale_y_continuous(name="%")+
  theme_ipsum()+
  theme(legend.position="none")+
  xlab("month")
ggsave(file=file.path(savedir , "box_rhum.mon.png"), bg="white", width=8, height=0.618*8)
# relative_humidity, by hour
ggplot(data=dt_run, aes(x=factor(Hour), y=relative_humidity, fill=factor(Hour)))+
  geom_boxplot()+
  scale_fill_viridis(discrete=TRUE, alpha=0.6, option="A")+
  scale_y_continuous(name="%")+
  theme_ipsum()+
  theme(legend.position="none")+
  xlab("hour")
ggsave(file=file.path(savedir , "box_rhum.hr.png"), bg="white", width=10, height=0.618*8)

## Correlation ------------
# Scatter plot: elec and temp
# (1) load minus MA
plot_xy_hr <- ggplot()+ 
  geom_point(data = dt_run, 
             mapping = aes(x=temperature, y=load_dt, color=factor(Hour)), alpha=0.3)+
  labs(x="Temperature C", y="Load (MW)", color="hour")+
  theme_classic()
#ggsave(file=file.path(savedir , "plot_temploaddt.hr.png"), width=10, height=0.618*10)

## Bootstrap --------------
tl <- floor(quantile(dt_run$temperature, prob=0.001, na.rm=TRUE))
th <- ceiling(quantile(dt_run$temperature, prob=0.999, na.rm=TRUE))
#yl <- floor(quantile(dt_run$load_dt, prob=0.001, na.rm=TRUE))
#yh <- ceiling(quantile(dt_run$load_dt, prob=0.999, na.rm=TRUE))
dt_run %>%
  filter(temperature>=tl & temperature<=th) %>%
  data_finalize(.) %>%
  mutate(logy = log(load), 
         s_t = ts*ntd, s2_t = ts2*ntd, sin1_t = sin1*ntd, cos1_t = cos1*ntd,
         s_p = ts*prcp, s2_p = ts2*prcp, sin1_p = sin1*prcp, cos1_p = cos1*prcp,
         s_h = ts*rhum, s2_h = ts2*rhum, sin1_h = sin1*rhum, cos1_h = cos1*rhum,
         s_w = ts*wsp, s2_w = ts2*wsp, sin1_w = sin1*wsp, cos1_w = cos1*wsp,
         s_c = ts*skc, s2_c = ts2*skc, sin1_c = sin1*skc, cos1_c = cos1*skc) -> dt_pool

## Run
# TRF ## Pooled results are necessary
# (1) y is load_dt
x_terms <- c("ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2")
dt_pool %>% 
  select(all_of(c("load_dt", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "load_dt", x_terms)
trf$create_sim_matrix(c(tl, th), c(tl, th), c(3, 2))
# First argument is the range of plotting curve and must be a subset of
# the second argument which must equal or contain the range of computing Fourier transform
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
trf_plot <- trf$plot(base_plot=plot_xy_hr)
ggsave(file=file.path(savedir , "trf_plot_loaddt.bt.png"), bg="white", width=10, height=0.618*10)

# TRF (2) -- regressing out weather covariates for regional comparison
x_terms <- c(
  "ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2", 
  "ntd", "prcp", "rhum", "wsp", "skc"
)
dt_pool %>% 
  select(all_of(c("logy", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "logy", x_terms)
trf$create_sim_matrix(c(tl, th), c(tl, th), c(3, 2))
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
df_trf_plot <- trf$aggregated_data
df_trf_plot <- df_trf_plot[, c("temperature", "CI_lb", "fv_median", "CI_ub")]
write.csv(df_trf_plot, file=file.path(savedir , "NC_trf_bt2.csv"), row.names=FALSE)

# CTRF # DO NOT DO POOLED CTRF.
# Coast -----------------------------------------------------
## Description --------------
region <- regions[3]
savedir  <- here(paste0("output_", region))
dt_run <- readRDS(here("script", "S1_modelanalysis", paste0(region, "_run.RDS")))
## Series: temp and elec
ggplot(data=dt_run, aes(x=Date))+
  geom_line(aes(y=load), color="orangered", linewidth=0.25, show.legend=FALSE)+
  geom_line(aes(y=load_ma), color="grey", linewidth=0.5, show.legend=FALSE)+
  scale_y_continuous(name="Megawatt")+
  labs(x="")
ggsave(file=file.path(savedir , "plot_series_load.png"), width=8.5, height=0.618*6.5)

## Bootstrap --------------
dt_run <- dt_run[dt_run$Ike==0, ] # Nobs=200874
# Add moving average and demeaned load after removing Ike
n <- sum(dt_run$Year==2002) + 1
dt_run$load_ma <- data.table::frollmean(dt_run$load, n, align="right", na.rm=TRUE)
avg2002 <- mean(dt_run$load[1:(n-1)])
dt_run$load_ma[1:(n-1)] <- avg2002
dt_run$load_dt <- dt_run$load - dt_run$load_ma
plot_xy_hr <- ggplot()+ 
  geom_point(data = dt_run, 
             mapping = aes(x=temperature, y=load_dt, color=factor(Hour)), alpha=0.3)+
  labs(x="Temperature C", y="Load (MW)", color="hour")+
  theme_classic()

tl <- floor(quantile(dt_run$temperature, prob=0.001, na.rm=TRUE))
th <- ceiling(quantile(dt_run$temperature, prob=0.999, na.rm=TRUE))
yl <- floor(quantile(dt_run$load_dt, prob=0.001, na.rm=TRUE))
yh <- ceiling(quantile(dt_run$load_dt, prob=0.999, na.rm=TRUE))
dt_run %>%
  filter(temperature>=tl & temperature<=th & load_dt>=yl) %>%
  data_finalize(.) %>%
  mutate(logy = log(load)) -> dt_pool # Nobs=200430

# TRF - (1) y is load_dt
x_terms <- c("ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2")
dt_pool %>% 
  select(all_of(c("load_dt", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "load_dt", x_terms)
trf$create_sim_matrix(c(tl, th), c(tl, th), c(3, 2))
# First argument is the range of plotting curve and must be a subset of
# the second argument which must equal or contain the range of computing Fourier transform
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
trf_plot <- trf$plot(base_plot=plot_xy_hr)
ggsave(file=file.path(savedir , "trf_plot_loaddt.bt.png"), bg="white", width=10, height=0.618*10)

# TRF (2) logy -- regressing out weather covariates for regional comparison
x_terms <- c(
  "ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2", 
  "ntd", "prcp", "rhum", "wsp", "skc"
)
dt_pool %>% 
  select(all_of(c("logy", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "logy", x_terms)
trf$create_sim_matrix(c(tl, th), c(tl, th), c(3, 2))
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
df_trf_plot <- trf$aggregated_data
df_trf_plot <- df_trf_plot[, c("temperature", "CI_lb", "fv_median", "CI_ub")]
write.csv(df_trf_plot, file=file.path(savedir , "Coast_trf_bt2.csv"), row.names=FALSE)

# South ----------------------------------------------------
## Description --------------
region <- regions[4]
savedir  <- here(paste0("output_", region))
dt_run <- readRDS(here("script", "S1_modelanalysis", paste0(region, "_run.RDS")))
## Series: temp and elec
ggplot(data=dt_run, aes(x=Date))+
  geom_line(aes(y=load), color="orangered", linewidth=0.25, show.legend=FALSE)+
  geom_line(aes(y=load_ma), color="grey", linewidth=0.5, show.legend=FALSE)+
  scale_y_continuous(name="Megawatt")+
  labs(x="")
ggsave(file=file.path(savedir , "plot_series_load.png"), width=8.5, height=0.618*6.5)

plot_xy_hr <- ggplot()+ 
  geom_point(data = dt_run, 
             mapping = aes(x=temperature, y=load_dt, color=factor(Hour)), alpha=0.3)+
  labs(x="Temperature C", y="Load (MW)", color="hour")+
  theme_classic()

## Bootstrap  --------------
tl <- floor(quantile(dt_run$temperature, prob=0.001, na.rm=TRUE))
th <- ceiling(quantile(dt_run$temperature, prob=0.999, na.rm=TRUE))
dt_run %>%
  filter(temperature>=tl) %>%
  data_finalize(.) %>%
  mutate(logy = log(load)) -> dt_pool # Nobs=201426
a <- min(dt_pool$temperature); b <- max(dt_pool$temperature)
## TRF - (1) y is load_dt
x_terms <- c("ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2")
dt_pool %>% 
  select(all_of(c("load_dt", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "load_dt", x_terms)
trf$create_sim_matrix(c(tl, th), c(a, b), c(3, 2)) # Caveat holds as above
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
trf_plot <- trf$plot(base_plot=plot_xy_hr)
ggsave(file=file.path(savedir , "trf_plot_loaddt.bt.png"), bg="white", width=10, height=0.618*10)

## TRF (2) logy
# regressing out weather covariates for regional comparison
x_terms <- c(
  "ts", "ts2", "ts3", "sin1", "cos1", "sin2", "cos2", 
  "ntd", "prcp", "rhum", "wsp", "skc"
)
dt_pool %>% 
  select(all_of(c("logy", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "logy", x_terms)
trf$create_sim_matrix(c(tl, th), c(a, b), c(3, 2))
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
df_trf_plot <- trf$aggregated_data
df_trf_plot <- df_trf_plot[, c("temperature", "CI_lb", "fv_median", "CI_ub")]
write.csv(df_trf_plot, file=file.path(savedir , "South_trf_bt2.csv"), row.names=FALSE)

# SC -------------------------------------------------------
region <- regions[2]
savedir  <- here(paste0("output_", region))
dt_run <- readRDS(here("script", "S1_modelanalysis", paste0(region, "_run.RDS")))
dt_run <- dt_run[dt_run$temperature<=45, ] # Data error (2025/04/20) Nobs=201610
plot_xy_hr <- ggplot()+ 
  geom_point(data = dt_run, 
             mapping = aes(x=temperature, y=load_dt, color=factor(Hour)), alpha=0.3)+
  labs(x="Temperature C", y="Load (MW)", color="hour")+
  theme_classic()

## Bootstrap  -----------
tl <- floor(quantile(dt_run$temperature, prob=0.001, na.rm=TRUE))
th <- ceiling(quantile(dt_run$temperature, prob=0.999, na.rm=TRUE))
yl <- floor(quantile(dt_run$load_dt, prob=0.001, na.rm=TRUE))
yh <- ceiling(quantile(dt_run$load_dt, prob=0.999, na.rm=TRUE))
dt_run %>%
  filter(temperature>=tl & temperature<=th & load_dt<=yh) %>%
  data_finalize(.) %>%
  mutate(logy = log(load)) -> dt_pool

## TRF - (1) y is load_dt
x_terms <- c("ts", "ts2", "sin1", "cos1", "sin2", "cos2", "sin3", "cos3")
dt_pool %>% 
  select(all_of(c("load_dt", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "load_dt", x_terms)
trf$create_sim_matrix(c(tl, th), c(tl, th), c(2, 3)) # Caveat holds as above
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
trf_plot <- trf$plot(base_plot=plot_xy_hr)
ggsave(file=file.path(savedir , "trf_plot_loaddt.bt.png"), bg="white", width=10, height=0.618*10)

## TRF (2) logy
# regressing out weather covariates for regional comparison
x_terms <- c(
  "ts", "ts2", "sin1", "cos1", "sin2", "cos2", "sin3", "cos3", 
  "ntd", "prcp", "rhum", "wsp", "skc"
)
dt_pool %>% 
  select(all_of(c("logy", x_terms))) %>%  # keep only needed variables
  na.omit() -> dt_b

trf <- TRF_Model$new(dt_b, "logy", x_terms)
trf$create_sim_matrix(c(tl, th), c(tl, th), c(2, 3))
system.time({
  trf$run_bootstrap(1000)
})
trf$aggregate_results()
df_trf_plot <- trf$aggregated_data
df_trf_plot <- df_trf_plot[, c("temperature", "CI_lb", "fv_median", "CI_ub")]
write.csv(df_trf_plot, file=file.path(savedir , "SC_trf_bt2.csv"), row.names=FALSE)
# End of script.