# Plot some functions shown in SM C

# 1. Critical value of RH --------------------------------------
# at given T and T_skin (K) for sweat to evaporate based on Clausius-Clapeyron.
# E.g. Let T_skin=308 K (34.85 C), c=100*exp(5423(1/T-1/308)) if T>=308

# Define
cRH_T_ <- function(T.vec, T_skin) {
    #' T.vec is an numeric vector of air temp
    #' T_skin is an number of skin temp
    #' Output is an vector of CV(T)
    ifelse(
        T.vec < T_skin, 
        100.5, 
        100*exp(5423*(1/T.vec - 1/T_skin))
    )
}

# Plot the graph and save
x <- seq(16, 42, 0.5); xk <- x+273.15
cvx1 <- cRH_T_(xk, T_skin=308)
cvx2 <- cRH_T_(xk, T_skin=292)
# Save
thisdir <- getwd()
png(
    filename = file.path(thisdir, "plot_cv.RH_T.png"), 
    width = 1000, # pixels
    height = 618, # pixels
    res = 150 # dpi
)
plot(
    x, cvx1, type="l", col="orangered", lwd=1.5,
    xlab = "Air temperature (C)",
    ylab = "RH (%)",
    main = "Critical RH",
    cex.lab = 1,
    cex.main = 1.2
)
lines(x, cvx2, col="royalblue", lwd=1.5)
# Show texts above curves
text(38, 84, "34.85", col="orangered", pos=3)
text(22, 82.5, "18.85", col="royalblue", pos=3)
dev.off()

# 2. N prop.to vapor pressure diff ----------------------------
library(ggplot2)

N_vp_ <- function(RH, T.vec, T_skin) {
    #' @param RH numeric vector, relative humidity \in [0, 100]
    #' @param T.vec numeric vector, air temperature in Kelvins
    #' @param T_skin number, skin temperature in Kelvins. usually chosen from [292, 308] 
    #' @return numeric vector of the \min(vp_diff, 100-RH)
    #' 
    dvp <- 100*exp(5423*(1/T.vec - 1/T_skin)) - RH
    dvp0 <- ifelse(dvp > 0, dvp, 0)
    return(dvp0)
}

# Simulation
x <- seq(24, 42, 0.5); xk <- x+273.15
y <- seq(65, 100, 0.5)
df_n <- expand.grid(temp_K=xk, rel_hum=y)
df_n$z <- N_vp_(df_n$rel_hum, df_n$temp_K, T_skin=308)
df_n$temp_C <- df_n$temp_K - 273.15

thisdir <- getwd()
ggplot(df_n, aes(x=temp_C, y=rel_hum, fill=z))+
    geom_tile()+
    scale_fill_viridis_c()+
    ggtitle("Evaporation strength")+
    labs(x="Air temperature (C)", y="RH (%)", fill=NULL)+
    annotate("text", x=40, y=90, label="Skin: 34.85", color="white")+
    theme_classic()+    
    theme(
        plot.title = element_text(hjust=0.5, size=15, face="bold"),
        plot.title.position = "plot",
        legend.position = "none"
    )
ggsave(file.path(thisdir, "hmap_evap_T_RH.png"), width=7, height=5)
# END
