
library(sf)
library(sp)
library(gstat)
library(ggplot2)

library(sf)
library(viridis)
library(lubridate) # for working with dates
library(patchwork)
library(grid)
library(gridExtra)


MPA_palmyra <- read_sf("input/shapefiles/Palmyra_Kingman.shp")

sample_biomass <- read.csv("input/sample_biomass_tot_13082024.csv",sep=';', dec = '.')
sample_position <- read.csv("input/sample_position_tot_13082024.csv",sep=';', dec = '.')

sample_biomass$date_time_local <- ifelse(nchar( sample_biomass$date_time_local) == 10,
                                         paste( sample_biomass$date_time_local, "00:00:00"),
                                         sample_biomass$date_time_local)

sample_biomass$date_time_local <- as.POSIXct(sample_biomass$date_time_local,format = "%Y-%m-%d %H:%M", tz = "Pacific/Midway")

#sample_biomass$date_time_local <- ymd_hms(sample_biomass$date_time_local,  tz = "Pacific/Midway")

sample_biomass$date <- as.POSIXct(format(sample_biomass$date_time_local, "%x") ,format="%x", tz = "Pacific/Midway")
sample_biomass$time <- as.POSIXct(format(sample_biomass$date_time_local, "%H:%M") ,format="%H:%M", tz = "Pacific/Midway")
sample_biomass$month <-month(sample_biomass$date_time_local)

year_1 <- which(sample_biomass$year == "2021" | sample_biomass$year == "2022" & sample_biomass$month %in% c(1:7))
year_3 <- which(sample_biomass$year == "2024" | sample_biomass$year == "2023" & sample_biomass$month %in% c(8:12))


sample_biomass$years <- 2
sample_biomass$years[year_1] <- 1
sample_biomass$years[year_3] <- 3


sample_biomass <- subset(sample_biomass, years %in% c(1,2))

sample_biomass <- subset(sample_biomass, Kiribati_EEZ == "0")
# 
sample_biomass <- subset(sample_biomass, date_time_local < "2024-07-01 00:00:00")



# # spatial format
# coordinates(sample_biomass_max) <- ~ Longitude + Latitude
# proj4string(sample_biomass_max) <- CRS("+proj=longlat +datum=WGS84")
# class(sample_biomass_max)

coordinates(sample_biomass) <- ~ Longitude + Latitude
proj4string(sample_biomass) <- CRS("+proj=longlat +datum=WGS84")
class(sample_biomass)



sample_biomass_2 <- sample_biomass

df <- data.frame(Biomass = sample_biomass_2$ESSum, Longitude = sample_biomass_2@coords[,1], Latitude = sample_biomass_2@coords[,2])
# Organise dataframe

n <- 100 # Number of points to create a grid for estimating fish biomass. This can be changed.
X <- seq(min(sample_biomass@coords[,1]), max(sample_biomass@coords[,1]), length.out = n) # Lon points for creating grid
Y  <- seq(min(sample_biomass@coords[,2]), max(sample_biomass@coords[,2]), length.out = n) # Lat points for creating grid
new_df <- data.frame(expand.grid(Longitude = X, Latitude = Y)) # Coordinates where you want to estimate the fish biomass
coordinates(new_df)  <- ~ Longitude + Latitude

# Inverse weighted estimates based on Lon and Lat and 
idwmodel <- idw(formula = Biomass ~ 1, 
                locations = ~ Longitude + Latitude, 
                data = df, 
                newdata = new_df, 
                idp = 1) #0.75

# idwmodel <- idw(ESSum ~Latitude + Longitude, sample_biomass,grid,
#                maxdist = Inf, idp = 4)


idwmodel_2 <- data.frame(idwmodel)
idwmodel_2$pred_mean <- idwmodel_2$var1.pred - mean(idwmodel_2$var1.pred)

MPA_palmyra_other <-  MPA_palmyra%>%
  dplyr::filter(ORGNAME != "PACIFIC REMOTE ISLANDS MARINE NATIONAL MONUMENT") #%>%
# dplyr::filter(OBJECTID != "6164")

IDW_1 <- ggplot() +
  #labs(fill="Year 2021-2022") +
  geom_raster(data = idwmodel_2, aes(x=Longitude, y=Latitude, fill=pred_mean))+
  scale_fill_viridis() +
  geom_contour(data = idwmodel_2, aes(x=Longitude, y=Latitude, z = pred_mean), color = "black")+
  geom_sf(data =  MPA_palmyra,
          color = 'red',
          alpha = 0.2,
          size = 0.5,
          fill=NA )+
  geom_sf(data =  MPA_palmyra_other,
          color = 'green',
          alpha = 0.2,
          size = 0.5,
          fill=NA )+
  geom_polygon(data = data.frame(x = c(-162.04, -160.07, -160.07), y = c(4.1, 6.84, 4.1)),
               aes(x = x, y = y), fill = "white", color ="white")+
  geom_point(data =position_kingman, aes(x = Longitude, y = Latitude) , size = 2, color = 'green') +
  scale_x_continuous("Longitude (°)",expand = c(0,0),labels = ~ .x,breaks = c(-161,-162,-163),limits = c(NA,-160.07) )+ #limits = c(-164.1392,-160.0922)
  scale_y_continuous("Latitude (°)",expand = c(0,0),labels = ~ .x,breaks = c(5,6,7))+ #limits = c(4.131465,8.108201)
  #  coord_cartesian(xlim = c(NA,-160.07))+
  labs(fill = "Relative biomass \ndensity (t)")+
  theme_classic()+
  theme(text=element_text(size = 18) )#,axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0)

IDW_025 <- IDW_025 + ggtitle("Idp = 0.25")
IDW_050 <- IDW_050 + ggtitle("Idp = 0.5")
IDW_075 <- IDW_075 + ggtitle("Idp = 0.75")
IDW_1 <- IDW_1 + ggtitle("Idp = 1")

# ggsave(filename = "output/figure_article/map_IDW_power_1.png", height = 8.5, width = 8.5, dpi = 300)

IDW_025+ IDW_050+ IDW_075+ IDW_1+
  plot_layout(ncol = 2, nrow = 2) + plot_annotation(tag_levels = "a")

ggsave(filename = "output/figure_article/IDW_comparison.png", height = 12, width = 12, dpi = 900)

