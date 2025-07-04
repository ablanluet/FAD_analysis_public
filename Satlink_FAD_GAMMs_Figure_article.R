
#  Figure GAM
#  ------------------------------------------------------------------------
#   Arthur Blanluet
#   04/07/2025
#  ------------------------------------------------------------------------
#   Input
#   GAM results
#  ------------------------------------------------------------------------
#   Output 
#   Results Figures
#  ------------------------------------------------------------------------

 # rm(list=ls())

library(tidyverse) # for general data wrangling and plotting
library(lubridate) # for working with dates
library(reshape2)
library(mgcv) # GAMMs
library(sf) # for vector data 
library(visreg)
library(mgcViz)
library(gratia)
library(viridis)  # better colors for everyone
library(raster) # for working with rasters
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(parallel)
library(ggpubr)




# load MPA polygon
MPA_palmyra <- read_sf("input/shapefiles/Palmyra_Kingman.shp")


position_kingman <- as.data.frame(cbind(c(-162.416667),c(6.383333)))
names(position_kingman) <- c("Longitude", "Latitude")

sample_biomass <- read.csv("input/sample_biomass_tot2_04102023.csv",sep=';', dec = '.')
sample_position <- read.csv("input/sample_position_tot2_04102023.csv",sep=';', dec = '.')


load(file="output/gam_mod_tweedie_tot_28112023.RData")



# Format for model all data

sample_biomass$date_time_local <- as.POSIXct(sample_biomass$date_time_local, tz = "Pacific/Midway")

sample_biomass$date <- as.POSIXct(format(sample_biomass$date_time_local, "%x") ,format="%x", tz = "Pacific/Midway")
sample_biomass$time <- as.POSIXct(format(sample_biomass$date_time_local, "%H:%M") ,format="%H:%M", tz = "Pacific/Midway")
sample_biomass$month <-month(sample_biomass$date_time_local)

year_1 <- which(sample_biomass$year == "2021" | sample_biomass$year == "2022" & sample_biomass$month %in% c(1:7))

sample_biomass$years <- 2
sample_biomass$years[year_1] <- 1

 #sample_biomass <- subset(sample_biomass, years == "1")

sample_biomass <- subset(sample_biomass, Kiribati_EEZ == "0")
# 
 sample_biomass <- subset(sample_biomass, date_time_local < "2023-07-01 00:00:00")
 
 

sample_biomass$first_time_position <- as.POSIXct(sample_biomass$first_time_position, tz = "Pacific/Midway")
#sample_biomass <- subset(sample_biomass, time_in_dataset > 1)


sample_biomass_vertical <- melt(sample_biomass, id.vars = c("FAD_Name", "Name", "date_time_local","Latitude","Longitude", "time_in_MPA","distance_to_land","date","time","month","current_velocity",
                                                            "thermocline_depth", "temperature", "sal","month","Chl","season","first_day_position","Speed",
                                                            "first_day_biomass","Moon_Illumination","distance_to_nearest_dFAD","presence_in_MPA","time_in_dataset",
                                                            "SST_front","Chl_front","years") #, 
                                ,measure.vars =c( "Layer1","Layer2", "Layer3", "Layer4","Layer5","Layer6","Layer7","Layer8","Layer9","Layer10")
                                , variable.name = "Layer", value.name = "Biomass")

sample_biomass_vertical$Name <- as.factor(sample_biomass_vertical$Name)
sample_biomass_vertical$FAD_Name <- as.factor(sample_biomass_vertical$FAD_Name)
sample_biomass_vertical$depth <- as.factor(sample_biomass_vertical$Layer)
sample_biomass_vertical$first_day_position <- as.factor(sample_biomass_vertical$first_day_position)
sample_biomass_vertical$first_day_biomass <- as.factor(sample_biomass_vertical$first_day_biomass)
sample_biomass_vertical$years <- as.factor(sample_biomass_vertical$years)

sample_biomass_vertical$time_day <- hour(sample_biomass_vertical$time) + minute(sample_biomass_vertical$time)/60


levels(sample_biomass_vertical$depth) <- c(1:10)*10
sample_biomass_vertical$depth <- as.numeric(as.character(sample_biomass_vertical$depth))
sample_biomass_vertical$month <- as.factor(sample_biomass_vertical$month)
sample_biomass_vertical$presence_in_MPA <- as.factor(sample_biomass_vertical$presence_in_MPA)
sample_biomass_vertical$season <- as.factor(sample_biomass_vertical$season)

### dumie data
# sample_dumie <- as.data.frame("Dumie", "Dumie", sample_biomass_vertical$date_time_local,  mean(sample_biomass_vertical$Latitude),  mean(sample_biomass_vertical$Longitude),
# mean(sample_biomass_vertical$time_in_MPA, na.rm=T),mean(sample_biomass_vertical$distance_to_land) , as.POSIXct(mean(sample_biomass_vertical$date)) ,as.POSIXct(mean(sample_biomass_vertical$time))  , '5',
# mean(sample_biomass_vertical$current_velocity, na.rm=T),mean(sample_biomass_vertical$thermocline_depth, na.rm=T),mean(sample_biomass_vertical$temperature, na.rm=T),
# mean(sample_biomass_vertical$sal, na.rm=T),5,mean(sample_biomass_vertical$Chl, na.rm=T),'fall', 'SE', 'quantile_3',mean(sample_biomass_vertical$Moon_Illumination, na.rm=T),
# mean(sample_biomass_vertical$distance_to_nearest_dFAD, na.rm=T), 'Layer10', 0, 100,5  )





sample_dumie <- sample_biomass_vertical[1,]

sample_dumie$FAD_Name <- "Dumie"# day time of max biomass
sample_dumie$Name <- "Dumie"
sample_dumie$date_time_local <- mean(sample_biomass_vertical$date_time_local)
sample_dumie$Latitude <- mean(sample_biomass_vertical$Latitude)
sample_dumie$Longitude <-  mean(sample_biomass_vertical$Longitude)
sample_dumie$time_in_MPA <- mean(sample_biomass_vertical$time_in_MPA, na.rm=T)
sample_dumie$distance_to_land <- mean(sample_biomass_vertical$distance_to_land)
sample_dumie$date <- mean(sample_biomass_vertical$date)
sample_dumie$time <- mean(sample_biomass_vertical$time)
sample_dumie$month <- 10
sample_dumie$current_velocity <- mean(sample_biomass_vertical$current_velocity, na.rm=T)
sample_dumie$thermocline_depth <- mean(sample_biomass_vertical$thermocline_depth, na.rm=T)
sample_dumie$temperature <- mean(sample_biomass_vertical$temperature, na.rm=T)
sample_dumie$sal <- mean(sample_biomass_vertical$sal, na.rm=T)
sample_dumie$month <- 10
sample_dumie$Chl  <- mean(sample_biomass_vertical$Chl, na.rm=T)
sample_dumie$season <- 'fall'
sample_dumie$first_day_position <- 'SE'
sample_dumie$first_day_biomass <-'quantile_1'
sample_dumie$Moon_Illumination <- mean(sample_biomass_vertical$Moon_Illumination, na.rm=T)
sample_dumie$distance_to_nearest_dFAD  <- mean(sample_biomass_vertical$distance_to_nearest_dFAD, na.rm=T)
sample_dumie$Layer <- 'Layer5'
sample_dumie$Biomass <- 0
sample_dumie$depth <- 110
sample_dumie$time_day <- mean(sample_biomass_vertical$time_day, na.rm=T)
sample_dumie$time_in_dataset <- mean(sample_biomass_vertical$time_in_dataset, na.rm=T)
sample_dumie$Speed <- mean(sample_biomass_vertical$Speed, na.rm=T)
sample_dumie$SST_front <- mean(sample_biomass_vertical$SST_front, na.rm=T)
sample_dumie$Chl_front <- mean(sample_biomass_vertical$Chl_front, na.rm=T)
sample_dumie$years <- '1'




sample_biomass_vertical <- rbind(sample_biomass_vertical, sample_dumie)


# 
##### Figure #######

time_in_MPA <- visreg(gam_mod_tweedie_season, "time_in_dataset", scale='response', #  by = "presence_in_MPA",, time_in_dataset
                     partial = F, gg = T, rug = F, overlay = TRUE)+
   geom_rug(sides="b",outside = FALSE, alpha = 1/2)+
   scale_y_continuous("Biomass (t)")+ #limits = c(NA,0.2),
   scale_x_continuous("Time (days)")+
   # geom_vline(xintercept=1)+
   # geom_vline(xintercept=2)+
  coord_cartesian(xlim =c(NA,30), ylim =c(0.05,0.13))+ #, ylim = c(NA, 0.15)
  theme(legend.position = c(0.9, 0.2),legend.key.size = unit(0.1, 'cm'))+
  theme_bw()#+
  # theme(text=element_text(size = 24))



# Speed <- visreg(gam_mod_tweedie_season, "Speed", scale='response', partial = F, rug =F, gg = T)+
#   # geom_rug(sides="b")+
#    scale_y_continuous("Biomass (t)")+
#    scale_x_continuous("Current speed (m/s)")+
#   theme_bw()#+
#theme(text=element_text(size = 24))

Chl <- visreg(gam_mod_tweedie_season, "Chl", scale='response', partial = F, rug =F, gg = T)+
  geom_rug(sides="b", alpha = 0.1)+
  scale_y_continuous("Biomass (t)")+
  scale_x_continuous(bquote('Chlorophyll (mg.'~ m^{-3}~')'))+
  coord_cartesian(ylim =c(0.079,0.141))+
  theme_bw()#+
  # theme(text=element_text(size = 24))

sal <- visreg(gam_mod_tweedie_season, "sal", scale='response', partial = F, rug = F, gg = T)+
   geom_rug(sides="b", alpha = 0.1)+
scale_x_continuous("Salinity")+ 
  scale_y_continuous("Biomass (t)")+
  coord_cartesian(ylim =c(0.05,0.13))+
  theme_bw()# +
  # theme(text=element_text(size = 24))
# sal
# ggsave(filename ="output/figure_rapport/sal.png", height = 8.5, width = 8.5, dpi = 300)

 
temperature <- visreg(gam_mod_tweedie_season, "temperature", scale='response', partial = F, rug = F, gg = T)+
  geom_rug(sides="b", alpha = 0.1)+
  scale_y_continuous("Biomass (t)")+
scale_x_continuous("Sea surface temperature (°C)")+
  coord_cartesian(ylim =c(0.08,0.18))+
 theme_bw() #+
  # theme(text=element_text(size = 24))
  # temperature
# 
  # ggsave(filename ="output/figure_rapport/temperature_0603.png", height = 8.5, width = 8.5, dpi = 300)


thermocline_depth <- visreg(gam_mod_tweedie_season, "thermocline_depth", scale='response', partial = F, rug = F, gg = T)+
  geom_rug(sides="b", alpha = 0.1)+
  scale_y_continuous("Biomass (t)")+
  scale_x_continuous("Thermocline depth (m)")+
  coord_cartesian(ylim =c(0.07,0.14))+
 theme_bw()#+
  # theme(text=element_text(size = 24))
# thermocline_depth
# 
 
current_velocity <- visreg(gam_mod_tweedie_season, "current_velocity", scale='response', partial = F, rug = F, gg = T)+
  geom_rug(sides="b", alpha = 0.1)+
scale_x_continuous(bquote('Current velocity (m. '~ s^{-1}~')'))+
  scale_y_continuous("Biomass (t)")+
  coord_cartesian(ylim =c(0.09,0.30))+
  theme_bw()#+
  # theme(text=element_text(size = 24))

# scale_x_continuous(bquote('Chlorophyll (mg.'~ m^{-3}~')'))+


Moon_Illumination <- visreg(gam_mod_tweedie_season, "Moon_Illumination", scale='response', partial = F, rug = F, gg = T)+
  geom_rug(sides="b", alpha = 0.1)+
  scale_y_continuous("Biomass (t)")+
  scale_x_continuous("Moon fraction")+
  coord_cartesian(ylim =c(0.078,0.151))+
   theme_bw()#+
 # theme(text=element_text(size = 24))
# 
# Moon_Illumination
# 
# ggsave(filename ="output/figure_rapport/Moon_Illumination_0603.png", height = 8.5, width = 8.5, dpi = 300)

distance_to_nearest_dFAD <- visreg(gam_mod_tweedie_season, "distance_to_nearest_dFAD", scale='response', partial = F, rug = F, gg = T)+
  geom_rug(sides="b", alpha = 0.1)+
  scale_x_continuous("Distance to nearest dFAD (km)")+
  scale_y_continuous("Biomass (t)")+
  # scale_y_continuous(limits = c(NA,1.5),"Biomass (t)")+
  coord_cartesian(xlim =c(NA,100), ylim = c(0.07, 0.125))+
theme_bw() #   +

distance_to_nearest_dFAD

distance_to_nearest_dFAD_2 <- visreg(gam_mod_tweedie_season, "distance_to_nearest_dFAD", scale='response', partial = F, rug = F)
  

# 
Years <-visreg(gam_mod_tweedie_season, "years", scale='response', partial = F, rug = F, gg = T)+
  scale_y_continuous("Biomass (t)")+
  # scale_y_continuous(limits = c(0,5))+
  # scale_x_discrete(labels=c("Layer1" = "0-14.2", "Layer2" = "14.2-25.4","Layer3" = "25.4-36.6","Layer4" = "36.6-47.8","Layer5" = "47.8-59",
  #                           "Layer6" = "59-70.2","Layer7" = "70.2-81.4","Layer8" = "81.4-92.6","Layer9" = "92.6-103.8","Layer10" = "103.8-115"), 'Layer depth (m)')+
  theme_bw()#+
# theme(text=element_text(size = 24))
# Years
# ggsave(filename ="output/figure_rapport/Years_comp.png", height = 8.5, width = 8.5, dpi = 300)

Years$scales$scales[[1]]$labels <- c("2021-2022","2022-2023")
Years$scales$scales[[1]]$name <- "Years"

layer <-visreg(gam_mod_tweedie_season, "Layer", scale='response', partial = F, rug = F, gg = T)+
  scale_y_continuous("Biomass (t)")+
  scale_x_discrete("MPA")+
  # scale_y_continuous(limits = c(0,5))+
  # scale_x_discrete(labels=c("Layer1" = "0-14.2", "Layer2" = "14.2-25.4","Layer3" = "25.4-36.6","Layer4" = "36.6-47.8","Layer5" = "47.8-59",
  #                           "Layer6" = "59-70.2","Layer7" = "70.2-81.4","Layer8" = "81.4-92.6","Layer9" = "92.6-103.8","Layer10" = "103.8-115"), 'Layer depth (m)')+
  theme_bw()
# theme(text=element_text(size = 24), axis.text.x = element_text(face="bold", angle=45,hjust=1))
 # layer
 # 

sample_biomass_vertical$presence_in_MPA  <- fct_relevel(sample_biomass_vertical$presence_in_MPA, "1","0")

Presence_MPA <-visreg(gam_mod_tweedie_season, "presence_in_MPA", scale='response', partial = F, rug = F, gg = F)+
  scale_y_continuous("Biomass (t)")+
  #scale_x_discrete("MPA")+
  # scale_y_continuous(limits = c(0,5))+
  # scale_x_discrete(labels=c("Layer1" = "0-14.2", "Layer2" = "14.2-25.4","Layer3" = "25.4-36.6","Layer4" = "36.6-47.8","Layer5" = "47.8-59",
  #                           "Layer6" = "59-70.2","Layer7" = "70.2-81.4","Layer8" = "81.4-92.6","Layer9" = "92.6-103.8","Layer10" = "103.8-115"), 'Layer depth (m)')+
  theme_bw()#+
 # 
 # 
 # ggsave(filename ="output/figure_rapport/layer_22042023.png", height = 8.5, width = 8.5, dpi = 300)


Presence_MPA$scales$scales[[1]]$labels <- c("Inside","Outside")
Presence_MPA$scales$scales[[1]]$name <- "MPA"


sample_biomass_vertical$season  <- fct_relevel(sample_biomass_vertical$season, "summer", "fall", "winter", "spring")


season <-visreg(gam_mod_tweedie_season, "season", scale='response', partial = F, rug = F, gg = T)+
  scale_y_continuous("Biomass (t)")+
  xlab("Season")+
  #scale_x_discrete(limits = c("winter","spring", "summer", "fall"))+
  theme_bw()#+
  # theme(text=element_text(size = 24))
 season
 #ggsave(filename ="output/figure_rapport/season_0603.png", height = 8.5, width = 8.5, dpi = 300)

 season$scales$scales[[1]]$labels <- c( "Jul-Sep", "Oct-Dec", "Jan-Mar", "Apr-Jun")
 
 random_effect_2 <- visreg(gam_mod_tweedie_season, "FAD_Name" , scale='response', rug = F)

 ggplot(random_effect_2$fit, aes(x=visregFit)) + geom_histogram(binwidth=0.05)+
   theme_bw()
 
 test <- random_effect_2$fit
 
 subset(test, visregFit > 0.5)
 # ISL+281009
 # SLX+347749
 # SLX+353752
 # SLX+364281
 # SLX+391031
 # 
 # SLX+382119 => fishing? 
 
random_effect <- visreg(gam_mod_tweedie_season, "FAD_Name" , scale='response', partial = F, rug = F, gg = T)+
  # scale_y_continuous(limits = c(0,0.04))+
  theme_bw()#+
  # theme(text=element_text(size = 24))


daytime <-visreg(gam_mod_tweedie_season, "time_day", overlay = TRUE, scale='response',
                 by = "Layer", partial = F, rug =F, gg = T)+   #  
  # geom_rug(sides="b")+
  scale_y_continuous("Biomass (t)")+
  scale_x_continuous(name = "Time of day (h)", breaks=seq(0, 24, 4))+
theme_bw()

daytime$scales$scales[[1]]$name <- "Depth Layer (m)"
daytime$scales$scales[[1]]$labels <- c("3-14.2", "14.2-25.4","25.4-36.6", "36.6-47.8","47.8-59", "59-70.2","70.2-81.4","81.4-92.6","92.6-103.8","103.8-115")
daytime$scales$scales[[2]]$name <- "Depth Layer (m)"
daytime$scales$scales[[2]]$labels <- c("3-14.2", "14.2-25.4","25.4-36.6", "36.6-47.8","47.8-59", "59-70.2","70.2-81.4","81.4-92.6","92.6-103.8","103.8-115")

# daytime$scales$scales[[1]]$name <- "Mean Layer depth (m)"
# daytime$scales$scales[[1]]$labels <- c("9.1", "20.3","31.5", "42.7","53.9", "65.1","76.3","87.5","98.7","109.9")
# daytime$scales$scales[[2]]$name <- "Mean Layer depth (m)"
# daytime$scales$scales[[2]]$labels  <- c("9.1", "20.3","31.5", "42.7","53.9", "65.1","76.3","87.5","98.7","109.9")

# daytime_2 <- daytime+
#   guides(fill = guide_legend(ncol=1))+
#   theme(legend.position = c(0.8, 0.75), 
#         legend.background = element_rect(fill = "white", color = "black"),
#         legend.key.size = unit(0.18, "cm"),
#         legend.title=element_text(size=5.5),
#         legend.text = element_text(size=4))

daytime_2 <- daytime+
   guides(fill = guide_legend(ncol=1))+
  theme(legend.position = c(0.8, 0.65), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.key.size = unit(0.2, "cm"),
        legend.title=element_text(size=6),
        legend.text = element_text(size=5.5))


#  theme(text=element_text(size = 24))  
leg <- get_legend(daytime_2) 
leg_2 <- as_ggplot(leg)

#   #legend.position = "none",
#   
# ggsave(filename = "output/daytime_layer_22042023.png", height = 8.5, width = 8.5, dpi = 300)

#time_in_MPA <- time_in_MPA +  theme(legend.position = "none")

daytime_3 <- daytime +  theme(legend.position = "none")

 
  Presence_MPA + Years + season +  #time_in_MPA+ #Front_SSL+ Front_chl+ #+presence_in_MPA+  #distance_to_land  + 
  daytime_2 + time_in_MPA + Moon_Illumination + distance_to_nearest_dFAD + 
  current_velocity + sal + Chl + temperature + thermocline_depth +  #first_day_position + first_day_biomass +   +#Speed+ #random_effect+
  plot_layout(ncol = 3, nrow = 4) + plot_annotation(tag_levels = "a")
  
ggsave(filename = "output/results_GAM_21032024_2.png", height = 8.5, width = 8.5, dpi = 900)

