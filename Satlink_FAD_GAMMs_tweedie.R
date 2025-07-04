#  GAM tweedie model
#  ------------------------------------------------------------------------
#   Arthur Blanluet
#   04/07/2025
#  ------------------------------------------------------------------------
#   Input
#   sample_biomass
#   sample_position
#   MPA polygone
#  ------------------------------------------------------------------------
#   Output 
#   tweedie_model_R.data (season or month)
#  ------------------------------------------------------------------------

start_time <- Sys.time()

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

library(tweedie)
library(statmod)

library(parallel)

# load MPA polygon

#MPA_palmyra <- read_sf(paste(path,"Palmyra_Kingman.shp",sep=""),sep=';', dec = '.')
# MPA_palmyra <- read_sf('E:/Save_Brisbane/carto/GIS/PalmyraAtoll_Kingman_FWS_NWR_Boundaries/Palmyra_Kingman.shp')
MPA_palmyra <- read_sf("input/shapefiles/Palmyra_Kingman.shp")

position_kingman <- as.data.frame(cbind(c(-162.416667),c(6.383333)))
names(position_kingman) <- c("Longitude", "Latitude")

sample_biomass <- read.csv("input/sample_biomass_tot_13082024.csv",sep=';', dec = '.')
sample_position <- read.csv("input/sample_position_tot_13082024.csv",sep=';', dec = '.')

# sample_biomass_max <- read.csv(paste(path,"sample_biomass_max.csv",sep=""),sep=';', dec = '.')


# Format for model all data



# Format for model all data

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
# 
 sample_biomass <- subset(sample_biomass, years %in% c(1,2))

sample_biomass <- subset(sample_biomass, Kiribati_EEZ == "0")
# 
 sample_biomass <- subset(sample_biomass, date_time_local < "2023-06-01 00:00:00")



sample_biomass$first_time_position <- as.POSIXct(sample_biomass$first_time_position, tz = "Pacific/Midway")
#sample_biomass <- subset(sample_biomass, time_in_dataset > 1)


sample_biomass_vertical <- melt(sample_biomass, id.vars = c("FAD_Name", "Name", "date_time_local","Latitude","Longitude", "time_in_MPA","distance_to_land","date","time","current_velocity",
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

sample_biomass_vertical$time_day <- hour(sample_biomass_vertical$date_time_local) + minute(sample_biomass_vertical$date_time_local)/60


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


large_distance <- which(sample_biomass_vertical$distance_to_nearest_dFAD > 100)
sample_biomass_vertical$distance_to_nearest_dFAD[large_distance] <- 100

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

# Format for model max data
# sample_biomass_max$date <- as.POSIXct(sample_biomass_max$date, tz = "Pacific/Midway")
# 
# sample_biomass_max$month <- as.factor(sample_biomass_max$month)
# sample_biomass_max$season <- as.factor(sample_biomass_max$season)
# sample_biomass_max$Name <- as.factor(sample_biomass_max$Name)
# sample_biomass_max$depth_group  <- as.factor(sample_biomass_max$depth_group )



# parallelisation 
# method 1: ,method = "fREML",discrete=TRUE,nthreads=c(2,1)
# method 2: 
# nc <- 2   ## cluster size, set for example portability
# if (detectCores()>1) { ## no point otherwise
#   cl <- makeCluster(nc) 
#   ## could also use makeForkCluster, but read warnings first!
# } else cl <- NULL
# , method = "REML",cluster=cl) # chunk size = 5000 ?

# model tweedie
# log link:
# bam(biomass ~ ...  ,family= tw(theta=NULL, link=log, a = 1.01, b=1.99)...
# poisson link:
# bam(log10(Biomass+1) ~ ...  ,family= tw(theta=NULL, link=power(0), a = 1.01, b=1.99)...


 nc <- 4   ## cluster size, set for example portability
 if (detectCores()>1) { ## no point otherwise
  cl <- makeCluster(nc) 
   ## could also use makeForkCluster, but read warnings first!
 } else cl <- NULL

gam_mod_tweedie_season <- bam( Biomass ~ 
                                s(time_in_dataset,k =6, bs = 'cr')+ # , by = season
                                #s(Chl, bs = 'cr',k =6)+ 
                                #s(sal,k =6,bs = 'cr')+  
                                #s(Speed,k =6,bs = 'cr')+#+
                                #s(temperature,k =6,bs = 'cr')# + s(distance_to_land,k =6,bs = 'cr')+
                                #s(thermocline_depth,k =6,bs = 'cr')  + 
                                #s(current_velocity,k =6,bs = 'cr')+  
                                s(time_day,k =6,by = Layer, bs = "cc") + #
                                s(FAD_Name, bs = 're') + # 
                                s(Moon_Illumination,k =6,bs = 'cr') + 
                                s(distance_to_nearest_dFAD,k =6,bs = 'cr')+
                                Layer + 
                                season + 
                                years +
                               # presence_in_MPA 
                               te(Longitude,Latitude) #+ 
                                #s(SST_front,k =6,bs = 'cr') #+ 
                                #s(Chl_front,k =6,bs = 'cr')
                               #+  #first_day_position+ first_day_biomass+
                                ,data = sample_biomass_vertical,
                                family= tw(theta=NULL, link=log, a = 1.01, b=1.99), 
                                method = "fREML",
                                cluster=cl,
                                chunk.size = 5000, 
                                control = gam.control(trace = TRUE)) 

 #   s(Longitude,Latitude) + 

summary(gam_mod_tweedie_season)

AIC(gam_mod_tweedie_season)

gam_mod_tweedie_season[["family"]][["family"]] # Tweedie(p=1.08) 


# concurvity(gam_mod_tweedie_season)
# concurvity(gam_mod_tweedie_season,full=FALSE)

save(gam_mod_tweedie_season, file="Output/gam_mod_tweedie_spatial_tensor_24022025_2021_2023.RData")
 
end_time <- Sys.time()
end_time - start_time



windows()
par(mfrow = c(2,2))
gam.check(gam_mod_tweedie_season)

dres <- resid(gam_mod_tweedie_season)
pres <- resid(gam_mod_tweedie_season, type = 'pearson')

qres1 <- qresid(gam_mod_tweedie_season)
qres2 <- qresid(gam_mod_tweedie_season)

windows()
par(mfrow = c(2,2))
stats::qqnorm(dres,main = "Deviance residuals")
stats::qqline(dres)

stats::qqnorm(pres,main = "Pearson residuals"); stats::qqline(pres)
stats::qqnorm(qres1,main = "Quantile residuals (set 1)");
stats::qqline(qres1)
stats::qqnorm(qres2,main = "Quantile residuals (set 2)");
stats::qqline(qres2)

