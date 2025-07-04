#  Processing of dFADs data :
#  split between biomass and position data
#  extrapolation of biomass coordinate position
#  addition of oceanographic data to the biomass data
#  ------------------------------------------------------------------------
#   Arthur Blanluet
#   04/07/2025
#  ------------------------------------------------------------------------
#   Input
#   PalmyraExport_"month".csv: raw biomass/position data from satlink for each month
#   MPA_square.shp: shapefile of the PRIMNM boundaries
#   oceanographic_raster.RData: oceanographic raster map for each day and oceanographic variable (too heavy for the GitHub, see the oceanographic_data.r)
#  ------------------------------------------------------------------------
#   Output 
#   sample_biomass 
#   sample_position
#  ------------------------------------------------------------------------



library(tidyverse) # for general data wrangling and plotting
#library(furrr) # for parallel operations on lists
library(lubridate) # for working with dates
library(sf) # for vector data 
library(raster) # for working with rasters
# library(maps) # additional helpful mapping packages
# library(maptools)
# library(rgeos)
library(reshape2)
library(sp) 

library(suncalc)

library(geosphere)



# load oceangraphic data
# load("input/oceanographic_raster.RData") 
# load("input/oceanographic_raster_2023_2024.RData") 



# load MPA polygon and depth raster

MPA_square <- read_sf('input/shapefiles/MPA_square.shp')

depth_raster <- raster('input/gebco_2020_depth_raster.tif')


## load FAD data 

sample <-  read.csv("input/exemple_public_data.csv",sep =';', dec = '.', na.strings = "NULL") 
sample$date_time <- as.POSIXct(sample_2$Timestamp, tz = "GMT", format = "%d/%m/%Y %H:%M")



sample <- sample[!duplicated(sample), ]


# remove message other than basic biomass and position message
sample$MD <- as.character(sample$MD)
sample <- subset(sample, !(MD %in% c("168", "162")))

# time format
#sample$StoredTime <- str_replace_all(sample$StoredTime, "/", "-")

sample$date_time_local <- as.POSIXct(format(sample$date_time, tz = "Etc/GMT-11"), tz = "Pacific/Midway")

sample$date <- format(sample$date_time_local, "%x")
sample$time <- format(sample$date_time_local, "%H:%M")

sample$month <- month(sample$date_time_local) 
sample$year <- year(sample$date_time_local) 

sample$season <- NA

summer_data <- which(sample$month %in% c(7,8,9))
sample$season[summer_data] <- "summer"

fall_data <- which(sample$month %in% c(10,11,12))
sample$season[fall_data] <- "fall"

winter_data <- which(sample$month %in% c(1,2,3))
sample$season[winter_data] <- "winter"

spring_data <- which(sample$month %in% c(4,5,6))
sample$season[spring_data] <- "spring"

# remove pbmatic data
#sample <- sample[-which(sample$Name %in% "SLX+368312_bis" & sample$date_time < "2022-02-02 00:00:00"),]


#### position data ###

sample_position <- subset(sample, MD == "161")

# sample_position_corrected <- subset(sample_position, Speed < 3) # remove aberrant speed (FAD still onboard) => don't need it

#### remove unique FADs ####
# alone_FAD <- sample_position$Name[ave(seq_along(sample_position$Name), sample_position$Name, FUN = length) == 1]

# sample <- subset(sample, !(Name %in% alone_FAD)) #sample <- subset(sample, Name != alone_FAD)

# sample_position <- subset(sample_position, !(Name %in% alone_FAD))

#### biomass data ####

sample_biomass <- subset(sample, MD == "174")
sample_biomass$non_tuna <- sample_biomass$Layer1+sample_biomass$Layer2 # fish from 0 to 20 m, i.e. non tuna fish
sample_biomass$small_tuna <- sample_biomass$Layer3+sample_biomass$Layer4+sample_biomass$Layer5+sample_biomass$Layer6+sample_biomass$Layer7 # fish from 20 to 80 m, i.e. small tuna
sample_biomass$large_tuna <- sample_biomass$Layer8+sample_biomass$Layer9+sample_biomass$Layer10 # fish from 80 to 120 m, i.e. large tuna

sample_biomass$sum_biomass <-sample_biomass$non_tuna+sample_biomass$small_tuna+sample_biomass$large_tuna 

# position interpolation

for(b in 1:length(unique(sample_position$Name))){
  # x <- 0
  
  sample_biomass_buoy <- subset(sample_biomass, Name == unique(sample_position$Name)[b])
  sample_position_buoy <- subset(sample_position, Name == unique(sample_position$Name)[b])
  sample_position_buoy_tmp <- subset(sample_position_buoy, select=c( "Longitude","Latitude"))
  
  
  sample_biomass_buoy_tmp <- subset(sample_biomass_buoy, date_time  < sample_position_buoy$date_time[1])
  if(dim(sample_biomass_buoy_tmp)[1]!= 0){
    # x <- 1
    distance <- mean(sample_position_buoy$Speed[1])*as.numeric(difftime(sample_position_buoy$date_time[1],sample_biomass_buoy_tmp$date_time , units = "hours"))
    angle <- sample_position_buoy$Drift[1]-180
    position <- destPoint(sample_position_buoy_tmp[1,], angle, distance*1852) # knots?
    speed <- sample_position_buoy$Speed[1]
    
    tmp_init <- which(sample_biomass$Name == unique(sample_position$Name)[b] & sample_biomass$date_time  < sample_position_buoy$date_time[1])
    sample_biomass$Latitude[tmp_init] <- position[,2]
    sample_biomass$Longitude[tmp_init] <- position[,1]
    sample_biomass$Speed[tmp_init] <- speed
  }
  
  for(l in 1:length(sample_position_buoy$date_time)){
    rm(speed)
    speed <- sample_position_buoy$Speed[l]
    
    if(l == length(sample_position_buoy$date_time)){
      sample_biomass_buoy_tmp <- subset(sample_biomass_buoy, date_time  >= sample_position_buoy$date_time[l])
      if(dim(sample_biomass_buoy_tmp)[1]!= 0){
        distance <- sample_position_buoy$Speed[l]*as.numeric(difftime(sample_biomass_buoy_tmp$date_time, sample_position_buoy$date_time[l], units = "hours"))
        angle <- sample_position_buoy$Drift[l]
        
        
        position <- destPoint(sample_position_buoy_tmp[l,], angle, distance*1852) # knots?
        
        tmp_last <- which(sample_biomass$Name == unique(sample_position$Name)[b] & sample_biomass$date_time  >= sample_position_buoy$date_time[l])
        sample_biomass$Latitude[tmp_last] <- position[,2]
        sample_biomass$Longitude[tmp_last] <- position[,1]
        sample_biomass$Speed[tmp_last] <- speed
      }
    }else{
      # extract the biomass data between two position data point
      sample_biomass_buoy_tmp <- subset(sample_biomass_buoy, date_time  >= sample_position_buoy$date_time[l] & date_time  < sample_position_buoy$date_time[l+1] )
      if(dim(sample_biomass_buoy_tmp)[1] != 0){
        
        # date difference to 0-1 value between the two position
        time_to_distance <- as.numeric(difftime(sample_biomass_buoy_tmp$date_time, sample_position_buoy$date_time[l], units = "hours"))/as.numeric(difftime(sample_position_buoy$date_time[l+1],sample_position_buoy$date_time[l], unit = "hours"))
        # distance vector between the biomass data lon/lat and the position data point
        distance <-  pointDistance(sample_position_buoy_tmp[l,], sample_position_buoy_tmp[l+1,], type='Euclidean', lonlat=TRUE)
        # angle between the two consecutive position point
        angle <- bearing(sample_position_buoy_tmp[l,], sample_position_buoy_tmp[l+1,])
        # matrix of position of biomass point between point l and l+1
        position <- destPoint(sample_position_buoy_tmp[l,], angle, distance*time_to_distance)
        
        
        tmp <- which(sample_biomass$Name == unique(sample_position$Name)[b] & sample_biomass$date_time  >= sample_position_buoy$date_time[l] & sample_biomass$date_time  < sample_position_buoy$date_time[l+1])
        sample_biomass$Latitude[tmp] <- position[,2]
        sample_biomass$Longitude[tmp] <- position[,1]
        sample_biomass$Speed[tmp] <- speed
        
      }
    }
  }
}

sample_biomass <- sample_biomass %>% 
  filter(!is.na(Latitude) & !is.na(Longitude))

sample_biomass <- sample_biomass %>% 
  filter(Speed < 3)

### duplicate dFAD ###



for(f in 1:length(unique(sample_biomass$FAD_Name))){
  rm(sample_biomass_buoy)
  rm(sample_position_buoy)
  sample_biomass_buoy <- subset(sample_biomass, FAD_Name == unique(sample_biomass$FAD_Name)[f])
  sample_position_buoy <- subset(sample_position, FAD_Name == unique(sample_biomass$FAD_Name)[f])
  sample_biomass_buoy$diff_time <- 0
  
  if(dim(sample_biomass_buoy)[1]> 1){
    
    for(t in 2:length(sample_biomass_buoy$date_time_local)){
      sample_biomass_buoy$diff_time[t] <-  difftime(sample_biomass_buoy$date_time_local[t], sample_biomass_buoy$date_time_local[t-1], units = "hours" )
    }
    
    dFAD_separation <- which(sample_biomass_buoy$diff_time > 5*24)
    
    if(length(dFAD_separation) > 0){
      for(a in 1:length(dFAD_separation)){
        sample_biomass_buoy$Name[dFAD_separation[a]:length(sample_biomass_buoy$Name)] <- paste(unique(sample_biomass$FAD_Name)[f], '_',a,sep='')
        time_po <- which(sample_position_buoy$date_time_local > sample_biomass_buoy$date_time_local[dFAD_separation[a]-1])
        sample_position_buoy$Name[time_po] <- paste(unique(sample_biomass$FAD_Name)[f], '_',a,sep='')
      }
    }
  }
  
  tmp <- which(sample_biomass$FAD_Name == unique(sample_biomass$FAD_Name)[f])
  sample_biomass$Name[tmp] <- sample_biomass_buoy$Name
  
  tmp_2 <- which(sample_position$FAD_Name == unique(sample_biomass$FAD_Name)[f])
  sample_position$Name[tmp_2] <- sample_position_buoy$Name
}



### time of first position ###
sample_biomass$first_time_position <- NA

for(f in 1:length(unique(sample_biomass$Name))){
  # sample_biomass_buoy <- subset(sample_biomass, Name == unique(sample_biomass$Name)[f])
  sample_position_buoy <- subset(sample_position, Name == unique(sample_biomass$Name)[f])
  
  min_time <- min(sample_position_buoy$date_time_local)
  
  tmp <- which(sample_biomass$Name == unique(sample_biomass$Name)[f])
  sample_biomass$first_time_position[tmp] <-  as.character(min_time)
  
}



### inside/outside of MPA ###


pnts_sf <- st_as_sf(sample_biomass, coords = c('Latitude', 'Longitude'), crs = st_crs(MPA_square))

limit_square_Longitude <- c(-163.1878, -161.2008, -161.2008, -161.4228, -163.1878, -163.1878)
limit_square_Latitude <- c(7.243881, 7.243882,5.339714, 5.026103,5.026103,7.243881)

sample_biomass$presence_in_MPA <- point.in.polygon(sample_biomass$Longitude, sample_biomass$Latitude, limit_square_Longitude, limit_square_Latitude, mode.checked=FALSE)

sample_biomass$presence_in_MPA <- as.character(sample_biomass$presence_in_MPA)


### large biomass/low biomass ###

max_biomass <- NA
for(b in 1:length(unique(sample_biomass$Name))){
  
  sample_biomass_buoy <- subset(sample_biomass, Name == unique(sample_biomass$Name)[b])
  sample_biomass_buoy_24h <- subset(sample_biomass_buoy, date_time_local <= sample_biomass_buoy$date_time_local[1]+3600*36)
  
  max_biomass[b] <- max(sample_biomass_buoy_24h$sum_biomass, na.rm = T )
  
}

quantile_biomass <- quantile(max_biomass, probs = seq(0, 1, 0.25), na.rm = T)

sample_biomass$first_day_biomass <- NA
max_biomass <- NA

for(b in 1:length(unique(sample_biomass$Name))){
  
  sample_biomass_buoy <- subset(sample_biomass, Name == unique(sample_biomass$Name)[b])
  sample_biomass_buoy_24h <- subset(sample_biomass_buoy, date_time_local <= sample_biomass_buoy$date_time_local[1]+3600*36)
  
  max_biomass[b] <- max(sample_biomass_buoy_24h$sum_biomass, na.rm = T )
  
  if(max_biomass[b] < quantile_biomass[2]){
    quantile_biomass_tmp <- "quantile_1"
  }else if(max_biomass[b] >= quantile_biomass[2] && max_biomass[b] < quantile_biomass[3]){
    quantile_biomass_tmp <- "quantile_2"
  }else if(max_biomass[b] >= quantile_biomass[3] && max_biomass[b] < quantile_biomass[4]){
    quantile_biomass_tmp <- "quantile_3" 
  }else if(max_biomass[b] >= quantile_biomass[4] ){
    quantile_biomass_tmp <- "quantile_4"
  }
  
  tmp <- which(sample_biomass$Name == unique(sample_biomass$Name)[b])
  sample_biomass$first_day_biomass[tmp] <- quantile_biomass_tmp
  
}


### position initial ###
point_1 <- c(-163.185,5.025)
point_2 <- c(-161.2,7.24)

# y1 = ax1 + b
# y2 = ax2 + b
# a = (y1-b)/x1
# b = y2 - ax2
# y1 = ax1 + y2 - ax2
# 
# a = y1-y2/x1-x2

a <- (point_1[2]-point_2[2])/(point_1[1]-point_2[1])
b <- point_2[2]-a*point_2[1]


sample_biomass$first_day_position <- NA

for(c in 1:length(unique(sample_position$Name))){
  
  sample_position_buoy <- subset(sample_position, Name == unique(sample_position$Name)[c])
  tmp_position = a*sample_position_buoy$Longitude[1] +b
  if(tmp_position > sample_position_buoy$Latitude[1]){
    first_day_position_tmp <- "SE"
  }else{
    first_day_position_tmp <- "NW"
  }
  tmp <- which(sample_biomass$Name == unique(sample_position$Name)[c])
  sample_biomass$first_day_position[tmp] <- first_day_position_tmp
}


### time spend in the MPA ###
# 

sample_biomass$time_in_dataset <- NA

for(b in 1:length(unique(sample_biomass$Name))){
  rm(sample_biomass_buoy)
  sample_biomass_buoy <- subset(sample_biomass, Name == unique(sample_biomass$Name)[b])
  
  
  
  time_in_dataset_buoy <- difftime(sample_biomass_buoy$date_time_local, min(sample_biomass_buoy$date_time_local), units = "days") # sample_biomass_buoy$date_time
  # print(unique(sample_biomass$Name)[b])
  # print(max(time_in_MPA_buoy))
  
  tmp <- which(sample_biomass$Name == unique(sample_biomass$Name)[b])
  sample_biomass$time_in_dataset[tmp] <- time_in_dataset_buoy
  
}


# time_in_MPA

sample_biomass$time_in_MPA <- NA

for(b in 1:length(unique(sample_biomass$Name))){
  rm(sample_biomass_buoy)
  sample_biomass_buoy <- subset(sample_biomass, Name == unique(sample_biomass$Name)[b])
  
  
  MPA_boundaries <- min(which(sample_biomass_buoy$presence_in_MPA == "1"))
  # when no data point in the MPA?
  
  
  
  time_in_MPA_buoy <- difftime(sample_biomass_buoy$date_time, sample_biomass_buoy$date_time[MPA_boundaries], units = "days") # sample_biomass_buoy$date_time
  # print(unique(sample_biomass$Name)[b])
  # print(max(time_in_MPA_buoy))
  
  tmp <- which(sample_biomass$Name == unique(sample_biomass$Name)[b])
  sample_biomass$time_in_MPA[tmp] <- time_in_MPA_buoy
  
}

# distance to closer dFAD

distance_to_dFAD_fct <- function(Longitude_sample, Latitude_sample, Longitude_dFAD, Latitude_dFAD) {
  distance_to_dFAD <- pointDistance(c(Longitude_sample,Latitude_sample), c(Longitude_dFAD,Latitude_dFAD), type='Euclidean', lonlat=TRUE)
  return(distance_to_dFAD)
}

sample_position$distance_to_nearest_dFAD <- NA 
sample_biomass$distance_to_nearest_dFAD <- NA 

for(b in 1:length(unique(sample_biomass$Name))){
  
  sample_position_buoy <- subset(sample_position, Name == unique(sample_position$Name)[b])
  
  
  distance_min <- NA
  
  for(c in 1:length(sample_position_buoy$date_time)){
    tmp_position_sample <- subset(sample_position, Name != unique(sample_position$Name)[b]& date_time < sample_position_buoy$date_time[c]+12*60*60& date_time > sample_position_buoy$date_time[c]-12*60*60)
    if(length(tmp_position_sample$date_time) != 0){
      distance_min[c] <-min(mapply(distance_to_dFAD_fct, Longitude_sample=tmp_position_sample$Longitude, Latitude_sample= tmp_position_sample$Latitude,
                                   Longitude_dFAD= sample_position_buoy$Longitude[c],Latitude_dFAD=sample_position_buoy$Latitude[c] )/1000) # km
    }else{
      distance_min[c] <- NA
    }
  }
  tmp <- which(sample_position$Name == unique(sample_position$Name)[b])
  sample_position$distance_to_nearest_dFAD[tmp] <- distance_min
}

for(b in 1:length(unique(sample_biomass$Name))){
  sample_biomass_buoy <- subset(sample_biomass, Name == unique(sample_biomass$Name)[b])
  sample_position_buoy <- subset(sample_position, Name == unique(sample_biomass$Name)[b])
  
  min_time_diff <- NA
  
  for(c in 1:length(sample_biomass_buoy$date_time)){
    min_time_diff[c] <-  which.min(abs(difftime(sample_position_buoy$date_time, sample_biomass_buoy$date_time[c], units = "days")))
  }
  tmp2 <- which(sample_biomass$Name == unique(sample_biomass$Name)[b])
  sample_biomass$distance_to_nearest_dFAD[tmp2] <- sample_position_buoy$distance_to_nearest_dFAD[min_time_diff]
  
}



### distance to the Island ###
position_palmyra = cbind(c(-162.078333),c(5.883611))
position_washington = cbind(c(-160.377778), c(4.683333))
position_kingman = cbind(c(-162.40),c(6.41))

distance_to_island_fct <- function(Longitude, Latitude) {
  distance_to_palmyra <- pointDistance(c(Longitude,Latitude), position_palmyra, type='Euclidean', lonlat=TRUE)
  distance_to_washington <- pointDistance(c(Longitude,Latitude), position_washington, type='Euclidean', lonlat=TRUE)
  distance_to_kingman <- pointDistance(c(Longitude,Latitude), position_kingman, type='Euclidean', lonlat=TRUE)
  distance_to_land_tmp <- min(distance_to_palmyra, distance_to_washington, distance_to_kingman)/1000 # km
  return(distance_to_land_tmp)
}

distance_to_land <- mapply( distance_to_island_fct, Longitude=sample_biomass$Longitude, Latitude= sample_biomass$Latitude)
sample_biomass <- cbind(sample_biomass,distance_to_land)

# remove data too close from the Island
to_remove <- which(sample_biomass$distance_to_land < 10)
sample_biomass <- sample_biomass[-to_remove,]

### Kiribati EEZ ###
sample_biomass$Kiribati_EEZ <- 0


limit_kiribati_Longitude <- c(-159.333, -163.064, -163.064, -159.333)
limit_kiribati_Latitude <- c(7.878, 2.661,0,0)

sample_biomass$Kiribati_EEZ <- point.in.polygon(sample_biomass$Longitude, sample_biomass$Latitude, limit_kiribati_Longitude, limit_kiribati_Latitude, mode.checked=FALSE)

sample_biomass$Kiribati_EEZ <- as.character(sample_biomass$Kiribati_EEZ)



### depth ###

sample_biomass_coord <- sample_biomass
coordinates(sample_biomass_coord)= ~ Longitude + Latitude

sea_floor_depth <- raster::extract(depth_raster, sample_biomass_coord)

sample_biomass <- cbind(sample_biomass,sea_floor_depth)

### Moon phase ###


MoonIllumination <- getMoonIllumination(sample_biomass$date_time_local, keep = c("fraction"))

sample_biomass$Moon_Illumination <- MoonIllumination[,2]

# 
# sunrise_fct <- function(Longitude, Latitude, date) {
#   sunrise_tmp <- getSunlightTimes( date = date, lat = Latitude, lon = Longitude,
#                                     keep = c("nauticalDawn"), tz = "Pacific/Midway")
#   sunrise <- as.character(sunrise_tmp$nauticalDawn)
#   return(sunrise)
# }
# 
# sunrise_time <- as.POSIXct(mapply( sunrise_fct, Longitude=sample_biomass$Longitude, Latitude= sample_biomass$Latitude, 
#                                           date= sample_biomass$date), tz = "Pacific/Midway")
# 
# 
# 
# sample_biomass <- cbind(sample_biomass,sunrise_time)
# 




### dawn ###

# sample_biomass$date <- as.Date(sample_biomass$date , format = "%d/%m/%y")
# 
# sunrise_fct <- function(Longitude, Latitude, date) {
#   sunrise_tmp <- getSunlightTimes( date = date, lat = Latitude, lon = Longitude,
#                                     keep = c("nauticalDawn"), tz = "Pacific/Midway")
#   sunrise <- as.character(sunrise_tmp$nauticalDawn)
#   return(sunrise)
# }
# 
# sunrise_time <- as.POSIXct(mapply( sunrise_fct, Longitude=sample_biomass$Longitude, Latitude= sample_biomass$Latitude, 
#                                           date= sample_biomass$date), tz = "Pacific/Midway")
# 
# 
# 
# sample_biomass <- cbind(sample_biomass,sunrise_time)
# 



### Oceanographic data ###


sample_biomass$current_velocity <- NA
sample_biomass$thermocline_depth <- NA
sample_biomass$temperature <- NA
sample_biomass$sal <- NA
sample_biomass$Chl <- NA
sample_biomass$SST_front <- NA
sample_biomass$Chl_front <- NA


for(d in 1:length(times_currents)){
  
  tmp <- which(sample_biomass$date_time_local >= as.POSIXct(times_currents[d]) & sample_biomass$date_time_local< as.POSIXct(times_currents[d+1]))
  
  
  if(length(tmp) != 0 ){
    sample_biomass$current_velocity[tmp] <- raster::extract(raster_currents[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    sample_biomass$sal[tmp] <- raster::extract(raster_sal[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    sample_biomass$thermocline_depth[tmp] <-  raster::extract(raster_thermocline[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    sample_biomass$temperature[tmp] <-  raster::extract(raster_temp[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    sample_biomass$SST_front[tmp] <-  raster::extract(raster_SST_front[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
  }
}

for(d in 1:length(times_chl)){
  
  tmp_2 <- which(sample_biomass$date_time_local >= as.POSIXct(times_chl[d]) & sample_biomass$date_time_local< as.POSIXct(times_chl[d+1]))
  
  
  if(length(tmp_2) != 0 ){
    sample_biomass$Chl[tmp_2] <-  raster::extract(raster_chl[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp_2],sample_biomass$Latitude[tmp_2])), method='simple')
    
    sample_biomass$Chl_front[tmp_2] <-  raster::extract(raster_Chl_front[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp_2],sample_biomass$Latitude[tmp_2])), method='simple')
  }
}

# 2023 2024
for(d in 1:length(times_current_2024)){
  
  tmp <- which(sample_biomass$date_time_local >= as.POSIXct(times_current_2024[d]) & sample_biomass$date_time_local< as.POSIXct(times_current_2024[d+1]))
  
  
  if(length(tmp) != 0 ){
    sample_biomass$current_velocity[tmp] <- raster::extract(raster_current_2024[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    #sample_biomass$sal[tmp] <- raster::extract(raster_sal[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    sample_biomass$thermocline_depth[tmp] <-  raster::extract(raster_thermocline_2024[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    sample_biomass$temperature[tmp] <-  raster::extract(raster_SST_2024[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
    
    sample_biomass$Chl[tmp] <-  raster::extract(raster_CHL_2024[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
  
    
    # sample_biomass$SST_front[tmp] <-  raster::extract(raster_SST_front[[d]], SpatialPoints(cbind(sample_biomass$Longitude[tmp],sample_biomass$Latitude[tmp])), method='simple')
  }
}





write.table(sample_biomass, paste('input/sample_biomass_tot','.csv', sep=''),sep=';', dec = '.',row.names=FALSE)
write.table(sample_position, paste('input/sample_position','.csv', sep=''),sep=';', dec = '.',row.names=FALSE)



