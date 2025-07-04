#  Oceanographic data :
#  
#  ------------------------------------------------------------------------
#   Arthur Blanluet
#   04/07/2025
#  ------------------------------------------------------------------------
#   Input
#   Copernicus data
#  ------------------------------------------------------------------------
#   Output 
#   oceanographic_raster
#  ------------------------------------------------------------------------

library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
# library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(abind)
library(grec) # front detection


path_save <- "F:/Save_Brisbane/script/GitHub/FAD_analysis_new/Input/"

#### thermocline #### 

raster_thermocline_2024 <- list()

nc_data <- nc_open("F:/Save_Brisbane/Oceanographic data/copernicus/cmems_mod_glo_phy_anfc_0.083deg_P1D-m_1723510091824.nc")

lon <- -ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")

times_thermocline_2024 <- as.character(as.POSIXct((t), origin = "1970-01-01", tz = "UTC"))

thermocline.array <- ncvar_get(nc_data, "mlotst")



for(f in 1:length(t)){

thermocline.array_t <- thermocline.array[,,f]

thermocline.array_t <- t(thermocline.array_t)

r <- raster(thermocline.array_t, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# t(ndvi.slice)

r <- flip(r, direction='y')

plot(r, main= paste(times_thermocline_2024[f],"thermocline"))

raster_thermocline_2024[[f]] <- r 

}



#### current speed #### 

raster_current_2024 <- list()

nc_data <- nc_open("F:/Save_Brisbane/Oceanographic data/copernicus/cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m_1723509447611.nc")

lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")

times_current_2024 <- as.character(as.POSIXct((t), origin = "1970-01-01", tz = "UTC"))

current.EW.array <- ncvar_get(nc_data, "uo")

current.NS.array <- ncvar_get(nc_data, "vo")

for(f in 1:length(t)){
  
  # current.EW.array_t <- abs(current.EW.array[,,f])
  # current.NS.array_t <- abs(current.NS.array[,,f])
  
  current.array_t <- sqrt(current.EW.array[,,f]^2 + current.NS.array[,,f]^2)
  
  
  #current.array_t <- pmax(current.EW.array_t, current.NS.array_t)
  
  current.array_t <- t(current.array_t)
  
  r <- raster(current.array_t, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # t(ndvi.slice)
  
  #r <- flip(r, direction='y')
  
  plot(r, main= paste(times_current_2024[f],"current max speed"))
  
  raster_current_2024[[f]] <- r 
  
}


#### CHL #### 

raster_CHL_2024 <- list()

nc_data <- nc_open("F:/Save_Brisbane/Oceanographic data/copernicus/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D_1723507406287.nc")

lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")

times_CHL_2024 <- as.character(as.POSIXct((t), origin = "1970-01-01", tz = "UTC"))

CHL.array <- ncvar_get(nc_data, "CHL")



for(f in 1:length(t)){
  
  CHL.array_t <- CHL.array[,,f]
  
  CHL.array_t <- t(CHL.array_t)
  
  r <- raster(CHL.array_t, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # t(ndvi.slice)
  
  r <- flip(r, direction='y')
  
  plot(r, main= paste(times_CHL_2024[f],"CHL"))
  
  raster_CHL_2024[[f]] <- r 
  
}


#### SST #### 

raster_SST_2024 <- list()

nc_data <- nc_open("F:/Save_Brisbane/Oceanographic data/copernicus/METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2_1723508941248.nc")

lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
t <- ncvar_get(nc_data, "time")

times_SST_2024 <- as.character(as.POSIXct((t), origin = "1970-01-01", tz = "UTC"))

SST.array <- ncvar_get(nc_data, "analysed_sst")-273.15



for(f in 1:length(t)){
  
  SST.array_t <- SST.array[,,f]
  
  SST.array_t <- t(SST.array_t)
  
  r <- raster(SST.array_t, xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  # t(ndvi.slice)
  
  r <- flip(r, direction='y')
  
  plot(r, main= paste(times_SST_2024[f],"SST"))
  
  raster_SST_2024[[f]] <- r 
  
}



save(raster_current_2024,raster_SST_2024,raster_CHL_2024,raster_thermocline_2024 ,times_current_2024, 
     times_SST_2024, times_CHL_2024, times_thermocline_2024, file=paste(path_save,"oceanographic_raster_2023_2024.RData", sep=""))

