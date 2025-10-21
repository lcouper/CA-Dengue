####################################################################
## CODE TO ESTIMATE TEMPERATURE-BASED DENGUE TRANSMISSION RISK IN CA
## USING CA-BCM MONTHLY TEMPERATURE PROJECTIONS FOR RCP 4.5
## WRITTEN BY: LISA COUPER
###################################################################

##### Load libraries, set working directory, source trait functions #####

setwd("Your_path_to_CADengue_DataFiles")

library(ncdf4)
library(raster)
library(scales)
library(tidyr)
library(sf)
library(maptools)
library(RNetCDF)
library(tigris)
library(RColorBrewer)
library(tigris)
library(data.table)

source("CADengue_Temperature_R0_Functions_Public.R")

##### Load CA-BCM data and estimate R0(T) ######

# Note, the CA-BCM data is publicly available here: https://ca.water.usgs.gov/projects/reg_hydro/basin-characterization-model.html
# We have upload a single set (for 2010), as an example

mtmaxdata = nc_open("CA_BCM_MIROC_rcp45_Monthly_tmx_2010.nc")
mtmindata = nc_open("CA_BCM_MIROC_rcp45_Monthly_tmn_2010.nc")

latm <- ncvar_get(mtmaxdata, "x")
lonm <- ncvar_get(mtmaxdata, "y")
time <- ncvar_get(mtmaxdata, "time")

# Note there is some weird indexing going on with the native files
# Double check you are working with the correct month
timeDate <- as.Date(time, origin = '2006-10-01')
#timeDesired <- "2023-09-01"
#timeIndex = which(timeDate == timeDesired)

indices <- c("2010-01-01" ,"2010-02-01", "2010-03-01", "2010-04-01", "2010-05-01", "2010-06-01", "2010-07-01", "2010-08-01", "2010-09-01")
#indices <- c("2010-10-01", "2012-11-01", "2012-12-01")

for (i in 1:length(indices)) {
  timeIndex = which(timeDate == indices[i])
  
  dfmax <- ncvar_get(nc = mtmaxdata, varid = "tmx", start = c(1,1,timeIndex), count = c(3486, 4477, 1))
  dfmin <- ncvar_get(nc = mtmindata, varid = "tmn", start = c(1,1,timeIndex), count = c(3486, 4477, 1))
  
  # Take average of monthly max and min (element-wise average)
  l <- list(dfmax, dfmin)
  dfarr <-array( unlist(l) , c(3486,4477,2) )
  dfavg <- apply( dfarr , 1:2 , mean )
  
  ##### Calculate R0 for the given month, rasterize, plot, and export ####
  
  dfR0 = R0multi(dfavg)
  # Convert NAs to 0
  dfR0[is.na(dfR0)] <- 0

  # Rescale matrix (based on max value across all months, years, scenarios)
  mn = 0; mx = 13.59043;
  dfR0r = (dfR0 - mn) / (13.59043 - mn)
  
  # convert df to raster 
  dfR0raster = RasterizeBCM(dfR0r)
  
  # Crop R0 raster to CA
  dfCrop <-crop(dfR0raster, CA, snap = 'near')
  dfmask <- mask(dfCrop, CA) # mask out other regions
  
  breakpoints <- seq(0, 1, 0.05)
  col2 <- rev(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd"))
  colors <- c("gray96", colorRampPalette(col2)(21))
  
  plot(dfmask, box = FALSE, axes = FALSE, breaks = breakpoints, col = colors, legend = F)
  plot(CA$geometry, add = TRUE)
  
  # Export raster if desired
  #writeRaster(dfmask, "2024-08-01.tif", format = "GTIFF")
}



