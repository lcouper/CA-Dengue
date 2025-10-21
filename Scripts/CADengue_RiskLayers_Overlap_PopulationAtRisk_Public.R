####################################################################
## CODE TO ESTIMATE REGIONS AT RISK BASED ON 'THRESHOLD'
## WRITTEN BY: LISA COUPER
####################################################################

##### Load libraries, set working directory, define time periods #####

setwd("Your_path_to_CADengue_DataFiles")

library(ncdf4)
library(terra)
library(raster)
library(scales)
library(tidyr)
library(ggplot2)
library(sf)
library(exactextractr)
library(maptools)
library(RNetCDF)
library(tigris)
library(RColorBrewer)
library(data.table)
library(dplyr)

# Note, in the manuscript we estimate the population at risk across multiple time periods
# (eg current, mid-century, end-of-century) and climate warming scenarios (eg RCP 4.5, 6.0, 8.5)
# Below, we provide the code for a single time period and scenario
# But the same script can be followed for estimation under other time periods / RCPs

#### Step 1. Define risk threshold #####

# The thresholds below were generated in the script 'CADengue_CalibratingRisk_Public.R'
# and stored in the file: "pooled_threshold_summary.csv"
# Below, we have just copied the values from this file, for ease

#RiskVals <- fread("~/Dropbox/CurrentProjects/CADengue/OtherDataFiles/pooled_threshold_summary.csv")
#threshold <- RiskVals$mean[RiskVals$metric == "overall" ]
#threshold_lower <- RiskVals$lwr95[RiskVals$metric == "overall" ]
#threshold_upper <- RiskVals$upr95[RiskVals$metric == "overall" ]

# Pick one for below 

# Ecological threshold:
threshold <- 0.4148545
threshold_lower <- 0.4078462
threshold_upper <- 0.4212455

# Eco-epidemiological threshold:
threshold <- 0.00166349
threshold_lower <- 0.001643019
threshold_upper <- 0.001683373

##### Step 2. Calculate regions above risk threshold  #####

# the code below first redefines the raster based on if pixel is above or below threshold
# then draws a 1000 m buffer around all pixels with value above this threshold
# this is then exported to qGIS to estimate populaiton size in each of these buffers

# No rasters for Jan, Feb, Mar, Nov, or Dec as risk is not > threshold for any pixels
 
# Define month names and matching two-digit codes
months <- c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")

for (m in months) {

  rast_path <- paste0("EcoRisk_", m, "2010to2020.tif")
  r <- rast(rast_path)
  
  make_layer <- function(r_layer, suffix) {
    r_thresh <- ifel(r_layer < get(paste0("threshold", suffix)), NA, r_layer)
    plot(r_thresh)
    p <- as.points(r_thresh)
    b <- buffer(p, 1000)
    a <- aggregate(b)
    out_path <- paste0("ShapefileForPopAtRisk_", m, suffix)
    writeVector(a, out_path)}
  
  make_layer(r, "") # mean layer
  make_layer(r, "_lower") # lower threshold layer
  make_layer(r, "_upper")} # upper threshold layer


# Note these raster pertain to *regions* where the risk threshold is exceeded
# To get the *population size* in these regions, we processed the shapefiles
# generated above in qGIS using scripts:
# Batch_ZonalStatistics_Script & Batch_Output_PopAtRisk_Script



