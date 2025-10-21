#################################################################
## FUNCTIONS FOR TEMPERATURE PROCESSING AND R0(T) ESTIMATION
## WRITTEN BY: LISA COUPER
## LAST EDITED: CURRENT
#################################################################

# Functions included in this script:

# 1. Converting K to C
# 2. Biting rate (T)
# 3. Probability of pathogen tranmission (T)
# 4. Probability of pathogen acquisition (T)
# 5. Adult mosquito life span (T)
# 6. Parasite development rate (T)
# 7. Eggs laid per female per day (T)
# 8. Probably egg to adult survival (T)
# 9. Mosquito egg-to-adult development rate (T)
# 10. Full R0(T) 
# 11. R0 model for matrix/dataframe input
# 12. Rasterize R0 or temperature for plotting (based on BCM data inputs)

# Note these thermal response curves were developed by Mordecai et al. 2017:
# https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005568 

#### 1. Convert Kelvin to Celsius ####
kelvin_to_celsius <- function(kelvin) {return(kelvin - 273.15)}

#### 2. Biting rate, a #####
# briere function. 
# units = 1/days
# For aedes aegypti, c = 0.000202, Tmin = 13.35, Tmax = 40.08
BitingRate = function(Temp)
{if (is.na(Temp)) {a <- NA}
 else if (Temp > 13.35 & Temp < 40.08) 
{a = 0.000202 * Temp * (Temp - 13.35) * sqrt((40.08 - Temp))}
  else {a= 0} # if temp entered is below Tmin or above Tmax, return 0
  return(a)} 

#### 3. Probability an infected mosquito becomes infectious ####
# b, Briere function
# For aedes aegytpi & DENV, c = 0.000849, Tmin = 17.05, Tmax = 35.83
ProbInfB = function(Temp)
{if (is.na(Temp)) {b <- NA}
  else if (Temp > 17.05 & Temp < 35.83)
{b = 0.000849 * Temp * (Temp - 17.05) * sqrt((35.83 - Temp))}
  else {b = 0}
  return(b)}

#### 4. Probability of getting infected from infectious blood meal ####
# c, Briere function
# For aedes aegypti & DENV, c = 0.000491, Tmin = 12.22, Tmax = 37.46
ProbInfC = function(Temp)
{if (is.na(Temp)) {c <- NA}
else if (Temp > 12.22 & Temp < 37.46)
{c = 0.000491 * Temp * (Temp - 12.22) * sqrt((37.46 - Temp))}
  else {c = 0}
  return(c)}

#### 5. Adult mosquito life span ####
# units: days, 1/mu, Quadratic
# For aedes aegypti, c = -0.148, Tmin = 9.16, Tmax = 37.73
Lifespan = function(Temp)
{if (is.na(Temp)) {lf <- NA}
  else if (Temp > 9.16 & Temp < 37.73)
{lf = (-0.148) * (Temp - 9.16) * (Temp - 37.73)}
  else {lf = 0}
  return(lf)}

#### 6. Parasite development rate (PDR) ####
# units: 1/days, Briere
# For DENV, c = 0.0000665, Tmin = 10.68, Tmax = 45.90
PDR = function(Temp)
{if (is.na(Temp)) {pdr <- NA}
else if (Temp > 10.68 & Temp < 45.90)
{pdr = 0.0000665 * Temp * (Temp - 10.68) * sqrt((45.90 - Temp))}
  else {pdr = 0}
  return(pdr)}

#### 7. Eggs laid per female per gonotrophic cycle, TFD ####
# Units :#/female, Briere
# For Aedes aegypti, c = 0.00856, Tmin = 14.58, Tmax = 34.61
TFD = function(Temp)
{if (is.na(Temp)) {tfd <- NA}
 else if (Temp > 14.58 & Temp < 34.61)
{tfd = 0.00856 * Temp * (Temp - 14.58) * sqrt((34.61 - Temp))}
  else {tfd = 0}
  return(tfd)}

#### 8. Probability egg-to-adult survival ####
# pEA, quadratic
# for Aedes aegypti, c = -0.00599, Tmin = 13.56, Tmax = 38.29
pEA = function(Temp)
{if (is.na(Temp)) {pea <- NA}
else if (Temp > 13.56 & Temp < 38.29)
{pea = (-0.00599) * (Temp - 13.56) * (Temp - 38.29)}
  else {pea = 0}
  return(pea)}

#### 9. Mosquito egg-to-adult development rate ####
# units: (1/days), MDR, Briere
# for Aedes aegypti, c = 0.0000786, Tmin = 11.36, Tmax = 39.17
MDR = function(Temp)
{if (is.na(Temp)) {mdr <- NA}
else if (Temp > 11.36 & Temp < 39.17)
{mdr = 0.0000786 * Temp * (Temp - 11.36) * sqrt((39.17 - Temp))}
  else {mdr = 0}
  return(mdr)}

#### 10. R0 model #####
R0 = function(Temp)
{if (Temp < 17.05)
{R0 = NA
return(R0)}
  else
  {R0 = sqrt((BitingRate(Temp)^2 * ProbInfB(Temp) * ProbInfC(Temp) * exp(1)^-(1/(Lifespan(Temp) * PDR(Temp))) * 
                (TFD(Temp)*BitingRate(Temp)) * pEA(Temp) * MDR(Temp))/ (1/Lifespan(Temp))^3)
  return(R0)}}

#### 11. R0 model for matrix/dataframe input #####
R0multi = function(z)  { # here 'Temp' is a matrix (i.e., temperature for single day across CA lat longs)
  R0mat = matrix(nrow = dim(z)[1], ncol = dim(z)[2])
    for (i in 1:dim(z)[1]) {
      for (j in 1:dim(z)[2]) {
  Temp = z[i,j]
  if (is.na(Temp)) {R0mat[i,j] <-NA}
  else if (Temp < 17.05) {R0mat[i,j] <- NA}
  else {R0mat[i,j] <- sqrt((BitingRate(Temp)^2 * ProbInfB(Temp) * ProbInfC(Temp) * exp(1)^-(1/(Lifespan(Temp) * PDR(Temp))) * 
                (TFD(Temp)*BitingRate(Temp)) * pEA(Temp) * MDR(Temp))/ (1/Lifespan(Temp))^3)}}}
  # Convert NAs to 0
  R0mat[is.na(R0mat)] <- 0
  return(R0mat)
  }


#### 12. Rasterize R0 or temperature for plotting (based on BCM data inputs) #####

RasterizeBCM = function(df) { # Here df is a matrix
  df = as.data.frame(df)
  colnames(df) = lonm
  rownames(df) = latm
  df$lon = rownames(df)
  df = df[,c(ncol(df),1:(ncol(df)-1))] # re-order so lon is listed first
  df_R0 <-gather(df, name, value, -c(lon))
  colnames(df_R0) = c("lat", "lon", "temp")
  # define projection used in CA BCM data
  dfR0_raster <- raster::rasterFromXYZ(df_R0, crs = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +no_defs +a=6378137") 
  return(dfR0_raster)
    }

# helpful notes on understanding crs syntax here: https://mgimond.github.io/Spatial/coordinate-systems-in-r.html
# info about the CA BCM crs here: https://ca.water.usgs.gov/projects/reg_hydro/projects/basin-characterization-model-v8/basin-characterization-model-v8-metadata.html#org


