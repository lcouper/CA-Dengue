########################################################################
## CODE TO ESTIMATE RISK AT INSTANCES OF LOCAL TRANSMISSION
## WRITTEN BY: LISA COUPER
########################################################################

##### Load libraries and set working directory #####

setwd("Your_path_to_CADengue_DataFiles")

library(ncdf4)
library(raster)
library(scales)
library(dismo)
library(dplyr)
library(stars)
library(tidyr)
library(sf)
library(maptools)
library(tigris)
library(RColorBrewer)
library(tigris)
library(exactextractr)
library(plotrix)
library(data.table)

# Data files needed:



# Instances of local transmission:
# Note months reflect suspected timing of transmission, not when the case was reported
# Pasadena, August 2023
# Long Beach, September 2023
# Baldwin Park, August 2024
# Panorama City, August 2024
# Escondido, September 2024
# El Monte, September 2024
# Hollywood Hills, September 2024
# San Bernardino, September 2024

##### Read in R0 rasters & Aedes SDMs, for relevant months in 2023/2024 #####

aug23r0 <- raster("2023-08-01_rcp45.tif")
aug23sdm <- raster("pred_aug2023.tif")
sep23r0 <- raster("2023-09-01_rcp45.tif")
sep23sdm <- raster("pred_sep2023.tif")
aug24r0 <- raster("2024-08-01_rcp45.tif")
aug24sdm <- raster("pred_aug2024.tif")
sep24r0 <- raster("2024-09-01_rcp45.tif")
sep24sdm <- raster("pred_sep2024.tif")
oct24r0 <- raster("2024-10-01_rcp45.tif")
oct24sdm <- raster("pred_sep2024.tif")
augimp <- raster("ImpAug_CensusTract_annual.tif")
sepimp <- raster("ImpSep_CensusTract_annual.tif")
augimpCounty <- raster("CasesCounty_Aug.tiff")
sepimpCounty <- raster("CasesCounty_Sep.tiff")

# transform CRS 
augimp <- projectRaster(augimp, crs = crs(aug23sdm))
sepimp <- projectRaster(sepimp, crs = crs(aug23sdm))
augimpCounty <- projectRaster(augimpCounty, crs = crs(aug23sdm))
sepimpCounty <- projectRaster(sepimpCounty, crs = crs(aug23sdm))

# ensure rasters have the same extent:
aug23sdm <- projectRaster(aug23sdm, to = aug23r0, method = "bilinear")
sep23sdm <- projectRaster(sep23sdm, to = sep23r0, method = "bilinear")
aug24sdm <- projectRaster(aug24sdm, to = aug24r0, method = "bilinear")
sep24sdm <- projectRaster(sep24sdm, to = sep24r0, method = "bilinear")
oct24sdm <- projectRaster(oct24sdm, to = oct24r0, method = "bilinear")
augimp <- projectRaster(augimp, to = aug23r0, method = "ngb")  
augimpCounty  <- projectRaster(augimpCounty, to = aug23r0, method = "ngb")  
sepimp <- projectRaster(sepimp, to = sep23r0, method = "ngb")  
sepimpCounty  <- projectRaster(sepimpCounty, to = sep23r0, method = "ngb")  

# Bring in CA shapefile and transform CRS
CA <- tigris::states() %>% subset(NAME == "California")
CA <- st_transform(CA, crs = st_crs(aug23sdm)$proj4string)

##### Pull GEOIDs for each location #########

options(tigris_use_cache = TRUE)
options(tigris_class = "sf")

# 1) Load California places and tracts (2020)
places_ca <- places("CA", year = 2020) %>% st_make_valid()
tracts_ca <- tracts(state = "CA", year = 2020, cb = TRUE) %>% st_make_valid()
tracts_ca <- st_transform(tracts_ca, st_crs(places_ca))

# 2) Helpers
get_place_geoids <- function(place_name) {
  poly <- places_ca %>% filter(NAME == place_name)
  if (nrow(poly) == 0) return(character(0))
  hits <- st_filter(tracts_ca, poly, .predicate = st_intersects)
  sort(unique(hits$GEOID))}

get_zcta_geoids <- function(zcta_codes) {
  z <- zctas(year = 2020) %>% st_make_valid()
  z <- z %>% filter(ZCTA5CE20 %in% zcta_codes) %>% st_transform(st_crs(tracts_ca))
  if (nrow(z) == 0) return(character(0))
  hits <- st_filter(tracts_ca, z, .predicate = st_intersects)
  sort(unique(hits$GEOID))}

# 3) Pull GEOIDs (separate vectors)

# Places (direct)
pasadena_ids      <- get_place_geoids("Pasadena")
longbeach_ids     <- get_place_geoids("Long Beach")
baldwinpark_ids   <- get_place_geoids("Baldwin Park")
escondido_ids     <- get_place_geoids("Escondido")
elmonte_ids       <- get_place_geoids("El Monte")
sanbernardino_ids <- get_place_geoids("San Bernardino")

# Panorama City (not a Census "place") → ZCTA proxy (main ZIP)
panoramacity_ids <- get_place_geoids("Panorama City")
if (length(panoramacity_ids) == 0) {
  panoramacity_ids <- get_zcta_geoids(c("91402"))  # primary ZIP for Panorama City
}

# Hollywood Hills (not a Census "place") → ZCTA proxy (union)
# 90068 covers most of Hollywood Hills; 90046 adds Hollywood Hills West/Laurel Canyon.
hollywoodhills_ids <- get_place_geoids("Hollywood Hills")
if (length(hollywoodhills_ids) == 0) {
  hollywoodhills_ids <- get_zcta_geoids(c("90068","90046"))}

# Counts in each
sapply(list(
  Pasadena = pasadena_ids,
  Long_Beach = longbeach_ids,
  Baldwin_Park = baldwinpark_ids,
  Panorama_City = panoramacity_ids,
  Escondido = escondido_ids,
  El_Monte = elmonte_ids,
  Hollywood_Hills = hollywoodhills_ids,
  San_Bernardino = sanbernardino_ids), length)

##### Estimate risk at times of observed local transmission ####
# using a boostrapping approach, drawing samples of risk estimates 
# from each census tract / month where local transmission was observed
# weighting the random draws by population density

catracts <- st_read("tl_2020_06_tract.shp")
catracts <- st_transform(catracts, crs(aug23sdm))

RiskVals <- data.frame(matrix(NA, nrow = 8, ncol = 13))
colnames(RiskVals) <- c("Location", "R0", "SDM", "Imp", "Eco_mean", "Eco_lower", "Eco_upper", 
                        "Overall_mean", "Overall_lower", "Overall_upper", "Overall_County_mean", "Overall_County_lower", "Overall_County_upper")
RiskVals$Location <- c("Pasadena", "Long Beach", "Baldwin Park", "Panorama City", 
                       "Escondido", "El Monte", "Hollywood Hills", "San Bernardino")

##### 1. Aug 2023; Pasadena #########

Pasadena <- catracts[catracts$GEOID %in% pasadena_ids,]
aug23eco <- aug23sdm * aug23r0
aug23overall <- aug23sdm * aug23r0 * augimp
aug23overallCounty <- aug23sdm * aug23r0 * augimpCounty

######### Bootstrap approach for risk estimation #########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(aug23r0, aug23sdm, augimp, augimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, Pasadena, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_pas       <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_pas           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_pas <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[1, 2] <- est_r0
RiskVals[1, 3] <- est_sdm
RiskVals[1, 4] <- est_imp

RiskVals[1, 5] <- est_eco
RiskVals[1, 8] <- est_overall
RiskVals[1,11] <- est_overallCounty

# 95% CIs from the population-weighted bootstrap
RiskVals[1, 6] <- quantile(boot_mean_eco_pas,           0.025, na.rm = TRUE)
RiskVals[1, 7] <- quantile(boot_mean_eco_pas,           0.975, na.rm = TRUE)

RiskVals[1, 9] <- quantile(boot_mean_overall_pas,       0.025, na.rm = TRUE)
RiskVals[1,10] <- quantile(boot_mean_overall_pas,       0.975, na.rm = TRUE)

RiskVals[1,12] <- quantile(boot_mean_overallCounty_pas, 0.025, na.rm = TRUE)
RiskVals[1,13] <- quantile(boot_mean_overallCounty_pas, 0.975, na.rm = TRUE)


##### 2. Sep 2023; Long Beach ########

LB <- catracts[catracts$GEOID %in% longbeach_ids,]
sep23eco <- sep23sdm * sep23r0
sep23overall <- sep23sdm * sep23r0 * sepimp
sep23overallCounty <- sep23sdm * sep23r0 * sepimpCounty

######## Bootstrap approach for risk estimation ########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(sep23r0, sep23sdm, sepimp, sepimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, LB, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_lb      <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_lb           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_lb <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[2, 2] <- est_r0
RiskVals[2, 3] <- est_sdm
RiskVals[2, 4] <- est_imp

RiskVals[2, 5] <- est_eco
RiskVals[2, 8] <- est_overall
RiskVals[2,11] <- est_overallCounty

RiskVals[2, 6] <- quantile(boot_mean_eco_lb,           0.025, na.rm = TRUE)
RiskVals[2, 7] <- quantile(boot_mean_eco_lb,           0.975, na.rm = TRUE)

RiskVals[2, 9] <- quantile(boot_mean_overall_lb,       0.025, na.rm = TRUE)
RiskVals[2,10] <- quantile(boot_mean_overall_lb,       0.975, na.rm = TRUE)

RiskVals[2,12] <- quantile(boot_mean_overallCounty_lb, 0.025, na.rm = TRUE)
RiskVals[2,13] <- quantile(boot_mean_overallCounty_lb, 0.975, na.rm = TRUE)


##### 3. Aug 2024; Baldwin Park ######

BP <- catracts[catracts$GEOID %in% baldwinpark_ids,]
aug24eco <- aug24sdm * aug24r0
aug24overall <- aug24sdm * aug24r0 * augimp
aug24overallCounty <- aug24sdm * aug24r0 * augimpCounty

######### Bootstrap approach for risk estimation #########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(aug24r0, aug24sdm, augimp, augimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, BP, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_bp      <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_bp           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_bp <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[3, 2] <- est_r0
RiskVals[3, 3] <- est_sdm
RiskVals[3, 4] <- est_imp

RiskVals[3, 5] <- est_eco
RiskVals[3, 8] <- est_overall
RiskVals[3,11] <- est_overallCounty
RiskVals[3, 6] <- quantile(boot_mean_eco_bp,           0.025, na.rm = TRUE)
RiskVals[3, 7] <- quantile(boot_mean_eco_bp,           0.975, na.rm = TRUE)
RiskVals[3, 9] <- quantile(boot_mean_overall_bp,       0.025, na.rm = TRUE)
RiskVals[3,10] <- quantile(boot_mean_overall_bp,       0.975, na.rm = TRUE)
RiskVals[3,12] <- quantile(boot_mean_overallCounty_bp, 0.025, na.rm = TRUE)
RiskVals[3,13] <- quantile(boot_mean_overallCounty_bp, 0.975, na.rm = TRUE)


##### 4. Aug 2024; Panorama City ######

PC <- catracts[catracts$GEOID %in% panoramacity_ids,]
aug24eco <- aug24sdm * aug24r0
aug24overall <- aug24sdm * aug24r0 * augimp
aug24overallCounty <- aug24sdm * aug24r0 * augimpCounty

######### Bootstrap approach for risk estimation #########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(aug24r0, aug24sdm, augimp, augimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, PC, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_pc      <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_pc           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_pc <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[4, 2] <- est_r0
RiskVals[4, 3] <- est_sdm
RiskVals[4, 4] <- est_imp

RiskVals[4, 5] <- est_eco
RiskVals[4, 8] <- est_overall
RiskVals[4,11] <- est_overallCounty
RiskVals[4, 6] <- quantile(boot_mean_eco_pc,           0.025, na.rm = TRUE)
RiskVals[4, 7] <- quantile(boot_mean_eco_pc,           0.975, na.rm = TRUE)
RiskVals[4, 9] <- quantile(boot_mean_overall_pc,       0.025, na.rm = TRUE)
RiskVals[4,10] <- quantile(boot_mean_overall_pc,       0.975, na.rm = TRUE)
RiskVals[4,12] <- quantile(boot_mean_overallCounty_pc, 0.025, na.rm = TRUE)
RiskVals[4,13] <- quantile(boot_mean_overallCounty_pc, 0.975, na.rm = TRUE)


##### 5. Sep 2024; Escondido  ########

Esc <- catracts[catracts$GEOID %in% escondido_ids,]
sep24eco <- sep24sdm * sep24r0
sep24overall <- sep24sdm * sep24r0 * sepimp
sep24overallCounty <- sep24sdm * sep24r0 * sepimpCounty

######### Bootstrap approach for risk estimation #########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(sep24r0, sep24sdm, sepimp, sepimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, Esc, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_esc      <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_esc           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_esc <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[5, 2] <- est_r0
RiskVals[5, 3] <- est_sdm
RiskVals[5, 4] <- est_imp

RiskVals[5, 5] <- est_eco
RiskVals[5, 8] <- est_overall
RiskVals[5,11] <- est_overallCounty
RiskVals[5, 6] <- quantile(boot_mean_eco_esc,           0.025, na.rm = TRUE)
RiskVals[5, 7] <- quantile(boot_mean_eco_esc,           0.975, na.rm = TRUE)
RiskVals[5, 9] <- quantile(boot_mean_overall_esc,       0.025, na.rm = TRUE)
RiskVals[5,10] <- quantile(boot_mean_overall_esc,       0.975, na.rm = TRUE)
RiskVals[5,12] <- quantile(boot_mean_overallCounty_esc, 0.025, na.rm = TRUE)
RiskVals[5,13] <- quantile(boot_mean_overallCounty_esc, 0.975, na.rm = TRUE)



##### 6. Sep 2024; El Monte ########

ElMonte <- catracts[catracts$GEOID %in% elmonte_ids,]
sep24eco <- sep24sdm * sep24r0
sep24overall <- sep24sdm * sep24r0 * sepimp
sep24overallCounty <- sep24sdm * sep24r0 * sepimpCounty

######### Bootstrap approach for risk estimation #########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(sep24r0, sep24sdm, sepimp, sepimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, ElMonte, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_elm      <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_elm           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_elm <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[6, 2] <- est_r0
RiskVals[6, 3] <- est_sdm
RiskVals[6, 4] <- est_imp

RiskVals[6, 5] <- est_eco
RiskVals[6, 8] <- est_overall
RiskVals[6,11] <- est_overallCounty
RiskVals[6, 6] <- quantile(boot_mean_eco_elm,           0.025, na.rm = TRUE)
RiskVals[6, 7] <- quantile(boot_mean_eco_elm,           0.975, na.rm = TRUE)
RiskVals[6, 9] <- quantile(boot_mean_overall_elm,       0.025, na.rm = TRUE)
RiskVals[6,10] <- quantile(boot_mean_overall_elm,       0.975, na.rm = TRUE)
RiskVals[6,12] <- quantile(boot_mean_overallCounty_elm, 0.025, na.rm = TRUE)
RiskVals[6,13] <- quantile(boot_mean_overallCounty_elm, 0.975, na.rm = TRUE)


##### 7. Sep 2024; Hollywood Hills ########

HH <- catracts[catracts$GEOID %in% hollywoodhills_ids,]
sep24eco <- sep24sdm * sep24r0
sep24overall <- sep24sdm * sep24r0 * sepimp
sep24overallCounty <- sep24sdm * sep24r0 * sepimpCounty

######### Bootstrap approach for risk estimation #########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(sep24r0, sep24sdm, sepimp, sepimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, HH, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_hh     <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_hh           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_hh <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[7, 2] <- est_r0
RiskVals[7, 3] <- est_sdm
RiskVals[7, 4] <- est_imp

RiskVals[7, 5] <- est_eco
RiskVals[7, 8] <- est_overall
RiskVals[7,11] <- est_overallCounty
RiskVals[7, 6] <- quantile(boot_mean_eco_hh,           0.025, na.rm = TRUE)
RiskVals[7, 7] <- quantile(boot_mean_eco_hh,           0.975, na.rm = TRUE)
RiskVals[7, 9] <- quantile(boot_mean_overall_hh,       0.025, na.rm = TRUE)
RiskVals[7,10] <- quantile(boot_mean_overall_hh,       0.975, na.rm = TRUE)
RiskVals[7,12] <- quantile(boot_mean_overallCounty_hh, 0.025, na.rm = TRUE)
RiskVals[7,13] <- quantile(boot_mean_overallCounty_hh, 0.975, na.rm = TRUE)


##### 8. Sep 2024; San Bernardino ########

SBern <- catracts[catracts$GEOID %in% sanbernardino_ids,]
sep24eco <- sep24sdm * sep24r0
sep24overall <- sep24sdm * sep24r0 * sepimp
sep24overallCounty <- sep24sdm * sep24r0 * sepimpCounty

######### Bootstrap approach for risk estimation #########

# Stack rasters and extract pixel values inside selected tracts
# Note we have already ensured rasters are of same extent
rs <- stack(sep24r0, sep24sdm, sepimp, sepimpCounty, pop); names(rs) <- c("r0","sdm","imp", "impCounty", "pop")
vals <- raster::extract(rs, SBern, df = TRUE)[,-1] %>%
  as_tibble() %>%
  drop_na(r0, sdm, imp, impCounty, pop) %>%
  filter(pop > 0)  # ignore zero-pop pixels for weighting

# Compute per-cell products
vals <- dplyr::mutate(vals, eco = r0 * sdm, overall = r0 * sdm * imp, overallCounty = r0 * sdm * impCounty)

# Population-weighted estimates
w <- vals$pop
wmean <- function(x, w) sum(w * x) / sum(w)

est_r0            <- wmean(vals$r0, w)
est_sdm           <- wmean(vals$sdm, w)
est_imp           <- wmean(vals$imp, w)
est_eco           <- wmean(vals$eco, w)
est_overall       <- wmean(vals$overall, w)
est_overallCounty <- wmean(vals$overallCounty, w)

# Population-weighted bootstrap of the mean
set.seed(1)
B <- 1000
n <- nrow(vals)
wp <- w / sum(w)  # sampling probabilities ~ population

boot_mean_overall_sb     <- replicate(B, mean(vals$overall[      sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_eco_sb           <- replicate(B, mean(vals$eco[          sample.int(n, n, replace=TRUE, prob=wp)]))
boot_mean_overallCounty_sb <- replicate(B, mean(vals$overallCounty[sample.int(n, n, replace=TRUE, prob=wp)]))

RiskVals[8, 2] <- est_r0
RiskVals[8, 3] <- est_sdm
RiskVals[8, 4] <- est_imp

RiskVals[8, 5] <- est_eco
RiskVals[8, 8] <- est_overall
RiskVals[8,11] <- est_overallCounty
RiskVals[8, 6] <- quantile(boot_mean_eco_sb,           0.025, na.rm = TRUE)
RiskVals[8, 7] <- quantile(boot_mean_eco_sb,           0.975, na.rm = TRUE)
RiskVals[8, 9] <- quantile(boot_mean_overall_sb,       0.025, na.rm = TRUE)
RiskVals[8,10] <- quantile(boot_mean_overall_sb,       0.975, na.rm = TRUE)
RiskVals[8,12] <- quantile(boot_mean_overallCounty_sb, 0.025, na.rm = TRUE)
RiskVals[8,13] <- quantile(boot_mean_overallCounty_sb, 0.975, na.rm = TRUE)




#### Pool bootstraps across locations/months ######

## Collect the boostrap vectors from each location-month-year

boots_overall <- list(
  Pasadena = boot_mean_overall_pas,
  LongBeach = boot_mean_overall_lb,
  BP = boot_mean_overall_bp,
  PC = boot_mean_overall_pc,
  Esc = boot_mean_overall_esc,
  ElMonte = boot_mean_overall_elm,
  HH    = boot_mean_overall_hh,
  SBern = boot_mean_overall_sb)

boots_eco <- list(
  Pasadena = boot_mean_eco_pas,
  LongBeach = boot_mean_eco_lb,
  BP = boot_mean_eco_bp,
  PC = boot_mean_eco_pc,
  Esc = boot_mean_eco_esc,
  ElMonte = boot_mean_eco_elm,
  HH    = boot_mean_eco_hh,
  SBern = boot_mean_eco_sb)

boots_overallCounty <- list(
  Pasadena = boot_mean_overallCounty_pas,
  LongBeach = boot_mean_overallCounty_lb,
  BP = boot_mean_overallCounty_bp,
  PC = boot_mean_overallCounty_pc,
  Esc = boot_mean_overallCounty_esc,
  ElMonte = boot_mean_overallCounty_elm,
  HH    = boot_mean_overallCounty_hh,
  SBern = boot_mean_overallCounty_sb)

# Pool
pooled_boot <- function(boot_list, Bpool = 5000, seed = 1) {
  set.seed(seed)
  pooled <- replicate(Bpool, {
    draws <- vapply(boot_list, function(v) sample(v, 1L), numeric(1))
    mean(draws)  })
  point <- mean(vapply(boot_list, mean, numeric(1)))  # average of instance means
  ci    <- quantile(pooled, c(.025, .975))
  list(point = point, lwr95 = ci[1], upr95 = ci[2], draws = pooled)}

overall_pooled        <- pooled_boot(boots_overall)
eco_pooled            <- pooled_boot(boots_eco)
overallCounty_pooled  <- pooled_boot(boots_overallCounty)


######## output #######

# 1. Values from each location-month-year
fwrite(RiskVals, "AllLocations_AllRiskValues_With95CI.csv")

# 2. Summary of pooled values
pooled_summary <- tibble(
  metric = c("eco", "overall", "overallCounty"),
  mean   = c(eco_pooled$point,        overall_pooled$point,        overallCounty_pooled$point),
  lwr95  = c(as.numeric(eco_pooled$lwr95),  as.numeric(overall_pooled$lwr95),  as.numeric(overallCounty_pooled$lwr95)),
  upr95  = c(as.numeric(eco_pooled$upr95),  as.numeric(overall_pooled$upr95),  as.numeric(overallCounty_pooled$upr95)),
  Bpool  = c(length(eco_pooled$draws), length(overall_pooled$draws), length(overallCounty_pooled$draws)))
fwrite(pooled_summary, "~/Dropbox/CurrentProjects/CADengue/OtherDataFiles/pooled_threshold_summary.csv")

# 3. Full bootstrap vectors from each location
fwrite(as.data.frame(boots_eco), "BootstrapsByLocation_EcoRisk.csv")
fwrite(as.data.frame(boots_overall), "BootstrapsByLocation_OverallRisk.csv")
fwrite(as.data.frame(boots_overallCounty), "BootstrapsByLocation_OverallCountyRisk.csv")

