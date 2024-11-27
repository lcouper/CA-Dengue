# Steps for calculating population "at risk"

#### Input data 
- CAPOP_2020_100m_TOTAL.tif (100m resolution CA population size estimates)
- Monthly ecological risk layers rasters

## Steps:

### In R:
1. Open 'CAdengue_RiskLayers_Overlap.R'
2. Use terra package to open a given monthly risk layer
3. Use if/else statement to set as "NA", any risk metric values below 0.119 (the 'cut-off', or the lowest value of ecological risk observed during actual local transmission)
4. Convert to points layer
5. Set a 1km buffer around all points
6. Aggregate buffers to create shapefile of all regions within 1 km buffer of a pixel with risk metric > 0.119
7. Export shapefile

```
library(terra)

# September
rast09<- rast("EcoRiskEstimates/Current/CurrentEcoRisk_Sep.tif")
r09 <- ifel(rast09 < 0.119, NA, rast09)
plot(r09)
p09 <- as.points(r09)
b09 <- buffer(p09, 1000) 
aa09 <- aggregate(b09)
writeVector(aa09, "CAPopulationAtRisk/ShapefileForPopAtRisk_Sep")
# used guidance from :https://stackoverflow.com/questions/66547138/rr-raster-manipulation-extract-values-from-buffer-area-without-overlap
```


### In qGIS:
1. Open 'qGIS_PopulationAtRiskEstimates.qgz' 
2. Import shapefile from step 7 above
3. Import CAPOP_2020_100m_TOTAL.tif layer
4. Use Zonal Statistics to calculate the **sum** of CA POP Within the buffers
5. Open attribute table to see value
6. To obtain CA population as a whole, use 'Raster Layer Statistics', which creates a temporarly layer showing the sum of all CA POP pixels (e.g. total population size)

<img width="499" alt="image" src="https://github.com/user-attachments/assets/22c7222d-2b2f-47e1-8700-a193bb7ed7f2">

