# Boostrap approach for calibrating risk values

We used a bootstrap approach to draw samples from tract-months in which local transmission was observed. We drew 1,000 samples from each tract-month, with probabilities weighted by population density within that tract based on CA POP estimates (100m resolution).   

First, we used qGIS to select from the CA POP raster only the census tracts with local transmission (as otherwise this raster is very large and slow to work with). We did so using the tract designation vector layer from the U.S. Census Bureau (cb_2018_06_tract_500k), and filtering for specific tracts based on 'NAME', e.g.: 

<img width="442" height="170" alt="image" src="https://github.com/user-attachments/assets/5cfa4de0-cb74-495d-bb60-78fe118e8779" />

We then clipped the CA POP raster to this filtered vector layer (using raster extraction), and ouput the resulting raster. Subsequent bootstrap steps were conducted in R in the file 'CADengue_CalibratingR0.R'




