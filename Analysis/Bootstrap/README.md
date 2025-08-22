# Boostrap approach for calibrating risk values

We used a bootstrap approach to draw samples from tract-months in which local transmission was observed. We drew 1,000 samples from each tract-month, with probabilities weighted by population density within that tract based on CA POP estimates (100m resolution). 
First, we used qGIS to select from the CA POP raster only the census tracts with local transmission (as otherwise this raster is very large and slow to work with). We did so using the tract designation vector layer from the U.S. Census Bureau (cb_2018_06_tract_500k), and filtering for specific tracts based on 'NAME', e.g.: 

<img width="454" height="166" alt="image" src="https://github.com/user-attachments/assets/68d06853-354c-4a95-b363-9c4a1fe6fc80" />

We then clipped the CA POP raster to this filtered vector layer (using raster extraction), and ouput the resulting raster. Subsequent bootstrap steps were conducted in R in the file 'CADengue_CalibratingR0.R'




