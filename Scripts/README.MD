# R Scripts used in analysis

1. <ins>CADengue_Temperature_R0_Functions_Public.R </ins>   This script includes functions for temperature-dependent trait value estimations, calculating overall R0, and creating temperature rasters.
2. <ins>CADengue_CalibratingRisk_Public.R </ins> This script estimates the risk metrics and their component parts (vector presence probaility, R0(T), travel-associated case rates) for each instance of local transmission (8 unique location-months). 
3. <ins>CADengue_R0_Estimation_CA_BCM4.5_Public </ins>   This script estimates R0(T) across California, using temperature data from CA BCM 
4. <ins>CADengue_RiskLayers_Overlap_PopulationAtRisk_Public </ins>    This script outputs rasters identifying the regions 'at risk' for dengue, based on being above the ecological or eco-epidemiological 'threshold' for local transmission. These rasters are then procssed in qGIS, along with human population density estimates from CA-POP ([Depsky et al. 2022](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0270746).)

Note: to access data files used in these scripts, please see [dengue_project_release](https://github.com/lcouper/CA-Dengue/releases/tag/denge_project_release)
