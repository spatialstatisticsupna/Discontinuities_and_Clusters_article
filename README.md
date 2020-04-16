# Estimating disease risks in the presence of local discontinuities and clusters
This repository contains the R code to fit with INLA the spatial and spatio-temporal models described in the work entitled _"A two-stage approach to estimate spatial and spatio-temporal disease risks in the presence of local discontinuities and clusters"_ (Adin et al., 2019).


## Table of contents

- [Data](#Data)
- [R code](#R-code)
- [References](#References)


# Data
Two motivating applications are described in this work: stomach cancer mortality data in Spanish provinces during the year 2013, and brain cancer incidence data in 27 administrative regions of Navarre and
the Basque Country during the period 2000-2008.

- [**StomachCancer_ESP.Rdata**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/data/StomachCancer_ESP.Rdata)
  
  This .Rdata contains the following objects:
	- **_Data_**: `data.frame` object with the number of observed and expected cases (_'obs'_ and _'exp'_ variables, respectively) and standardized mortality ratio (_'SMR'_) for each province (_'prov'_) and time period (_'year'_) for stomach cancer mortality data.
	- **_Carto.ESP_**: `SpatialPolygonsDataFrame` object with the cartography of the 47 continental Spanish provinces.
	- **_W_**: Spatial adjacency matrix of the Spanish provinces.


- [**BrainCancer_MUN.Rdata**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/data/BrainCancer_MUN.Rdata)
  
  This .Rdata contains the following objects:
	- **_Data.INCI_** and **_Data.MORT_**: `data.frame` objects with the number of observed and expected cases (_'obs'_ and _'exp'_ variables, respectively) and standardized incidence/mortality ratio (_'SMR'_) for each administrative region (_'region'_) and time period (_'year'_) for brain cancer incidence and mortality data.
	- **_Carto.COM_**: `SpatialPolygonsDataFrame` object with the cartography of the 27 administrative regions of Navarre and the Basque Country.
	- **_W_**: Spatial adjacency matrix of the administrative regions.


# R code
R code to fit with INLA (http://www.r-inla.org/) the two-stage spatial and spatio-temporal models described in Adin et al. (2019). All the R files are written by the authors of the paper.

**1. Main files**

- [**Example1_Spatial_CountData.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/Example1_Spatial_CountData.R)

  R code to fit two-stage spatial cluster models using Spanish stomach mortality data (stored in `StomachCancer_ESP.Rdata` file). It reproduces the results obtained in Section 5.1. of the present work.
  
- [**Example2_SpatioTemporal_CountData.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/Example2_SpatioTemporal_CountData.R)

  R code to fit two-stage spatio-temporal cluster models using brain cancer incidence data in the regions of Navarre and Basque Country (stored in `BrainCancer_MUN.Rdata` file). It reproduces the results obtained in Section 5.2. of the present work.
  
**2. Algorithms**

- [**AHC_function.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/AHC_function.R)

  Agglomerative hierarchical clustering algorithm for spatial data (Anderson et al., 2014).

- [**AHC_clusteringST_function.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/AHC_clusteringST_function.R)

  Agglomerative hierarchical clustering algorithm for spatio-temporal data.

- [**Model1.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/Model1.R) and [**Model2.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/Model2.R)

  R code to fit the spatial models described in Equations (1) and (2) of the present work for each of the cluster configuration candidate.
  
 - [**TLmodel1a.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/TLmodel1a.R), [**TLmodel1b.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/TLmodel1b.R), [**TLmodel2a.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/TLmodel2a.R) and  [**TLmodel2b.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/TLmodel2b.R)
 
   R code to fit the spatio-temporal models with purely spatial cluster structures described in Section 3.1 of the present work for each cluster configuration candidate.
  
 - [**OptionI.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/OptionI.R), [**OptionII.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/OptionII.R), [**OptionIII.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/OptionIII.R) and [**OptionIV.R**](https://github.com/spatialstatisticsupna/Discontinuities_and_Clusters_article/blob/master/R/OptionIV.R)
 
   R code to fit the spatio-temporal models with space-time cluster structures described in Section 3.2 of the present work for each cluster configuration candidate.
  
  
# References
[Adin, A., Lee, D., Goicoa, T., and Ugarte, M.D. (2019). A two-stage approach to estimate spatial and spatio-temporal disease risks in the presence of local discontinuities and clusters. _Statistical Methods in Medical Research_, __28(9)__, 2595-2613.](https://doi.org/10.1177/0962280218767975)

[Anderson, C., Lee, D., and Dean, N. (2014). Identifying clusters in Bayesian disease mapping. _Biostatistics_, __15__, 457–469.](https://doi.org/10.1093/biostatistics/kxu005)
