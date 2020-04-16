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



# References
[Adin, A., Lee, D., Goicoa, T., and Ugarte, M.D. (2019). A two-stage approach to estimate spatial and spatio-temporal disease risks in the presence of local discontinuities and clusters. _Statistical Methods in Medical Research_, __28(9)__, 2595-2613.](https://doi.org/10.1177/0962280218767975)
