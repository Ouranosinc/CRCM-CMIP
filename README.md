# CRCM-CMIP
[![DOI](https://zenodo.org/badge/790831449.svg)](https://zenodo.org/doi/10.5281/zenodo.11061924)


### Data reference
Paquin, D., C. McCray, C. B. Gauthier, M. Giguère, O. Asselin, P .Bourgault, M.-P. Labonté and D. Matte. The CRCM5-CMIP6 Ouranos’ ensemble : A dynamically-downscaled ensemble of CMIP6 simulations over North America. to be submitted to Scientific Data


[Ouranos](https://www.ouranos.ca/en) : Canadian Regional Climate Model – version 5

**Martynov et al. 2013, Separovic et al. 2013**

Based on GEM 3.3.3.1

#### Configuration
NAM-11 CORDEX North American domain at 0.11° 695x668 grid points including a 20-point sponge (and halo) zone surrounding the domain, 5-minute time steps, 
xlat1=28.525 xlon2=145.955. 56 vertical levels and a top at 10 hPa. 17 surface levels and a bottom at 15 m.

#### Spectral Nudging
A spectral nudging is applied to the horizontal wind component with a half-response wavelength of 1177km and a 
relaxation time of 13.34 h. The nudging strength is set to zero from the surface to a height of 500 hPa and increases 
linearly onward to the top of the model’s simulated atmosphere (10 hPa)

### Parameterization
#### Atmosphere
Precipitation: modified Sundqvist  (1998); precipitation partition Bourgouin  (2000) ; Implicit vertical diffusion. 
Shallow convection: Kuo (1965) transient shallow, Non‐cloudy boundary layer formulation. 
Deep convection: Kain-Fritsch (1990); 
Radiation: Li & Barker (2005)

#### Surface
CLASS3.5c (Verseghy, 1993)

Lake model: FLake

#### Ocean
Prescribed SST & sea ice fraction

#### Aerosols
Prescribed

### Data Access
Due to its large size, the full dataset can't yet be shared publicly.

A subset of the variables are stored on Ouranos' THREDDS server.

- Annual files : https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/catalog/birdhouse/disk2/ouranos/CORDEX/catalog.html
- Aggregated datasets : https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/catalog/datasets/simulations/RCM-CMIP6/catalog.html

Other variables can be provided upon request by writing to simulations_ouranos@ouranos.ca.

### Acknowlegments
Developed by the [ESCER Centre](https://escer.uqam.ca/) at UQAM (Université du Québec à Montréal) with the collaboration
of Environment and Climate Change Canada (ECCC). **CRCM5; Martynov et al. 2013, Separovic et al. 2013**

The CRCM5 data has been generated and supplied by Ouranos.

CRCM5 computations were made on the supercomputers beluga and narval managed by Calcul Québec and the [Digital Research 
Alliance of Canada](https://alliancecan.ca/en). The operation of this supercomputer received financial support from 
Innovation, Science and Economic Development Canada and the Ministère de l’Économie et de l’Innovation du Québec.

### Some references for CRCM5
Asselin, M. Leduc, D. Paquin, K. Winger, A. Di Luca, M. Bukovsky, B. Music, and M. Giguère (2022). On the 
Intercontinental Transferability of Regional Climate Model Response to Severe Forestation.  MDPI's Climate 
https://doi.org/10.3390/cli10100138 

Bresson, E., R. Laprise, D. Paquin, J. M. Thériault, R. de Elia, 2017: Evaluating CRCM5 ability to simulate mixed 
precipitation. Atmosphere-Ocean. 55(2); 79-93. http://dx.doi.org/10.1080/07055900.2017.1310084 

Leduc, M., A. Mailhot, A. Frigon, J.-L. Martel, R. Ludwig, G.B. Brietzke, M. Giguère, F. Brissette, R. Turcotte, M. 
Braun, (2019) ClimEx project: a 50-member ensemble of climate change projections at 12-km resolution over Europe and 
northeastern North America with the Canadian Regional Climate Model (CRCM5). Journal of Applied Meteorology and 
Climatology. doi: 10.1175/JAMC-D-18-0021.1

Martynov A, R Laprise, L Sushama, K Winger, L Separovic, B Dugas. 2013. Reanalysis-driven climate simulation over CORDEX
North America domain using the Canadian Regional Climate Model, version 5: model performance evaluation. Clim Dyn 
41:2973-3005. DOI 10.1007/s00382-013-1778-9

Martynov A, L Sushama, R Laprise, K Winger, B Dugas. 2012. Interactive lakes in the Canadian regional climate model 
version 5: the role of lakes in the regional climate of North America. Tellus A 64, 016226. DOI: 
10.3402/tellusa.v64i0.16226.

Martynov A, L Sushama, R Laprise. 2010. Simulation of temperate freezing lakes by one-dimensional lake models: 
performance assessment for interactive coupling with regional climate models. Boreal Env Res 15:143-164.

Matte, D., Thériault, J. M., & Laprise, R. (2019). Mixed precipitation occurrences over southern Québec, Canada, under 
warmer climate conditions using a regional climate model. Climate Dynamics, 53(1), 1125–1141. 
https://doi.org/10.1007/s00382-018-4231-2

McCray, C. D., D. Paquin, J. M. Thériault, É. Bresson (2022). A multi-algorithm analysis of projected changes to 
freezing rain over North America in an ensemble of regional climate model simulations. Journal of Geophysical Research -
Atmospheres https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JD036935

McCray, D. C., J. M. Thériault, D. Paquin, É. Bresson, 2022. Quantifying the impact of precipitation-type algorithm 
selection on the representation of freezing rain in an ensemble of regional climate model simulations. Journal of 
Applied Meteorology and Climatology. 
https://journals.ametsoc.org/view/journals/apme/aop/JAMC-D-21-0202.1/JAMC-D-21-0202.1.xml 

McCray, C.D., G. Schmidt, D. Paquin, M. Leduc, Z. Bi, M. Radiyat, C. Silverman, M. Spitz, B. Brettschneider (2023). 
Changing Nature of High-Impact Snowfall Events in Eastern North America. Journal of Geophysical Research: Atmospheres. 
https://doi.org/10.1029/2023JD038804

Mironov D, E Heise, E Kourzeneva, B Ritter, N Schneider, A Terzhevik. 2010. Implementation of the lake parameterisation 
scheme FLake into the numerical weather prediction model COSMO. Boreal Env Res 15:218-230.

Mittermeier, M., E. Bresson, D. Paquin, R. Ludwig, 2021 A deep learning approach for the identification of long-duration
mixed precipitation in Montréal (Canada). Atmosphere-Ocean. https://doi.org/10.1080/07055900.2021.1992341

Riette S, D Caya. 2002. Sensitivity of short simulations to the various parameters in the new CRCM spectral nudging. – 
In: RITCHIE, H. (Ed.): Research activities in Atmospheric and Oceanic Modeling, WMO/TD No. 1105, Report No. 32: 
7.39–7.40.

Pérez Bello, A., A. Mailhot and D. Paquin, 2021 The response of daily and sub-daily extreme precipitations to changes in
surface and dew point temperatures. Journal of Geophysical Research – Atmospheres http://dx.doi.org/10.1029/2021JD034972

Pérez Bello, A., A. Mailhot, D. Paquin and D. Paquin-Ricard (2022). Temperature-precipitation scaling rates: a rainfall 
event-based perspective. Journal of Geophysical Research – Atmospheres. 
https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JD037873

Separovic L, A Alexandru, R Laprise, A Martynov, L Sushama, K Winger, K Tete, M Valin. 2013. Present climate and climate
change over North America as simulated by the fifth-generation Canadian regional climate model. Clim Dyn 41:3167-3201. 
DOI 10.1007/s00382-013-1737-5.

St-Pierre, M., J. Thériault and D. Paquin, 2019. Influence of the model spatial resolution on atmospheric conditions
leading to freezing rain in regional climate simulations. Atmosphere-Ocean, 
https://doi.org/10.1080/07055900.2019.1583088.
