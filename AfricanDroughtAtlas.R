##################################################
## -------------------------------------------------
## Project: AfricanDroughtAtlas
## Script purpose: 
## Updated: 09-05-2023
## Date Created: 02-08-2018
## Author 1: J. Nuñez
## Author 2: H. Maureira
## Author 3: P. Rojas
## Email: hmaureria@cazalac.org
## -------------------------------------------------
## Notes: 
## - CRU data from: http://data.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_3.21/data/pre/cru_ts3.21.1921.1930.pre.dat.gz
## - CHIRPS data from IRIDL: http://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/
##################################################

# Libraries ----------------------
library(raster);library(rgdal);library(countrycode);library(rts);require(ncdf4);library(HelpersMG)
library(plyr);library(latticeExtra);library(reshape2);library(gdata);library(corrplot);library(sqldf)
library(zoo);library(Kendall);library(zyp);library(car);library(gtools);# Para mantener orden texto-numero en id_station
library(rgeos);library(lmom);library(lmomRFA);library(sp);library(rrcov);library(nsRFA);library(ModelMap)
library(maptools);library(stringr);library(rasterVis);library(hydroGOF);library(randomForest);library(progress);#Check proper installation of SAG-GIS and RSAGA
library(gtools);library(here);library(chron);library(lattice);library(RColorBrewer);
library(sf);library(circular);library(reshape);library(deldir)
#library(RSAGA);
# Config ----------------------
#Option 1: Replace with the corresponding three ISO letters. In this case, Botswana should be 
#replaced with 'BWA'.
#Option 2: The file name of the shape without the extension, located in the 'shape' folder. 
#For example, if the file is named 'file.shape', the input should be 'file'.
country="DZA"


#  Optinal Config =======================

#workdir = "C:/Users/pablo/OneDrive/Escritorio/CAZALAC/ADA3/ADA/"
workdir <- here()
setwd(workdir)

# Load Functions -------------------
source('DroughAtlasFunctions.R')

# I. DATABASE CONSTRUCTION ------------
# Block i.a. database construction from CRU 3.21 or CHIRPS
# Choose one of the two model options based on country size    
# Option 1: "CRU"  CRU 3.21
# Option 2: "CHIRPS" from http://iridl.ldeo.columbia.edu resolution 0.25'
# Clip_method: Rectangle or Shape

database_creation(model="CHIRPS", country = country, clip_method="Shape" )

# II. VARIABLES AND INDICES  CALCULATION --------------

output_2 <- step2_variable_calculation(country = country)

# III. REGIONALIZATION USING MINIMUM DISTANCE SEARCHING  ALGORITHM -----------------
# Minimum Distance Searching Algorithm

output_3 <- regionalization(output_2$BaseDatosEstaciones,output_2$BaseRegistrosPr)

# IV: EXPLORATORY DATA ANALYSIS (OPTIONAL) #                                     #

exploratory(BaseDatosEstacionesClust = output_3$BaseDatosEstacionesClust, country = country, VarInter="CumSumDec")

# V. L-MOMENTS BASED REGIONAL FREQUENCY ANALYSIS ------------------------------

output5 <- regional_frequency_analysis(output_3$ClustLevels,output_3$NombreClusters,output_3$VarInter, output_3$BaseDatosEstacionesClust, output_2$BaseRegistrosPr, output_2$z,output_3$Hc, country)

# VI. FREQUENCY/QUANTIL/RETURN PERIOD ESTIMATION ---------------------------------
output6 <- period_estimation(output5$BaseSummaryStatistics, output5$BaseDatosEstacionesClust, output5$BaseBaseRegiones, output5$ResumeTable, output5$BaseMediaCompleta, output5$BaseProporCeros, output5$Basep0bias, output5$Baserfitdist, output5$Baserfitpara)

# VII: FREQUENCY/RETURN PERIOD/QUANTIL MAPPING ------------------
period_mapping(output5$ResumeTable, output6$BaseModelMapCor, country, output_2$Boundaries)

