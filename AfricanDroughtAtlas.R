##################################################
## -------------------------------------------------
## Project: AfricanDroughtAtlas
## Script purpose: 
## Date Created: 09-05-2023
## Baed on CRU: 02-08-2018
## Author 1: J. Nuñez
## Author 2: H. Maureria
## Author 3: P. Rojas
## Email: hmaureria@cazalac.org
## -------------------------------------------------
## Notes: 
##        
##        -THIS TRY IS WITH REG ADAPT WITHOUT MODIFIERS
##################################################

# Libraries ----------------------
library(raster);library(rgdal);library(countrycode);library(rts);require(ncdf4);library(HelpersMG)
library(plyr);library(latticeExtra);library(reshape2);library(gdata);library(corrplot);library(sqldf)
library(zoo);library(Kendall);library(zyp);library(car);library(gtools);# Para mantener orden texto-numero en id_station
library(rgeos);library(lmom);library(lmomRFA);library(sp);library(rrcov);library(nsRFA);library(ModelMap)
library(maptools);library(stringr);library(rasterVis);library(hydroGOF);library(randomForest);library(progress);#Check proper installation of SAG-GIS and RSAGA
library(RSAGA);library(gtools);library(here);library(chron);library(lattice);library(RColorBrewer);
library(sf)


# Config ----------------------
#Option 1: Replace with the corresponding three ISO letters. In this case, Botswana should be 
#replaced with 'BWA'.
#Option 2: The file name of the shape without the extension, located in the 'shape' folder. 
#For example, if the file is named 'file.shape', the input should be 'file'.
country="BWA"


#  Optinal Config =======================

#workdir = "C:/Users/pablo/OneDrive/Escritorio/CAZALAC/ADA3/ADA/"
workdir <- here()
setwd(workdir)

# Load Functions -------------------
source('DroughAtlasFunctions.R')

# Test -------------------
rsaga.env()
rsaga.get.version()


# I. DATABASE CONSTRUCTION ------------
# Block i.a. database construction from CRU 3.21 or CHIRPS
# Choose one of the two model options based on country size    
# Option 1: "CRU"  CRU 3.21
# Option 2: "CHIRPS" from http://iridl.ldeo.columbia.edu resolution 0.25'

database_creation(model="CRU", country = country )

# II. VARIABLES AND INDICES  CALCULATION --------------

output_2 <- step2_variable_calculation(country = country)

# III. REGIONALIZATION USING MINIMUM DISTANCE SEARCHING  ALGORITHM -----------------
# Minimum Distance Searching Algorithm

output_3 <- regionalization(output_2$BaseDatosEstaciones,output_2$BaseRegistrosPr)

# V. L-MOMENTS BASED REGIONAL FREQUENCY ANALYSIS ------------------------------

output4 <- regional_frequency_analysis(output_3$ClustLevels,output_3$NombreClusters,output_3$VarInter, output_3$BaseDatosEstacionesClust, output_2$BaseRegistrosPr, output_2$z,output_3$Hc, country)

# VI. FREQUENCY/QUANTIL/RETURN PERIOD ESTIMATION ---------------------------------
output5 <- period_estimation(output4$BaseSummaryStatistics, output4$BaseDatosEstacionesClust, output4$BaseBaseRegiones, output4$ResumeTable, output4$BaseMediaCompleta, output4$BaseProporCeros, output4$Basep0bias, output4$Baserfitdist, output4$Baserfitpara)

# VII: FREQUENCY/RETURN PERIOD/QUANTIL MAPPING ------------------
period_mapping(output4$ResumeTable, output5$BaseModelMapCor, country, output_2$Boundaries)


#END

#################################################################################################################
#                                           BLOCK IV: EXPLORATORY DATA ANALYSIS    (OPTIONAL)                   #                                     #
# #################################################################################################################
# library("sqldf"); library("circular");library("lmom");library("reshape");library("zoo");library("reshape")
# #Several plots showing main charateristics of the stations in the specified region
# region.of.interest=3 #Some of the available regions in ClustReg_adapt or a number between 1 and ClustLevels
# #VarInter<-"CumSumDec"
# ListadoEstaciones<-sqldf(paste("select id_station from BaseDatosEstacionesClust where ClustReg_adapt=='",region.of.interest,"'",sep=""))
# LEst<-levels(factor(ListadoEstaciones["id_station"][,]))
# Listado<-as.character(LEst)
# 
# #.......Save plots to a PDF
# pdf(paste0("EDA_",country,"_Region_",region.of.interest,".pdf"), width = 16 , height = 10, title = paste0("EDA Plots for ",country,": Region",region.of.interest))
# par(mfrow=c(3,2))
# for (m in 1:length(Listado)){
#   #Fig1. Boxplot
#   Boxplotmensual<-sqldf(paste("select Jan,Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec from BaseCompletaAdapt where id_station=='",Listado[m],"'",sep=""))
#   boxplot(Boxplotmensual,   main=Listado[m],pch=20,cex=1, cex.main=1, cex.axis=0.8, cex.lab=0.1,  las=1, xlab="Month", ylab="Precipitation[mm]")
#   
#   #Fig2. Barplot
#   BP<-matrix(colMeans(Boxplotmensual,na.rm=T),nrow=1)
#   barplot(BP,names.arg=c("E","F","M","A","M","J","J","A","S","O","N","D"),main=Listado[m], cex.names=0.8,cex.main=1,  xlab="Mes", ylab="Precipitation[mm]")
#   
#   Var<-sqldf(paste("select ",VarInter," from BaseCompletaAdapt where id_station=='",Listado[m],"'",sep=""))
#   #Fig3. Histogram
#   hist(as.numeric(Var[,1]),prob=T,xlab="Precipitation [mm]",ylab="Density",main=paste(Listado[m]," en ",VarInter),cex.main=1)
#   lines(density(as.numeric(Var[,1]),na.rm=T),col=451,lwd=2)
#   
#   # Fig4. Extreme value plot
#   evplot(as.numeric(Var[,1]),xlab="Gumbel Reduced Variable",ylab="Quantil",main=paste(Listado[m]," en ",VarInter),cex.main=1)
#   evdistq(quanor, pelnor(samlmu(as.numeric(Var[,1]))))
#   
#   # Fig5. L-moment ratio diagram
#   lmrd(samlmu(as.numeric(Var[,1])),xlim=c(-0.3,0.6),pch=19,col="red",cex=1.5,legend.lmrd=list(cex=1),main=paste(Listado[m]," en ",VarInter),cex.main=1,cex.lab=1)
#   
#   #Fig6.Time series plot 
#   TSP<-sqldf(paste("select Year, Jan,Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec from BaseCompletaAdapt where id_station=='",Listado[m],"'",sep=""))
#   dfm <- melt(TSP, id.var="Year", variable_name="Mes")
#   names(dfm) <- c("Year", "Mes", "Precipitacion")
#   pp_num <- unclass(dfm$Mes)
#   pp_frac <- as.numeric( (pp_num-0.5)/12   )
#   yr_frac <- as.numeric(dfm$Year) + pp_frac
#   dfm <- data.frame(yr_frac,  dfm)
#   dfm <- dfm[order(dfm$yr_frac), ]
#   dfmzoo<-zoo(dfm$Precipitacion,as.yearmon(dfm$yr_frac))
#   plot(dfmzoo,type="h",xlab="Year",ylab="Precipitation[mm]",main=Listado[m],cex.main=1,cex.lab=1)
#   
# }
# dev.off()
# # ..............................................END OF BLOCK IV .................................................
