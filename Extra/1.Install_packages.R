#Install Drought Atlas required packages
requiredPackages=c("raster","rgdal","countrycode","rts","ncdf4","HelpersMG",
                   "plyr","latticeExtra","reshape2","gdata","corrplot","sqldf",
                   "zoo","Kendall","zyp","car","gtools",
                   "rgeos","lmom","lmomRFA","sp","rrcov","nsRFA","ModelMap",
                   "maptools","stringr","rasterVis","hydroGOF","randomForest","progress",
                   "gtools","here","chron","lattice","RColorBrewer",
                   "sf","circular","reshape","deldir","shiny","leaflet","shinybusy","shinyjs","geodata")

for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
