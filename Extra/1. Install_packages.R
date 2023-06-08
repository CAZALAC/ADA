#Install Drought Atlas required packages
TotalADAPackages=c("raster","rgdal","countrycode","rts","ncdf4","HelpersMG",
                   "plyr","latticeExtra","reshape2","gdata","corrplot","sqldf",
                   "zoo","Kendall","zyp","car","gtools",
                   "rgeos","lmom","lmomRFA","sp","rrcov","nsRFA","ModelMap",
                   "maptools","stringr","rasterVis","hydroGOF","randomForest","progress",
                   "gtools","here","chron","lattice","RColorBrewer",
                   "sf","circular","reshape","deldir")

usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
}
usePackage(TotalADAPackages)
