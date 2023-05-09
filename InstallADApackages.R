#Install Drought Atlas required packages
TotalADAPackages=c("raster", "rgdal", "countrycode", "rts", "ncdf4","HelpersMG", "latticeExtra", "reshape2", 
                   "gdata", "gtools", "plyr", "sqldf", "zoo", "Kendall","zyp", "car", "rgeos", "lmom" , "lmomRFA",
                   "sp", "rrcov", "nsRFA", "circular", "reshape" , "ModelMap", "maptools", "RSAGA" ,"stringr", 
                   "rasterVis", "hydroGOF" , "randomForest","progress", "deldir", "ggplot2", "edarf", "chron",
                   "lattice", "homtest")


usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
}
usePackage(TotalADAPackages)
