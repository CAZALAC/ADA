workdir <- here()
setwd(workdir)
print(workdir)
# Open CRU nc file --------------------
nc <- nc_open("cru_ts3.24.01.1901.2015.pre.dat.nc")
print(nc)

# Extract variables ====================
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
time <- ncvar_get(nc, "time")
time.index=as.Date(as.numeric(time),origin="1900-01-01")
pre = list()
pre$x = ncvar_get(nc, "lon")
pre$y = ncvar_get(nc, "lat")
pre$z = ncvar_get(nc, "pre", start=c(300, 90, 829), count=c(213, 183, 552))#1970 a 2015, lon=-30.25 to 75.75 lat=-45.25 to 45.75
# Just for Africa
xmin=lon[300];xmax=lon[300+213];ymin=lat[90];ymax=lat[90+183]
#New version
#xmin=min(lon);xmax=max(lon);ymin=min(lat);ymax=max(lat)

# Create rasterstack with precipitation layers ==================
pre.layers<-raster::stack()
for (i in 1:dim(pre$z)[3]){
  r.pre=flip(raster(t(pre$z[,,i])),direction='y')
  extent(r.pre)=c(xmin,xmax,ymin,ymax)
  crs(r.pre) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  pre.layers<-stack(pre.layers,r.pre)
}
plot(pre.layers)


pre.monthly.AFR = pre.layers

plot(pre.monthly.AFR[[1]])
class(pre.monthly.AFR)
#    writeRaster(pre.monthly.AFR, 'testo.tif')
# Rasterstack time series creation with monthly data  ================
d= time.index[seq(829,length.out=dim(pre$z)[3])]
nlayers(pre.layers)==length(d)
nlayers(pre.monthly.AFR)==length(d)
pre.l.ts<-rts(pre.monthly.AFR,d)

#Rasterstack time series creation with annual data ===================
pre.anual <- apply.yearly(pre.l.ts, sum)
plot(mean(pre.anual@raster))
#pre.anual.2=mask(pre.anual[[1]],polygons(Boundaries))
#pre.anual.2=crop(pre.anual.2,polygons(Boundaries))

country="demodata"
# DataBase creation for a given country -------------------------
#Create working directory and set working directory for a given country
dir.create(paste0("./",country))# Only the first time and if it has not been created
setwd(paste0("./",country))

#writeOGR(Boundaries, dsn = '.', layer = country, driver = "ESRI Shapefile",overwrite_layer=TRUE)
## S4 method for signature 'Spatial'
#shapefile(Boundaries, filename=country, overwrite=TRUE)
#plot(Boundaries)
#rm(Boundaries)
#.....................
# manual correction of country's shapefile in SAGA GIS and then save it to the folder
#.....................
country.rts = pre.monthly.AFR 

country.rts.monthly=rts(country.rts,d)
country.rts.anual <- apply.yearly(country.rts.monthly, sum)
plot(mean(country.rts.anual@raster))

#Here a maximum number of stations is defined ======================================

useful.pixels=length(na.omit(getValues(country.rts[[1]])))

rs=sampleRandom(country.rts[[1]], cells=TRUE,xy=TRUE,size=useful.pixels,sp=TRUE)

points(rs,col="red")
plot(rs)
#writeOGR(obj=rs, dsn=getwd(),layer="randomsample", driver="ESRI Shapefile",overwrite_layer=TRUE)
shapefile(rs, filename="randomsample", overwrite=TRUE)
