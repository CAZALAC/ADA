workdir <- here()
setwd(workdir)
print(workdir)
country = "DEMODATA2"
http://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.daily-improved/.global/.0p05/.prcp/T/%28Jan%201981%29/%28Dec%202016%29/RANGE/weeklytomonthly/T/253.5/VALUE/Y/-49.975/49.975/RANGEEDGES/data.nc
#Create working directory and set working directory
dir.create(paste0("./",country,sep=""))
setwd(paste0("./",country,sep=""))


#If download fail, then copy paste url in a web browser and download the data.nc file manually
#After download, put the file into the country's folder.
# Download netCDF file

dfile=c("data.nc")      
print(url)
#Additional time for slow connections and if DL takes time to process

datos_ncdf <- nc_open(paste0(getwd(),"/data.nc"))
#Extraer longitud
lon <- ncvar_get(datos_ncdf, "X")
#Extraer latitud
lat <- ncvar_get(datos_ncdf, "Y")

#Extraer tiempo (esto puede ser confuso, pero como tenemos las fechas, no habr?a problema en modificarlo) (1 jan 1981) - (31 Dic 2016)
tpo <- ncvar_get(datos_ncdf, "T")
tpo2 <- seq(from = as.Date("1981-01-01"), to = as.Date("2016-12-31"), by = "month")

#get the name of the variable
var_name_ndf = names(datos_ncdf$var)[1]

# get precipitation
#prcp_array <- ncvar_get(datos_ncdf, "precipitation") #prcp aparece en float
prcp_array <- ncvar_get(datos_ncdf, var_name_ndf) #prcp aparece en float

dlname <- ncatt_get(datos_ncdf,var_name_ndf,"long_name")
dunits <- ncatt_get(datos_ncdf,var_name_ndf,"units")
fillvalue <- ncatt_get(datos_ncdf,var_name_ndf,"_FillValue")
dim(prcp_array)


# create dataframe -- reshape data
# matrix (nlon*nlat rows by 2 cols) of lons and lats
# De esta forma, se obtiene la cantidad de estaciones virtuales en el CUADRADO.
lonlat <- as.matrix(expand.grid(lon,lat))
dim(lonlat)

#Transformar el array en un rasterbrick (no es necesario parece)
#array_brick <- brick(prcp_array, xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))

# Se selecciona mapa del pa?s para seleccionar puntos de interes ==========


#Se usan los l?mites del mapa para recortar la base de datos


#Crear shapefile de puntos a partir de lonlat (corresponde a todas las estaciones virtuales, o sea, el rect?ngulo)
est_cuadro <- data.frame(lonlat)
est_cuadro$ID <- paste("id_", 1:nrow(est_cuadro), sep = "")
#est_cuadro <- SpatialPointsDataFrame(coords = lonlatDF[,1:2], proj4string = proj4string(mapa), bbox = lonlatDF)
coordinates(est_cuadro) <- ~ Var1 + Var2

#Se edita el archivo lonlat para seleccionar aquellos puntos que est?n dentro del l?mite pa?s
est_selec <- est_cuadro #perfecto, se guarda
plot(est_selec)




#Se obtiene un archivo (dataframe) lonlat con los puntos seleccionados
est_selecDF <- data.frame(est_selec)

#Se crea un data frame que almacenar? en cada lista, los valores por d?a
tpo2 <- seq(from = as.Date("1981-01-01"), to = as.Date("2016-12-31"), by = "month")

DF_est <- data.frame(date = tpo2, prcp_chirps = NA)

#Se crea una lista cuyos nombres corresponde a lat long (que en el fondo son las estaciones)
lista_est <- list()

for (i in 1:nrow(est_selecDF)){ #Continuar desde ac?...
  
  #print(paste0("station ",i," of ",nrow(est_selecDF)))
  
  lista_est[[i]] <- DF_est
  names(lista_est)[i] <- as.character(est_selecDF[i,3])
  
}

#........................................
# Extraer valores desde el prcpc_array usando "est_selecDF"
# almacenar en la lista

l = 1
for (j in 1:dim(prcp_array)[3]){ #ac? no es array brick, OJO
  print(paste0("month ",j," of ",dim(prcp_array)[3]))
  
  raster_uni <- prcp_array[,,j]
  raster_uni <- raster(t(raster_uni), xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
  raster_uni <- flip(raster_uni, "y")
  
  for (k in 1:nrow(est_selecDF)){
    
    valor_fecha <- extract(raster_uni, matrix(c(est_selecDF[k,1], est_selecDF[k,2]), ncol = 2))
    lista_est[[k]][l,2] <- valor_fecha
    
  }
  
  l = 1+l
  
} # ac? finalizado debiera estar ordenados los valores por fecha.


#...........................................
# Exportar datos como txt
dir.create(paste0(getwd(),"/est_procesadas/"))

#fix by pablo, elimina columnas enteras que tengan na
lista_est <- Filter(function(x)!all(is.na(x$prcp_chirps)), lista_est)



proj4string(lista_est)<- CRS("+proj=longlat +datum=WGS84")
#the same format that nc file
raster::shapefile(lista_est, "randomsample2.shp", overwrite=TRUE)
setwd(workdir)
















for (m in 1:length(lista_est)){
  print(paste0("station ", m, " of ",length(lista_est)))
  
  export <- lista_est[[m]]
  write.table(export, paste0(getwd(),"/est_procesadas/", names(lista_est)[m], "_latam.txt"), sep =";", row.names = FALSE)
}

# Se exporta la metadata (Estaciones virtuales)
metaexp <- est_selecDF[,c(3,1,2)]
#si se encuentra dentro de las filtradas sin nan
metaexp[metaexp$ID %in% names(lista_est),]

colnames(metaexp) <- c("station_id", "lon", "lat")

write.table(metaexp, paste0(getwd(),"/est_procesadas/","metadata.txt"), sep = ";", row.names = FALSE)

#.............................................................
#Convertir al formato CAZALAC BaseDatosEstaciones y BasedatosRegistros

EstacionesLong=list()
latam.list=list.files(path=paste0(getwd(),"/est_procesadas") ,
                      pattern = "\\latam.txt$")
latam.list=substr(latam.list, 1, nchar(latam.list)-10)
for (j in 1:length(latam.list)){
  EstacionesLong[[j]]=read.table(paste0(getwd(),"/est_procesadas/",latam.list[j],"_latam.txt"),sep=";",header=TRUE)
}

BaseDatosRegistros=data.frame(id_station=rep(latam.list[1],36),
                              Year=1981:2016,
                              matrix(EstacionesLong[[1]][,2],ncol=12,byrow=T))

for (k in 2:length(latam.list)){
  bd=data.frame(id_station=rep(latam.list[k],36),
                Year=1981:2016,
                matrix(EstacionesLong[[k]][,2],ncol=12,byrow=T))
  BaseDatosRegistros=rbind(BaseDatosRegistros,bd)
}


# CORRECCION DE CEROS. TODOS LOS CEROS SON REEMPLAZADOS POR UN VALOR RANDOM UNIF ENTRE 0 Y 1



#BaseDatosEstaciones. Defino el tama?o total a 500
BaseDatosEstaciones=read.csv(paste0(getwd(),"/est_procesadas/metadata.txt"),sep=";",header=TRUE)
colnames(BaseDatosEstaciones)[1]="id_station"
row.names(BaseDatosEstaciones)=NULL
BaseDatosEstaciones = BaseDatosEstaciones[BaseDatosEstaciones$id_station %in% names(lista_est),]


newshapefile <- BaseDatosEstaciones
coordinates(newshapefile)=~lon+lat
proj4string(newshapefile)<- CRS("+proj=longlat +datum=WGS84")
#the same format that nc file
newshapefile@data["x"] <- BaseDatosEstaciones$lon
newshapefile@data["y"] <- BaseDatosEstaciones$lat
newshapefile@data["cell"] <- 1
raster::shapefile(newshapefile, "randomsample.shp", overwrite=TRUE)
setwd(workdir)