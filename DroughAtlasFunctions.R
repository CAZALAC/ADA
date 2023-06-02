#DroughtAtlasFunctions
#AfricanDroughtAtlas2018v2 17/10/2018

#Funcion 1. Parameters estimation for the Gaucho distribution model
pelgau=function (x) 
{
  rLm=regavlmom(x)
  lambda1 <- rLm[1]
  lambda2 <- rLm[2]
  tau3 <- rLm[3]
  tau4 <- rLm[4]
  sumquad.tau3tau4 = function(k.h, t3.t4) {
    k <- k.h[1]
    h <- 0.5
    t3 <- t3.t4[1]
    t4 <- t3.t4[2]
    error <- FALSE
    if (((k < -1) && (h >= 0)) || ((h < 0) && ((k <= -1) || 
                                               (k >= -1/h)))) {
      stop("L-moments are defined if h>=0 and k>-1, or if h<0 and -1<k<-1/h")
      error = TRUE
    }
    g <- c(0, 0, 0, 0)
    if (h == 0) {
      tau3 <- 2 * (1 - 3^(-k))/(1 - 2^(-k)) - 3
      tau4 <- (5 * (1 - 4^(-k)) - 10 * (1 - 3^(-k)) + 6 * 
                 (1 - 2^(-k)))/(1 - 2^(-k))
    }
    else {
      for (r in 1:4) {
        if (h > 0) {
          g[r] <- (r * gamma(1 + k) * gamma(r/0.5))/(0.5^(1 + 
                                                            k) * gamma(1 + k + r/0.5))
        }
        else {
          g[r] = (r * gamma(1 + k) * gamma(-k - r/h))/((-h)^(1 + 
                                                               k) * gamma(1 - r/h))
        }
      }
      tau3 <- (-g[1] + 3 * g[2] - 2 * g[3])/(g[1] - g[2])
      tau4 <- -(-g[1] + 6 * g[2] - 10 * g[3] + 5 * g[4])/(g[1] - 
                                                            g[2])
    }
    if (error == FALSE) {
      output <- (t3 - tau3)^2 + (t4 - tau4)^2
    }
    else {
      output <- -1
    }
    return(output)
  }
  xi.alpha = function(lambda1, lambda2, k, h) {
    if (((k < -1) && (h >= 0)) || ((h < 0) && ((k <= -1) || 
                                               (k >= -1/h)))) {
      stop("L-moments are defined if h>=0 and k>-1, or if h<0 and -1<k<-1/h")
    }
    g <- c(0, 0)
    if (h == 0) {
      alpha <- (lambda2 * k)/((1 - 2^(-k)) * gamma(1 + k))
      xi <- lambda1 - alpha * (1 - gamma(1 + k))/k
    }
    else {
      for (r in 1:2) {
        if (h > 0) {
          g[r] <- (r * gamma(1 + k) * gamma(r/0.5))/(0.5^(1 + 
                                                            k) * gamma(1 + k + r/0.5))
        }
        else {
          g[r] = (r * gamma(1 + k) * gamma(-k - r/h))/((-h)^(1 + 
                                                               k) * gamma(1 - r/h))
        }
      }
      alpha <- (lambda2 * k)/(g[1] - g[2])
      xi <- lambda1 - alpha * (1 - g[1])/k
    }
    output <- list(xi = xi, alpha = alpha)
    return(output)
  }
  quanti <- length(tau3)
  k <- rep(NA, quanti)
  h <- rep(NA, quanti)
  xi <- rep(NA, quanti)
  alpha <- rep(NA, quanti)
  for (i in 1:quanti) {
    minimo <- optim(c(1,0.5), sumquad.tau3tau4, t3.t4 = c(tau3[i], 
                                                          tau4[i]))
    if (minimo$value != -1) {
      k[i] <- minimo$par[1]
      h[i] <- 0.5
      pp <- xi.alpha(lambda1[i], lambda2[i], k[i], 0.5)
      xi[i] <- pp$xi
      alpha[i] <- pp$alpha
    }
  }
  #output <- list(xi = xi, alpha = alpha, k = k, h = h)
  #return(output)
  
  
  index=x[,3]
  names(index)=x[,1]
  outputgau<-list(dist="gau",para=c(xi=xi,alpha=alpha,k=k,h=h),qfunc="kappah0.5",rmom=rLm,index=index,class="rfd")
  return(outputgau)
}

#Function extracted from package homtest, that is disabled from CRAN
#https://rdrr.io/cran/homtest/src/R/Lmoments.R
Lmoments <- function(x) {
  
  camp <- sort(x)
  n <- length(camp)
  
  nn <- rep(n-1,n)
  pp <- seq(0,n-1)
  p1 <- pp/nn
  p2 <- p1 * (pp-1)/(nn-1)
  p3 <- p2 * (pp-2)/(nn-2)
  
  b0 <- sum(camp)/n
  b1 <- sum(p1*camp)/n
  b2 <- sum(p2*camp)/n
  b3 <- sum(p3*camp)/n
  
  l1 <- b0
  l2 <- 2*b1-b0
  lcv <- 2*b1/b0-1
  lca <- 2*(3*b2-b0)/(2*b1-b0)-3
  lkur <- 5*(2*(2*b3-3*b2)+b0)/(2*b1-b0)+6
  
  Lmom <- c(l1,l2,lcv,lca,lkur)
  names(Lmom) <- c("l1","l2","lcv","lca","lkur")
  
  return(Lmom)
}


#Function 2. Function to compute ARF-LM statistics as in rgtest from lmomRFA
H1ZDi=function (y,Nsim) 
{
  #library(homtest) disabled so it can run. Lmoments function is enabled manually
  #if (length(x) != length(cod)) {
  #  stop("x and cod must have the same length")
  #}
  cod<-rep(1:length(y),lengths(y))
  fac <- factor(cod)
  x=unlist(y)
  ni <- tapply(x, fac, length)
  k <- nlevels(fac)
  Lm <- sapply(split(x, fac), Lmoments)
  rLm <- regionalLmoments(x,fac)
  ti <- as.numeric(Lm[3, ])
  t3i <- as.numeric(Lm[4, ])
  t4i <- as.numeric(Lm[5, ])
  lambda1Reg <- as.numeric(rLm[1])
  lambda2Reg <- as.numeric(rLm[2])
  tauReg <- as.numeric(rLm[3])
  tau3Reg <- as.numeric(rLm[4])
  tau4Reg <- as.numeric(rLm[5])
  V1 <- (sum(ni * (ti - tauReg)^2)/sum(ni))^0.5
  V2 <- sum(ni * ((ti - tauReg)^2 + (t3i - tau3Reg)^2)^0.5)/sum(ni)
  V3 <- sum(ni * ((t3i - tau3Reg)^2 + (t4i - tau4Reg)^2)^0.5)/sum(ni)
  bestdist=ifelse(tau4Reg<0.16667*tau3Reg^0+0.83333*tau3Reg^2,"kap","glo")
  #parkappa <- pelkap(c(lambda1Reg, lambda2Reg, tau3Reg, tau4Reg)) Cambiado por lo de abajo para enfrentar kappa
  if(bestdist=="kap"){
    pardist=pelkap(c(lambda1Reg, lambda2Reg, tau3Reg, tau4Reg))} else{
      pardist=pelglo(c(lambda1Reg, lambda2Reg, tau3Reg, tau4Reg))}
  
  #xi <- parkappa[1]#parkappa$xi
  #alfa <- parkappa[2]#parkappa$alfa
  #kappa <- parkappa[3]#parkappa$kappa
  #hacca <- parkappa[4]#parkappa$hacca
  V1s <- rep(NA, Nsim)
  V2s <- rep(NA, Nsim)
  V3s <- rep(NA, Nsim)
  bias<-rep(NA,Nsim)
  output<-list()
  
  for (i in 1:Nsim) {
    ti.sim <- rep(NA, k)
    t3i.sim <- rep(NA, k)
    t4i.sim <- rep(NA, k)
    for (j in 1:k) {
      if(bestdist=="kap"){
        campione <- rand.kappa(ni[j], pardist[1], pardist[2], pardist[3], pardist[4])} else {
          campione <- rand.genlogis(ni[j], pardist[1], pardist[2], pardist[3])}
      
      campione.ad <- campione/mean(campione)
      lmom <- Lmoments(campione.ad)
      ti.sim[j] <- lmom[3]
      t3i.sim[j] <- lmom[4]
      t4i.sim[j] <- lmom[5]
    }
    tauReg.sim <- sum(ni * ti.sim)/sum(ni)
    tau3Reg.sim <- sum(ni * t3i.sim)/sum(ni)
    tau4Reg.sim <- sum(ni * t4i.sim)/sum(ni)
    V1s[i] <- (sum(ni * (ti.sim - tauReg.sim)^2)/sum(ni))^0.5
    V2s[i] <- sum(ni * ((ti.sim - tauReg.sim)^2 + (t3i.sim - tau3Reg.sim)^2)^0.5)/sum(ni)
    V3s[i] <- sum(ni * ((t3i.sim - tau3Reg.sim)^2 + (t4i.sim - tau4Reg.sim)^2)^0.5)/sum(ni)
    bias[i]=tau4Reg.sim-tau4Reg
  }
  muV1 <- mean(V1s)
  stdV1 <- sd(V1s)
  muV2 <- mean(V2s)
  stdV2 <- sd(V2s)
  muV3 <- mean(V3s)
  stdV3 <- sd(V3s)
  
  B4=sum(bias)/Nsim
  std4=((sum(bias^2)-Nsim*B4^2)/(Nsim-1))^0.5
  
  
  if (k<=2) {
    H1=0.0
    H2=0.0
    H3=0.0
  } else {
    H1 <- (V1 - muV1)/stdV1
    H2 <- (V2 - muV2)/stdV2
    H3 <- (V3 - muV3)/stdV3}
  H<-c(H1,H2,H3)
  names(H)<-c("H1","H2","H3")
  #ZDIST calculus
  ZGLO=(lmrglo(pelglo(rLm[c(1,2,4,5)]),nmom=4)[4]-tau4Reg+B4)/std4
  ZGEV=(lmrgev(pelgev(rLm[c(1,2,4,5)]),nmom=4)[4]-tau4Reg+B4)/std4
  ZGNO=(lmrgno(pelgno(rLm[c(1,2,4,5)]),nmom=4)[4]-tau4Reg+B4)/std4
  ZPE3=(lmrpe3(pelpe3(rLm[c(1,2,4,5)]),nmom=4)[4]-tau4Reg+B4)/std4
  ZGPA=(lmrgpa(pelgpa(rLm[c(1,2,4,5)]),nmom=4)[4]-tau4Reg+B4)/std4
  ZGAU=(lmrkap(pelgau(regsamlmu(y))$para,nmom=4)[4]-tau4Reg+B4)/std4
  
  names(ZGLO)="glo"
  names(ZGEV)="gev"
  names(ZGNO)="gno"
  names(ZPE3)="pe3"
  names(ZGPA)="gpa"
  names(ZGAU)="gau"
  #Di calculus
  # .......Calculo de Di Cl?sico y Robusto...# Se calcula una medida alternativa de discordancia
  if (k<=3) {
    Di<-rep(1,k)
    
  } else {
    Dimatrix<-data.frame(t(Lm))[,3:5]
    Dimcd<-Cov(Dimatrix)
    Di<- sqrt(getDistance(Dimcd))}
  
  if (k<=6) {
    rDi<-rep(1,k)
  } else {
    rDimatrix<-data.frame(t(Lm))[,3:5]
    rDimcd<-CovMcd(rDimatrix)
    rDi<- sqrt(getDistance(rDimcd))}
  names(Di)=names(y)
  names(rDi)=names(y)
  
  #output
  output[[1]]<- H
  output[[2]]<-c(ZGLO,ZGEV,ZGNO,ZPE3,ZGPA,ZGAU)
  output[[3]]<-Di
  output[[4]]<-rDi
  names(output)[[1]] <- "H"
  names(output)[[2]] <- "Z"
  names(output)[[3]] <- "Di"
  names(output)[[4]] <-"rDi"
  
  #if (k<=3) {
  # return(output)
  #} else {
  #par(mfrow=c(2,1))
  #par(cex = 0.6)
  #par(mar = c(4, 4, 1, 1), oma = c(3, 3, 3, 3))
  #plot(Dimcd, which = "dist",main="Classic Di Plot")
  #plot(rDimcd, which = "dist",main="Robust Di Plot")
  return(output)}
#}

# Function 3. Function to compute Voronoi maps

voronoipolygons <- function(x) {
  require(deldir)
  require(sp)
  if (.hasSlot(x, 'coords')) {
    crds <- x@coords  
  } else crds <- x
  z <- deldir(crds[,1], crds[,2])
  w <- tile.list(z)
  polys <- vector(mode='list', length=length(w))
  for (i in seq(along=polys)) {
    pcrds <- cbind(w[[i]]$x, w[[i]]$y)
    pcrds <- rbind(pcrds, pcrds[1,])
    polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP <- SpatialPolygons(polys)
  voronoi <- SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1],
                                                          y=crds[,2], row.names=sapply(slot(SP, 'polygons'), 
                                                                                       function(x) slot(x, 'ID'))))
}

# Function 4. Function for random sample

randomSample = function(df,n) { 
  return (df[sample(nrow(df), n),])}


# Function 5. Function to create Drought Atlas Maps

Mapping.Function <- function(x, filename,model,folder.mapas) {
  out <- raster(x[[1]])
  out <- writeStart(out, filename, overwrite=TRUE)
  for (r in 1:nrow(out)) {
    v <- data.frame(getValues(x,r))
    v.pred <- rep(NA, length = nrow(v))
    nonPredict <- apply(((v == -999) | is.na(v)), 1, any)
    if(all(nonPredict)==TRUE) v.pred[!nonPredict]<-NA else v.pred[!nonPredict] <- predict(model, newdata=v[!nonPredict,],OOB=TRUE,type="response")
    out <- writeValues(out, v.pred, r)
    print(paste("Printing row=",r,sep=""))}
  out <- writeStop(out)
  writeRaster(out,filename=paste0(getwd(),folder.mapas,"/",filename,".sdat"),format="SAGA", overwrite=TRUE,NAflag=-999)}

# Function 6. Function to extract slope and P value from a linear regression model
lm.coeficients <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  s<-summary(modelobject)$coefficients[2]
  attributes(s) <- NULL
  attributes(p) <- NULL
  salida<-list(slope=s,Pvalue=p)
  return(salida)}


##################################################
## -------------------------------------------------
## Project: AfricanDroughtAtlas
## Script purpose: Functions to create Stations and Records DataBases from CRU and CHIRPS
## Date Created: 09-05-2023
## Author 1: J. Nuñez
## Author 2: H. Maureria
## Author 3: P. Rojas
## Email: hmaureria@cazalac.org
## -------------------------------------------------
## Notes:
## - CRU data from: http://data.ceda.ac.uk/badc/cru/data/cru_ts/cru_ts_3.21/data/pre/cru_ts3.21.1921.1930.pre.dat.gz
##       
##################################################

#server down, check
#library(geodata)
#debido a que se llama mucho
get_country_shape  <- function(country){
  sf_objet <- sf::st_as_sf(geodata::gadm(country=country, level=0, path=tempdir()))
  #legacy
  sf_objet <- as(sf_objet, "Spatial")
  return(sf_objet)
}

#Functions =====================================
country_list <- function(country){
  ISO.codes=read.csv("CountryISOCodes.csv",sep=";")
  Afr.country.list=as.character(ISO.codes$ThreeLetter)
  Afr.country.name=countrycode(Afr.country.list, "iso3c","country.name")
  return(Afr.country.name)
}
#Shape or country boundaries
shape_country <- function(country){
  # List of countries according to ISO 
  ISO.codes=read.csv("CountryISOCodes.csv",sep=";")
  Afr.country.list=as.character(ISO.codes$ThreeLetter)
  Afr.country.name=countrycode(Afr.country.list, "iso3c","country.name")
  
  # The input and return of variables can be improved 
  if(country %in% Afr.country.list){
    Boundaries  <- get_country_shape(country=country)
    #gadm(country=country, level=1, path=tempdir())
  }else{
    shape=paste0("Shape/",country,".shp",sep="")
    BoundariesAFR <- raster::shapefile(shape)
    #el script funciona con lonlat
    try(BoundariesAFR <- spTransform(BoundariesAFR, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"), silent = TRUE)
    #por si el archivo viene sin prj
    try(raster::projection(BoundariesAFR) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
    Boundaries <- BoundariesAFR
  }
  Bound.Pol = extent(Boundaries)
  if(Bound.Pol[1]>0) xmin=Bound.Pol[1]*0.5 else xmin=Bound.Pol[1]*1.5
  if(Bound.Pol[2]>0) xmax=Bound.Pol[2]*1.5 else xmax=Bound.Pol[2]*0.5
  if(Bound.Pol[3]>0) ymin=Bound.Pol[3]*0.5 else ymin=Bound.Pol[3]*1.5
  if(Bound.Pol[4]>0) ymax=Bound.Pol[4]*1.5 else ymax=Bound.Pol[4]*0.5
  x_coord <- c(xmin, xmax, xmax, xmin)
  y_coord <- c(ymin, ymin, ymax, ymax)
  xym <- cbind(x_coord, y_coord)
  p <- Polygon(xym)
  ps <- Polygons(list(p),1)
  sps <- SpatialPolygons(list(ps))
  try(proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") , silent = TRUE)
  return(list("Boundaries" = Boundaries, "sps" = sps))
}
#Zero correction
zero_correction <- function(BaseDatosRegistros){
  for (i in 3:14){
    baseline.index <- which(BaseDatosRegistros[,i] == 0)
    noise <- runif(length (  baseline.index ))
    BaseDatosRegistros[baseline.index,i] <- BaseDatosRegistros[baseline.index,i] + noise
  }
  return(BaseDatosRegistros)
}

#Funcion para muestreo de estaciones cuando hay mas de 3500
randomSample = function(df,n) { 
  return (df[sample(nrow(df), n),])}


# Config =================
#Option 1: Replace with the corresponding three ISO letters. In this case, Botswana should be 
#replaced with 'BWA'.
#Option 2: The file name of the shape without the extension, located in the 'shape' folder. 
#For example, if the file is named 'file.shape', the input should be 'file'.

database_creation <- function(model = "CRU", country = "BWA") {
  # BLOCK I.A. DATABASE CONSTRUCTION FROM CRU 3.21 ------------
  # debug: model = "CRU"; country = "Huasco"
  if(model=="CRU"){
    # Optional Config Setup, Country Code, Shape and raster creation =================
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
    #xmin=lon[300];xmax=lon[300+213];ymin=lat[90];ymax=lat[90+183]
    #New version
    xmin=min(lon);xmax=max(lon);ymin=min(lat);ymax=max(lat)
    
    # Create rasterstack with precipitation layers ==================
    pre.layers<-raster::stack()
    for (i in 1:dim(pre$z)[3]){
      r.pre=flip(raster(t(pre$z[,,i])),direction='y')
      extent(r.pre)=c(xmin,xmax,ymin,ymax)
      crs(r.pre) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
      pre.layers<-stack(pre.layers,r.pre)
    }
    plot(pre.layers)

    #Crop and mask rasterstack with continent boundaries ==========
    dummyreturn <- shape_country(country)
    Boundaries <- dummyreturn$Boundaries
    sps <- dummyreturn$sps
    rm(dummyreturn)
    #BoundariesAFR=readOGR("AfricaDA.shp") #from http://www.maplibrary.org/library/stacks/Africa/index.htm
    pre.monthly.AFR=mask(pre.layers,polygons(sps))
    pre.monthly.AFR=crop(pre.layers,polygons(sps))
 
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
    
    
    # DataBase creation for a given country -------------------------
    #Create working directory and set working directory for a given country
    dir.create(paste0("./",country))# Only the first time and if it has not been created
    setwd(paste0("./",country))

    #writeOGR(Boundaries, dsn = '.', layer = country, driver = "ESRI Shapefile",overwrite_layer=TRUE)
    ## S4 method for signature 'Spatial'
    shapefile(Boundaries, filename=country, overwrite=TRUE)
    #plot(Boundaries)
    #rm(Boundaries)
    #.....................
    # manual correction of country's shapefile in SAGA GIS and then save it to the folder
    #.....................
    
    plot(Boundaries)
    
    #Stations Database Generacion (BaseDatosEstaciones) for a given country ============
    #country.rts=crop(pre.monthly.AFR,Boundaries)
    #porque asi evito pedazos del area sin estaciones en mallas gruesas
    country.rts=crop(pre.monthly.AFR,sps)
    country.rts=mask(country.rts,sps)
    country.rts.monthly=rts(country.rts,d)
    country.rts.anual <- apply.yearly(country.rts.monthly, sum)
    plot(mean(country.rts.anual@raster))

    #Here a maximum number of stations is defined ======================================
    
    useful.pixels=length(na.omit(getValues(country.rts[[1]])))
    if(useful.pixels<1000) size=useful.pixels else size=1000
    
    rs=sampleRandom(country.rts[[1]], cells=TRUE,xy=TRUE,size=size,sp=TRUE)
    points(rs,col="red")
    plot(rs)
    #writeOGR(obj=rs, dsn=getwd(),layer="randomsample", driver="ESRI Shapefile",overwrite_layer=TRUE)
    shapefile(rs, filename="randomsample", overwrite=TRUE)
    
    write.csv(rs,"randomsample.csv")
    
    BaseDatosEstaciones=data.frame(id_station=paste0("station_",1:size),
                                   lon=data.frame(rs)[,2],
                                   lat=data.frame(rs)[,3])
    
    #Here the Stations DataBase ("BaseDatosEstaciones.csv") is saved into the working directory
    write.csv(BaseDatosEstaciones,"BaseDatosEstaciones.csv",row.names = FALSE)
    #write.csv(BaseDatosEstaciones,"BaseDatosEstacionesBackup.csv",row.names = FALSE)
    
    #Records DataBase Generation ("BaseDatosRegistros.csv") for a given country
    ppDB=raster::extract(country.rts,data.frame(rs)[,2:3])
    #ppDB.DF <- data.frame(matrix(unlist(ppDB), nrow=size, byrow=F),stringsAsFactors=FALSE)
    ppDB.DF2 <- t(data.frame(matrix(unlist(ppDB), nrow=size, byrow=F),stringsAsFactors=FALSE))
    ppDB.DF2.ts=ts(ppDB.DF2,start=c(1970,1),end=c(2015,12),frequency = 12)
    
    BaseDatosRegistros=data.frame(matrix(ppDB.DF2.ts[,1],ncol=12,byrow=TRUE))
    #If there is more than one station, the binding is performed. 
    #Otherwise, it is omitted to avoid an out of bounds error.
    if(dim(BaseDatosEstaciones)[1]>1){
      for(i in 2:dim(BaseDatosEstaciones)[1]){
        temp.df=data.frame(matrix(ppDB.DF2.ts[,i],ncol=12,byrow=TRUE))
        BaseDatosRegistros=rbind(BaseDatosRegistros,temp.df)
      }
    } else {
      message("Warning: Only 1 station selected in CRU")
    }
    colnames(BaseDatosRegistros)[1:12]=month.abb
    BaseDatosRegistros$id_station=rep(BaseDatosEstaciones$id_station,each=46)
    BaseDatosRegistros$Year=rep(1970:2015,dim(BaseDatosEstaciones)[1])
    BaseDatosRegistros=BaseDatosRegistros[,c(13,14,1:12)]
    
    
    # Zero values are corrected
    BaseDatosRegistros <- zero_correction(BaseDatosRegistros)
    
    #Here the Records DataBase ("BaseDatosRegistros.csv") is saved into the working directory
    write.csv(BaseDatosRegistros,"BaseDatosRegistros.csv",row.names = FALSE)
    #write.csv(BaseDatosRegistros,"BaseDatosRegistrosBackup.csv",row.names = FALSE)
    #................................................. END ...........................................................
    setwd(workdir)
    
  }
  else if(model=="CHIRPS"){
    
    # BLOCK I.B. DATABASE CONSTRUCTION FROM CHIRPS ------------------
    #AFRICAN DROUGHT ATLAS CHIRPS
    
    # Optional Setup, Country Code, Shape and raster creation ========
    workdir = here()
    setwd(workdir)
    
    # Listado de paises segun codigo ISO
    ISO.codes <- read.csv("CountryISOCodes.csv",sep=";")
    Afr.country.list <- as.character(ISO.codes$ThreeLetter)
    Afr.country.name <- countrycode(Afr.country.list, "iso3c","country.name")
    
    
    shape_return <- shape_country(country)
    Boundaries <- shape_return$Boundaries
    sps <- shape_return$sps
    rm(shape_return)
    
    #BoundariesAFR = raster::shapefile("AfricaDA.shp")
    #raster::projection(BoundariesAFR)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
    
    #Create working directory and set working directory
    dir.create(paste0(workdir,country))
    setwd(paste0(getwd(),"/",country))
    
    shapefile(Boundaries, filename=country, overwrite=TRUE)
    #writeOGR(Boundaries, dsn = '.', layer = country, driver = "ESRI Shapefile",overwrite_layer=TRUE)
    
    #rm(Boundaries)
    #save.image(paste0("C:/Users/jnune/Documents/AtlasSequiaALCCRU/",country,"/",country,".RData"))
    ##                CORREGIR MAPA MANUALMENTE
    #.......................................................................
    
    
    #Boundaries=readOGR(".",country)   #Obtenido de http://www.maplibrary.org/library/stacks/Africa/index.htm
    plot(Boundaries)
    Bound.Pol=extent(Boundaries)
    if(Bound.Pol[1]>0) xmin=Bound.Pol[1]*0.95 else xmin=Bound.Pol[1]*1.05
    if(Bound.Pol[2]>0) xmax=Bound.Pol[2]*1.05 else xmax=Bound.Pol[2]*0.95
    if(Bound.Pol[3]>0) ymin=Bound.Pol[3]*0.95 else ymin=Bound.Pol[3]*1.05
    if(Bound.Pol[4]>0) ymax=Bound.Pol[4]*1.05 else ymax=Bound.Pol[4]*0.95
    xmin=round(xmin,4)
    xmax=round(xmax,4)
    ymin=round(ymin,4)
    ymax=round(ymax,4)
    
    
    #If download fail, then copy paste url in a web browser and download the data.nc file manually
    #After download, put the file into the country's folder.
    # Download netCDF file
    url=paste0('http://iridl.ldeo.columbia.edu/SOURCES/.UCSB/.CHIRPS/.v2p0/.daily-improved/.global/.0p25/.prcp/X/',xmin,'/',xmax,'/RANGEEDGES/Y/',ymin,'/',ymax,'/RANGEEDGES/T/(Jan%201981)/(Dec%202016)/RANGE/weeklytomonthly/data.nc')
    dfile=c("data.nc")
    print(url)
    #Additional time for slow connections and if DL takes time to process
    options(timeout = max(6000, getOption("timeout")))
    
    tryCatch(download.file(url, mode="wb",destfile=paste0(getwd(),"/data.nc")), 
             error = function(e) print(paste(url, 'The data cannot be downloaded from the datalibrary.')))  
    
    # Lectura de datos mediante netCDF ===========
    # 06/07/2017
    # H.Maureira
    
    #Importar datos (Esto debe ser descargado de manera manual desde IRI)
    #datos_ncdf <- nc_open("./input/data_10N_20S.nc")
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
    
    #Quick look
    image(lon,lat,prcp_array[,,1], col=rev(brewer.pal(10,"RdBu")))
    
    # create dataframe -- reshape data
    # matrix (nlon*nlat rows by 2 cols) of lons and lats
    # De esta forma, se obtiene la cantidad de estaciones virtuales en el CUADRADO.
    lonlat <- as.matrix(expand.grid(lon,lat))
    dim(lonlat)
    
    #Transformar el array en un rasterbrick (no es necesario parece)
    #array_brick <- brick(prcp_array, xmn = min(lon), xmx = max(lon), ymn = min(lat), ymx = max(lat))
    
    # Se selecciona mapa del pa?s para seleccionar puntos de interes ==========
    
    
    #Se usan los l?mites del mapa para recortar la base de datos
    mapa<-Boundaries
    
    #Crear shapefile de puntos a partir de lonlat (corresponde a todas las estaciones virtuales, o sea, el rect?ngulo)
    est_cuadro <- data.frame(lonlat)
    est_cuadro$ID <- paste("id_", 1:nrow(est_cuadro), sep = "")
    #est_cuadro <- SpatialPointsDataFrame(coords = lonlatDF[,1:2], proj4string = proj4string(mapa), bbox = lonlatDF)
    coordinates(est_cuadro) <- ~ Var1 + Var2
    proj4string(est_cuadro) <- proj4string(mapa)
    
    #Se edita el archivo lonlat para seleccionar aquellos puntos que est?n dentro del l?mite pa?s
    est_selec <- est_cuadro[mapa,] #perfecto, se guarda
    
    #Se obtiene un archivo (dataframe) lonlat con los puntos seleccionados
    est_selecDF <- data.frame(est_selec)
    
    #Se crea un data frame que almacenar? en cada lista, los valores por d?a
    DF_est <- data.frame(date = tpo2, prcp_chirps = NA)
    
    #Se crea una lista cuyos nombres corresponde a lat long (que en el fondo son las estaciones)
    lista_est <- list()
    for (i in 1:nrow(est_selecDF)){ #Continuar desde ac?...
      
      print(paste0("estacion ",i," de ",nrow(est_selecDF)))
      
      lista_est[[i]] <- DF_est
      names(lista_est)[i] <- as.character(est_selecDF[i,3])
      
    }
    
    #........................................
    # Extraer valores desde el prcpc_array usando "est_selecDF"
    # almacenar en la lista
    
    l = 1
    for (j in 1:dim(prcp_array)[3]){ #ac? no es array brick, OJO
      print(paste0("mes ",j," de ",dim(prcp_array)[3]))
      
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
    
    for (m in 1:length(lista_est)){
      print(paste0("estacion ", m, "de ",length(lista_est)))
      
      export <- lista_est[[m]]
      write.table(export, paste0(getwd(),"/est_procesadas/", names(lista_est)[m], "_latam.txt"), sep =";", row.names = FALSE)
    }
    
    
    # Se exporta la metadata (Estaciones virtuales)
    metaexp <- est_selecDF[,c(3,1,2)]
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
    
    colnames(BaseDatosRegistros)[3:14]=month.abb
    
    # CORRECCION DE CEROS. TODOS LOS CEROS SON REEMPLAZADOS POR UN VALOR RANDOM UNIF ENTRE 0 Y 1
    
    BaseDatosRegistros <- zero_correction(BaseDatosRegistros)
    
    #row.names(BaseRegistrosPr)=NULL
    write.csv(BaseDatosRegistros,"BaseDatosRegistros.csv",row.names=FALSE)
    #write.csv(BaseDatosRegistros,"BaseDarosRegistrosBackup.csv",sep=",",row.names=FALSE)
    
    #BaseDatosEstaciones. Defino el tama?o total a 500
    BaseDatosEstaciones=read.csv(paste0(getwd(),"/est_procesadas/metadata.txt"),sep=";",header=TRUE)
    colnames(BaseDatosEstaciones)[1]="id_station"
    row.names(BaseDatosEstaciones)=NULL
    if (dim(BaseDatosEstaciones)[1]>500){
      BaseDatosEstaciones=randomSample(BaseDatosEstaciones,500)
    } else {
      BaseDatosEstaciones=BaseDatosEstaciones
    }
    
    write.csv(BaseDatosEstaciones,"BaseDatosEstaciones.csv", row.names=FALSE)
    #write.csv(BaseDatosEstaciones,"BaseDatosEstacionesBackup.csv", row.names=FALSE)
    #..........................................................END ...........................................
    setwd(workdir)
  }
  else{
    message("Error selecting the model. Only the options CRU and CHIRPS are available.")
  }
}






