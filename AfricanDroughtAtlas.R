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
#####
#import libraries
library(raster);library(rgdal);library(countrycode);library(rts);require(ncdf4);library(HelpersMG)
library(plyr);library(latticeExtra);library(reshape2);library(gdata);library(corrplot);library(sqldf)
library(zoo);library(Kendall);library(zyp);library(car);library(gtools);# Para mantener orden texto-numero en id_station
library(rgeos);library(lmom);library(lmomRFA);library(sp);library(rrcov);library(nsRFA);library(ModelMap)
library(maptools);library(stringr);library(raster);library(rasterVis);library(hydroGOF);library(randomForest);library(progress);#Check proper installation of SAG-GIS and RSAGA
library(RSAGA);library(gtools)

#####
#config
workdir = "C:/Users/pablo/OneDrive/Escritorio/CAZALAC/ADA/"
country="DJI" #PLEASE REPLACE WITH CORRESPONDING THREE ISO LETTER . IN THIS CASE BOTZWANA=BWA
#THREE LETTER CODES CAN BE FOUND IN "CountryISOCodes.csv"

#####
#Function
setwd(workdir)
source('DroughAtlasFunctions.R')

#####
#Check
rsaga.env()
rsaga.get.version()
setwd(paste0(workdir, country)) #changed slashes
getwd()

###################################################################################################################
#                                    BLOCK I.A. DATABASE CONSTRUCTION FROM CRU 3.21 or CHIRPS                     #
#                                    CHOOSE ONE OF THE TWO ALTERNATIVES BASED ON COUNTRY SIZE                     #
###################################################################################################################



##################################################################################################################
#                                         BLOCK II. VARIABLES AND INDICES  CALCULATION                           #
##################################################################################################################



Boundaries=readOGR(".",country) #Obtained from http://www.maplibrary.org/library/stacks/Africa/index.htm
projection(Boundaries)="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
plot(Boundaries)

#BaseDatosEstaciones with 3 columns: factor, num, num
BaseDatosEstaciones <- read.csv("BaseDatosEstaciones.csv",sep=",")
head(BaseDatosEstaciones,3)
BaseDatosRegistros <- read.csv("BaseDatosRegistros.csv",sep=",")
head(BaseDatosRegistros,3)


  #BaseDatosRegistros with 14 columns:factor, int, 12 num
  BaseRegistrosPr=BaseDatosRegistros
  
  # Monthly accumulated precipitation values 
  RegionDF<-BaseRegistrosPr
  RegionDF[,1] <- factor(RegionDF[,1]) #Defino solamente los niveles del DF, ya que el subset mantiene todos los anteriores.
  
  #SIMULACION
  #RegionDF <- RegionDFSim
  #RegionDF[,1] <- factor(RegionDF[,1])
  
  #Almacenamiento de Estaciones pertenecientes a la region, en una lista
  #(Paso Necesario, ya que se cambiara el formato a Long)
  #---------------------------------------------------------------------
  
  EstacionesWide <- list()
  #levels(reorder(RegionDF[,1]))
  for (i in mixedsort(levels(RegionDF[,1]))){
    EstacionesWider <- subset(RegionDF, id_station == paste(i))
    EstacionesWide[[i]] <- EstacionesWider
    rm(EstacionesWider)
  }
  
  #Almacenamiento de Estaciones en formato Long y en formato necesario para c?lculo de H1
  #--------------------------------------------------------------------------------------
  EstacionesLong <- list()
  #DatosH1 <- list()
  
  for (j in 1:length(EstacionesWide)){
    dfm <- EstacionesWide[[j]]
    #De wide a long 
    # Conveertir data.frame de formato ancho a largo, se resta 1 para eliminar columna "Region"
    dfm <- melt(dfm[,2:(ncol(dfm))], id.var="Year", variable_name="Month") #Dejar ANTES COMO YEAR #melt(dfm[,2:(ncol(dfm)-1)]
    names(dfm) <- c("Year", "Month", "Lluvia")
    # dfm$Mointh is factor, Convert to month number & calc yr_frac
    lluvia_mo_num <- unclass(dfm$Month)
    lluvia_mo_frac <- as.numeric(lluvia_mo_num/12 )
    yr_frac <- dfm$Year + lluvia_mo_frac
    # build consolidated data.frame
    dfm <- data.frame(yr_frac, dfm)
    dfm <- dfm[order(dfm$yr_frac), ]
    #dfm <- dfm[!is.na(dfm$Lluvia),]
    dfm <-dfm[2:4]
    dfm[,2] <- as.character(dfm[,2])
    
    EstacionesLong[[j]] <- dfm
    #DatosH1[[j]] <- dfm[,3]
  }
  
  #Sumatoria de cada Estacion, que sera almacenado en otra lista
  #---------------------------------------------------------
  k <- 12 #This is a critical values because defines the time step to be analyzed latter
  
  #Creacion de la lista donde se va a almacenar la suma.
  EstacionesSuma <- list()
  
  #Sumatorias, en base al valor de k, que se almacena en EstacionesSuma.
  for (l in 1:length(EstacionesLong)){#Selecciono la Estacion en Long
    
    datospp <- EstacionesLong[[l]][,3] #almaceno los valores de pp que seran sumados despues
    
    EstacionSumar <- EstacionesLong[[l]] #se escoge la estacion en formato long
    EstacionSumar[,3] <- NA #Se asigna columna de valores de pp como NA.
    colnames(EstacionSumar)[3] <- names(EstacionesWide)[[l]]    #"Suma" 
    
    for (m in 1:nrow(EstacionSumar)){ #aca es donde se empieza a determinar si la sumatoria se hace o no.
      FilaFinal <- m
      FilaInicial <- (FilaFinal-(k-1))
      
      if (FilaInicial<=0){
        EstacionSumar[m,3]<-NA
      } else{
        EstacionSumar[m,3] <- sum(datospp[FilaInicial:FilaFinal])
      }
    }
    EstacionesSuma[[l]] <- EstacionSumar #finalmente, se almacena la estacion sumada en una lista.
  }
  
  #Same names from "EstacionesWide" are assigned to "EstacionesSuma"
  names(EstacionesSuma) <- names(EstacionesWide)
  
  ### JNunez request (bind BaseRegistrosPr w/ monthly sums) ---
  BaseSums <- list()
  for (u in 1:12){
    
    BaseSums[[u]] <- sapply(EstacionesSuma, function(x) subset(x, Month == month.abb[u], select = colnames(x)[3]))
    BaseSums[[u]] <- unlist(BaseSums[[u]])
    names(BaseSums[[u]]) <- NULL
    
    names(BaseSums)[u] <- paste(month.abb[u], "_", k, "_months", sep = "")
    
  }
  
  BaseRegistrosPr_sumas <- data.frame(matrix(unlist(BaseSums), nrow = nrow(BaseRegistrosPr), byrow=F))
  colnames(BaseRegistrosPr_sumas) <- paste(month.abb, "_", k, sep = "")
  
  ## bind
  BaseRegistrosPr <- cbind(BaseRegistrosPr, BaseRegistrosPr_sumas)
  CumSum=data.frame(matrix(ncol = 12,nrow=dim(BaseRegistrosPr)[1]))
  for (i in 1:dim(BaseRegistrosPr)[1]){
    CumSum[i,1:12]=cumsum(as.numeric(BaseRegistrosPr[i,3:14]))
  }
  colnames(CumSum)=paste0("CumSum",month.abb)
  
  BaseRegistrosPr=cbind(BaseRegistrosPr,CumSum)
  write.csv(BaseRegistrosPr,"BaseRegistrosPr.csv",row.names=FALSE)
  
  #.....................................................................................................................
  #                   Calculation of several Indices like Julian Mean Day, Seasonality Index and Modified Fourier Index
  
  JanDec=matrix(BaseRegistrosPr$CumSumDec)
  DiaJulianoAng<-matrix(seq(15,345,30)*2*pi/365,nrow=length(BaseRegistrosPr[[1]]),ncol=12,byrow=TRUE)
  Prec<-BaseRegistrosPr[match(month.abb,names(BaseRegistrosPr))]
  x<-Prec*cos(DiaJulianoAng)/JanDec
  y<-Prec*sin(DiaJulianoAng)/JanDec
  xcos<-matrix(rowMeans(x),nrow=length(BaseRegistrosPr[[1]]),ncol=1)
  ysin<-matrix(rowMeans(y),nrow=length(BaseRegistrosPr[[1]]),ncol=1)
  xcossum<-matrix(rowSums(x),nrow=length(BaseRegistrosPr[[1]]),ncol=1)
  ysinsum<-matrix(rowSums(y),nrow=length(BaseRegistrosPr[[1]]),ncol=1)
  angulo<-atan(ysin/xcos)
  angulo_corregido<-matrix(0,nrow=length(BaseRegistrosPr[[1]]),ncol=1)
  for (k in 1:length(BaseRegistrosPr[[1]])) {
    if (sum(is.na(xcos[k]>0))|sum(is.na(ysin[k]>0))) angulo_corregido[k]<-NA else
      if (xcos[k]>0&ysin[k]>0)  angulo_corregido[k]<-angulo[k] else if (xcos[k]<0&ysin[k]<0) angulo_corregido[k]<-angulo[k]+pi  else if (ysin[k]>0&xcos[k]<0) angulo_corregido[k]<-angulo[k]+pi else   angulo_corregido[k]<-angulo[k]+2*pi 
  }
  JMD<-(angulo_corregido*365)/(2*pi)
  SI<-sqrt(xcossum^2+ysinsum^2)
  IFM<-matrix(rowSums((Prec*Prec)/JanDec),nrow=length(BaseRegistrosPr[[1]]),ncol=1) # Esta linea realiza el calculo del Indice de Fourier Modificado
  
  # C.5.3. Aggregate record indices and add to Stations's DataBase
  BaseDatosIntermedia<-cbind(BaseRegistrosPr,SI,JMD,IFM)# Uno las columnas con los Indices calculados para cada registro. Contiene los valores anuales.
  write.csv(BaseDatosIntermedia,"BaseIntermedia.csv")
  rm(BaseDatosIntermedia)
  BaseDatosIntermedia <- read.csv("BaseIntermedia.csv",sep=",")
  BaseDatosIntermedia$id_station=as.factor(BaseDatosIntermedia$id_station) #se habilita esta parte. HMC
  
  # C.5.4. Derivated variables
  # C.5.4.1.  Mean Annual Precipitation
  PMA_por_Estacion<-as.matrix(tapply(BaseDatosIntermedia$CumSumDec,BaseDatosIntermedia$id_station,mean,na.rm=TRUE))
  
  # C.5.4.2. Seasonality Index (SI)
  SI_por_Estacion<-as.matrix(tapply(BaseDatosIntermedia$SI,BaseDatosIntermedia$id_station,mean,na.rm=TRUE))
  
  # C.5.4.3. JUlian Mean Day (JMD)
  JMD_por_Estacion<-as.matrix(tapply(BaseDatosIntermedia$JMD,BaseDatosIntermedia$id_station,mean,na.rm=TRUE))
  
  # C.5.4.4. Modified Fourier Index (IFM)
  IFM_por_Estacion<-as.matrix(tapply(BaseDatosIntermedia$IFM,BaseDatosIntermedia$id_station,mean,na.rm=TRUE))
  
  # C.5.4.5. Record Length (LongRec)
  Longitud<-function(x) (length(x)-sum(is.na(x)))
  LR_por_Estacion<-as.matrix(tapply(BaseDatosIntermedia$CumSumDec,BaseDatosIntermedia$id_station,Longitud))
  
  # C.5.4.6. Firt year of record (FirstYear)
  PrimerAnio_por_Estacion<-as.matrix(tapply(BaseDatosIntermedia$Year,BaseDatosIntermedia$id_station,min))
  
  # C.5.4.7. Last year of record (LastYear)
  UltimoAnio_por_Estacion<-as.matrix(tapply(BaseDatosIntermedia$Year,BaseDatosIntermedia$id_station,max))
  
  
  # C.5.5. Add indices to Stations DataBase and remove objects from memory
  id_station<-levels(BaseDatosIntermedia$id_station)
  BaseDatosIndices<-data.frame(id_station=id_station,MAP=PMA_por_Estacion, SI_Medio=SI_por_Estacion,
                               JMD_Medio=JMD_por_Estacion, IFM_Medio=IFM_por_Estacion, RL_Station=LR_por_Estacion, FirstYear=PrimerAnio_por_Estacion,
                               LastYear=UltimoAnio_por_Estacion) #Ok, está arreglado habilitando la opción de usar factores.
  
  write.csv(BaseDatosIndices,"BaseDatosIndices.csv",row.names=FALSE)
  rm(BaseDatosIndices)
  BaseDatosIndices<-read.csv("BaseDatosIndices.csv",header=T,sep=",")



#Please check correct import. Dimensions and structure or both data.frames
#Please check that names under id_station are similar in both data.frames

BaseDatosEstaciones<-join(BaseDatosEstaciones,BaseDatosIndices,by="id_station") # Ac? uno las dos bases de datos a nivel de estaciones solamente, no de todos los registros. Contiene los valores medios

#Here an updated version of Stations DataBase is saved
write.csv(BaseDatosEstaciones,"BaseDatosEstaciones.csv",row.names=FALSE)
remove(BaseDatosEstaciones)
BaseDatosEstaciones<-read.csv("BaseDatosEstaciones.csv",header=T,sep=",")
head(BaseDatosEstaciones,5)

#...................................................................................................................
#                                     Trend, serial correlation, Mann-Kendall and Sen slope calculation




  
  RegionTotalS<-sqldf(paste("SELECT id_station, Year, CumSumDec FROM BaseDatosEstaciones join BaseDatosIntermedia USING(id_station)",sep=""))
  #RegionTotalS<-sqldf("select id_station, Year, JanDec from BaseCompleta where RL_Station>=15")
  RegionTotalS_dat<-RegionTotalS[c("CumSumDec","Year")][,]
  RegionTotalS_fac<-factor(RegionTotalS["id_station"][,])
  RegTotalS<-split(RegionTotalS_dat,RegionTotalS_fac)#
  estaciones<-length(RegTotalS)
  
  #....................................................................................
  for (est in 1:estaciones){
    
    Y<-RegTotalS[[est]]$CumSumDec/mean(RegTotalS[[est]]$CumSumDec,na.rm=T)
    X<-RegTotalS[[est]]$Year
    XY<-cbind(X,Y)
    orderX=order(X)
    XY<-data.frame(XY[orderX,])
    z <- read.zoo(XY)
    gz <- zoo( x = NULL ,seq(min(time(z)), max(time(z))))
    XYzoo<-merge(z, gz)
    XYdf<-data.frame(X=index(XYzoo),Y=XYzoo,row.names=NULL)
    
    mk<-as.numeric(MannKendall(as.ts(XYzoo)))[2]# Paquete Kendall
    ss<-zyp.sen(Y~X,data=XYdf)[[1]][2] # Paquete zyp
    
    fit=lm(Y~X,data=XYdf)
    DW<-durbinWatsonTest(fit)#paquete car
    
    lm.coeficients <- function (modelobject) {
      if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
      f <- summary(modelobject)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      s<-summary(modelobject)$coefficients[2]
      attributes(s) <- NULL
      attributes(p) <- NULL
      salida<-list(slope=s,Pvalue=p)
      return(salida)
    }
    RegTotalS[[est]]$Pendiente<-lm.coeficients(fit)$slope
    RegTotalS[[est]]$Pvalue<-lm.coeficients(fit)$Pvalue
    RegTotalS[[est]]$Autocorr<-DW[[1]]
    RegTotalS[[est]]$DW_Pvalue<-DW[[3]]
    RegTotalS[[est]]$Sen<-ss
    RegTotalS[[est]]$MKPvalue<-mk
  }
  par(mfrow=c(1,2))
  RegTotalS2<-ldply(RegTotalS, data.frame)
  colnames(RegTotalS2)[1]<-'id_station'
  Pendiente_por_Estacion<-as.matrix(tapply(RegTotalS2[[4]],RegTotalS2[[1]],mean))
  hist(Pendiente_por_Estacion)
  boxplot(Pendiente_por_Estacion)
  Pvalue_por_Estacion<-as.matrix(tapply(RegTotalS2[[5]],RegTotalS2[[1]],mean))
  hist(Pvalue_por_Estacion)
  boxplot(Pvalue_por_Estacion)
  Autocor_por_Estacion<-as.matrix(tapply(RegTotalS2[[6]],RegTotalS2[[1]],mean))
  hist(Autocor_por_Estacion)
  boxplot(Autocor_por_Estacion)
  DW_por_Estacion<-as.matrix(tapply(RegTotalS2[[7]],RegTotalS2[[1]],mean))
  hist(DW_por_Estacion)
  boxplot(DW_por_Estacion)
  Sen_por_Estacion<-as.matrix(tapply(RegTotalS2[[8]],RegTotalS2[[1]],mean))
  hist(Sen_por_Estacion)
  boxplot(Sen_por_Estacion)
  MKPvalue_por_Estacion<-as.matrix(tapply(RegTotalS2[[9]],RegTotalS2[[1]],mean))
  hist(MKPvalue_por_Estacion)
  boxplot(MKPvalue_por_Estacion)
  par(mfrow=c(1,1))
  
  id_estaciones<-levels(as.factor(RegTotalS2[,1]))
  BaseDatosTendenciaDW<-cbind(id_estaciones,Pendiente_por_Estacion,Pvalue_por_Estacion,Autocor_por_Estacion,
                              DW_por_Estacion,Sen_por_Estacion,MKPvalue_por_Estacion)
  colnames(BaseDatosTendenciaDW)[c(1,2,3,4,5,6,7)]<-c('id_station','Pendiente','Pvalue','Autocorr','DW_Pvalue',
                                                      'Sen','MK')
  write.csv(BaseDatosTendenciaDW, "BaseDatosTendenciaDW.csv", row.names=FALSE)
  rm(BaseDatosTendenciaDW) # Aca debo remover y volver a cargar la  Base de Datos en el sistema
  BaseDatosTendenciaDW <- read.csv("BaseDatosTendenciaDW.csv",sep=",",header=T)# Y vuelvo a cargar la Base de Datos
  
  # Here an updated version of Stations DataBase is saved
  BaseDatosEstaciones<-join(BaseDatosEstaciones,BaseDatosTendenciaDW,by="id_station") # Ac? uno las dos bases de datos a nivel de estaciones solamente.
  write.csv(BaseDatosEstaciones,  "BaseDatosEstaciones.csv", row.names=FALSE)
  rm(BaseDatosEstaciones)
  BaseDatosEstaciones<-read.csv("BaseDatosEstaciones.csv",header=T, sep=",")
  head(BaseDatosEstaciones,5) #Hasta el final del Block II, el script se ejecuta de forma normal, salvo pequeños ajustes.
  #.................................................END OF BLOCK II..........................................................
  



###############################################################################################################
#                           BLOCK III. REGIONALIZATION USING MINIMUM DISTANCE SEARCHING  ALGORITHM            #
###############################################################################################################

# Minimum Distance Searching Algorithm

BaseCompletaAdapt=merge(BaseDatosEstaciones,BaseRegistrosPr, by.x = "id_station", by.y = "id_station") # Ac? uno la Base de Datos de Estaciones con Indices Medios y  la Base de Datos de Registros
write.csv(BaseCompletaAdapt, "BaseCompletaAdapt.csv", row.names=FALSE)
rm(BaseCompletaAdapt) # Aca debo remover y volver a cargar la  Base de Datos en el sistema
#BaseCompletaAdapt <- read.csv("BaseCompletaAdapt.csv",header=T, sep=",")# Y vuelvo a cargar la Base de Datos. #CAZ: Probablemente busca factores en col id
BaseCompletaAdapt <- read.csv("BaseCompletaAdapt.csv",header=T, sep=",", stringsAsFactors = TRUE) #adaptación CAZ

VarInter="CumSumDec" ################Indicate the name of the variable to be analyzed
Hc=2 # Maximum Acepted Heterogeneity Measure
Hf=0.95 # H1 factor
H1_critical=Hc*Hf
spBEPr_adapt=BaseDatosEstaciones
coordinates(spBEPr_adapt) <- ~lon+lat
d_adapt <- gDistance(spBEPr_adapt, byid=T)
row.names(d_adapt)=as.character(BaseDatosEstaciones$id_station)
colnames(d_adapt)=as.character(BaseDatosEstaciones$id_station)

RegStations=data.frame(id_station=NA,ClustReg_adapt=NA,Criteria=NA)
counter=1
Skips=3
MaxRS=21 #Maximum number of stations per region allowed

#cuando d_adapt[1,] es 0, genera un error de dimensiones, por lo que se rompe el ciclo abajo
while(TRUE){ #Acá aparece el error de library "homtest": Solucionado
  REG_adapt=list()# Almacena el nombre de las estaciones disponibles en cada iteracion
  BaseReg_adapt=list()#Almacena los registros de las estaciones para calcular H1
  H1s=0
  #H1sorig=0
  
  for (i in 1:length(d_adapt[1,])){
    print(i)
    crit=4# Predeterminado. La region se cierra cuando se alcanza el maximo.
    REG_adapt[[i]]=row.names(d_adapt)[order(d_adapt[1,],decreasing=F)[i]]
    Region_adapt<-sqldf(paste("select id_station,",VarInter," from BaseCompletaAdapt where id_station=='",REG_adapt[[i]],"'",sep=""))
    Region_adapt_dat<-Region_adapt[VarInter][,]
    Region_adapt_fac<-factor(Region_adapt["id_station"][,])
    Reg_adapt<-split(Region_adapt_dat,Region_adapt_fac)# Con esto separo los registros segun la estacion
    Reg_adapt=lapply(Reg_adapt, function(x) {x[x!=0]})#Agregado para analizar solo sin ceros
    if(length(Reg_adapt[[1]])==0){Reg_adapt[[1]][1:5]=rand.GEV(5, 0.6986,0.2860,-0.3294)} else 
      if(length(Reg_adapt[[1]])==1){Reg_adapt[[1]][2:5]=rand.GEV(4, 0.6986,0.2860,-0.3294)} else
        if(length(Reg_adapt[[1]])==2){Reg_adapt[[1]][3:5]=rand.GEV(3, 0.6986,0.2860,-0.3294)} else
          if(length(Reg_adapt[[1]])==3){Reg_adapt[[1]][4:5]=rand.GEV(2, 0.6986,0.2860,-0.3294)} else
            Reg_adapt=Reg_adapt
    BaseReg_adapt[i]=Reg_adapt
    names(BaseReg_adapt)[i]=names(Reg_adapt)
    #SumSt<-regsamlmu(BaseReg_adapt)# Calculate sample L-moments de la regi?n de testeo
    #SumStRegData<-as.regdata(SumSt)# Change format
    #Rlm<-regavlmom(SumStRegData)#Calculate regional L-moments
    #HW_H1orig<-regtst(SumStRegData, nsim=500)
    HW_H1<-H1ZDi(BaseReg_adapt, Nsim=1000)#H1 calculado con nueva funci?n
    #H1sorig[i+1]=HW_H1orig$H[1]
    H1s[i+1]=HW_H1$H[1]
    #H1Record[i+1]=HW_H1$H[1]# INVENTAR UN H1 PARA i-1<2.5. Ir almacenando los H1s y 
    #Preguntar si el actual y el anterior fueron homogeneos. Solamente cuando
    #dos consecutivas tengan heterogeneo, ahi parar. Eso permite incrementar el tama?o de region
    #para evitar rechazos alternados, se pueden almacenar la cantidad de rechazos y para cuando sean 2
    if(is.na(HW_H1$H[1]>H1_critical)==TRUE) {crit=1;break}# Si se produce na entonces break por criterio 1
    if((HW_H1$H[1]>H1_critical & H1s[i]<H1_critical & sum(H1s>H1_critical)>=Skips)==TRUE) {crit=2;break}# Si se cumplen los skip maximos
    if((HW_H1$H[1]>H1_critical & H1s[i]>H1_critical)==TRUE) {crit=3; break}# Si actual y anterior fueron criticos. Dos cr?ticos consecutivos
    if(length(REG_adapt)>MaxRS) {break}# Cierro region cuando es muy grande
  }
  if(crit==1) { RSt=unlist(REG_adapt)[-i]} else
    if(crit==2) {RSt=unlist(REG_adapt)[-i]} else 
      if (crit==3) {RSt=unlist(REG_adapt)[c(-i,-(i-1))]} else {RSt=unlist(REG_adapt)}# Corregir aca para que en la ?ltima iteraci?n no elimine a la estaci?n
  #Puede ser un d_adaptint que se evalue si sigue teniendo valores mayores a cero o no
  RegSttemp=data.frame(id_station=RSt,ClustReg_adapt=counter,Criteria=crit)
  RegStations=rbind(RegStations,RegSttemp)
  #ya no queda nada en el df
  if(!(dim(d_adapt)[1] == 0)){
   break; 
  }
  d_adapt=d_adapt[-grep(paste0(RSt,"$",collapse="|"),row.names(d_adapt)),-grep(paste0(RSt,"$",collapse="|"),colnames(d_adapt))]# Secreto es adicionar signo $
  counter=sum(counter,1)
} #Al parecer, el error puede bypasearse, dado que el enfoque de while no está bien implementado.
#
RegStations=(RegStations[-1,])#CAZ: Perfecto, ya que primera fila es NA
#RegStations[dim(BaseEstacionesPrIII)[1],1]=unlist(REG_adapt)[i]# Con esto corrijo la p?rdida de la ?ltima estaci?n
#RegStations[dim(BaseEstacionesPrIII)[1],2]=counter# COn esto corrijo la p?rdida de la ?ltima estaci?n
RegStations$id_station=as.factor(RegStations$id_station)
RegStations$ClustReg_adapt=as.integer(RegStations$ClustReg_adapt)
pie(table(RegStations$Criteria))

BaseDatosEstacionesClust=join(BaseDatosEstaciones,RegStations,by="id_station")
write.csv(BaseDatosEstacionesClust,"BaseDatosEstacionesClust.csv",row.names=FALSE)
rm(BaseDatosEstacionesClust)
BaseDatosEstacionesClust=read.csv("BaseDatosEstacionesClust.csv",sep=",") #Al parecer, stringsasfactors antes era true POR DEFECTO.
BaseDatosEstacionesClust=read.csv("BaseDatosEstacionesClust.csv",sep=",", stringsAsFactors = TRUE)
head(BaseDatosEstacionesClust,5)
# D.1.8. Fix database names
NombreClusters=colnames(BaseDatosEstacionesClust[grep('^Clu',names(BaseDatosEstacionesClust))]) #acá existe un error al parecer
#NombreClusters <- as.numeric(levels(as.factor(BaseDatosEstacionesClust$ClustReg_adapt))) #ESTA ES UNA ADAPTACIÓN DE CAZALAC, Parece que Jorge estaba ok


# Calculo de la cantidad de clusters en cada Tratamiento y si coinciden con el rango de clusters.
ClustLevels=length(levels(as.factor(BaseDatosEstacionesClust$ClustReg_adapt)))
#...........................................END OF BLOCK III.....................................................


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


#################################################################################################################
#                                           BLOCK V: L-MOMENTS BASED REGIONAL FREQUENCY ANALYSIS                #
#################################################################################################################

# E.1. Create homogeneous regions
#Incluir linea con ifelse que evalue la dimension de NombresCluster y asigne MaxClusters
#VarInter="CumSumDec"
MaxClusters=ClustLevels
#apply(BaseDatosEstaciones[,grep('^Clu',names(BaseDatosEstaciones))],2,FUN=max)#Mas de un cluster
#MaxClusters=max(BaseEstacionesPrIII$ClustReg_adapt)#En caso uso solamente ClustRegAdapt
#BaseDatosIntermedia <- read.csv("BaseIntermedia.csv",sep=";")
#BaseCompletaIII=merge(BaseEstacionesPrIII,BaseDatosIntermedia, by.x = "id_station", by.y = "id_station") # Ac? uno la Base de Datos de Estaciones con Indices Medios y  la Base de Datos de Registros
#write.table(BaseCompletaIII, paste("BaseCompletaIII",ISO_Country,".csv",sep=""),sep=";", row.names=FALSE)
#rm(BaseCompletaIII) # Aca debo remover y volver a cargar la  Base de Datos en el sistema
#BaseCompletaIII <- read.table(paste("BaseCompletaIII",ISO_Country,".csv",sep=""),header=T, sep=";")# Y vuelvo a cargar la Base de Datos

BaseResumenCluster=data.frame(H1=NA,DistOptima=NA,abs_z_value=NA,
                              NDiscSites=NA,Sitios=NA,Cluster=NA)#Almaceno todas las TablasResumen
BaseSummaryStatistics=list()
Baserfitdist=list()
Baserfitpara=list()
BaseBaseRegiones=list()
BaseProporCeros=list()
Basep0bias=list()
BaseMediaCompleta=list()

for (k in 1:length(NombreClusters)){# Start for loop
  print(k)
  BaseRegiones=list()
  ProporCeros=list()
  p0bias=list()
  MediaCompleta=list()
  NombreRegiones=character()
  for(l in 1:MaxClusters[k]){
    Region<-sqldf(paste("SELECT id_station, ", VarInter," FROM BaseDatosEstacionesClust join BaseRegistrosPr USING(id_station) where ",NombreClusters[k],"==",l,sep="")) #NombreClusters k es ClustReg_adapt
    #Region<-sqldf(paste("select id_station,JanDec from BaseCompletaIII where ",NombreClusters[k],"==",l,sep=""))
    #Region_dat<-Region["JanDec"][,]
    #Region_fac<-factor(Region["id_station"][,])
    Region_dat<-Region[,2]
    Region_fac<-factor(Region[,1])
    Reg<-split(Region_dat,Region_fac)# Con esto separo los registros segun la estacion
    Reg<-lapply(Reg, function(x) if (sum(x==0)>length(x)-4) {x[1:4]=rand.GEV(4, 0.58,0.29,-0.46);x} else {x})
    ProporCeros[[l]]=lapply(Reg, function(x) {sum(x==0)/(length(x)+1)})
    p0bias[[l]]=lapply(Reg, function(x) {(sum(x==0)+1)/(2*(length(x)+1))})
    MediaCompleta[[l]]=lapply(Reg, mean)
    RegsinCero=lapply(Reg, function(x) {x[x!=0]})
    BaseRegiones[[l]]=RegsinCero
    
    NombreRegiones[l]=paste("Region ",NombreClusters[k],l,sep=".")
  }
  
  # E.2. Test missing stations
  N=c()
  for (h in 1:length(BaseRegiones)){
    N[h]=length(BaseRegiones[[h]])}
  TotalEstaciones<-sum(N)
  
  print(paste("The total number of analyzed stations is: ",TotalEstaciones,sep=""))
  print (paste(" from ",length(levels(BaseDatosEstacionesClust$id_station)),"  available stations in BaseCompletaIII",sep=""))
  
  # E.3. Regionalization
  Regiones<-length(BaseRegiones)
  
  MaxEstpR<-function(x) { 
    N<-c()
    for (i in 1:length(x)){
      N[i]=length(x[[i]])} 
    Nmax=max(N)
    return(Nmax)}

  MaxEst<-MaxEstpR(BaseRegiones)
  ResultadosSummaryStatistics<-array(0,dim=c(MaxEst,7,Regiones))#MaxEst=Cantidad de estaciones de la region que mas estaciones tiene
  ResultadosSummaryStatisticsRegData<-array(0,dim=c(MaxEst,7,Regiones))#MaxEst=Cantidad de estaciones de la regi?n que m?s estaciones tiene
  ResultadosRlmoments<-array(0,dim=c(5,Regiones))#5= L-momentos regionales
  ResultadosARFD<-array(0,dim=c(MaxEst,Regiones))#MaxEst=Cantidad estaciones m?xima por regi?n
  ResultadosARFDRobust<-array(0,dim=c(MaxEst,Regiones))#MaxEst=Cantidad estaciones m?xima por regi?n
  ResultadosARFH<-array(0,dim=c(3,Regiones))# 3= Indices de homogeneidad H1,H2,H3
  row.names(ResultadosARFH)<-c("H1","H2","H3")
  ResultadosARFZ<-array(0,dim=c(6,Regiones))# 6= numero de modelos de probabilidad a calcular su Bondad de Ajuste (glo, gev, gno, pe3, gpa)
  row.names(ResultadosARFZ)<-c("glo", "gev", "gno", "pe3", "gpa","gau")
  Resultadosrfitdist<-array(0,dim=c(1,Regiones))# 1=Un solo ajuste por regi?n
  Resultadosrfitpara<-array(0,dim=c(5,Regiones))#5= cantidad de parametros de la Wakeby
  ResultadosRegionalQuantiles<-array(0,dim=c(19,Regiones))# 19=M?xima cantidad de cuantiles a calcular
  ResultadosRMAP<-array(0,dim=c(1,Regiones))# 1= Un solo valor de precipitacion media anual por region
  SitiosenRegion<-as.data.frame(array(0,dim=c(MaxEst,Regiones))) #Para almacenar el nombre de los sitios que quedan en cada regi?n
  TablaResumen_Regiones<-as.data.frame(array(0,dim=c(Regiones,5)))
  names(TablaResumen_Regiones)<-c("H1","DistOptima","abs_z_value","NDiscSites","Sitios")
  rownames(TablaResumen_Regiones)<-NombreRegiones
  
  #.......Save plots to a PDF
  pdf(paste0("LMRFA_",country,".pdf"), width = 16 , height = 10, title = paste0("L-moment Ratio Diagram for Region ",z))
  par(mfrow=c(1,1))
  for (z in 1:Regiones) { 
    par(mfrow=c(1,1))
    SummaryStatistics<-regsamlmu(BaseRegiones[[z]])# Calculate sample L-moments #EN BASEreGIONES APARECE LE PROBLEMA
    SummaryStatisticsRegData<-as.regdata(SummaryStatistics)# Change format
    Rlmoments<-regavlmom(SummaryStatisticsRegData)#Calculate regional L-moments
    #ARF<-regtst(SummaryStatisticsRegData, nsim=1000)# Calculate region statistics
    ARF<-H1ZDi(BaseRegiones[[z]],Nsim=1000)#Elaborada especialmente por mi basado en homtest
    a<-length(BaseRegiones[[z]])
    ResultadosSummaryStatistics[1:a,1:7,z]<-as.matrix(SummaryStatistics) #Se supone que aca almaceno todos los L-momentos
    ResultadosRlmoments[1:5,z]<-Rlmoments
    ResultadosARFD[1:a,z]<-ARF$Di # Se almacenan las medidas de discordancia
    ResultadosARFDRobust[1:a,z]<-ARF$rDi
    ResultadosARFH[1:3,z]<-ARF$H # Se almacenan las medidas de homogeneidad
    ResultadosARFZ[1:6,z]<-ARF$Z # Se almacenan las medidas de bondad de ajuste
    SitiosenRegion[1:a,z]<-SummaryStatistics[1] #Genero un dataframe con los sitios que almacena cada regi?n
    
    #  SELECCION Y AJUSTE DEL MODELO DE DISTRIBUCION DE PROBABILIDAD
    #Hay que modificar regfit para poder incluir gaucho
    #Si ARF$Z es gau, entonces calcular internamente, si no, con
    DistOpt<-names(which.min(abs(ARF$Z)))               
    
    if (DistOpt=="gau"){
    rfit<-pelgau(SummaryStatistics)
    RegionalQuantiles<-quakap(seq(0.05, 0.95, by=0.05), rfit$para)#Quantile calculation
    names(RegionalQuantiles)=seq(0.05, 0.95, by=0.05)
    Resultadosrfitdist[z]<-rfit$dist # Se identifica la distribucion utilizada
    Resultadosrfitpara[1:3,z]<-rfit$para[1:3]
    ResultadosRegionalQuantiles[1:19,z]<-RegionalQuantiles # Para cada regi?n "z", almaceno sus resultados
    } else {
    rfit<-regfit(SummaryStatisticsRegData, DistOpt) # Aca se selecciona autom?ticamente la distribuci?n con menor Z
    RegionalQuantiles<-regquant(seq(0.05, 0.95, by=0.05), rfit)#Quantile calculation
    Resultadosrfitdist[z]<-rfit$dist # Se identifica la distribucion utilizada
    Resultadosrfitpara[1:3,z]<-rfit$para # Se presentan los parametros de la distribuci?n ajustada
    ResultadosRegionalQuantiles[1:19,z]<-RegionalQuantiles # Para cada regi?n "z", almaceno sus resultados
    }
    ResultadosRMAP[z]<-weighted.mean(SummaryStatisticsRegData[[3]],SummaryStatisticsRegData[[2]]) # Se calcula la precipitaci?n media cada Regi?n
    TablaResumen_Regiones[z,1:5]<-c(round(ARF$H[1],digits=1),DistOpt,min(abs(ARF$Z)),length(which(ARF$D>3)),length(ARF$D))
    TablaResumen_Regiones[["Cluster"]]=NombreClusters[k]#Acá se debiese agregar la l, ya que define cada cluster. NO SE AGREGA LA L
    #TablaResumen_Regiones[["Cluster"]][z]=paste0(NombreClusters[k],"_",z)#Adaptación CAZALAC. ESTE ARCHIVO NO CUENT ACON ESTA ADAPTACIÓN
    lmrd(SummaryStatisticsRegData,xlim=c(-0.3,0.6),pch="*",legend.lmrd=list(cex=0.5),main=paste("Region ",NombreRegiones[z],sep=""),cex.main=0.8)
    lmrdpoints(Rlmoments, type="p", pch=19, cex=1.2, col="red" )
    #rgc <- regqfunc(rfit2)# Calcular Curva Crecimiento Regional
    #rgc(seq(0.05, 0.95, by=0.05))
    #curve(rgc, 0.01, 0.99, xlab="Probabilidad de no-excedencia, F", ylab="Curva Regional", main=paste("Region ",NombreRegiones[z],sep=""),cex.main=0.8)
    
  } # Cierra el ciclo for
  dev.off()
  
  # EVALUO RESULTADOS
  BaseResumenCluster=rbind(BaseResumenCluster,TablaResumen_Regiones)
  BaseSummaryStatistics[[k]]=ResultadosSummaryStatistics
  Baserfitdist[[k]]=Resultadosrfitdist
  Baserfitpara[[k]]=Resultadosrfitpara
  BaseBaseRegiones[[k]]=BaseRegiones
  BaseMediaCompleta[[k]]=MediaCompleta
  BaseProporCeros[[k]]=ProporCeros
  Basep0bias[[k]]=p0bias
} #Close for loop #Al parecer, se puede bypasear el error

BaseResumenCluster=BaseResumenCluster[-1,]
BaseResumenCluster$H1<- as.numeric(as.character(BaseResumenCluster$H1)) 
BaseResumenCluster$abs_z_value<- as.numeric(as.character(BaseResumenCluster$abs_z_value))
BaseResumenCluster$NDiscSites<- as.numeric(as.character(BaseResumenCluster$NDiscSites))
BaseResumenCluster$Sitios<- as.numeric(as.character(BaseResumenCluster$Sitios))
BaseResumenCluster$Cluster<- as.factor(BaseResumenCluster$Cluster)
write.csv(BaseResumenCluster,"BaseResumenCluster.csv",row.names =FALSE)
rm(BaseResumenCluster)
#BaseResumenCluster=read.csv("BaseResumenCluster.csv",sep=",") #Ojo con esta parte, quizás hay que ver el tema de los factores
BaseResumenCluster=read.csv("BaseResumenCluster.csv",sep=",", stringsAsFactors = TRUE) #Adapt Cazalac
View(BaseResumenCluster)

#............................................................................................................
#creacion de un DF, con dimensiones iguales a la cantidad de clusters donde se van 
#a almacenar las sumas de regiones que cumplen con criterios de eficiencia.
rl=function(x){length(x[x<=Hc])/length(x)*100}
Sumas <-data.frame(matrix(ncol=length(levels(BaseResumenCluster[,6])), nrow = 1))
colnames(Sumas) <- levels(BaseResumenCluster[,6])#ES IMPORTANTE QUE LA COLUMNA CLUSTER FIGURE COMO FACTOR, caso contrario, habria que hacer un "as.factor"

#loop que realiza las sumatorias y las almacena en el data frame creado.
for (i in 1:length(levels(BaseResumenCluster[,6]))){
  #Primer subset referido a cluster
  suma <- subset(BaseResumenCluster, Cluster == levels(BaseResumenCluster[,6])[i])
  
  #Segundo subset referido a H
  suma <- subset(suma, H1 <=Hc)
  
  #Sumatoria de columna "Estaciones"
  suma <- colSums(as.matrix(suma[,5]))
  
  #Se almacena el valor en la columna que corresponda (En el data frame creado anteriormente)
  Sumas[,which(colnames(Sumas) == levels(BaseResumenCluster[,6])[i])] <- suma/length(BaseDatosEstacionesClust[,1])*100 #con esto me seguro que encuentre la letra indicada
  #Sumas$(noquote(levels(datos[,3])[i])) <- suma
}

#FIGURAS 3-4-5 (MoranKernel.R)

ResumeTable=data.frame(Cluster=aggregate(H1~Cluster,data=BaseResumenCluster,rl)[1],
                       RE=aggregate(H1~Cluster,data=BaseResumenCluster,rl)[2],
                       TotalSites=aggregate(Sitios~Cluster,data=BaseResumenCluster,sum)[2],
                       AvrgSpR=aggregate(Sitios~Cluster,data=BaseResumenCluster,mean)[2],
                       AE=t(Sumas),row.names=NULL)# Derepente puede haber una regi?n extra?a

colnames(ResumeTable)=c("Cluster","RE","Total.Sites","SitesperRegion","AE")
ColNom <- data.frame(Trial=1:length(NombreClusters), Cluster=NombreClusters)
ResumeTable <- ResumeTable
#ResumeTable <- merge(ColNom,ResumeTable,by="Cluster") #sacado cazalac
#ResumeTable <- ResumeTable[order(ResumeTable$Trial),] #sacado cazalac
View(ResumeTable)
write.csv(ResumeTable,"ResumeTable.csv",row.names=FALSE)

#.............................................END OF BLOCK V ....................................................

##################################################################################################################
#                                           BLOCK VI: FREQUENCY/QUANTIL/RETURN PERIOD ESTIMATION                 #
##################################################################################################################
# These are the L-moments of the Variable of Interest (not necessarily annual precipitation)
#These are the L-moments of the non zero records
lmom.df=data.frame(BaseSummaryStatistics[[1]][,,1])
#fix debido a que si es igual a 1 BaseSummaryStatistics queda fuera de rango
if(dim(BaseSummaryStatistics[[1]])[3] > 1){
for (a in 2:dim(BaseSummaryStatistics[[1]])[3]){
  lmom.df=rbind(lmom.df,data.frame(BaseSummaryStatistics[[1]][,,a]))
}
}
colnames(lmom.df)=c("id_station","n","mean","L_CV","L_Skewness","L_Kurtosis","t_5")
lmom.df$id_station=as.factor(lmom.df$id_station)
lmom.df$n=as.numeric(as.character(lmom.df$n))
lmom.df$mean=as.numeric(as.character(lmom.df$mean))
lmom.df$L_CV=as.numeric(as.character(lmom.df$L_CV))
lmom.df$L_Skewness=as.numeric(as.character(lmom.df$L_Skewness))
lmom.df$L_Kurtosis=as.numeric(as.character(lmom.df$L_Kurtosis))
lmom.df$t_5=as.numeric(as.character(lmom.df$t_5))
lmom.df=lmom.df[which(rowSums(lmom.df[,2:7]) > 0),]
par(mfrow=c(1,1))
lmrd(distributions = "GLO GEV GPA GNO PE3 WAK.LB ALL.LB",cex=0.6)
lmrdpoints(lmom.df$L_Skewness,lmom.df$L_Kurtosis,type="p",pch=19,col="red")

BaseDatosEstacionesClustLmom=join(BaseDatosEstacionesClust,lmom.df,by="id_station")
write.csv(BaseDatosEstacionesClustLmom,"BaseDatosEstacionesClustLmom.csv",row.names=FALSE)
#.................................................................................................................

# Frecuency and quantile values of interest
CuantilInteres<-c(0.5,0.6, 0.7,0.8,0.9,1)
ProbInteres=c(1/100,1/90,1/80,1/70,1/60,1/50,1/40,1/30,1/20,1/10,1/5,1/2)
nameProbInteres=c("1_in_100yr","1_in_90yr","1_in_80yr","1_in_70yr","1_in_60yr","1_in_50yr","1_in_40yr","1_in_30yr",
                  "1_in_20yr","1_in_10yr","1_in_5yr","1_in_2yr")

SClst=which.max(ResumeTable$AE)

#......... Frequency Estimation
FRECUENCIAS<-list()
for (j in 1:length(CuantilInteres)){
  FQEFrequencias<-list()
 
  for (zz in 1:length(BaseBaseRegiones[[SClst]])){
    FrequencyEstimation<-data.frame(id_station=BaseSummaryStatistics[[SClst]][1:length(BaseBaseRegiones[[SClst]][[zz]]),1,zz],
                                    MediaSinCero=BaseSummaryStatistics[[SClst]][1:length(BaseBaseRegiones[[SClst]][[zz]]),3,zz],
                                    MediaConCero=unlist(BaseMediaCompleta[[SClst]][[zz]]))
    #FrequencyEstimation<- data.frame(id_station=ResultadosSummaryStatistics[1:length(BaseRegiones[[zz]]),1,zz],Media=ResultadosSummaryStatistics[1:length(BaseRegiones[[zz]]),3,zz])
    FrequencyEstimation$MediaSinCero<-as.numeric(as.character(FrequencyEstimation$MediaSinCero))# ?as.numeric(levels(f))[f]?
    FrequencyEstimation$MediaConCero<-as.numeric(as.character(FrequencyEstimation$MediaConCero))
    FrequencyEstimation$propCero<-unlist(BaseProporCeros[[SClst]][zz])
    FrequencyEstimation$biasCero<-unlist(Basep0bias[[SClst]][zz])
    FrequencyEstimation$Dist<-Baserfitdist[[SClst]][zz]
    #FrequencyEstimation=FrequencyEstimation[mixedorder(FrequencyEstimation$id_station),]
    #cuant<-Int.pp[row.names(Int.pp) %in% FrequencyEstimation$id_station,]#Se debe garantizar el orden
    #cuant=cuant[mixedorder(cuant$id_station),]
    #FrequencyEstimation$acumulado=cuant$int.pp
    #FrequencyEstimation$qSinCero<-cuant[,2]/FrequencyEstimation$MediaSinCero
    #FrequencyEstimation$qConCero<-cuant[,2]/FrequencyEstimation$MediaConCero
    #FrequencyEstimation$desvioSinCero<-(FrequencyEstimation$qSinCero-1)*100
    #FrequencyEstimation$desvioConCero<-(FrequencyEstimation$qConCero-1)*100
    switch(Baserfitdist[[SClst]][zz],
           "glo"= FrequencyEstimation$EstFreq<-cdfglo(CuantilInteres[j], para = Baserfitpara[[SClst]][1:3,zz]), #**** Valores de X van en el nombre, donde X es el Cuantil de Interes
           "gev"= FrequencyEstimation$EstFreq<-cdfgev(CuantilInteres[j], para = Baserfitpara[[SClst]][1:3,zz]), #****
           "gpa"= FrequencyEstimation$EstFreq<-cdfgpa(CuantilInteres[j], para = Baserfitpara[[SClst]][1:3,zz]),  #***
           "gno"= FrequencyEstimation$EstFreq<-cdfgno(CuantilInteres[j], para = Baserfitpara[[SClst]][1:3,zz]),  #***
           "pe3"= FrequencyEstimation$EstFreq<-cdfpe3(CuantilInteres[j], para = Baserfitpara[[SClst]][1:3,zz]),
           "gau"= FrequencyEstimation$EstFreq<-cdfkap(CuantilInteres[j], para = c(Baserfitpara[[SClst]][1:3,zz],0.5)))  #***
    #CORREGIR EN TOOOOODOS LOS SCRIPT EXISTENTES
    FrequencyEstimation$Px=ifelse(CuantilInteres[j]>0,
                                  FrequencyEstimation$propCero+(1-FrequencyEstimation$propCero)*FrequencyEstimation$EstFreq,
                                  FrequencyEstimation$biasCero)
    FrequencyEstimation$PR=ifelse (FrequencyEstimation$Px>0.5,
                                   1/(1-FrequencyEstimation$Px),
                                   1/FrequencyEstimation$Px)
    #FrequencyEstimation$IPR=ifelse (FrequencyEstimation$qConCero>=1,
    #                               1/(1-FrequencyEstimation$Px),
    #                               -1/FrequencyEstimation$Px)
    #FrequencyEstimation$Z=qnorm(FrequencyEstimation$Px)
    
    names(FrequencyEstimation)[7]<-paste0("EstFreq",CuantilInteres[j])
    names(FrequencyEstimation)[8]<-paste0("Px",CuantilInteres[j])
    names(FrequencyEstimation)[9]<-paste0("RP",CuantilInteres[j])
    #names(FrequencyEstimation)[6]<-"EstFreq"
    #names(FrequencyEstimation)[7]<-"PR"
    
    FQEFrequencias[[zz]]<-FrequencyEstimation
  }
  
  FQEFrecuenciasDataFrame<-FQEFrequencias[[1]]
  for (xxx in 2:length(FQEFrequencias)){
    FQEFrecuenciasDataFrame<-rbind(FQEFrecuenciasDataFrame,FQEFrequencias[[xxx]])}
  FRECUENCIAS[[j]]<-FQEFrecuenciasDataFrame
}

 
#......... Quantile Estimation
CUANTILES<-list()
FQECuantiles<-list()
for (l in 1:length(ProbInteres)){
  FQECuantiles<-list()
  for (yy in 1:length(BaseBaseRegiones[[SClst]])){
    QuantilEstimation<-data.frame(id_station=BaseSummaryStatistics[[SClst]][1:length(BaseBaseRegiones[[SClst]][[yy]]),1,yy],Media=BaseSummaryStatistics[[SClst]][1:length(BaseBaseRegiones[[SClst]][[yy]]),3,yy])
    #QuantilEstimation<- data.frame(id_station=ResultadosSummaryStatistics[1:length(BaseRegiones[[yy]]),1,yy],Media=ResultadosSummaryStatistics[1:length(BaseRegiones[[yy]]),3,yy])
    QuantilEstimation$Media<-as.numeric(as.character(QuantilEstimation$Media))# ?as.numeric(levels(f))[f]?
    QuantilEstimation$Dist<-Baserfitdist[[SClst]][yy]
    
    
    switch(Baserfitdist[[SClst]][yy],
           "glo"= QuantilEstimation$EstQuant<-quaglo(ProbInteres[l], para = Baserfitpara[[SClst]][1:3,yy]), #*** Valores de f van el el nombre, donde f es la probabilidad de interes
           "gev"= QuantilEstimation$EstQuant<-quagev(ProbInteres[l], para = Baserfitpara[[SClst]][1:3,yy]),#***
           "gpa"= QuantilEstimation$EstQuant<-quagpa(ProbInteres[l], para = Baserfitpara[[SClst]][1:3,yy]),#***
           "gno"= QuantilEstimation$EstQuant<-quagno(ProbInteres[l], para = Baserfitpara[[SClst]][1:3,yy]),#***
           "pe3"= QuantilEstimation$EstQuant<-quape3(ProbInteres[l], para = Baserfitpara[[SClst]][1:3,yy]),
           "gau"= QuantilEstimation$EstQuant<-quakap(ProbInteres[l], para = c(Baserfitpara[[SClst]][1:3,zz],0.5)))  #***

    
    QuantilEstimation$EstQuant<-QuantilEstimation$Media*QuantilEstimation$EstQuant #***
    QuantilEstimation$Balance<-(QuantilEstimation$EstQuant-QuantilEstimation$Media)/QuantilEstimation$Media*100
    names(QuantilEstimation)[4]<-paste0("EstQuant_",nameProbInteres[l])
    names(QuantilEstimation)[5]<-paste0("Balance_",nameProbInteres[l])
    FQECuantiles[[yy]]<-QuantilEstimation
  }
  
  FQECuantilesDataFrame<-FQECuantiles[[1]]
  for (vvv in 2:length(FQECuantiles)){
    FQECuantilesDataFrame<-rbind(FQECuantilesDataFrame,FQECuantiles[[vvv]])
  }
  CUANTILES[[l]]<-FQECuantilesDataFrame
}  

# F.4. Dataframe construction# 
DFFrecuencias<-FRECUENCIAS[[1]]
for (zzz in 2:length(FRECUENCIAS)){
  DFFrecuencias<-cbind(DFFrecuencias,FRECUENCIAS[[zzz]][7:9])
}


DFCuantiles<-CUANTILES[[1]]
for (zzz in 2:length(CUANTILES)){
  DFCuantiles<-cbind(DFCuantiles,CUANTILES[[zzz]][4:5])
}

DFCuantiles_DFFrecuencias<-cbind(DFFrecuencias,DFCuantiles[4:length(DFCuantiles)])

#................... Creation of the DataBase "BaseModelMap" for mapping
BaseModelMap<-merge(BaseDatosEstacionesClust,DFCuantiles_DFFrecuencias,by.x="id_station",by.y="id_station")
colnames(BaseModelMap)=gsub(".", "", names(BaseModelMap), fixed = TRUE)
#colnames(BaseModelMap)=gsub("=", "", names(BaseModelMap), fixed = TRUE)
BaseModelMapCor=BaseModelMap

# Remove Inf values
PRCols=grep('^RP',names(BaseModelMapCor))
BaseModelMapCor[PRCols][BaseModelMapCor[PRCols]==Inf]=NA
for (i in 1:length(PRCols)){
  BaseModelMapCor[PRCols[i]][is.na(BaseModelMapCor[PRCols[i]])]=max(BaseModelMapCor[PRCols[i]],na.rm=T)
}

# Transform negative quantiles to zero
FQCols=grep('^EstQuant',names(BaseModelMapCor))
for (i in 1:length(FQCols)){
  BaseModelMapCor[FQCols[i]][BaseModelMapCor[FQCols[i]]<0]=0
}

write.csv(BaseModelMapCor, file = "BaseModelMapCor.csv",row.names=FALSE)
remove(BaseModelMapCor) # Aca debo remover y volver a cargar la  Base de Datos en el sistema
BaseModelMapCor <- read.csv(file="BaseModelMapCor.csv",sep=",",head=T)# Y vuelvo a cargar la Base de Datos
View(BaseModelMapCor)
# At this stage, the Final DataBase "BaseModelMapCor", with al necessary variables for mapping is already available
#The user can produce preliminary maps with conventional interpolation methods
#....................................................END F BLOCK VI...........................................

########################################################################################################
#                                BLOCK VII: FREQUENCY/RETURN PERIOD/QUANTIL MAPPING                    #
########################################################################################################

# Prepare Thiessen polygons mask shapefile and plotting

# Calculo Periodo de Retorno m?ximo posible de calcular 5t- rule Jakob et al, 1999
#RecLength=quantile(BaseModelMapCorII$RL_Station,0.5)
RecLength=36
sitespr=round(tail(ResumeTable$SitesperRegion,1),0)
max_RP=as.numeric(RecLength*sitespr/5)
minFreq=1/max_RP

BaseModelMapMask=BaseModelMapCor[!duplicated(BaseModelMapCor[,c("lon","lat")]),]#Eliminar duplicados de coordenadas. Cuidado con columnas que indican coordenadas
#Maximum Return Period to be estimated (5000)
PRCols=grep('^RP',names(BaseModelMapMask))
for (i in 1:length(PRCols)){
  BaseModelMapMask[PRCols[i]][BaseModelMapMask[PRCols[i]]>5000]=5000
}

# Minimum frequency to be estimated (0.001)
FQCols=grep('^EstFreq',names(BaseModelMapMask))
for (i in 1:length(FQCols)){
  BaseModelMapMask[FQCols[i]][BaseModelMapMask[FQCols[i]]<0.001]=0.001
}

BaseMPMaskPlot=BaseModelMapMask[c(grep('lon',names(BaseModelMapMask)),
                                  grep('lat',names(BaseModelMapMask)),
                                  grep('ClustReg_adapt',names(BaseModelMapMask)),#Adicionado para mapeo de clusters
                                  grep(paste0('Dist',"$",collapse="|"),names(BaseModelMapMask)),#Adicionado para mapeo de clusters
                                  grep('EstFreq',names(BaseModelMapMask)),
                                  grep('^RP',names(BaseModelMapMask)),
                                  grep('EstQuant',names(BaseModelMapMask)))]
coordinates(BaseMPMaskPlot)= c("lon","lat")
writeSpatialShape(BaseMPMaskPlot, "BaseMPMaskPlot")

#OPcion Thiessen con Voronoi dentro de R
Thiessen <- voronoipolygons(BaseMPMaskPlot)
Thiessen<-spCbind(Thiessen, data.frame(BaseMPMaskPlot))
Thiessen <- Thiessen[,!(names(Thiessen) %in% c("optional"))]
Thiessen <- Thiessen[,-c(1,2)]
plot(Thiessen)
writeSpatialShape(Thiessen, "Thiessen")
#rsaga.geoprocessor("shapes_points",16,list(POINTS="BaseMPMaskPlot.shp",POLYGONS="Thiessen",FRAME=4))
getinfo.shape(paste0(country,".shp"))# Este se edita y guarda aparte con SAGA
getinfo.shape("Thiessen.shp") 
rsaga.geoprocessor("shapes_polygons",11,list(CLIP=paste0(country,".shp"),
                                             S_INPUT="Thiessen.shp",
                                             S_OUTPUT="CutThiessen.shp",
                                             M_INPUT="Thiessen.shp",
                                             MULTIPLE=0),display.command=TRUE)
#Plot 
CutThiessen <- readShapeSpatial("CutThiessen.shp")
CutThiessen <- CutThiessen[,!(names(CutThiessen) %in% c("SP_ID"))]
names(CutThiessen)[3:length(names(CutThiessen))]=names(BaseMPMaskPlot)
plot(CutThiessen)


#.......Variables for mapping
MapVarEstFreq=names(BaseModelMapCor)[grep('^EstFreq',names(BaseModelMapCor))]
MapVarRetPeriod=names(BaseModelMapCor)[grep('^RP',names(BaseModelMapCor))]
MapVarEstQuant=names(BaseModelMapCor)[grep('^^EstQuant',names(BaseModelMapCor))]
MapVarBalance=names(BaseModelMapCor)[grep('^^Balance',names(BaseModelMapCor))]

# ..................................Quantile maps (precipitation maps) for a given return period
folder.mapas="/Mapas"
dir.create(paste0(getwd(),folder.mapas))

#....... Creation of list objects for saving results
CalGOFDB=list()
ValidGOFDB=list()
VarImpDB=list()
folder <- getwd() # Defino path a la carpeta donde se van a almacenar los resultados
#...Progress bar for rasterbrick
pb <- progress_bar$new(
  format = "(:spin) [:bar] :percent",
  total = 21, clear = FALSE, width = 60)

Predictor<-paste0("AFRP",1:22)
# Incorporaci?n a BaseModelMapCor de los valores de los predictores
Prtrs=data.frame(matrix(NA,nrow=dim(BaseModelMapCor)[1],ncol=22))
colnames(Prtrs)=paste0("AFRP",1:22)
coords =  SpatialPoints(BaseModelMapCor[, c("lon", "lat")],proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#extract stations coordinates
plot(CutThiessen,col="blue")
points(coords,col="red")

#Rasterbrick creation with the maps of all predictors
raster_list<-Predictor   
rast.list<-list()
temp.r.l=raster(paste0("../../ADAFolderPredictors/",raster_list[1],".img"))
temp.r.l=crop(temp.r.l,polygons(Boundaries))
temp.r.l<-mask(temp.r.l,polygons(Boundaries))
plot(temp.r.l,main="P1")
rast.list[[1]]=temp.r.l
Prtrs[,1]=extract(temp.r.l, coords, method='bilinear')
for(l in 2:length(raster_list)){
  country.raster<-raster(paste0("../../ADAFolderPredictors/",raster_list[l],".img"))
  country.raster<-crop(country.raster,polygons(Boundaries))
  country.raster<-mask(country.raster,polygons(Boundaries))
    if(compareRaster(country.raster,temp.r.l)==TRUE){
    country.raster=country.raster} else {
    country.raster<-resample(country.raster,temp.r.l)  
    }
  print(compareRaster(country.raster,temp.r.l))
  Prtrs[,l]=extract(country.raster, coords, method='bilinear')
  rast.list[[l]] <-country.raster
  plot(country.raster,main=paste0("P",l))
  pb$tick()
  Sys.sleep(1 / 21)
}

print(paste0(length(Predictor), " layers used from ",22," layers available"))
brQ=brick(rast.list)
names(brQ) <- raster_list #CAZALAC mod so predict function can be executed (best solution ever)
BaseModelMapCor2=na.omit(cbind(BaseModelMapCor,Prtrs))
#BaseModelMapCor2=randomSample(BaseModelMapCor2,500)# Activar en caso muy grande
write.csv(BaseModelMapCor2,"BaseModelMapCor2.csv",row.names=FALSE)
View(BaseModelMapCor2)

#....Modelling and Mapping of precipitation quantiles
#.......Save calibration, validation and variable importance plots to a PDF
pdf(paste0("Quantiles_",country,"_.pdf"), width = 16 , height = 10)
par(mfrow=c(2,2))
for(i in 1:length(MapVarEstQuant)){
  Predictant<- MapVarEstQuant[i] # Defino el nombre del archivo de salida
  qdatafn <- BaseModelMapCor2[BaseModelMapCor2[,MapVarEstQuant[i]]>=0,] #Datos con cuantiles no negativos
  Mask <- CutThiessen[CutThiessen[[MapVarEstQuant[i]]]>=0,] #Datos con cuantiles no negativos
  
  qdata.trainfn <- paste0(Predictant,"_BaseModelMapCor_TRAIN.csv") # Este va a ser el archivo (o data.frame) que se va a usar para entrenamiento
  qdata.testfn <-  paste0(Predictant,"_BaseModelMapCor_TEST.csv") # Este va a ser el archivo (o data.frame) que se va a usar como test de validaci?n del modelo
  
  get.test( proportion.test=0.3,
            qdatafn=qdatafn,
            seed=,
            folder=folder,
            qdata.trainfn=qdata.trainfn,
            qdata.testfn=qdata.testfn) # Con este comando  selecciono aleatoriamente un set de casos de validaci?n
  
  #Randomforest model calibration
  RFModelcalib=read.csv(paste0(Predictant,"_BaseModelMapCor_TRAIN.csv"))
  RFModelcalib=RFModelcalib[,c(Predictant,Predictor)]
  f <- paste(names(RFModelcalib)[1], "~", paste(names(RFModelcalib)[-1], collapse=" + "))
  mcal<-randomForest(as.formula(f),data=RFModelcalib,mtry=10,ntree=1000)
  calib=data.frame(Obs=RFModelcalib[1],Simul=predict(mcal))
  names(calib)=c("Obs","Simul")
  write.csv(calib,file=paste0(folder,"/calibrationTablefor",Predictant,".csv"),row.names=F)
  partylmcalib=lm(Simul~Obs,data=calib)
  summary(partylmcalib)
  plot(calib[,1],calib[,2],xlab="Observed",ylab="Simulated",col="gray",pch=19,main=paste0("Calibration for ",Predictant))  
  abline(partylmcalib,col="black",lwd=1)  
  
  varImpPlot(mcal,main=paste0("Variable Importance for calibrated ",Predictant))

  CalGOFDB[[i]]=gof(calib$Sim,calib$Obs,digits=3) # Diagn?stico del modelo de calibraci?n
  VarImpDB[[i]]=importance(mcal)
  
  #Randomforest model validation
  RFModelvalid=read.csv(paste0(folder,"/",Predictant,"_BaseModelMapCor_TEST.csv"))
  RFModelvalid=RFModelvalid[,c(Predictant,Predictor)]
  valid=data.frame(Obs=RFModelvalid[1],Simul=predict(mcal,newdata=RFModelvalid[,-1]))
  names(valid)=c("Obs","Simul")
  write.csv(valid,file=paste0(folder,"/ValidationTablefor",Predictant,".csv"),row.names=F)
  partylmvalid=lm(Simul~Obs,data=valid)
  summary(partylmvalid)
  plot(valid[,1],valid[,2],xlab="Observed",ylab="Simulated",col="gray",pch=19,main=paste0("Validation for ",Predictant))  
  abline(partylmvalid,col="black")  
  ValidGOFDB[[i]]=gof(valid$Sim,valid$Obs,digits=3)# Diagn?stico del modelo de validaci?n
  
  # Mapping using calibrated randomforest model
  Mapping.Function(brQ,Predictant,mcal,folder.mapas) #Acá aparece un error
  FinalQMap<-mask(raster(paste0(getwd(),folder.mapas,"/",Predictant,".sdat")),polygons(Mask))
  projection(FinalQMap)="+proj=longlat +ellps=WGS84 + datum=WGS84"
  writeRaster(FinalQMap,filename=paste0(getwd(),folder.mapas,"/",Predictant,".tif"),format="GTiff", overwrite=TRUE,NAflag=-999)
  plot(FinalQMap,main=Predictant)
}
dev.off()

CalGOFmeasure=matrix(unlist(CalGOFDB),nrow=20);colnames(CalGOFmeasure)=MapVarEstQuant;rownames(CalGOFmeasure)=rownames(CalGOFDB[[i]])
ValGOFmeasure=matrix(unlist(ValidGOFDB),nrow=20);colnames(ValGOFmeasure)=MapVarEstQuant;rownames(ValGOFmeasure)=rownames(ValidGOFDB[[i]])
VarImpmeasure=matrix(unlist(VarImpDB),nrow=i);colnames(VarImpmeasure)=Predictor;rownames(VarImpmeasure)=MapVarEstQuant
VarImpmeasure=t(VarImpmeasure)

write.csv(CalGOFmeasure,"CalGOFmeasure.csv", row.names=TRUE)
write.csv(ValGOFmeasure,"ValGOFmeasure.csv", row.names=TRUE)
write.csv(VarImpmeasure,"VarImpmeasure.csv", row.names=TRUE)


# ..................................Return Period maps for a given proportion of mean precipitation

CalGOFDBRP=list()
ValidGOFDBRP=list()
VarImpDBRP=list()
folder <- getwd()

#Modelling and mapping of Return Periods
#.......Save calibration, validation and variable importance plots to a PDF
pdf(paste0("ReturnPeriods_",country,"_.pdf"), width = 16 , height = 10)
par(mfrow=c(2,2))
for(i in 1:length(MapVarRetPeriod)){
  # Creaci?n de la basedatos qdatafnRP
  qdatafnRP <- BaseModelMapCor2[BaseModelMapCor2[,MapVarRetPeriod[i]]<50000,]
  #M?scara  
  Mask <- CutThiessen[CutThiessen[[MapVarRetPeriod[i]]]<50000,]
  FacCorMedia=aggregate(MediaConCero~ClustReg_adapt, data=qdatafnRP,mean)
  FacCorDS=aggregate(MediaConCero~ClustReg_adapt, data=qdatafnRP,sd)
  FacCor=data.frame(ClustReg_adapt=as.character(FacCorMedia$ClustReg_adapt),FCmedia=FacCorMedia$Media,FCds=FacCorDS$Media)
  RPs=qdatafnRP[grep('^RP',names(qdatafnRP))]
  RPs=qdatafnRP[MapVarRetPeriod[i]]
  RPs2=cbind(qdatafnRP$id_station,qdatafnRP$ClustReg_adapt,qdatafnRP$MediaConCero,RPs)
  colnames(RPs2)[1:3]=c("id_station","ClustReg_adapt","Media")
  RPs3=merge(RPs2,FacCor,by="ClustReg_adapt")
  RPs3[MapVarRetPeriod[i]]=(((RPs3$Media-RPs3$FCmedia)/RPs3$FCds)*0.05*RPs3[MapVarRetPeriod[i]]+RPs3[MapVarRetPeriod[i]])
  RPs3=RPs3[order(RPs3$id_station),]
  for (j in 1:dim(RPs3)[1]){
    if (is.na(RPs3[j,"FCds"])==FALSE) RPs3[j,MapVarRetPeriod[i]]=RPs3[j,MapVarRetPeriod[i]] else RPs3[j,MapVarRetPeriod[i]]=RPs2[j,MapVarRetPeriod[i]]
  }
  qdatafnRP[MapVarRetPeriod[i]]=RPs3[MapVarRetPeriod[i]]
  
  #Identificacion de predictando
  Predictant<- MapVarRetPeriod[i] # Defino el nombre del archivo de salida
  
  #Generacion muestras de calibracion y validacion
  qdata.trainfnRP <- "BaseModelMapCor_TRAINRP.csv" # Este va a ser el archivo (o data.frame) que se va a usar para entrenamiento
  qdata.testfnRP <- "BaseModelMapCor_TESTRP.csv" # Este va a ser el archivo (o data.frame) que se va a usar como test de validaci?n del modelo
  get.test( proportion.test=0.3,
            qdatafn=qdatafnRP,
            seed=,
            folder=folder,
            qdata.trainfn=qdata.trainfnRP,
            qdata.testfn=qdata.testfnRP) # Con este comando  selecciono aleatoriamente un set de casos de validaci?n
  
  #Elaboracion, visualizacion y almacenamiento del modelo de calibracion
  RFModelcalibRP=read.csv("BaseModelMapCor_TRAINRP.csv")
  RFModelcalibRP=RFModelcalibRP[,c(Predictant,Predictor)]
  f <- paste(names(RFModelcalibRP)[1], "~", paste(names(RFModelcalibRP)[-1], collapse=" + "))
  mcalRP<-randomForest(as.formula(f),data=RFModelcalibRP,mtry=10,ntree=1000)
  calibRP=data.frame(Obs=RFModelcalibRP[1],Simul=predict(mcalRP))
  names(calibRP)=c("Obs","Simul")
  partylmcalibRP=lm(Simul~Obs,data=calibRP)
  summary(partylmcalibRP)
 
  plot(calibRP[,1],calibRP[,2],xlab="Observed",ylab="Simulated",col="gray",pch=19,main=paste0("Calibration RP for ",Predictant))  
  abline(partylmcalibRP,col="black",lwd=1)  
  
  varImpPlot(mcalRP,main=paste0("Variable Importance for calibrated ",Predictant))
  
  CalGOFDBRP[[i]]=gof(calibRP$Sim,calibRP$Obs,digits=3) # Diagn?stico del modelo de calibraci?n
  VarImpDBRP[[i]]=importance(mcalRP)
  
  #Elaboraci?n del modelo de validacion
  RFModelvalidRP=read.csv("BaseModelMapCor_TESTRP.csv")
  RFModelvalidRP=RFModelvalidRP[,c(Predictant,Predictor)]
  validRP=data.frame(Obs=RFModelvalidRP[1],Simul=predict(mcalRP,newdata=RFModelvalidRP[,-1]))
  names(validRP)=c("Obs","Simul")
  write.csv(validRP,file=paste0("ValidationRPTablefor",Predictant,".csv"),row.names=F)
  partylmvalidRP=lm(Simul~Obs,data=validRP)
  summary(partylmvalidRP)
  plot(validRP[,1],validRP[,2],xlab="Observed",ylab="Simulated",col="gray",pch=19,main=paste0("Validation RP for ",Predictant))  
  abline(partylmcalibRP,col="black",lwd=1)
 
  
  # Diagn?stico del modelo de validaci?n.Sin embargo, solamente vale pseudo R2
  ValidGOFDBRP[[i]]=gof(validRP$Sim,validRP$Obs,digits=3)# Diagn?stico del modelo de validaci?n
  
  # Mapeo usando modelo de calibraci?n
  Mapping.Function(brQ,Predictant,mcalRP,folder.mapas)
  FinalRPMap<-mask(raster(paste0(getwd(),folder.mapas,"/",Predictant,".sdat")),polygons(Mask))
  projection(FinalRPMap)="+proj=longlat +ellps=WGS84 + datum=WGS84"
  writeRaster(FinalRPMap,filename=paste0(getwd(),folder.mapas,"/",Predictant,".tif"),format="GTiff", overwrite=TRUE,NAflag=-999)
  plot(FinalRPMap,main=Predictant)
}# Finalizado Generaci?n de Mapas de Periodo de Retorno
dev.off()

#Genero, guardo y visualizo las tablas de diagnostico de calibracion y validacion
CalGOFRPmeasure=matrix(unlist(CalGOFDBRP),nrow=20);colnames(CalGOFRPmeasure)=MapVarRetPeriod;rownames(CalGOFRPmeasure)=rownames(CalGOFDBRP[[i]])
ValGOFRPmeasure=matrix(unlist(ValidGOFDBRP),nrow=20);colnames(ValGOFRPmeasure)=MapVarRetPeriod;rownames(ValGOFRPmeasure)=rownames(ValidGOFDBRP[[i]])
VarImpRPmeasure=matrix(unlist(VarImpDBRP),nrow=i);colnames(VarImpRPmeasure)=Predictor;rownames(VarImpRPmeasure)=MapVarRetPeriod
VarImpRPmeasure=t(VarImpRPmeasure)

write.csv(CalGOFRPmeasure,"CalGOFRPmeasure.csv", row.names=TRUE)
write.csv(ValGOFRPmeasure,"ValGOFRPmeasure.csv", row.names=TRUE)
write.csv(VarImpRPmeasure,"VarImpRPmeasure.csv", row.names=TRUE)
