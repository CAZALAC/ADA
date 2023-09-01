#Identificacion de predictando
Predictant<- "Mediaconcero" # Defino el nombre del archivo de salida

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
RFModelcalibRP=RFModelcalibRP[,c("MediaConCero",Predictor)]
f <- paste(names(RFModelcalibRP)[1], "~", paste(names(RFModelcalibRP)[-1], collapse=" + "))
mcalRP<-randomForest(as.formula(f),data=RFModelcalibRP,mtry=10,ntree=1000)
calibRP=data.frame(Obs=RFModelcalibRP[1],Simul=predict(mcalRP))
names(calibRP)=c("Obs","Simul")
partylmcalibRP=lm(Simul~Obs,data=calibRP)
summary(partylmcalibRP)



Mapping.Function(brQ,Predictant,mcalRP,folder.mapas)

FinalQMap<-mask(raster(paste0(getwd(),folder.mapas,"/",Predictant,".sdat")),polygons(Mask))
projection(FinalQMap)="+proj=longlat +ellps=WGS84 + datum=WGS84"
terra::writeRaster(FinalQMap,filename=paste0(getwd(),folder.mapas,"/",Predictant,".tif"),filetype="GTiff", overwrite=TRUE,NAflag=-999)
plot(FinalQMap,main=Predictant)