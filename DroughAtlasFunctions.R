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

