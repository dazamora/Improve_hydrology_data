####################################################################
## This script improve data in hydrological time series of precipitation,
## temperature and flow through outlier detection (spatial-temporal) and fill missing data
###################################################################


# Version 1.0.0
# Copyright David Zamora
# e-mail: dazamoraa@unal.edu.co

#Total Path
root<-getwd()

#-----Inputs-----

#-----Outputs-----



#------Load Library------
library(magicaxis)
library(waterData)
library(ggplot2)
library(gstat)
library(sp)
library(maptools)

#------Load Data------
listdata<-list.files(paste(root,"/DATA/",sep=""))
# In situ observe
Prec<-as.data.frame(read.table(paste(root,"/DATA/",grep("P",listdata,value = TRUE),sep=""),header=T))
Tmax<-read.table(paste(root,"/DATA/",grep("T_Max",listdata,value = TRUE),sep=""),header=T)
Tmin<-read.table(paste(root,"/DATA/",grep("T_Min",listdata,value = TRUE),sep=""),header=T)
Tmean<-read.table(paste(root,"/DATA/",grep("T_Mean",listdata,value = TRUE),sep=""),header=T)
Flow<-read.table(paste(root,"/DATA/",grep("Flow",listdata,value = TRUE),sep=""),header=T)

#Satellite data 
Prec.ERA<-as.data.frame(read.table(paste(root,"/DATA/",grep("ERA",listdata,value = TRUE),sep=""),header=T))
Prec.TRMM<-as.data.frame(read.table(paste(root,"/DATA/",grep("TRMM",listdata,value = TRUE),sep=""),header=T))

# Spatial Coordenates
Coord.Station<-read.csv(paste(root,"/DATA/",grep("Weather",listdata,value = TRUE),sep=""),header=T,sep = ";")

# What is time period to model?
period<-c("1/01/1985","31/12/2009")
ind.period<-c(which(Prec[,1]==period[1]),which(Prec[,1]==period[2]))
ind.period.sat<-c(which(Prec.ERA[,1]==period[1]),which(Prec.ERA[,1]==period[2]))
H.days<-seq(from=as.Date("01-01-1985",format="%d-%m-%Y"),to=as.Date("31-12-2009",format="%d-%m-%Y"),by="day")

Prec<-Prec[ind.period[1]:ind.period[2],]
Tmax<-Tmax[ind.period[1]:ind.period[2],]
Tmin<-Tmin[ind.period[1]:ind.period[2],]
Tmean<-Tmean[ind.period[1]:ind.period[2],]
Flow<-Flow[ind.period[1]:ind.period[2],]
Prec.ERA<-Prec.ERA[ind.period.sat[1]:ind.period.sat[2],]

P.partial<-matrix(NA,dim(Prec.ERA)[1],dim(Prec.TRMM)[2])
P.partial[(dim(Prec.ERA)[1]-dim(Prec.TRMM)[1]+1):dim(Prec.ERA)[1],2:17]<-as.matrix(Prec.TRMM[,2:17])
P.partial[,1]<-Prec.ERA[,1]
Prec.TRMM<-P.partial

#-----Fill missing data with satellite information----
M.Total.Prec<-list(OBS=Prec,ERA=Prec.ERA,TRMM=Prec.TRMM)
M.Median.Prec<-c()

for(k in 2:dim(Prec)[2]){
  M.tempo<-cbind(M.Total.Prec[[1]][,k],M.Total.Prec[[2]][,k],M.Total.Prec[[3]][,k])
  Median.Prec<-apply(M.tempo,1,median,na.rm=T)
  M.Median.Prec<-cbind(Median.Prec,M.Median.Prec)
}


#------Draw Time Series------
par(mfrow=c(8,2))
colr<-rainbow(18)

for(i in 2:17){
  par(mgp=c(2,0.5,0))
  par(mar=c(2.8,3.2,1,0.5)+0.1)
  
  plot(M.Median.Prec[,i-1],type="l",col=colr[i],xlab="",ylab="",las=2)
  ind.out<-which(rowSums(Matrix.Outliers)>10)

}


#----Aggregate Time Series and outlier detection-----

mo <- strftime(H.days, "%m") # Clasification by months
Ma.Prec <- data.frame(mo, M.Median.Prec)

Mon<-unique(mo)
Matrix.Outliers<-matrix(0,dim(M.Median.Prec)[1],dim(M.Median.Prec)[2])

for(i in 1:dim(M.Median.Prec)[2]){
  Prec.agg <- aggregate(M.Median.Prec[,i]~ mo, Ma.Prec,boxplot, na.rm=T,plot=FALSE)
  for(j in 1:12){
    ind<-which(Mon[j]==mo)
    outliers<-ind[is.element(el=M.Median.Prec[ind,i],set = Prec.agg$`M.Median.Prec[, i]`[j,4]$out)]
    Matrix.Outliers[outliers,i]<-1
    rm(outliers)
    }
}

#-----Interpolate by IDW-----

x.range <- as.numeric(c(-75.75,-74.75))  # min/max longitude of the interpolation area
y.range <- as.numeric(c(4,4.75))  # min/max latitude of the interpolation area
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.01),
                   y = seq(from = y.range[1], to = y.range[2], by = 0.01))
coordinates(grd) <- ~x + y
gridded(grd) <- TRUE
# plot(grd, cex = 1.5, col = "grey")
# points(Coord.Station, pch = 1, col = "red", cex = 1)

Coord.Station1<-Coord.Station
event.rain<-which(rowSums(Matrix.Outliers)>=1)
event.rain.out<-which(rowSums(Matrix.Outliers)<=10 & rowSums(Matrix.Outliers)>0)


hist.out<-hist(rowSums(Matrix.Outliers[event.rain,]),breaks = c(0:16),labels = TRUE,xaxt="n",
     las=2,col = c(rep(gray(0.5),10),rep("red",6)),main="",xlab="Número de  estaciones outliers por evento",
     ylab ="Frecuencia", ylim = c(0,1200))

axis(1,at=hist.out$mids,labels = 1:16)


Matrix.Outliers2<-Matrix.Outliers
M.Median.Prec2<-M.Median.Prec

for(j in 1:length(event.rain.out)){
  
  print(paste("Event ",j,sep="="))
  out.up1<-which(Matrix.Outliers[event.rain.out[j],]==1)
  
  # Validated data - without outliers
  Coord.Station.V<-Coord.Station[-out.up1,]
  Coord.Station.V$Data<-as.numeric(M.Median.Prec[event.rain.out[j],-out.up1])
  Coord.Station1<-Coord.Station.V
  Coord.Station1$y<-Coord.Station1$Latitud
  Coord.Station1$x<-Coord.Station1$Longitud
  coordinates(Coord.Station1) = ~x + y
  Mod.idw1<- idw(formula = Data ~ 1, locations = Coord.Station1, 
                 newdata = grd,debug.level=0)
  idw.output1 = as.data.frame(Mod.idw1)  # output is defined as a data table
  names(idw.output1)[1:3] <- c("long", "lat", "var1.pred")
  # ggplot() + geom_tile(data = idw.output, aes(x = long, y = lat, fill = var1.pred)) 
  
  pred.rain1<-c()
  
  for(l in 1:(dim(Coord.Station.V)[1]+length(out.up1))){
    Coord.End<-rbind(Coord.Station.V[,1:3],Coord.Station[out.up1,1:3])
    ind.lat1<-which(round(Coord.End$Latitud[l],2)==idw.output1$lat)
    ind.lon1<-which(round(Coord.End$Longitud[l],2)==idw.output1$long)
    ind.sample1<-which(is.element(ind.lat1,ind.lon1)==TRUE)
    pred.rain1[l]<-idw.output1$var1.pred[ind.lat1[ind.sample1]]
  }
  
  for(k in 1:length(out.up1)){
    
    # Coord.Station2 are validated plus each outliers step by step
    Coord.Station.Out<-Coord.Station[out.up1[k],]
    Coord.Station.Out$Data<-as.numeric(M.Median.Prec[event.rain.out[j],out.up1[k]])
    Coord.Station2<-rbind(Coord.Station.V,Coord.Station.Out)
    Coord.Station2$y<-Coord.Station2$Latitud
    Coord.Station2$x<-Coord.Station2$Longitud
    coordinates(Coord.Station2) = ~x + y
    Mod.idw2<- idw(formula = Data ~ 1, locations = Coord.Station2, 
                   newdata = grd,debug.level=0)
    
    idw.output2 = as.data.frame(Mod.idw2)  # output is defined as a data table
    names(idw.output2)[1:3] <- c("long", "lat", "var1.pred")
    
    # ggplot() + geom_tile(data = idw.output, aes(x = long, y = lat, fill = var1.pred)) 
    
    pred.rain2<-c()
    
    for(l in 1:dim(Coord.Station2)[1]){
      ind.lat<-which(round(Coord.Station2$Latitud[l],2)==idw.output$lat)
      ind.lon<-which(round(Coord.Station2$Longitud[l],2)==idw.output$long)
      ind.sample<-which(is.element(ind.lat,ind.lon)==TRUE)
      pred.rain2[l]<-idw.output$var1.pred[ind.lat[ind.sample]]
    }
    
    set.valid<-pred.rain1[c(1:dim(Coord.Station.V)[1],dim(Coord.Station.V)[1]+k)]
    set.out<-pred.rain2
    
    final.out<-ifelse(wilcox.test(set.valid,set.out,paired = TRUE)$p.value<0.05,2,-1)
    
    if(final.out==2){
        Matrix.Outliers2[event.rain.out[j],out.up1[k]]<-final.out
        M.Median.Prec[event.rain.out[j],out.up1[k]]<-pred.rain1[dim(Coord.Station.V)[1]+1]
    }else{
      
    }
    print(paste("outlier ",k,sep="="))
  }
  
}

event.rain.def<-which(rowSums(Matrix.Outliers)>10)

for(i in 1:length(event.rain.def)){
  M.Median.Prec[event.rain.def[i],Matrix.Outliers[event.rain.def[i],]==1]<-NA
}

missdata<-importDVs("05054000", sdate="1985-01-01", edate="2009-12-31")

Final.Prec<-matrix(0,dim(M.Median.Prec)[1],dim(M.Median.Prec)[2])

for(l in 1:dim(M.Median.Prec)[2]){
  missdata[[2]]<-M.Median.Prec[,l]
  filldata<-fillMiss(missdata,block=30, pmiss=50)
  Final.Prec[,l]<-filldata$val
}

write.csv(cbind(H.days,Final.Prec),file=paste(root,"RESULTS","Precipitacion_prepo.csv",sep="/"))

#-----Fill Missing Data of Temperature-----

missdata[[2]]<-Tmean$CAJAMARCA_.21215100.
filldata<-fillMiss(missdata,block=70, pmiss=50)


plot(missdata$val,type="l",col=gray(0.75),lwd=0.6)
salida<-which(is.na(Tmean$CAJAMARCA_.21215100.)==TRUE)
lines(salida,rep(15,length(salida)),type="h",col="red")
lines(filldata$val,col="blue")

Matrix.Temp<-data.frame(Tme=filldata$val,Tmn=Tmin$CAJAMARCA_.21215100.,Tmx=Tmax$CAJAMARCA_.21215100.)
ind.na.Tmn<-which(is.na(Matrix.Temp$Tmn)==TRUE)
ind.na.Tmx<-which(is.na(Matrix.Temp$Tmx)==TRUE)

# Model for minimun temperatures
mod1<-nls(Tmn~A*Tme+B,data = Matrix.Temp[-ind.na.Tmn,],start = list(A=0.9,B=1))
Matrix.Temp[ind.na.Tmn,]<-coef(mod1)[1]*Matrix.Temp[ind.na.Tmn,1]+coef(mod1)[2]

# Model for maximum temperatures
mod2<-nls(Tmx~A*Tme+B,data = Matrix.Temp[-ind.na.Tmx,],start = list(A=0.9,B=1))
Matrix.Temp[ind.na.Tmx,]<-coef(mod2)[1]*Matrix.Temp[ind.na.Tmx,1]+coef(mod2)[2]

write.csv(cbind(H.days,Matrix.Temp),file=paste(root,"RESULTS","Temperatura_prepro.csv",sep="/"))

#-----Fill missing Data of Stream Flow-----

par(mfrow=c(5,2))

for(k in 2:10){
  par(mgp=c(2,0.5,0))
  par(mar=c(2.8,3.2,1,0.5)+0.1)
  plot(Flow[,k],col=rainbow(k)[k-1],type="l",main=names(Flow)[k],cex.main=0.8)
  salida<-which(is.na(Flow[,k])==TRUE)
  print(length(salida))
  lines(salida,rep(min(Flow[,k],na.rm=T)*1.1,length(salida)),col="black",type="h")
}


plot(rescale(Flow[,2],to=c(-1,1),range(Flow[,2:10],na.rm=T)),type="l",xlim = c(8000,9000),ylim=c(-1,0.7))
for(i in c(3,5,4,7,10,9,8))
lines(rescale(Flow[,i],to=c(-1,1),range(Flow[,2:10],na.rm=T)),col=i)


Flow.scale<-apply(Flow[,2:10],2,function(x) rescale(x,to=c(-1,1),range(x,na.rm=T)))



