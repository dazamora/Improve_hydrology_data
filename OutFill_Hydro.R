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
library(hydroGOF)
library(scales)

#------Load Data------
listdata<-list.files(paste(root,"/DATA/",sep=""))
# In situ observe
Prec<-as.data.frame(read.table(paste(root,"/DATA/",grep("Precipitacion",listdata,value = TRUE),sep=""),header=T))
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

png(filename =paste(root,"FIGURES","Temporal_Outliers.png",sep="/"), width = wth,height = hth1, 
    pointsize = 10, res = reso)

hist.out<-hist(rowSums(Matrix.Outliers[event.rain,]),breaks = c(0:16),labels = TRUE,xaxt="n",
     las=2,col = c(rep(gray(0.5),10),rep("red",6)),main="",xlab="Número de  estaciones outliers por evento",
     ylab ="Frecuencia", ylim = c(0,1200))

axis(1,at=hist.out$mids,labels = 1:16)

legend("topright",c("Evento atípico \n en la cuenca","Eventos atípicos en \n 10 o menos estaciones"),
       pch = 22 ,pt.bg=c("red",gray(0.5)),col="black",bty = "n",cex=0.8)

dev.off()

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
  #ggplot() + geom_tile(data = idw.output1, aes(x = long, y = lat, fill = var1.pred))
  
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
    
    # ggplot() + geom_tile(data = idw.output2, aes(x = long, y = lat, fill = var1.pred))
    
    pred.rain2<-c()
    
    for(l in 1:dim(Coord.Station2)[1]){
      ind.lat<-which(round(Coord.Station2$Latitud[l],2)==idw.output2$lat)
      ind.lon<-which(round(Coord.Station2$Longitud[l],2)==idw.output2$long)
      ind.sample<-which(is.element(ind.lat,ind.lon)==TRUE)
      pred.rain2[l]<-idw.output2$var1.pred[ind.lat[ind.sample]]
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

# Percentage of outliers spatial-temporal in precipitation time series
temporal.outlier<-apply(Matrix.Outliers[event.rain,],2,function(x)length(which(x==1)))
spatial.outlier<-apply(Matrix.Outliers2,2,function(x) length(which(x==2)))

OUT.PREC<-rbind(temporal.outlier,spatial.outlier)
rownames(OUT.PREC)<-c("Temporal","Spatial")
colnames(OUT.PREC)<-names(Prec)[2:17]

write.csv(t(OUT.PREC),file=paste(root,"RESULTS","Outliers_Precipitacion.csv",sep="/"))


#-----Fill outliers of all basin in precipitation time series-----

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

flow.stations<-c("Date","Payande","Carmen","PTE Carretera","Yuldaima","San Vicente",
                 "Cocora","PTE Luisa","Chuzo","PTE Bolívar")
colnames(Flow)<-flow.stations

at.year<-seq(min(H.days),max(H.days),by="year")

png(filename =paste(root,"FIGURES","Fill missing flow data.png",sep="/"), width = wth,height = hth1, 
    pointsize = 10, res = reso)

par(mfrow=c(5,2))

for(k in 2:10){
  par(mgp=c(1.1,0.5,0))
  par(mar=c(2.8,3.2,1,0.5)+0.1)
  plot(H.days,Flow[,k],col=rainbow(k)[k-1],type="l",main=names(Flow)[k],cex.main=0.8,las=2,cex.lab=0.7,xaxt="n",
       xlab="Tiempo [dia]",ylab=expression(paste("Caudal [m"^3,"/s]",sep="")),lwd=0.3,cex.lab=0.5,tck=-0.03,
       cex.axis=0.5,bty="n",lwd.ticks=0.5)
  salida<-which(is.na(Flow[,k])==TRUE)
  print(length(salida))
  lines(H.days[salida],salida,rep(max(Flow[,k],na.rm=T)*1.1,length(salida)),col=alpha("gray88",alpha = 0.1),
        type="h",lwd=0.5)
  axis.Date(side = 1, at = at.year, format="%Y",labels = TRUE, 
            tck=-0.03,las=2,cex.axis=0.5,col="gray45",lwd=0.5)
  box(lwd=0.5)
}

dev.off()

# Multivariate Flow Correlation Function  
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use = "na.or.complete"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.2)
}
# Figure about linear correlation between flow time series
png(filename =paste(root,"FIGURES","Correlation_Flows.png",sep="/"), width = wth,height = hth1, 
    pointsize = 10, res = reso)
par(mgp=c(2,0.5,0))
pairs(Flow[,2:10], lower.panel = panel.smooth, upper.panel = panel.cor,pch=16, 
      cex=0.3,col=gray(0.5),las=2,cex.axis=0.5,tck = -0.1)
dev.off()

# Fill missing data in Payande
missdata[[2]]<-Flow[,2]
filldata<-fillMiss(missdata,block=100, pmiss=50)
Flow[,2]<-filldata$val

# Fill missing data in Cocora

size.flow<-which(is.na(Flow[,7])==FALSE)

ind.cal<-sample(size.flow,round(length(size.flow)*0.6,2),replace = F)
ind.val<-size.flow[-ind.cal]
  
# mod3<-nls(Flow[ind,7]~A*Flow[ind,3]+B*Flow[ind,2]+C,start = list(A=0.1,B=0.5,C=0.2))
mod3<-nls(Flow[ind.cal,7]~A*Flow[ind.cal,3]+B,start = list(A=0.1,B=0.5))

Cocora.miss<-coef(mod3)[1]*Flow[-size.flow,3]+coef(mod3)[2]

Flow[-size.flow,7]<-Cocora.miss

# Fill missing data in Cocora PTE Bolivar
size.flow<-which(is.na(Flow[,10])==FALSE)



SSE<-c()
Coef.mod4<-c()
ind.cal<-list()

for(j in 1:1000){
  ind.cal[[j]]<-sample(size.flow,round(length(size.flow)*0.6,2),replace = F)
  ind.val<-size.flow[-ind.cal[[j]]]
  
  ind<-ind.cal[[j]]
  # mod4<-nls(Flow[ind.cal,10]~A*Flow[ind.cal,3]+B*Flow[ind.cal,4]+C,start = list(A=0.1,B=0.5,C=0.2))
  mod4<-nls(Flow[ind,10]~A*Flow[ind,3]+B,start = list(A=0.1,B=0.5))
  Coef.mod4<-rbind(Coef.mod4,coef(mod4))
  
  sim<-coef(mod4)[1]*Flow[ind.cal[[j]],3]+coef(mod4)[2]
  obs<-Flow[ind,10]
  SSE[j]<-NSE(sim,obs)
  
}

Coef.mod4<-Coef.mod4[which.max(SSE),]

Bolivar.miss<-Coef.mod4[1]*Flow[-size.flow,3]+Coef.mod4[2]

Flow[-size.flow,10]<-Bolivar.miss

# Save flow time series
write.csv(Flow,file=paste(root,"RESULTS","Caudales.csv",sep="/"))

#------Draw outliers in precipitation time series------
png(filename =paste(root,"FIGURES","Outliers_Precipitation_2.png",sep="/"), width = wth,height = hth1, 
    pointsize = 10, res = reso)
par(mfrow=c(2,1))

for(k in 3:4){
  
  i<-c(13,6,9,17)[k]
  par(mgp=c(1,0.4,0))
  par(mar=c(2.8,3.2,1,0.5)+0.1)
  plot(H.days,M.Median.Prec2[,i-1],type="l",col=alpha("gray84",alpha=0.7),las=2,xaxt="n",yaxt="n",lwd=0.3,
       main=c("Estación Cajamarca","Estación Pastales","Estación El Palmar","Estación Chicoral")[k],
       cex.main=0.7,ylab="Precipitación [mm/dia]", xlab="Tiempo [dia]",cex.lab=0.65,tck=-0.02,cex.axis=0.5,
       bty="n")
  
  ind.out.event<-which(Matrix.Outliers[,i-1]==1)
  ind.out.basin<-which(rowSums(Matrix.Outliers)>10 & Matrix.Outliers[,i-1]==1)
  ind.out.temspa<-which(Matrix.Outliers2[,i-1]==2)
  
  points(H.days[ind.out.event],M.Median.Prec2[ind.out.event,i-1],col="red",pch=2,cex=0.55,lwd=0.35)
  points(H.days[ind.out.basin],M.Median.Prec2[ind.out.basin,i-1],col="black",pch=16,cex=0.35,lwd=0.35)
  points(H.days[ind.out.temspa],M.Median.Prec2[ind.out.temspa,i-1],col="blue",pch=0,cex=0.45,lwd=0.35)
  axis.Date(side = 1, at = at.year, format="%Y",labels = TRUE, 
            tck=-0.02,las=2,cex.axis=0.5,lwd=0.5)
  axis(2,tck=-0.02,cex.axis=0.5,lwd=0.5,las=2)
  box(lwd=0.5)
  
  if(k==1){
    legend("top",c("Outlier Temporal-Estación","Outlier Temporal-Cuenca","Outlier Espacio-Temporal"),
           pch = c(2,16,0),col=c("red","black","blue"),cex=0.4,box.lwd = 0.5,pt.lwd =0.35)
  }
}

dev.off()







