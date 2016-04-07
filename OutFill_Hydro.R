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

#------Draw Time Series------
par(mfrow=c(8,2))
colr<-rainbow(18)
for(i in 2:17){
  par(mgp=c(2,0.5,0))
  par(mar=c(2.8,3.2,1,0.5)+0.1)
  
  plot(Prec[,i],type="l",col=colr[i],xlab="",ylab="")
  
}

#-----Fill missing data with satellite information----





#----Aggregate Time Series and outlier detection-----


mo <- strftime(H.days, "%m") # Clasification by months
Ma.Prec <- data.frame(mo, Prec[,2:17])
Mon<-unique(mo)

Matrix.Outliers<-matrix(0,dim(Prec)[1],dim(Prec)[2]-1)

for(i in 1:(dim(Prec)[2]-1)){
  Prec.agg <- aggregate(Ma.Prec[[i+1]]~ mo, Ma.Prec,boxplot, na.rm=T,plot=FALSE)
  for(j in 1:12){
    ind<-which(Mon[j]==mo)
    outliers<-ind[is.element(el=Prec[ind,i],set = Prec.agg$`Ma.Prec[[i + 1]]`[j,4]$out)]
    Matrix.Outliers[outliers,i]<-1
    rm(outliers)
    }
}








