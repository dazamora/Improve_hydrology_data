# Realiza los graficos de la interpolación IDW

library(rgeos)

#Coello.shape<-readShapePoly("C:/Users/DAVID ZAMORA/Desktop/Area_1.shp")


tempp<-as.data.frame(Coord.Station1[,2:4])
outlier<-as.data.frame(Coord.Station[out.up1,1:3])
tempp2<-as.data.frame(Coord.Station2[,2:4])


ggplot() + geom_tile(data = idw.output2, aes(x = long, y = lat, fill = var1.pred))+
  scale_fill_distiller(palette = "Spectral")+
  #geom_path(data = est_contour, aes(long, lat, group = group), colour = "grey") + 
  geom_point(data = tempp2, aes(x = Longitud, y = Latitud), shape = 21, 
             colour = "blue") + labs(fill = "Precipitacion", title = "Campo de precipitación de la cuenca Coello")+
geom_point(data = outlier, aes(x = Longitud, y = Latitud), shape = 21, 
           colour = "red") + labs(fill = "Precipitacion", title = "Campo de precipitación de la cuenca Coello")

ggsave("prediccion_con_outlier.png")
