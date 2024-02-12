
rm(list=ls())
library(raster)
library(sf)

# Se lee la tabla provista para esta prueba
tabla=read.delim('~/Downloads/0243969-200613084148143.csv', header=T)

# Selección de columnas útiles para el resto del código
tabla=data.frame(Sp=tabla$species, long=tabla$decimalLongitude, lat=tabla$decimalLatitude)

tabla=na.omit(tabla) #remover registros para los que no se conoce la especies, o para los cuales no hay coordenada completas

especies_lista=unique(tabla$Sp) #lista de especies en la base de datos
especies_lista=especies_lista[-which(especies_lista=='')]

# Cargar shapefile de América
America=st_read('~/America_Border/America_Border.shp', 'America_Border')
America <- America[1] 
America=as_Spatial(America)

# Creación de una raster vacío que será usado en la rasterización de los datos
plantilla=raster(ncol=102, nrow=108, xmn=-120, xmx=-30, ymn=-55, ymx=30)
plantilla[]=runif(ncell(plantilla))

# Creación de una lista donde se guardará el raster para cada especie, y una matríz de comunidad vacía.
lista_mapas<-list()
comunidad=as.data.frame(plantilla)

# El loop trabaja de la siguiente forma. Pasa especies por especie, hace un subset de los puntos por cada especies, los cuales se rasterizan de acuerdo a la raster plantilla vacío. Cada raster es posteriormente transformado a data.frame para ser unido a la matríz de comunidades.

for (i in 1:length(especies_lista))

{
	especie=especies_lista[i]
	
	coor_especie=subset(tabla, Sp==especie)
	coordinates(coor_especie)=c('long', 'lat')
	
	map_especie=rasterize(coor_especie, plantilla, field=1)
	names(map_especie)<-especie
	
	lista_mapas[[i]]<-map_especie
	com=as.data.frame(map_especie)
	comunidad<-cbind(comunidad, com)
	
}

# Se hace la suma por pixel para el cálculo de riqueza de especies
comunidad$layer<-NULL
riqueza=rowSums(comunidad, na.rm=T)
Riqueza_mapa=plantilla
values(Riqueza_mapa)<-riqueza

Riqueza_mapa[which(values(Riqueza_mapa)<1)]<- NA

# Generación del mapa
plot(Riqueza_mapa)
plot(America, add=T)

#####################
# Beta
#####################

library(BAT)
library(picante)

comunidad[is.na(comunidad)]<-0

# Generación de matríz de comunidades vacía para el cálculo de Beta
Beta_Diversity=plantilla
Beta_Diversity[]<-NA


## IMPORTANTE. Este código solo calcula betadiversidad para pixeles con más de res especies. Esto lo decidí porque calcular un beta con menos especies puede no tener sentido biológico.

Cells=which(values(Riqueza_mapa)>3)
  
  
  ## En términos generales, el loop para por cada pixel con más de tres especies, selecciona sus pixeles vecinos en las 8 direcciones, y hace un calculo de beta entre los 8 para. Finalmente se calcula un promedio entre los 8, cuyo valor regresa al pixel focal.
  for (h in 1:length(Cells))
    
  {
    
    print(paste(((round(h/length(Cells), 2))*100), "%", sep=" "))
    
    vecinos=adjacent(plantilla, Cells[h], directions=8)
    
    Vec=c()
    
    for (p in 1:nrow(vecinos))
      
    {
      AA=vecinos[p,]
      
      Vec=c(Vec, betadiver(comunidad[AA,], method="e")[1])
      #Vec1=c(Vec1, betadiver(community_2[AA,], method="I")[1])
      #Vec2=c(Vec2, betadiver(community_2[AA,], method="m")[1])
      #Vec3=c(Vec3, betadiver(community_2[AA,], method="19")[1])
      #Vec=c(Vec, vegdist(community_2[AA,which(colSums(community_2[AA,], na.rm = T)>0)], method="jaccard", binary = TRUE)[1]) #or method='chao'
      
      
    }
    
    Beta_Diversity[Cells[h]]<-mean(Vec)
    
   }

# Generación del mapa
plot(Beta_Diversity)
plot(America, add=T)


