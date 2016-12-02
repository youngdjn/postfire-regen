setwd("C:/Users/DYoung/Dropbox/Research projects/Post-fire regen shared")

library(sp)
library(raster)
library(reshape2)
library(rgdal)
library(rgeos)
library(maptools)

plot <- read.csv("Data/Regen DB tables/Plot_data.csv")

# remove plots with no coordinates
plot <- plot[!is.na(plot$Northing),]
plot <- plot[(plot$Easting < 835000) & (plot$Easting > 150000) & (plot$Northing > 2500000) & (plot$Northing < 6000000),]
plot.coords <- data.frame(x=plot$Easting,y=plot$Northing)


utm10 <- CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")
plots <- SpatialPointsDataFrame(coords=plot.coords,data=plot,proj4string=utm10)

fires <- unique(plots$Fire)

plot.firesev <- data.frame(Regen_Plot=NULL,fire.sev=NULL)
counter <- 1
for(fire in fires) {
  
  plots.fire <- plots[plots$Fire == fire,]
  raster.dir <- paste("C:/Users/DYoung/Documents/UC Davis/GIS/MTBS fire severity/",fire,"/",sep="")
  
  # get the name of the DNBR6 layer
  file.list <- try(list.files(raster.dir))
  
  
  raster.file <- file.list[grep("dnbr6.?.tif$",file.list)]
  
  if(length(raster.file) != 1) {
    cat("More than or less than one dnbr6 layer for Fire",fire,"\n")
    counter <- counter + 1
    next()
  }
  
  sev.layer <- raster(paste(raster.dir,raster.file,sep=""))
  
  firesev <- extract(sev.layer,plots.fire)
  
  sev.layer.lowsev.only <- sev.layer
  sev.layer.lowsev.only[sev.layer.lowsev.only %in% c(4)] <- NA
  sev.layer.lowsev.only[!is.na(sev.layer.lowsev.only)] <- 0
  dist.to.low.rast <- distance(sev.layer.lowsev.only)
  
  plots.fire <- spTransform(plots.fire,CRS(projection(sev.layer)))
  

  dist.to.low <- extract(dist.to.low.rast,plots.fire,method="bilinear")
  regen.plot <- plots.fire$Regen_Plot
  
  fire.vals <- data.frame(Regen_Plot=regen.plot,firesev=firesev,dist.to.low=dist.to.low)
  plot.firesev <- rbind(plot.firesev,fire.vals)
  
  cat("Finished severity analysis of fire ",fire," (",counter," of ",length(fires),")\n",sep="")
  counter <- counter + 1
}

write.csv(plot.firesev,"Data/plot_fire_data.csv",row.names=FALSE)

