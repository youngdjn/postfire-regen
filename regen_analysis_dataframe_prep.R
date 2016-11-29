setwd("C:/Users/DYoung/Dropbox/Research projects/Post-fire regen shared")

library("reshape")
library(ggplot2)
library(sp)
library(rgdal)
library(raster)

source("R/regen_utils_2_zscore.R")



### These only need to be run once per database export
##source("R/read_topowx_prism4_seasonal.R") # read in monthly temp and precip and save monthly values and water year averages (and normals) for each plot, save to Data/plot_climate.csv and Data/plot_climate_water_year_summary.csv
##source("R/extract_fire_sev.R") # extract fire severity and distance to low severity for each plot, save to Data/plot_fire_data.csv
##source("R/summarize_regen.R") # summarizes seedling and spling counts by age and species for each plot, and for each species and each plot finds the nearest seed tree (if any listed) # normal for it to generate many sap reps 

### Read in the summarized data files
plot.climate <- read.csv("Data/plot_climate_water_year_summary.csv",stringsAsFactors=FALSE)
plot.fire.data <- read.csv("Data/plot_fire_data.csv",stringsAsFactors=FALSE)
regen <- read.csv("Data/regen_seedtr_plot_sp.csv",stringsAsFactors=FALSE)
plot.seed.tree.any <- read.csv("Data/regen_seedtr_plot_any.csv")


ppt.cols <- grep("ppt",names(plot.climate))
ppt.df <- plot.climate[,ppt.cols]
ppt.df <- ppt.df[,15:29]
ppt.row <- ppt.df[plot.climate$Regen_Plot == "AMR1300134",]
ppt.row

ppt.df.new <- data.frame(plot=names(ppt.row),precip=t(ppt.row))

plot(ppt.df.new)

### Read in raw data files (direct from DB export)
fire.years <- read.csv("Data/Regen DB tables/fire_years.csv",header=TRUE,stringsAsFactors=FALSE)
plot <- read.csv("Data/Regen DB tables/Plot_data.csv",header=TRUE,stringsAsFactors=FALSE)

## get the survey year (only) out of the survey date
date.split <- strsplit(plot$Date,"/")
yearplus <- sapply(date.split, '[[', 3)
year.split <- strsplit(yearplus," ")
plot$Year <- as.numeric(sapply(year.split,'[[',1))



### merge the plot-level data into a single DF
plot.1a <- merge(plot,plot.fire.data,by="Regen_Plot",all.x=TRUE)
plot.1b <- merge(plot.1a,plot.seed.tree.any,by="Regen_Plot",all.x=TRUE)
plot.2 <- merge(plot.1b,fire.years,by="Fire",all.x=TRUE)
plot.2$survey.years.post <- plot.2$Year - plot.2$fire.year

### thin to only plots that were surveyed 5 years post-fire
#survey.years.post.fire <- 5
#plot.3 <- plot.2[plot.2$survey.years.post == survey.years.post.fire,]
#plot.3 <- plot.3[!is.na(plot.3$Regen_Plot),]

plot.3 <- plot.2[!is.na(plot.2$Regen_Plot),]


### get summarized climate data for each plot
plot.3.clim <- summarize.clim(plot.3,plot.climate,years.clim=1:3)
plot.3.clim2 <- summarize.clim(plot.3,plot.climate,years.clim=4:5)
names(plot.3.clim2) <- paste(names(plot.3.clim2),".late",sep="")

### get summarized regen data for each plot
plot.3.regen <- summarize.regen.ind(plot.3,regen,sp=c("ABCO","PSME","PIPO"),years.regen=4:5,all.sp=TRUE)
#plot.3.regen <- summarize.regen(plot.3,regen,sp=c("ABCO","PSME","PIPO","PILA","PIJE"),years.regen=4:5)

### merge the plot, climate, and regen data
plot.4 <- merge(plot.3,plot.3.clim,by="Regen_Plot",all.x=TRUE)
plot.4.1 <- merge(plot.4,plot.3.clim2,by.x="Regen_Plot",by.y="Regen_Plot.late",all.x=TRUE)
plot.5 <- merge(plot.4.1,plot.3.regen,by="Regen_Plot",all.x=TRUE)

d <- plot.5

d <- d[d$Regen_Plot != "SHR09000015",]


#### extract calveg data (in development--takes too long to open shapefiles) ####

# # explort csv of plots to open in QGIS
# plots_for_calveg <- d[,c("Regen_Plot","Northing","Easting")]
# 
# UTM10N <- CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# albers.proj <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
# plots.for.calveg <- SpatialPointsDataFrame(cbind(d$Easting,d$Northing),data=d[,c("Regen_Plot","Fire")],proj4string=UTM10N)
# plots.for.calveg <- spTransform(plots.for.calveg,albers.proj)
# writeOGR(plots.for.calveg,"plots_for_calveg","plots_for_calveg",driver="ESRI Shapefile",overwrite_layer=TRUE)
# 
# 
# write.csv(plots_for_calveg,"plots_for_calveg.csv",row.names=FALSE)



#### extract elevation ####
update.elev <- FALSE

if(update.elev) {
  dem <- raster("C:/Users/DYoung/Documents/UC Davis/GIS/CA abiotic layers/DEM/Camerged9_SNclip_nd.tif")
  UTM10N <- CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  d.unique <- d[d$species=="ALL",]
  plots.for.calveg <- SpatialPointsDataFrame(cbind(d.unique$Easting,d.unique$Northing),data=d.unique,proj4string=UTM10N)
  plots.for.calveg$elev.m <- extract(dem,plots.for.calveg,method="bilinear")
  write.csv(plots.for.calveg,"Data/plots_elev.csv")
}

plots.elev <- read.csv("Data/plots_elev.csv")
plots.elev <- plots.elev[,c("Regen_Plot","elev.m")]

d <- merge(d,plots.elev,all.x=TRUE)


write.csv(d,"Data/regen_clim_full.csv",row.names=FALSE)


