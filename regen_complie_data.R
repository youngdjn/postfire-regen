setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(raster)
library(ncdf4)

source("regen_compile_data_functions.R")


#### 0. Operations that apply to all subsequent steps below ####

# Load plot locs
plot <- read.csv("data_survey/Plot_data.csv")

plot <- plot[plot$Fire != "ZACA",] # remove Zaca fire

plot <- plot[!is.na(plot$Northing),]
plot <- plot[(plot$Easting < 1105000) & (plot$Easting > 150000) & (plot$Northing > 2500000) & (plot$Northing < 6000000),] #Easting 110 used to be 83
plot.coords <- data.frame(x=plot$Easting,y=plot$Northing)
utm10 <- CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
plots <- SpatialPointsDataFrame(coords=data.frame(x=plot$Easting,y=plot$Northing),data=plot,proj4string=utm10)



#### 1. Extract monthly climate data for each plot ####

CPUs <- 3
years <- 1985:2014

topowx.dir <- "~/UC Davis/GIS/CA abiotic layers/TopoWx/"
tmin.dir <- "tmin_annual/"
tmax.dir <- "tmax_annual/"
tmax.normal.dir <- "tmax_normal/"
tmin.normal.dir <- "tmin_normal/"


prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"
ppt.dir <- "ppt_annual/"
ppt.normal.dir <- "ppt_normal/"


mo.num <- 1:12
mo.chr <- sprintf("%02d",mo.num)

### Open climate layers and store each month as a list element ###

### TopoWx tmin and tmax monthly-by-yearly ###

tmax.yr.mo <- list()
tmin.yr.mo <- list()
for(year in years) {
  
  tmax.file <- paste(topowx.dir,tmax.dir,"tmax_",year,".nc",sep="")
  tmax.year <- brick(tmax.file)
  
  tmin.file <- paste(topowx.dir,tmin.dir,"tmin_",year,".nc",sep="")
  tmin.year <- brick(tmin.file)
  
  for(month in 1:12) {
    
    tmax.layer.index <- paste("tmax",year,mo.chr[month],sep=".")
    tmin.layer.index <- paste("tmin",year,mo.chr[month],sep=".")
    
    tmax.yr.mo[tmax.layer.index] <- tmax.year[[month]]
    tmin.yr.mo[tmin.layer.index] <- tmin.year[[month]]
  }
}


### PRISM precip monthly-by-yearly ###

ppt.yr.mo <- list()

for(year in years) {
  
  ppt.yr.files <- paste(prism.dir,ppt.dir,"PRISM_ppt_stable_4kmM3_",as.character(year),mo.chr,"_bil.bil",sep="")
  ppt.yr <- stack(ppt.yr.files)
  
  for(month in 1:12) {
    
    layer.index <- paste("ppt",year,mo.chr[month],sep=".")
    ppt.yr.mo[layer.index] <- ppt.yr[[month]]
  }
}


### TopoWx tmin and tmax monthly normal ###

tmax.normal.mo <- list()
tmin.normal.mo <- list()

tmax.file <- paste(topowx.dir,tmax.normal.dir,"normals_tmax.nc",sep="")
tmax.normal <- brick(tmax.file,varname="tmax_normal")

tmin.file <- paste(topowx.dir,tmin.normal.dir,"normals_tmin.nc",sep="")
tmin.normal <- brick(tmin.file,varname="tmin_normal")

for(month in 1:12) {
  
  tmax.layer.index <- paste("tmax.normal",mo.chr[month],sep=".")
  tmin.layer.index <- paste("tmin.normal",mo.chr[month],sep=".")
  
  tmax.normal.mo[tmax.layer.index] <- tmax.normal[[month]]
  tmin.normal.mo[tmin.layer.index] <- tmin.normal[[month]]
}


### PRISM precip monthly normals ###

prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"

ppt.normal.mo <- list()

ppt.files <- paste(prism.dir,ppt.normal.dir,"PRISM_ppt_30yr_normal_800mM2_",mo.chr,"_bil.bil",sep="")
ppt.normal <- stack(ppt.files)

for(month in 1:12) {
  layer.index <- paste("ppt.normal",mo.chr[month],sep=".")
  ppt.normal.mo[layer.index] <- ppt.normal[[month]]
}


### stack everything of the same resolution and extent into one stack

temp.800m.rast <- unlist(list(tmin.yr.mo,tmax.yr.mo,tmin.normal.mo,tmax.normal.mo))
ppt.4km.rast <- unlist(list(ppt.yr.mo))
ppt.800m.rast <- unlist(list(ppt.normal.mo))


#### Extract climate values for points ####

# Load points
plot <- read.csv("data_survey/Plot_data.csv")
plot <- plot[!is.na(plot$Northing),]
plot <- plot[(plot$Easting < 1105000) & (plot$Easting > 150000) & (plot$Northing > 2500000) & (plot$Northing < 6000000),] #Easting 110 used to be 83
plot.coords <- data.frame(x=plot$Easting,y=plot$Northing)
utm10 <- CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
plots <- SpatialPointsDataFrame(coords=data.frame(x=plot$Easting,y=plot$Northing),data=plot,proj4string=utm10)

## Optional to write plots as GIS
# library(sp)
# library(rgdal)
# writeOGR(plots,dsn="plots_welch",layer=1,driver="ESRI Shapefile")


ppt.4km <- sapply(X = ppt.4km.rast,FUN = extract.single, plots = plots) # takes ~1 min
ppt.800m <- sapply(X = ppt.800m.rast,FUN = extract.single, plots = plots)
temp.800m <- sapply(X = temp.800m.rast,FUN = extract.single, plots = plots) # takes ~10 min

plots.clim <- data.frame(Regen_Plot=plots$Regen_Plot,temp.800m,ppt.4km,ppt.800m)
plots.clim <- plots.clim[complete.cases(plots.clim),]
write.csv(plots.clim,"../data_intermediate_processing_local/plot_climate_monthly.csv",row.names=FALSE)



#### 2. Summarize monthly climate by water year ####

plots.clim <- read.csv("../data_intermediate_processing_local/plot_climate_monthly.csv",header=TRUE)

first.year <- years[1]
last.year <- years[length(years)]

nyears <- last.year - first.year + 1
if(nyears < 3) {stop("You must supply at least 3 years.")}
n.water.years <- nyears - 1
inter.years <- years[2:(nyears-1)]
inter.years.rep <- rep(inter.years,each=12)

first.year.mos <- c("10","11","12")
last.year.mos <- c("01","02","03","04","05","06","07","08","09")
inter.mos <- rep(mo.chr,nyears-2)

first.year.mo.yr <- paste(first.year,first.year.mos,sep=".")
last.year.mo.yr <- paste(last.year,last.year.mos,sep=".")
inter.mos.mo.yr <- paste(inter.years.rep,inter.mos,sep=".")

all.mo.yr <- c(first.year.mo.yr,inter.mos.mo.yr,last.year.mo.yr)
all.mo.yr.tmin <- paste("tmin",all.mo.yr,sep=".")
all.mo.yr.tmax <- paste("tmax",all.mo.yr,sep=".")
all.mo.yr.ppt <- paste("ppt",all.mo.yr,sep=".")

tmin.col.nums <- sapply(all.mo.yr.tmin,function(x) {grep(x,names(plots.clim))})
tmax.col.nums <- sapply(all.mo.yr.tmax,function(x) {grep(x,names(plots.clim))})
ppt.col.nums <- sapply(all.mo.yr.ppt,function(x) {grep(x,names(plots.clim))})

tmin.col.nums.mat <- matrix(tmin.col.nums,ncol=12,byrow=TRUE)
tmax.col.nums.mat <- matrix(tmax.col.nums,ncol=12,byrow=TRUE)
ppt.col.nums.mat <- matrix(ppt.col.nums,ncol=12,byrow=TRUE)

water.year.summary.annual <- matrix(nrow=nrow(plots.clim),ncol=0)
for (i in 1:n.water.years) {
  
  water.year.ending <- years[i]+1
  
  tmin.col.nums.year <- tmin.col.nums.mat[i,]
  tmin.cols.year <- plots.clim[,tmin.col.nums.year]
  
  tmax.col.nums.year <- tmax.col.nums.mat[i,]
  tmax.cols.year <- plots.clim[,tmax.col.nums.year]
  
  tmean.cols.year <- (tmin.cols.year + tmax.cols.year) / 2
  
  ppt.col.nums.year <- ppt.col.nums.mat[i,]
  ppt.cols.year <- plots.clim[,ppt.col.nums.year]
  
  ### snow: for each plot, sum ppt from sept (col 1) to may (col 9) for months when avg temp was <=0
  snow <- ppt.cols.year
  snow[tmean.cols.year > 0] <- 0 # no snow if temp is > 0
  
  ### rain: for each plot, sum ppt from jan (col 5) to aug (col 12) for months when avg temp was >0
  rain <- ppt.cols.year
  rain[tmean.cols.year <= 0] <- 0

  tmin.avg.year <- apply(tmin.cols.year,1,mean)
  tmax.avg.year <- apply(tmax.cols.year,1,mean)
  
  ppt.tot.year <- apply(ppt.cols.year,1,sum)
  
  snow.tot.year <- apply(snow[c(1:8,12)],1,sum)
  rain.tot.year <- apply(rain[4:11],1,sum)
  
  tmean.avg.year <- (tmin.avg.year+tmax.avg.year)/2
  
  # calculate JJA average tmin, tmax, tmean
  # calculate DJF averate tmin, tmax, tmean
  tmean.avg.JJA <- apply(tmean.cols.year[,9:11],1,mean)
  tmean.avg.DJF <- apply(tmean.cols.year[,3:5],1,mean)

  ### merge into DF
  colname.prefixes <- c("tmin","tmax","tmean","tmean.JJA","tmean.DJF","ppt","snow","rain")
  colnames <- paste(colname.prefixes,water.year.ending,sep=".")
  
  clim.year <- cbind(tmin.avg.year,tmax.avg.year,tmean.avg.year,tmean.avg.JJA,tmean.avg.DJF,ppt.tot.year,snow.tot.year,rain.tot.year)
  colnames(clim.year) <- colnames
  
  water.year.summary.annual <- cbind(water.year.summary.annual,clim.year)
}

### summarize normal values annually (in the case of normals, calendar and water year are identical) ###

tmin.col.nums <- grep("tmin.normal",names(plots.clim))
tmax.col.nums <- grep("tmax.normal",names(plots.clim))
ppt.col.nums <- grep("ppt.normal",names(plots.clim))

tmin.cols <- plots.clim[,tmin.col.nums]
tmax.cols <- plots.clim[,tmax.col.nums]
ppt.cols <- plots.clim[,ppt.col.nums]

tmean.cols <- (tmin.cols + tmax.cols)/2

### snow: for each plot, sum ppt from nov (col 3) to may (col 9)
snow <- ppt.cols
snow[tmean.cols > 0] <- 0 # no snow if temp is > 0

### rain: for each plot, sum ppt from jan (col 5) to aug (col 12) for months when avg temp was >0
rain <- ppt.cols
rain[tmean.cols <= 0] <- 0

tmin.avg.normal.annual <- apply(tmin.cols,1,mean)
tmax.avg.normal.annual <- apply(tmax.cols,1,mean)
tmean.avg.normal.annual <- (tmin.avg.normal.annual + tmax.avg.normal.annual) / 2

## tmean JJA and DJF
tmean.avg.normal.DJF <- apply(tmean.cols[,c(12,1,2)],1,mean)
tmean.avg.normal.JJA <- apply(tmean.cols[,6:8],1,mean)

ppt.tot.normal.annual <- apply(ppt.cols,1,sum)
snow.tot.normal.annual <- apply(snow[,c(11,12,1,2,3,4,5)],1,sum)
rain.tot.normal.annual <- apply(rain[,c(3,4,5,6,7,8,9)],1,sum)

water.year.summary.normal <- cbind(tmin.avg.normal.annual,tmax.avg.normal.annual,tmean.avg.normal.annual,tmean.avg.normal.JJA,tmean.avg.normal.DJF,ppt.tot.normal.annual,snow.tot.normal.annual,rain.tot.normal.annual)
colnames(water.year.summary.normal) <- c("tmin.normal.ann","tmax.normal.ann","tmean.normal.ann","tmean.normal.JJA","tmean.normal.DJF","ppt.normal.ann","snow.normal.ann","rain.normal.ann")
water.year.summary <- data.frame(Regen_Plot=plots.clim$Regen_Plot,water.year.summary.annual,water.year.summary.normal)

write.csv(water.year.summary,"data_intermediate_processing/plot_climate_water_year.csv",row.names=FALSE)








#### 3. Extract fire severity for each plot ####

fire.severity.rasters.dir <- "C:/Users/dyoung/Documents/UC Davis/GIS/MTBS fire severity/"


fires <- unique(plots$Fire)

plot.firesev <- data.frame(Regen_Plot=NULL,fire.sev=NULL)
counter <- 1
for(fire in fires) {
  
  plots.fire <- plots[plots$Fire == fire,]
  fire.raster.dir <- paste(fire.severity.rasters.dir,fire,"/",sep="")
  
  # get the name of the DNBR6 layer
  file.list <- try(list.files(fire.raster.dir))
  
  
  raster.file <- file.list[grep("dnbr6.?.tif$",file.list)]
  
  if(length(raster.file) != 1) {
    cat("More than or less than one dnbr6 layer for Fire",fire,"\n")
    counter <- counter + 1
    next()
  }
  
  sev.layer <- raster(paste(fire.raster.dir,raster.file,sep=""))
  
  firesev <- extract(sev.layer,plots.fire)
  
  sev.layer.lowsev.only <- sev.layer
  sev.layer.lowsev.only[sev.layer.lowsev.only %in% c(4)] <- NA # exclude high only so we compute "distance to low or moderate"
  sev.layer.lowsev.only[!is.na(sev.layer.lowsev.only)] <- 0
  dist.to.low.rast <- distance(sev.layer.lowsev.only) #! make sure distance units are correct
  
  plots.fire <- spTransform(plots.fire,CRS(projection(sev.layer)))
  
  
  dist.to.low <- extract(dist.to.low.rast,plots.fire,method="bilinear")
  regen.plot <- plots.fire$Regen_Plot
  
  fire.vals <- data.frame(Regen_Plot=regen.plot,firesev=firesev,dist.to.low=dist.to.low) #! notes that value of 5 here means increased greenness
  plot.firesev <- rbind(plot.firesev,fire.vals)
  
  cat("Finished severity analysis of fire ",fire," (",counter," of ",length(fires),")\n",sep="")
  counter <- counter + 1
}

write.csv(plot.firesev,"data_intermediate_processing/plot_fire_data.csv",row.names=FALSE)


