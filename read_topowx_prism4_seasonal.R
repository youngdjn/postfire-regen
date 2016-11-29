setwd("C:/Users/DYoung/Dropbox/Research projects/Post-fire regen shared")

library(raster)
library(ncdf4)

CPUs <- 3

years <- 1985:2014

topowx.dir <- "~/UC Davis/GIS/CA abiotic layers/TopoWx/"
tmin.dir <- "tmin_annual/"
tmax.dir <- "tmax_annual/"

mo.num <- 1:12
mo.chr <- sprintf("%02d",mo.num)

prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"
ppt.dir <- "ppt_annual/"


### TopoWx tmin and tmax monthly-yearly ###

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

#tmax.yr.mo.stack <- stack(tmax.yr.mo)
#tmin.yr.mo.stack <- stack(tmin.yr.mo)



### PRISM precip monthly-yearly ###

prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"
ppt.dir <- "ppt_annual/"

ppt.yr.mo <- list()

for(year in years) {
  
    ppt.yr.files <- paste(prism.dir,ppt.dir,"PRISM_ppt_stable_4kmM3_",as.character(year),mo.chr,"_bil.bil",sep="")
    ppt.yr <- stack(ppt.yr.files)
    
    for(month in 1:12) {
      
      layer.index <- paste("ppt",year,mo.chr[month],sep=".")
      ppt.yr.mo[layer.index] <- ppt.yr[[month]]
    }
}
#ppt.yr.mo.stack <- stack(ppt.yr.mo)


### TopoWx tmin and tmax monthly normal ###

tmax.normal.mo <- list()
tmin.normal.mo <- list()

tmax.file <- paste(topowx.dir,"tmax_normal/normals_tmax.nc",sep="")
tmax.normal <- brick(tmax.file,varname="tmax_normal")

tmin.file <- paste(topowx.dir,"tmin_normal/normals_tmin.nc",sep="")
tmin.normal <- brick(tmin.file,varname="tmin_normal")

for(month in 1:12) {
  
  tmax.layer.index <- paste("tmax.normal",mo.chr[month],sep=".")
  tmin.layer.index <- paste("tmin.normal",mo.chr[month],sep=".")
  
  tmax.normal.mo[tmax.layer.index] <- tmax.normal[[month]]
  tmin.normal.mo[tmin.layer.index] <- tmin.normal[[month]]
}

#tmax.normal.mo.stack <- stack(tmax.normal.mo)
#tmin.normal.mo.stack <- stack(tmin.normal.mo)


### PRISM precip monthly normals ###

prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"

ppt.normal.mo <- list()
  
ppt.files <- paste(prism.dir,"ppt_normal/PRISM_ppt_30yr_normal_800mM2_",mo.chr,"_bil.bil",sep="")
ppt.normal <- stack(ppt.files)

for(month in 1:12) {
  layer.index <- paste("ppt.normal",mo.chr[month],sep=".")
  ppt.normal.mo[layer.index] <- ppt.normal[[month]]
}
#ppt.normal.mo.stack <- stack(ppt.normal.mo)
  





### stack everything of the same resolution and extent into one stack

temp.800m.rast <- unlist(list(tmin.yr.mo,tmax.yr.mo,tmin.normal.mo,tmax.normal.mo))
ppt.4km.rast <- unlist(list(ppt.yr.mo))
ppt.800m.rast <- unlist(list(ppt.normal.mo))


#### Extract climate values for points ####

# Load points
plot <- read.csv("Data/Regen DB tables/Plot_data.csv")
plot <- plot[!is.na(plot$Northing),]
plot <- plot[(plot$Easting < 1105000) & (plot$Easting > 150000) & (plot$Northing > 2500000) & (plot$Northing < 6000000),] #Easting 110 used to be 83
plot.coords <- data.frame(x=plot$Easting,y=plot$Northing)
utm10 <- CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
plots <- SpatialPointsDataFrame(coords=data.frame(x=plot$Easting,y=plot$Northing),data=plot,proj4string=utm10)

# library(sp)
# library(rgdal)
# writeOGR(plots,dsn="plots_welch",layer=1,driver="ESRI Shapefile")


extract.single <- function(layer, plots) {
  vals <- extract(layer,plots,method="bilinear")
  return(vals)
}

ppt.4km <- sapply(X = ppt.4km.rast,FUN = extract.single, plots = plots)
ppt.800m <- sapply(X = ppt.800m.rast,FUN = extract.single, plots = plots)
temp.800m <- sapply(X = temp.800m.rast,FUN = extract.single, plots = plots)

plots.clim <- data.frame(Regen_Plot=plots$Regen_Plot,temp.800m,ppt.4km,ppt.800m)
plots.clim <- plots.clim[complete.cases(plots.clim),]
write.csv(plots.clim,"Data/plot_climate.csv",row.names=FALSE)



#### Summarize annual values by water year
plots.clim <- read.csv("Data/plot_climate.csv",header=TRUE)

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

#### summarize normal values annually ####
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

write.csv(water.year.summary,"Data/plot_climate_water_year_summary.csv",row.names=FALSE)




                            