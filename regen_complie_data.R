setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(raster)
library(ncdf4)
library(sp)
library(reshape2)
library(plyr)
library(rgdal)

source("regen_compile_data_functions.R")


#### -1. Read in and merge the various databases and excel sheets with plot data ####

db.dirs <- c("../data_survey/Cameron/","../data_survey/Jared/","../data_survey/John/","../data_survey/Michelle/")



for(i in 1:length(db.dirs)) {
  
  db.dir <- db.dirs[i]
  
  table.file <- "Plot_data.txt"
  table.loc <- paste0(db.dir,table.file)
  plot <- read.csv(table.loc,header=TRUE,stringsAsFactors=FALSE)

  table.file <- "surviving_trees.txt"
  table.loc <- paste0(db.dir,table.file)
  surviving.trees <- read.csv(table.loc,header=TRUE,stringsAsFactors=FALSE)
  
  if(i != 4) { # from Michelle we don't need these because hers were all control plots
    table.file <- "Resprouts.txt"
    table.loc <- paste0(db.dir,table.file)
    resprout <- read.csv(table.loc,header=TRUE,stringsAsFactors=FALSE)
    
    table.file <- "sapling_regen.txt"
    table.loc <- paste0(db.dir,table.file)
    sap <- read.csv(table.loc,header=TRUE,stringsAsFactors=FALSE)
    
    table.file <- "tree_regen.txt"
    table.loc <- paste0(db.dir,table.file)
    seedl <- read.csv(table.loc,header=TRUE,stringsAsFactors=FALSE)
  }

  if(i == 1) {
    
    plot.comb <- plot
    resprout.comb <- resprout
    sap.comb <- sap
    seedl.comb <- seedl
    surviving.trees.comb <- surviving.trees

  } else {
    
    plot.comb <- rbind.fill(plot.comb,plot)
    resprout.comb <- rbind.fill(resprout.comb,resprout)
    sap.comb <- rbind.fill(sap.comb,sap)
    seedl.comb <- rbind.fill(seedl.comb,seedl)
    surviving.trees.comb <- rbind.fill(surviving.trees.comb,surviving.trees)

  }
  
}


plot.comb$seed_tree_distance_general <- min(c(plot.comb$seed_tree_distance_conifer, plot.comb$seed_tree_distance_hardwood, plot.comb$seed.tree.distance.general), na.rm=TRUE)

# Find when only one quadrant was surveyed for seedlings, and multiply counts by 4 when true
seedl.count.columns <- grep("yr$",names(seedl.comb))
all.quadrants <- (seedl.comb$quadrants %in% c("ALL","4",""," ")) | is.na(seedl.comb$quadrants) # all quadrants were surveyed (T/F)
seedl.comb[!all.quadrants,seedl.count.columns] <- 4*seedl.comb[!all.quadrants,seedl.count.columns]


#populate seedling unk_yr column with values for CADE and all hardwoods
hardwoods <- c("QUKE","QUCH2","ARME","LIDE3","CHCH","QUGA4","ACMA","CEMO2","CONU4","POTR5","QUBE5","QUJO3","QUWI","UMCA")
not.ageable <- c("CADE27",hardwoods)
seedl.comb$Count_total <- rowSums(seedl.comb[,seedl.count.columns])
seedl.comb[seedl.comb$Species %in% not.ageable,"unk_yr"] <- (seedl.comb$unk_yr + seedl.comb$Count_total)[seedl.comb$Species %in% not.ageable]





#! when a 2016 plot had a shrub that was actually a hardwood, fake it so it counts as a hardwood (make fake hardwood column that is counted when calculating hardwood presence but not count)
#! exclude 2016 plots from hardwood count: not 0 but NA, not part of analysis


#! Resume here: compare tables with Welch, then open Welch tables and merge in






#! fix species names to have numbers: once tables merged, sort by species name to identify where.







#### 0.5 Operations that apply to all eteps below #####



# Load plot locs
plot <- read.csv("../data_survey/Plot_data.csv",stringsAsFactors=FALSE)

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


### Extract climate values for points ###
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
  
  fire.vals <- data.frame(Regen_Plot=regen.plot,firesev=firesev,dist.to.low=dist.to.low) #! note that value of 5 here means increased greenness
  plot.firesev <- rbind(plot.firesev,fire.vals)
  
  cat("Finished severity analysis of fire ",fire," (",counter," of ",length(fires),")\n",sep="")
  counter <- counter + 1
}

write.csv(plot.firesev,"../data_intermediate_processing_local/plot_fire_data.csv",row.names=FALSE)




#### 4. Summarize regen and surviving trees by species (and age for regen) for each plot ####

sap <- read.csv("../data_survey/sapling_regen.csv")
shrub <- read.csv("../data_survey/shrub_regen.csv")
seedl <- read.csv("../data_survey/tree_regen.csv")
resprout <- read.csv("../data_survey/Resprouts.csv")
surviving.trees <- read.csv("../data_survey/surviving_trees.csv",header=TRUE,stringsAsFactors=FALSE)


# specify only non-planted seedlings (assuming that if it's blank, it means it's seeded)
seedl <- seedl[which((seedl$seed_veg_plant != "P") | is.na(seedl$seed_veg_plant)),]

### check regenerating tree data for duplicates
seedl.rep <- aggregate(seedl[,10],by=list(species=seedl$Species,Regen_Plot=seedl$Regen_Plot),FUN=length)
seedl.rep.rows <- seedl.rep[seedl.rep$x>1,]
if(nrow(seedl.rep.rows)>0) {warning("Multiple seedling entries for some species-plots combinations. Duplicates listed in 'seedl.rep.rows'.")}

sap.rep <- aggregate(sap[,10],by=list(species=sap$Species,Regen_Plot=sap$Regen_Plot),FUN=length)
sap.rep.rows <- sap.rep[sap.rep$x>1,]
if(nrow(sap.rep.rows)>0) {warning("Multiple sapling entries for some species-plot combinations. Duplicates listed in 'sap.rep.rows'.")}

### aggregate seedling table by plot and species
#!!!! need to include unknown age here. This is once we add 2016 data which include that column
#! Can do that by making it a separate column that is only counted when computing TOTAL seedlings
#! Also, make all counts for CADE and all hardwoods fall in the unknown column.
seedl.ag <- aggregate(seedl[,5:16],by=list(species=seedl$Species,Regen_Plot=seedl$Regen_Plot),FUN=sum,na.rm=TRUE)
count.yrs <- paste("count.",0:11,"yr",sep="")
names(seedl.ag) <- c("species","Regen_Plot",count.yrs)

### aggregate sapling table by plot and species
sap$tot <- rowSums(sap[,6:14],na.rm=TRUE)
sap$X10yr[sap$tot == 0 | (is.na(sap$tot))] <- 1 # if the species had a row for the sapling but no age, assume it was 10yr; some saplings are not assigned an age; just interpret each row as a presence
sap.ag <- aggregate(sap[,c(6:14,20)],by=list(species=sap$Species,Regen_Plot=sap$Regen_Plot),FUN=sum,na.rm=TRUE)
count.yrs <- paste("count.",3:11,"yr",sep="")
names(sap.ag) <- c("species","Regen_Plot",count.yrs,"count.tot")
sap.ag <- sap.ag[,1:ncol(sap.ag)-1]

### aggregate resprout table by plot and species
resprout$COUNT.TOTAL <- ifelse(is.na(resprout$COUNT.TOTAL),1,resprout$COUNT.TOTAL) # sometimes it has a count, sometimes it has one record per individual (stem?)
resprout.ag <- aggregate(resprout$COUNT.TOTAL,by=list(species=resprout$Species,Regen_Plot=resprout$Regen_Plot),FUN=sum,na.rm=TRUE)
count.yrs <- paste("count.",5,"yr",sep="")
names(resprout.ag) <- c("species","Regen_Plot",count.yrs)

### merge seedling and sapling tables and resprout table
regen <- rbind.fill(seedl.ag,sap.ag)
regen <- rbind.fill(regen,resprout.ag) # add in resprout table (comment out here if desired)
regen.ag <- aggregate(regen[,3:14],by=list(species=regen$species,Regen_Plot=regen$Regen_Plot),FUN=sum,na.rm=TRUE)

### aggregate surviving tree table by plot and species
#! optionally add a cutoff here so we con't consider small trees
surviving.trees$ba <- (surviving.trees$DBH/2)^2*3.14
surviving.trees$count <- 1
surviving.trees.ag <- aggregate(surviving.trees[,c("count","ba")],by=list(species=surviving.trees$Species,Regen_Plot=surviving.trees$Regen_Plot),FUN=sum)
names(surviving.trees.ag) <- c("species","Regen_Plot","surviving.trees.count","surviving.trees.ba")

### merge regen and surviving tree DFs
regen.ag <- merge(regen.ag,surviving.trees.ag,by=c("Regen_Plot","species"),all=TRUE)


regen.ag$surviving.trees.count[is.na(regen.ag$surviving.trees.count)] <- 0
regen.ag$surviving.trees.ba[is.na(regen.ag$surviving.trees.ba)] <- 0

write.csv(regen.ag,"../data_intermediate_processing_local/tree_summarized_sp.csv",row.names=FALSE)





#### 5. Integrate all plot and species data ####

### Read in raw data files (direct from DB export)
fire.years <- read.csv("../data_fire/fire_years.csv",header=TRUE,stringsAsFactors=FALSE)
seed.tree <- read.csv("../data_survey/seed_tree.csv",header=TRUE,stringsAsFactors=FALSE)

### Read in the summarized data files
plot.climate <- read.csv("../data_intermediate_processing_local/plot_climate_water_year.csv",stringsAsFactors=FALSE)
plot.fire.data <- read.csv("../data_intermediate_processing_local/plot_fire_data.csv",stringsAsFactors=FALSE)
plot.tree.sp <- read.csv("../data_intermediate_processing_local/tree_summarized_sp.csv",stringsAsFactors=FALSE)

# if there is no 12yr column, add it to prevent bugs in summing tree ages
if(!("count.12yr" %in% names(plot.tree.sp))) {
  plot.tree.sp$count.12yr <- 0
}



## to plot ppt over time at any given plot
# ppt.cols <- grep("ppt",names(plot.climate))
# ppt.df <- plot.climate[,ppt.cols]
# ppt.df <- ppt.df[,15:29]
# ppt.row <- ppt.df[plot.climate$Regen_Plot == "AMR1300134",]
# ppt.df.new <- data.frame(plot=names(ppt.row),precip=t(ppt.row)) # start a new climate DF based on the old one
# plot(ppt.df.new)

## get the survey year (only) out of the survey date
date.split <- strsplit(plot$Date,"/")
yearplus <- sapply(date.split, '[[', 3)
year.split <- strsplit(yearplus," ")
plot$Year <- as.numeric(sapply(year.split,'[[',1))


## extract elevation for each plot
dem <- raster("C:/Users/DYoung/Documents/UC Davis/GIS/CA abiotic layers/DEM/Camerged9_SNclip_nd.tif")
plots$elev.m <- extract(dem,plots,method="bilinear")

## extract march solar rad for each plot
rad.march <- stack("C:/Users/DYoung/Documents/UC Davis/GIS/CA abiotic layers/solar rad/rsun nldas adjusted/glob_rad_monthly_dobr1_int.tif")[[3]]
plots$rad.march <- extract(rad.march,plots,method="bilinear")

#compile DF of geospatial data that was extracted for each plot
plots.extracted <- plots[,c("Regen_Plot","elev.m","rad.march")]


## merge the plot-level data into a single DF
plot.1 <- merge(plot,plot.fire.data,by="Regen_Plot",all.x=TRUE)
plot.2 <- merge(plot.1,fire.years,by="Fire",all.x=TRUE)
plot.2$survey.years.post <- plot.2$Year - plot.2$fire.year

plot.3 <- merge(plot.2,plots.extracted,all.x=TRUE)

### get summarized climate data for each plot
##!! If want to look at weather beyond 4 years post-fire (for those fires that had more than 4 years), will need to make this relative to number of years post-fire
plot.3.clim <- summarize.clim(plot.3,plot.climate,years.clim=1:3) #first three years after fire
plot.3.clim2 <- summarize.clim(plot.3,plot.climate,years.clim=3:4) # years 3-4 after fire
names(plot.3.clim2) <- paste(names(plot.3.clim2),".late",sep="")

### get summarized regen data for each plot
###!!! Need to account for unknown age here (will be a new column, created above; only add it to the "all" regen totals, not to young or old)
plot.3.regen.old <- summarize.regen.ind(plot.3,plot.tree.sp,sp=c("ABCO","PSME","PIPO"),regen.ages="old",all.sp=TRUE)
plot.3.regen.young <- summarize.regen.ind(plot.3,plot.tree.sp,sp=c("ABCO","PSME","PIPO"),regen.ages="young",all.sp=TRUE)[,1:3] #only take the regen data (because funct also outupts adults data but we get that from the first call, the previous line)
plot.3.regen.all <- summarize.regen.ind(plot.3,plot.tree.sp,sp=c("ABCO","PSME","PIPO"),regen.ages="all",all.sp=TRUE)[,1:3] #only take the regen data (because funct also outupts adults data but we get that from the first call)
names(plot.3.regen.old)[3] <- "regen.count.old"
names(plot.3.regen.young)[3] <- "regen.count.young"
names(plot.3.regen.all)[3] <- "regen.count.all"


plot.3.regen.pre <- merge(plot.3.regen.young,plot.3.regen.old)
plot.3.regen <- merge(plot.3.regen.pre,plot.3.regen.all)


### species-level seed tree distance
## if a species within a plot has multiple seed trees listed, get the shortest distance
seed.tree.sp <- aggregate(seed.tree$Dist_m,by=list(seed.tree$Regen_Plot,seed.tree$Species),FUN=min)
names(seed.tree.sp) <- c("Regen_Plot","species","seed.tree.sp")
plot.3.regen <- merge(plot.3.regen,seed.tree.sp,all.x=TRUE)
# species table ready for export


### add plot-level seed tree distance (shortest distance among all seed trees recorded for the plot)
seed.tree.any <- aggregate(seed.tree$Dist_m,by=list(seed.tree$Regen_Plot),FUN=min)
names(seed.tree.any) <- c("Regen_Plot","seed.tree.any")

### merge the plot, climate, and plot-level seed tree data
plot.clim <- merge(plot.3,plot.3.clim,by="Regen_Plot",all.x=TRUE)
plot.clim.seedtree <- merge(plot.clim,seed.tree.any,all.x=TRUE)

plot.clim.seedtree$FORB <- plot.clim.seedtree$FORBE
plot.clim.seedtree <- remove.vars(plot.clim.seedtree,"FORBE")


### write plot-level and species-level output files
write.csv(plot.clim.seedtree,"data_intermediate/plot_level.csv",row.names=FALSE)
write.csv(plot.3.regen,"data_intermediate/speciesXplot_level.csv",row.names=FALSE)

