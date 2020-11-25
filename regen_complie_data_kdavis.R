library(raster)
library(sp)
library(reshape2)
library(plyr)
library(rgdal)
library(sp)
library(data.table)

source("regen_compile_data_functions.R")



#### -2. Prepare Michelle (Power revisit) spreadsheet so it has same format as the DB export spreadsheets ####

## Merge in Michelle BA and percent cover to the Michelle plot sheet
plot.michelle.partial <- read.csv("data_survey/Michelle/Plot_data_partial.txt",header=TRUE, stringsAsFactors=FALSE)
#lifeform.michelle <- read.csv("../data_survey/Michelle/lifeform.txt",header=TRUE, stringsAsFactors=FALSE)
input.michelle <- read.csv("data_survey/Michelle/input.txt",header=TRUE, stringsAsFactors=FALSE)

## simplify input to one line per plot
input.min <- aggregate(input.michelle[,c("BAF","seed_tree_distance_conifer","seed_tree_distance_hardwood")],by=list(input.michelle$Regen_Plot),FUN=min,na.rm=TRUE)
input.min[input.min==Inf] <- NA # when there was no seed tree observed
input.sum <- aggregate(input.michelle[,c("BA_live_count","BA_dead_count")],by=list(input.michelle$Regen_Plot),FUN=sum)

names(input.sum)[1] <- names(input.min)[1] <- "Regen_Plot"

input.comp <- merge(input.min,input.sum,by="Regen_Plot")

## merge input and lifeform into main Michelle plot table
#plot.michelle.pre <- merge(plot.michelle.partial,lifeform.michelle,by="Regen_Plot",all.x=TRUE)
plot.michelle <- merge(plot.michelle.partial,input.comp,by="Regen_Plot",all.x=TRUE)

write.csv(plot.michelle,"data_survey/Michelle/Plot_data.txt",row.names=FALSE)


#### -1. Read in and merge the various databases and excel sheets with plot data ####


db.dirs <- c("data_survey/Cameron/","data_survey/Jared/","data_survey/John/","data_survey/Michelle/")



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
    surviving.trees.comb <- rbind.fill(surviving.trees.comb,surviving.trees)
    
    if(i == 4) next()
    
    resprout.comb <- rbind.fill(resprout.comb,resprout)
    sap.comb <- rbind.fill(sap.comb,sap)
    seedl.comb <- rbind.fill(seedl.comb,seedl)

  }
  
}



# 


# 
# ## Add Richter data
# 
# plot.richter <- read.csv("../data_survey/Clark/richter_plot_data.csv",header=TRUE, stringsAsFactors=FALSE)
# resprout.richter <- read.csv("../data_survey/Clark/richter_resprout.csv",header=TRUE, stringsAsFactors=FALSE)
# sapling.richter <- read.csv("../data_survey/Clark/richter_saplings.csv",header=TRUE, stringsAsFactors=FALSE)
# seedling.richter <- read.csv("../data_survey/Clark/richter_seedlings.csv",header=TRUE, stringsAsFactors=FALSE)
# 
# plot.richter$Plot.ID <- paste0("richter.",plot.richter$Plot.ID)
# plot.richter <- plot.richter[,c(1:3,6,7,13,14,19,20,21)]
# names(plot.richter)[1:3] <- c("Regen_Plot","Longitude","Latitude")
# plot.richter$Fire <- "POWER"
# 
# seedling.richter$Plot_ID <- paste0("richter.",seedling.richter$Plot_ID)
# seedling.richter <- seedling.richter[,-1]
# names(seedling.richter) <- c("Regen_Plot","Species","seed_veg_plant","X0yr","Ct_1yr","Ct_2yr","Ct_3yr","Ct_4yr","Ct_5yr","X6yr","X7yr","X8yr","X9yr","X10yr","X11yr","X12yr","X13yr","X14yr","X15yr","X16yr","X17yr","X18yr","X19yr","X20.yr",paste0("x",21:26,"yr"),"tallest_ht_cm","tallest_age","tallest_lastyr_cm")
# seedling.richter <- seedling.richter[!is.na(seedling.richter$Species) & seedling.richter$Species != "",]
# seedling.richter <- seedling.richter[,which(!is.na(names(seedling.richter)))]
# 
# sapling.richter$Plot_ID <- paste0("richter.",sapling.richter$Plot_ID)
# sapling.data.richter <- sapling.richter[,c(-1,-7:-10)]
# names(sapling.data.richter) <- c("Regen_Plot","Species","seed_veg_plant","DBH..cm.","age")
# sapling.tallest.richter <- sapling.richter[,c(2,7:10)]
# sapling.tallest.richter <- sapling.tallest.richter[!is.na(sapling.tallest.richter$TI_Species) & sapling.tallest.richter$TI_Species != "",]
# names(sapling.tallest.richter) <- c("Regen_Plot","Species","tallest_ht_cm","age","tallest_lastyr_cm")
# #only one height per species per age
# sapling.tallest.richter <- as.data.table(sapling.tallest.richter)
# sapling.tallest.richter <- sapling.tallest.richter[,list(tallest_ht_cm = max(tallest_ht_cm),
#                                                          tallest_lastyr_cm = max(tallest_lastyr_cm)),
#                                                    by=list(Regen_Plot,Species,age)]
# sapling.data.richter <- merge(sapling.data.richter,sapling.tallest.richter,by=c("Regen_Plot","Species","age"),all.X=TRUE)
# 
# resprout.richter <- resprout.richter[,c(2,24:29)]
# resprout.richter <- resprout.richter[resprout.richter$Resprout_Species != "" & !is.na(resprout.richter$Resprout_Species),]
# resprout.richter$Plot_ID <- paste0("richter.",resprout.richter$Plot_ID)
# names(resprout.richter) <- c("Regen_Plot","Species","._sprouts","DBH..cm.","tallest_ht_cm","tallest_lastyr_cm","talst_age")
# 
# 
# ##merge in all Richter tables
# plot.comb <- rbind.fill(plot.comb,plot.richter)
# resprout.comb <- rbind.fill(resprout.comb,resprout.richter)
# sap.comb <- rbind.fill(sap.comb,sapling.data.richter)
# seedl.comb <- rbind.fill(seedl.comb,seedling.richter)
# 





# Different methods for computing seed tree; convert into the value that can be obtained from all methods (closest tree regardless of identity)
# when seed tree distance is 0, >200, or 200, or NA (all indicating not found), make it 999
plot.comb$seed_tree_distance_conifer[plot.comb$seed_tree_distance_conifer %in% c("0",">200","200")] <- 999
plot.comb$seed_tree_distance_hardwood[plot.comb$seed_tree_distance_hardwood %in% c("0",">200","200")] <- 999
plot.comb$seed_tree_distance_general[plot.comb$seed_tree_distance_general %in% c("0",">200","200")] <- 999
# plot.comb$seed_tree_distance_conifer[is.na(plot.comb$seed_tree_distance_conifer)] <- 999
# plot.comb$seed_tree_distance_hardwood[is.na(plot.comb$seed_tree_distance_hardwood)] <- 999
# plot.comb$seed_tree_distance_general[is.na(plot.comb$seed_tree_distance_general)] <- 999


#! in fugure, consider computing seed tree distance separately for hardwood and conifer, just note that some of the early 2016 plots don't have separate values.

# make it numeric
plot.comb$seed_tree_distance_conifer <- as.numeric(plot.comb$seed_tree_distance_conifer)
plot.comb$seed_tree_distance_hardwood <- as.numeric(plot.comb$seed_tree_distance_hardwood)
plot.comb$seed_tree_distance_general <- as.numeric(plot.comb$seed_tree_distance_general)

#compute general seed tree distance as the minimum of hardwood and conifer
plot.comb$seed_tree_distance_general <- pmin(plot.comb$seed_tree_distance_conifer, plot.comb$seed_tree_distance_hardwood, plot.comb$seed_tree_distance_general, na.rm=TRUE)

# set survey year
plot.comb$Date <- "8/10/2016 0:00:00" # it's all 2016; the day doesn't matter
# except the new plots surveyed in 2017
surveyed.2017 <- c("CHI2220","CHI2428","CHI2216","CHI2475")
plot.comb[plot.comb$Regen_Plot %in% surveyed.2017,"Date"] <- "6/21/2017 0:00:00"

##!! NEW: set plot size


# # Find when only one quadrant was surveyed for seedlings, and multiply counts by 4 when true
# seedl.count.columns <- grep("yr$",names(seedl.comb))
# all.quadrants <- (seedl.comb$quadrants %in% c("ALL","4",""," ")) | is.na(seedl.comb$quadrants) # all quadrants were surveyed (T/F)
# seedl.comb[!all.quadrants,seedl.count.columns] <- 4*seedl.comb[!all.quadrants,seedl.count.columns]


# compute Fire column and Year.of.Fire based on first three plot letters
plot.comb$fire.prefix <- sapply(plot.comb$Regen_Plot,substr,start=1,stop=3)

fire.years <- data.frame(
  fire.prefix = c("BAG","PEA","PIT","STR","CHI","CUB","POW","ric"), #ric is from richter data (POWER only)
  Year.of.Fire = c(2012,2012,2008,2004,2012,2008,2004,2004),
  Fire = c("BAGLEY","PEAK","BTU LIGHTENING","STRAYLOR","CHIPS","CUB","POWER","POWER")
  )

#remove the fire name and year columns
plot.comb <- plot.comb[,! (names(plot.comb) %in% c("Year.of.Fire","Fire"))]
#merge in the fire name and year
plot.comb <- merge(plot.comb,fire.years,all.x=TRUE)


# swapped coords for CUB0409
plot.comb[which(plot.comb$Regen_Plot == "CUB0409"),c("Longitude","Latitude")] <- plot.comb[which(plot.comb$Regen_Plot == "CUB0409"),c("Latitude","Longitude")]




#there are a few plots with no coordinates written by field crew--add them
plot.comb[plot.comb$Regen_Plot == "CUB0294",c("Longitude","Latitude")] <- c(121.45444,40.21851)
plot.comb[plot.comb$Regen_Plot == "CUB0407",c("Longitude","Latitude")] <- c(121.43747,40.17117)


## Change lat/long to easting/northing (for all but power)
plot.comb.lat <- plot.comb[which(plot.comb$Latitude < 1000),]
plot.comb.utm <- plot.comb[which(plot.comb$Latitude >= 1000),]

geo <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
utm <- CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

plot.comb.lat$Longitude <- -plot.comb.lat$Longitude

plot.comb.lat.geo <- SpatialPoints(plot.comb.lat[,c("Longitude","Latitude")], proj4string=geo)
plot.comb.lat.utm <- spTransform(plot.comb.lat.geo,utm)
plot.comb.lat[,c("Easting","Northing")] <- coordinates(plot.comb.lat.utm)

plot.comb.utm[,c("Northing","Easting")] <- plot.comb.utm[,c("Latitude","Longitude")]

plot.comb <- rbind.fill(plot.comb.lat,plot.comb.utm)

#remove lat/long
plot.comb <- plot.comb[,!(names(plot.comb) %in% c("Latitude","Longitude"))]


### We now have all 2016 plots in plot.comb
plot.comb.notes <- plot.comb[,c("Regen_Plot","NOTES")]


### Change undecertain or incorrect species names based on photos taken by crew

plot.comb[which(toupper(plot.comb$dominant_shrub_1) == "SYMO?"),"dominant_shrub_1"] <- "SYMO"

plot.comb[which(toupper(plot.comb$dominant_shrub_1) == "MANZANITA?"),"dominant_shrub_1"] <- "LIDE3"
plot.comb[which(toupper(plot.comb$dominant_shrub_2) == "MANZANITA?"),"dominant_shrub_2"] <- "LIDE3"
plot.comb[which(toupper(plot.comb$dominant_shrub_3) == "MANZANITA?"),"dominant_shrub_3"] <- "LIDE3"

plot.comb[which(toupper(plot.comb$dominant_shrub_1) == "SYRO?"),"dominant_shrub_1"] <- "SYRO"
plot.comb[which(toupper(plot.comb$dominant_shrub_2) == "SYRO?"),"dominant_shrub_2"] <- "SYRO"
plot.comb[which(toupper(plot.comb$dominant_shrub_3) == "SYRO?"),"dominant_shrub_3"] <- "SYRO"

plot.comb[which(toupper(plot.comb$dominant_shrub_1) == "PREM?"),"dominant_shrub_1"] <- "SALIX"
plot.comb[which(toupper(plot.comb$dominant_shrub_2) == "PREM?"),"dominant_shrub_2"] <- "SALIX"
plot.comb[which(toupper(plot.comb$dominant_shrub_3) == "PREM?"),"dominant_shrub_3"] <- "SALIX"

should.be.ceve <- c("CEIN","CEIN3","CENI")
plot.comb[which(toupper(plot.comb$dominant_shrub_1) %in% should.be.ceve & plot.comb$Fire == "CUB"),"dominant_shrub_1"] <- "CEVE"
plot.comb[which(toupper(plot.comb$dominant_shrub_2) %in% should.be.ceve & plot.comb$Fire == "CUB"),"dominant_shrub_2"] <- "CEVE"
plot.comb[which(toupper(plot.comb$dominant_shrub_3) %in% should.be.ceve & plot.comb$Fire == "CUB"),"dominant_shrub_3"] <- "CEVE"

sap.comb[which(sap.comb$Species == "COSE?"),"Species"] <- "CONU"
resprout.comb[which(resprout.comb$Species == "COSE?"),"Species"] <- "CONU"
seedl.comb[which(seedl.comb$Species == "COSE?"),"Species"] <- "CONU"

plot.comb[which(plot.comb$Regen_Plot == "CUB0404"),"dominant_shrub_1"] <- "ARNE"

resprout.comb[which(resprout.comb$Species=="QUKE" & resprout.comb$Regen_Plot == "BAG1360"),"Species"] <- "QUCH2"

plot.comb[which(toupper(plot.comb$dominant_shrub_1) == "SCRUB OAK?"),"dominant_shrub_1"] <- "LIDE3"
plot.comb[which(toupper(plot.comb$dominant_shrub_2) == "SCRUB OAK?"),"dominant_shrub_2"] <- "LIDE3"
plot.comb[which(toupper(plot.comb$dominant_shrub_3) == "sCRUB OAK?"),"dominant_shrub_3"] <- "LIDE3"

plot.comb[which(toupper(plot.comb$dominant_shrub_1) == "BEAQ?"),"dominant_shrub_1"] <- "QUCH2"
plot.comb[which(toupper(plot.comb$dominant_shrub_2) == "BEAQ?"),"dominant_shrub_2"] <- "QUCH2"
plot.comb[which(toupper(plot.comb$dominant_shrub_3) == "BEAQ?"),"dominant_shrub_3"] <- "QUCH2"

plot.comb[which(plot.comb$Regen_Plot == "PIT0207"),"dominant_shrub_1"] <- "ARME"

plot.comb[which(toupper(plot.comb$dominant_shrub_1) == "ARVI?"),"dominant_shrub_1"] <- "QUVA"
plot.comb[which(toupper(plot.comb$dominant_shrub_2) == "ARVI?"),"dominant_shrub_2"] <- "QUVA"
plot.comb[which(toupper(plot.comb$dominant_shrub_3) == "ARVI?"),"dominant_shrub_3"] <- "QUVA"

plot.comb[which(plot.comb$Regen_Plot == "CUB0404"),"dominant_shrub_1"] <- "ARNE"

# LIDE3 as "shrub" in BAG1719,1683 and QUCH2 in BAG1537 is irrelevant even though LIDE3,QUCH2 is not a shrub because shrub cover was low to ground so it would not have been measured as a tree.




plot.welch <- read.csv("data_survey/Welch/Plot_data.txt",stringsAsFactors=FALSE)
sap.welch <- read.csv("data_survey/Welch/sapling_regen.txt",stringsAsFactors=FALSE)
seedl.welch <- read.csv("data_survey/Welch/tree_regen.txt",stringsAsFactors=FALSE)
resprout.welch <- read.csv("data_survey/Welch/Resprouts.txt",stringsAsFactors=FALSE)
surviving.trees.welch <- read.csv("data_survey/Welch/surviving_trees.txt",stringsAsFactors=FALSE)
seed.tree.welch <- read.csv("data_survey/Welch/seed_tree.txt",stringsAsFactors=FALSE)
shrub.welch <- read.csv("data_survey/Welch/shrub_regen.txt",stringsAsFactors=FALSE)

## From Welch plots, get modal height of the dominant shrub species and store it in the plot table

shrub.plots <- unique(shrub.welch$Regen_Plot)
for(shrub.plot in shrub.plots) {
  
  plot.shrubs <- shrub.welch[shrub.welch$Regen_Plot == shrub.plot,]
  max.cov <- max(plot.shrubs$Cover)
  max.cov.row <- row.names(plot.shrubs)[which(plot.shrubs$Cover == max.cov)[1]]
  dom.ht <- plot.shrubs[max.cov.row,]$modal_ht_cm
  
  plot.welch[plot.welch$Regen_Plot == shrub.plot,"dominant_shrub_ht_cm"] <- dom.ht
  
}




## Make seed tree column for welch
# take the nearest seed tree per plot
seed.tree.welch.nearest <- aggregate(seed.tree.welch$Dist_m,by=list(seed.tree.welch$Regen_Plot),FUN=min, na.rm=TRUE)
names(seed.tree.welch.nearest) <- c("Regen_Plot","seed_tree_distance_general")
seed.tree.welch.nearest[seed.tree.welch.nearest == Inf] <- NA
# merge into the welch plot table
plot.welch <- merge(plot.welch,seed.tree.welch.nearest,all.x=TRUE)




## Merge Welch plots with 2016 plots

plot <- rbind.fill(plot.welch,plot.comb)
sapling <- rbind.fill(sap.welch,sap.comb)
seedl <- rbind.fill(seedl.welch,seedl.comb)
resprout <- rbind.fill(resprout.welch,resprout.comb)
surviving.trees <- rbind.fill(surviving.trees.welch,surviving.trees.comb)
seed.tree <- seed.tree.welch


# Fix species names to use USDA PLANTS codes. Also fix an incorrectly-named species
from <- c("JU","AB","CADE","CONU","LIDE","QUCH","NODE","COSE","PRVI?","ACGI","CADE ","PIP0")
to <- c("JUNIPERUS","ABIES","CADE27","CONU4","LIDE3","QUCH2","LIDE3","CONU","PRVI","ACMA","CADE27","PIPO")


surviving.trees$Species <- mapvalues(surviving.trees$Species,from=from,to=to)
sapling$Species <- mapvalues(sapling$Species,from=from,to=to)
seedl$Species <- mapvalues(seedl$Species,from=from,to=to)
resprout$Species <- mapvalues(resprout$Species,from=from,to=to)
seed.tree$Species <- mapvalues(seed.tree$Species,from=from,to=to)


## there is a surviving tree called "CAPIM?" in PIT0164. No notes. Don't know what it is. Must drop plot.
plot <- plot[plot$Regen_Plot != "PIT0164",]


## add two resprouts that were incorrectly considered "shrubs" by the crew
resprout.add <- data.frame(ID=c(10001,10002),Regen_Plot=c("BAG1382","PIT0207"),COUNT.TOTAL=c(1000,1000),Species=c("LIDE3","ARME"),X5yr=c(1000,1000),._sprouts=c(1000,1000),DBH..cm.=c(1000,1000),talst_age=c(5,5),tallest_ht_cm=c(137,NA),tallest_lastyr_cm=c(NA,NA))
resprout <- rbind.fill(resprout,resprout.add)



# Calculate live BA
plot$BA.Live1 <- plot$BA_live_count * plot$BAF


#populate seedling unk_yr column with values for all hardwoods (but NOT CADE)
seedl.count.columns <- grep("yr$",names(seedl))
hardwoods<- c("QUKE","QUCH2","ARME","LIDE3","CHCH","QUGA4","ACMA","CEMO2","CONU4","POTR5","QUBE5","QUJO3","QUWI","UMCA")
#not.ageable <- c("CADE27",hardwoods)
not.ageable <- c(hardwoods)
seedl$Count_total <- rowSums(seedl[,seedl.count.columns],na.rm=TRUE) #this includes unk_yr (incase some of the seedlings were considered ageable and others not)
seedl[seedl$Species %in% not.ageable,"unk_yr"] <- seedl$Count_total[seedl$Species %in% not.ageable]
seedl.known.age.count.cols <- seedl.count.columns[1:(length(seedl.count.columns)-1)]

seedl.count.columns.names <- grep("yr$",names(seedl),value=TRUE)
seedl.count.columns.names <- seedl.count.columns.names[!(seedl.count.columns.names %in% "unk_yr")]

#for non-ageable seedlings, where ages were just put into unk_yr, set their known age columns to 0
seedl[seedl$Species %in% not.ageable,seedl.count.columns.names] <- 0





# save these data frames
write.csv(plot,"data_survey/Compiled/Plot_data_Davis.csv",row.names=FALSE)
write.csv(sapling,"data_survey/Compiled/sapling_regen_Davis.csv",row.names=FALSE)
write.csv(seedl,"data_survey/Compiled/tree_regen_Davis.csv",row.names=FALSE)
write.csv(resprout,"data_survey/Compiled/Resprouts_Davis.csv",row.names=FALSE)
write.csv(surviving.trees,"data_survey/Compiled/surviving_trees_Davis.csv",row.names=FALSE)
write.csv(seed.tree,"data_survey/Compiled/seed_tree_Davis.csv",row.names=FALSE)




#### 0.5 Operations that apply to all steps below #####



# Load plot locs
plot <- read.csv("data_survey/Compiled/Plot_data_Davis.csv",stringsAsFactors=FALSE)

plot <- plot[plot$Fire != "ZACA",] # remove Zaca fire

plot <- plot[!is.na(plot$Northing),]
plot <- plot[(plot$Easting < 1105000) & (plot$Easting > 150000) & (plot$Northing > 2500000) & (plot$Northing < 6000000),] #Easting 110 used to be 83
plot.coords <- data.frame(x=plot$Easting,y=plot$Northing)
utm10 <- CRS("+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
plots <- SpatialPointsDataFrame(coords=data.frame(x=plot$Easting,y=plot$Northing),data=plot,proj4string=utm10)



##!! exclude this plot because the location is wrong: SHR0900015


#write plots to shapefile
#writeOGR(plots, getwd(),"plots_shapefile",driver="ESRI Shapefile",overwrite=TRUE)

# 
# #### 1. Extract monthly climate data for each plot ####
# 
# CPUs <- 3
# years <- 1985:2015
# 
# topowx.dir <- "~/UC Davis/GIS/CA abiotic layers/TopoWx/"
# tmin.dir <- "tmin_annual/"
# tmax.dir <- "tmax_annual/"
# tmax.normal.dir <- "tmax_normal/"
# tmin.normal.dir <- "tmin_normal/"
# 
# 
# prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"
# ppt.dir <- "ppt_annual/"
# ppt.normal.dir <- "ppt_normal/"
# 
# 
# mo.num <- 1:12
# mo.chr <- sprintf("%02d",mo.num)
# 
# ### Open climate layers and store each month as a list element ###
# 
# ### TopoWx tmin and tmax monthly-by-yearly ###
# 
# tmax.yr.mo <- list()
# tmin.yr.mo <- list()
# for(year in years) {
#   
#   tmax.file <- paste(topowx.dir,tmax.dir,"tmax_",year,".nc",sep="")
#   tmax.year <- brick(tmax.file)
#   
#   tmin.file <- paste(topowx.dir,tmin.dir,"tmin_",year,".nc",sep="")
#   tmin.year <- brick(tmin.file)
#   
#   for(month in 1:12) {
#     
#     tmax.layer.index <- paste("tmax",year,mo.chr[month],sep=".")
#     tmin.layer.index <- paste("tmin",year,mo.chr[month],sep=".")
#     
#     tmax.yr.mo[tmax.layer.index] <- tmax.year[[month]]
#     tmin.yr.mo[tmin.layer.index] <- tmin.year[[month]]
#   }
# }
# 
# 
# ### PRISM precip monthly-by-yearly ###
# 
# ppt.yr.mo <- list()
# 
# for(year in years) {
#   
#   ppt.yr.files <- paste(prism.dir,ppt.dir,"PRISM_ppt_stable_4kmM3_",as.character(year),mo.chr,"_bil.bil",sep="")
#   ppt.yr <- stack(ppt.yr.files)
#   
#   for(month in 1:12) {
#     
#     layer.index <- paste("ppt",year,mo.chr[month],sep=".")
#     ppt.yr.mo[layer.index] <- ppt.yr[[month]]
#   }
# }
# 
# 
# ### MODIS snow cover monthly-by-yearly ###
# 
# snow.yr.mo <- list()
# 
# snow.brick <- brick("C:/Users/DYoung/Documents/UC Davis/GIS/Snowpack/n5eil01u.ecs.nsidc.org/MOST/snow_monthly_regen/snow_monthly_regen.grd")
# 
# for(year in years) {
#   for(month in 1:12) {
#     layer.index <- paste("snow",year,mo.chr[month],sep=".")
#     snow.yr.mo[[layer.index]] <- snow.brick[[layer.index]]
#   }
# }
# 
# ## get monthly normal out of this
# 
# snow.normal.mo <- list()
# for(month in mo.chr) {
#   
#   mo.search <- paste0(month,"$")
#   month.indices <- grep(mo.search,names(snow.yr.mo))
#   month.layers <- snow.yr.mo[month.indices]
#   layer.name <- paste0("snow.normal.",month)
#   snow.normal.mo[[layer.name]] <- mean(brick(month.layers),na.rm=TRUE)
# 
# }
# 
# 
# 
# 
# 
# ### TopoWx tmin and tmax monthly normal ###
# 
# tmax.normal.mo <- list()
# tmin.normal.mo <- list()
# 
# tmax.file <- paste(topowx.dir,tmax.normal.dir,"normals_tmax.nc",sep="")
# tmax.normal <- brick(tmax.file,varname="tmax_normal")
# 
# tmin.file <- paste(topowx.dir,tmin.normal.dir,"normals_tmin.nc",sep="")
# tmin.normal <- brick(tmin.file,varname="tmin_normal")
# 
# for(month in 1:12) {
#   
#   tmax.layer.index <- paste("tmax.normal",mo.chr[month],sep=".")
#   tmin.layer.index <- paste("tmin.normal",mo.chr[month],sep=".")
#   
#   tmax.normal.mo[tmax.layer.index] <- tmax.normal[[month]]
#   tmin.normal.mo[tmin.layer.index] <- tmin.normal[[month]]
# }
# 
# 
# ### PRISM precip monthly normals ###
# 
# prism.dir <- "~/UC Davis/GIS/CA abiotic layers/PRISM/"
# 
# ppt.normal.mo <- list()
# 
# ppt.files <- paste(prism.dir,ppt.normal.dir,"PRISM_ppt_30yr_normal_800mM2_",mo.chr,"_bil.bil",sep="")
# ppt.normal <- stack(ppt.files)
# 
# for(month in 1:12) {
#   layer.index <- paste("ppt.normal",mo.chr[month],sep=".")
#   ppt.normal.mo[layer.index] <- ppt.normal[[month]]
# }
# 
# 
# ### stack everything of the same resolution and extent into one stack
# 
# temp.800m.rast <- unlist(list(tmin.yr.mo,tmax.yr.mo,tmin.normal.mo,tmax.normal.mo))
# ppt.4km.rast <- unlist(list(ppt.yr.mo))
# ppt.800m.rast <- unlist(list(ppt.normal.mo))
# snow.5km.rast <- unlist(list(snow.yr.mo,snow.normal.mo))
# 
# 
# 
# ### Extract climate values for points ###
# ppt.4km <- sapply(X = ppt.4km.rast,FUN = extract.single, plots = plots) # takes ~1 min
# ppt.800m <- sapply(X = ppt.800m.rast,FUN = extract.single, plots = plots)
# temp.800m <- sapply(X = temp.800m.rast,FUN = extract.single, plots = plots) # takes ~10 min
# snow.5km <- sapply(X = snow.5km.rast,FUN = extract.single, plots = plots) # takes ~10 sec
# 
# 
# plots.clim <- data.frame(Regen_Plot=plots$Regen_Plot,temp.800m,ppt.4km,ppt.800m,snow.5km)
# #plots.clim <- plots.clim[complete.cases(plots.clim),]
# 
# 
# write.csv(plots.clim,"../data_intermediate_processing_local/plot_climate_monthly.csv",row.names=FALSE)
# 
# 
# 
# #### 2. Summarize monthly climate by water year ####
# 
# plots.clim <- read.csv("../data_intermediate_processing_local/plot_climate_monthly.csv",header=TRUE)
# 
# first.year <- years[1]
# last.year <- years[length(years)]
# 
# nyears <- last.year - first.year + 1
# if(nyears < 3) {stop("You must supply at least 3 years.")}
# n.water.years <- nyears - 1
# inter.years <- years[2:(nyears-1)]
# inter.years.rep <- rep(inter.years,each=12)
# 
# first.year.mos <- c("10","11","12")
# last.year.mos <- c("01","02","03","04","05","06","07","08","09")
# inter.mos <- rep(mo.chr,nyears-2)
# 
# first.year.mo.yr <- paste(first.year,first.year.mos,sep=".")
# last.year.mo.yr <- paste(last.year,last.year.mos,sep=".")
# inter.mos.mo.yr <- paste(inter.years.rep,inter.mos,sep=".")
# 
# all.mo.yr <- c(first.year.mo.yr,inter.mos.mo.yr,last.year.mo.yr)
# all.mo.yr.tmin <- paste("tmin",all.mo.yr,sep=".")
# all.mo.yr.tmax <- paste("tmax",all.mo.yr,sep=".")
# all.mo.yr.ppt <- paste("ppt",all.mo.yr,sep=".")
# all.mo.yr.snow <- paste("snow",all.mo.yr,sep=".")
# 
# tmin.col.nums <- sapply(all.mo.yr.tmin,function(x) {grep(x,names(plots.clim))})
# tmax.col.nums <- sapply(all.mo.yr.tmax,function(x) {grep(x,names(plots.clim))})
# ppt.col.nums <- sapply(all.mo.yr.ppt,function(x) {grep(x,names(plots.clim))})
# snow.col.nums <- sapply(all.mo.yr.snow,function(x) {grep(x,names(plots.clim))})
# 
# tmin.col.nums.mat <- matrix(tmin.col.nums,ncol=12,byrow=TRUE)
# tmax.col.nums.mat <- matrix(tmax.col.nums,ncol=12,byrow=TRUE)
# ppt.col.nums.mat <- matrix(ppt.col.nums,ncol=12,byrow=TRUE)
# snow.col.nums.mat <- matrix(snow.col.nums,ncol=12,byrow=TRUE)
# 
# 
# water.year.summary.annual <- matrix(nrow=nrow(plots.clim),ncol=0)
# for (i in 1:n.water.years) {
#   
#   water.year.ending <- years[i]+1
#   
#   tmin.col.nums.year <- tmin.col.nums.mat[i,]
#   tmin.cols.year <- plots.clim[,tmin.col.nums.year]
#   
#   tmax.col.nums.year <- tmax.col.nums.mat[i,]
#   tmax.cols.year <- plots.clim[,tmax.col.nums.year]
#   
#   tmean.cols.year <- (tmin.cols.year + tmax.cols.year) / 2
#   
#   ppt.col.nums.year <- ppt.col.nums.mat[i,]
#   ppt.cols.year <- plots.clim[,ppt.col.nums.year]
# 
#   snow.col.nums.year <- snow.col.nums.mat[i,]
#   snow.cols.year <- plots.clim[,snow.col.nums.year]  
#   snow <- snow.cols.year
#   
#   # ### snow: for each plot, sum ppt from sept (col 1) to may (col 9) for months when avg temp was <=0
#   # snow <- ppt.cols.year
#   # snow[tmean.cols.year > 0] <- 0 # no snow if temp is > 0
#   # 
#   # ### rain: for each plot, sum ppt from jan (col 5) to aug (col 12) for months when avg temp was >0
#   # rain <- ppt.cols.year
#   # rain[tmean.cols.year <= 0] <- 0
# 
#   tmin.avg.year <- apply(tmin.cols.year,1,mean)
#   tmax.avg.year <- apply(tmax.cols.year,1,mean)
#   
#   ppt.tot.year <- apply(ppt.cols.year,1,sum)
#   
#   snow.tot.year <- apply(snow[c(4:11)],1,sum)
#   #rain.tot.year <- apply(rain[4:11],1,sum)
#   
#   tmean.avg.year <- (tmin.avg.year+tmax.avg.year)/2
#   
#   # calculate JJA average tmin, tmax, tmean
#   # calculate DJF averate tmin, tmax, tmean
#   tmean.avg.JJA <- apply(tmean.cols.year[,9:11],1,mean)
#   tmean.avg.DJF <- apply(tmean.cols.year[,3:5],1,mean)
# 
#   ### merge into DF
#   colname.prefixes <- c("tmin","tmax","tmean","tmean.JJA","tmean.DJF","ppt","snow") #,"rain")
#   colnames <- paste(colname.prefixes,water.year.ending,sep=".")
#   
#   clim.year <- cbind(tmin.avg.year,tmax.avg.year,tmean.avg.year,tmean.avg.JJA,tmean.avg.DJF,ppt.tot.year,snow.tot.year) #,rain.tot.year)
#   colnames(clim.year) <- colnames
#   
#   water.year.summary.annual <- cbind(water.year.summary.annual,clim.year)
# }
# 
# ### summarize normal values annually (in the case of normals, calendar and water year are identical) ###
# 
# tmin.col.nums <- grep("tmin.normal",names(plots.clim))
# tmax.col.nums <- grep("tmax.normal",names(plots.clim))
# ppt.col.nums <- grep("ppt.normal",names(plots.clim))
# snow.col.nums <- grep("snow.normal",names(plots.clim))
# 
# tmin.cols <- plots.clim[,tmin.col.nums]
# tmax.cols <- plots.clim[,tmax.col.nums]
# ppt.cols <- plots.clim[,ppt.col.nums]
# snow.cols <- plots.clim[,snow.col.nums]
# 
# tmean.cols <- (tmin.cols + tmax.cols)/2
# 
# ### snow: for each plot, sum snow from nov (col 3) to may (col 9)
# snow <- snow.cols
# 
# 
# # ### rain: for each plot, sum ppt from jan (col 5) to aug (col 12) for months when avg temp was >0
# # rain <- ppt.cols
# # rain[tmean.cols <= 0] <- 0
# 
# tmin.avg.normal.annual <- apply(tmin.cols,1,mean)
# tmax.avg.normal.annual <- apply(tmax.cols,1,mean)
# tmean.avg.normal.annual <- (tmin.avg.normal.annual + tmax.avg.normal.annual) / 2
# 
# ## tmean JJA and DJF
# tmean.avg.normal.DJF <- apply(tmean.cols[,c(12,1,2)],1,mean)
# tmean.avg.normal.JJA <- apply(tmean.cols[,6:8],1,mean)
# 
# ppt.tot.normal.annual <- apply(ppt.cols,1,sum)
# snow.tot.normal.annual <- apply(snow[,1:7],1,sum)
# # rain.tot.normal.annual <- apply(rain[,c(3,4,5,6,7,8,9)],1,sum)
# 
# water.year.summary.normal <- cbind(tmin.avg.normal.annual,tmax.avg.normal.annual,tmean.avg.normal.annual,tmean.avg.normal.JJA,tmean.avg.normal.DJF,ppt.tot.normal.annual,snow.tot.normal.annual) #,rain.tot.normal.annual)
# colnames(water.year.summary.normal) <- c("tmin.normal.ann","tmax.normal.ann","tmean.normal.ann","tmean.normal.JJA","tmean.normal.DJF","ppt.normal.ann","snow.normal.ann") #,"rain.normal.ann")
# water.year.summary <- data.frame(Regen_Plot=plots.clim$Regen_Plot,water.year.summary.annual,water.year.summary.normal)
# 
# 
# 
# 
# ### add water balance values ###
# 
# library(reshape)
# 
# plots.wb <- read.csv("../data_water_balance/plots_wb.csv",header=TRUE)
# water.years <- years[-1]
# 
# plots.wb <- plots.wb[plots.wb$year.ending %in% water.years,]
# 
# plots.wb.def <- cast(plots.wb,ID ~ year.ending,value="def")
# names(plots.wb.def)[-1] <- paste0("def.",names(plots.wb.def[-1]))
# 
# plots.wb.aet <- cast(plots.wb,ID ~ year.ending,value="aet")
# names(plots.wb.aet)[-1] <- paste0("aet.",names(plots.wb.aet[-1]))
# 
# 
# water.year.summary <- merge(water.year.summary,plots.wb.def,by.x="Regen_Plot",by.y="ID",all.x=TRUE)
# water.year.summary <- merge(water.year.summary,plots.wb.aet,by.x="Regen_Plot",by.y="ID",all.x=TRUE)
# 
# write.csv(water.year.summary,"../data_intermediate_processing_local/plot_climate_water_year.csv",row.names=FALSE)
# 
# 
# 

# 
# #### 3. Extract fire severity for each plot ####
# 
# fire.severity.rasters.dir <- "C:/Users/dyoung/Documents/UC Davis/GIS/MTBS fire severity/"
# 
# 
# fires <- unique(plots$Fire)
# 
# plot.firesev <- data.frame(Regen_Plot=NULL,fire.sev=NULL)
# counter <- 1
# for(fire in fires) {
#   
#   #disable this
#   if(counter > 2) {
#     cat("Fire severity analysis disabled.\n")
#     next()
#   }
#   
#   plots.fire <- plots[plots$Fire == fire,]
#   fire.raster.dir <- paste(fire.severity.rasters.dir,fire,"/",sep="")
#   
#   # get the name of the DNBR6 layer
#   file.list <- try(list.files(fire.raster.dir))
#   
#   
#   raster.file <- file.list[grep("dnbr6.?.tif$",file.list)]
#   
#   if(length(raster.file) != 1) {
#     cat("More than or less than one dnbr6 layer for Fire",fire,"\n")
#     counter <- counter + 1
#     next()
#   }
#   
#   sev.layer <- raster(paste(fire.raster.dir,raster.file,sep=""))
#   
#   firesev <- extract(sev.layer,plots.fire)
#   
#   sev.layer.lowsev.only <- sev.layer
#   sev.layer.lowsev.only[sev.layer.lowsev.only %in% c(4)] <- NA # exclude high only so we compute "distance to low or moderate"
#   sev.layer.lowsev.only[!is.na(sev.layer.lowsev.only)] <- 0
#   dist.to.low.rast <- distance(sev.layer.lowsev.only) #! make sure distance units are correct
#   
#   plots.fire <- spTransform(plots.fire,CRS(projection(sev.layer)))
#   
#   
#   dist.to.low <- extract(dist.to.low.rast,plots.fire,method="bilinear")
#   regen.plot <- plots.fire$Regen_Plot
#   
#   fire.vals <- data.frame(Regen_Plot=regen.plot,firesev=firesev,dist.to.low=dist.to.low) #! note that value of 5 here means increased greenness
#   plot.firesev <- rbind(plot.firesev,fire.vals)
#   
#   cat("Finished severity analysis of fire ",fire," (",counter," of ",length(fires),")\n",sep="")
#   counter <- counter + 1
# }
# 
# write.csv(plot.firesev,"../data_intermediate_processing_local/plot_fire_data.csv",row.names=FALSE)
# 



#### 4. Summarize regen and surviving trees by species (and age for regen) for each plot ####

sap <- read.csv("data_survey/Compiled/sapling_regen_Davis.csv")
#shrub <- read.csv("../data_survey/shrub_regen.csv")
seedl <- read.csv("data_survey/Compiled/tree_regen_Davis.csv")
resprout <- read.csv("data_survey/Compiled/Resprouts_Davis.csv")
surviving.trees <- read.csv("data_survey/Compiled/surviving_trees_Davis.csv",header=TRUE,stringsAsFactors=FALSE)

sap.plt <- merge(sap,plot[,c("Regen_Plot","FIRE_SEV")],by="Regen_Plot")


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


## make a column for subsampled
seedl$subsampled = ifelse(toupper(seedl$quadrants) %in% c("ALL","4","") | is.na(seedl$quadrants), FALSE, TRUE)

#get rid of the extra 11yr column
seedl$X11.yr <- seedl$X11.yr + seedl$X11yr
seedl <- seedl[,!(names(seedl) %in% "X11yr")]

#put unk_yr at the end
unkyr <- seedl$unk_yr
seedl <- seedl[,!(names(seedl) %in% "unk_yr")]
seedl$unk_yr <- unkyr

seedl.count.columns <- grep("yr$",names(seedl))
seedl.ag <- aggregate(seedl[,seedl.count.columns],by=list(species=seedl$Species,Regen_Plot=seedl$Regen_Plot),FUN=sum,na.rm=TRUE)
count.yrs <- c(paste("count.",0:(length(seedl.count.columns)-2),"yr",sep=""),"unk_yr")
names(seedl.ag) <- c("species","Regen_Plot",count.yrs)





### aggregate sapling table by plot and species
sap$tot <- rowSums(sap[,6:14],na.rm=TRUE)

#if the species had a row for the sapling but not under an age column, see if it had an age assigned in the age columns. if so, put "1" under the appropriate age column, otherwise assume it was 10 yr (this was only hardwoods and all hardwoods are going to get re-aged anyway because it's impossible to know)
saps.no.age <- which(sap$tot == 0 | (is.na(sap$tot)))

for(i in saps.no.age) {
  
  sap.row <- sap[i,]
  age <- max(sap.row$tallest_age,sap.row$age,na.rm=TRUE)
  
  if(age != -Inf) { #there was an age specified
    
    search <- paste0("X",age,"\\.?yr$")
    age.col <- grep(search,names(sap.row))
    sap[i,age.col] <- 1 # set the age (if no matching age col, does not set the age)
    sap[i,"tot"] <- 1

  } else { # there was no age specified, so set age 10 to "1". In highsev plots, this only applies to hardwoods, which are all going to get reassigned to the age of the fire anyway
    
    sap[i,"X10yr"] <- 1
    
  }
}

sap.ag <- aggregate(sap[,c(6:14)],by=list(species=sap$Species,Regen_Plot=sap$Regen_Plot),FUN=sum,na.rm=TRUE)
count.yrs <- paste("count.",3:11,"yr",sep="")
names(sap.ag) <- c("species","Regen_Plot",count.yrs)


##For hardwood saps, set age to unk_yr

#populate seedling unk_yr column with values for all hardwoods (but NOT CADE)
sap.count.columns <- grep("yr$",names(sap.ag))
hardwoods<- c("QUKE","QUCH2","ARME","LIDE3","CHCH","QUGA4","ACMA","CEMO2","CONU4","POTR5","QUBE5","QUJO3","QUWI","UMCA")
#not.ageable <- c("CADE27",hardwoods)
not.ageable <- c(hardwoods)
sap.ag$tot <- rowSums(sap.ag[,sap.count.columns],na.rm=TRUE) #this includes unk_yr (incase some of the seedlings were considered ageable and others not)
sap.ag[sap.ag$species %in% not.ageable,"unk_yr"] <- sap.ag$tot[sap.ag$species %in% not.ageable]

#for non-ageable seedlings, where ages were just put into unk_yr, set their known age columns to 0
sap.ag[sap.ag$species %in% not.ageable,sap.count.columns] <- 0


### aggregate resprout table by plot and species
resprout$COUNT.TOTAL <- ifelse(is.na(resprout$COUNT.TOTAL),1,resprout$COUNT.TOTAL) # sometimes it has a count, sometimes it has one record per individual (stem?)
resprout.ag <- aggregate(resprout$COUNT.TOTAL,by=list(species=resprout$Species,Regen_Plot=resprout$Regen_Plot),FUN=sum,na.rm=TRUE)
count.yrs <- paste("count.",4,"yr",sep="") # use this because it will always be counted as "old" and "all" and not as "young" which is what we need.
names(resprout.ag) <- c("species","Regen_Plot",count.yrs)

### merge seedling and sapling tables and resprout table
regen <- rbind.fill(seedl.ag,sap.ag)
regen <- rbind.fill(regen,resprout.ag) # add in resprout table (comment out here if desired)
regen.ag <- aggregate(regen[,3:ncol(regen)],by=list(species=regen$species,Regen_Plot=regen$Regen_Plot),FUN=sum,na.rm=TRUE)
# for hardwoods, F is for seedlings/saplings and 5yr is for resprouts

## Pull in subsample status
library(dplyr)
regen.ag = left_join(regen.ag,seedl %>% select(Regen_Plot, species = Species,subsampled))


### aggregate surviving tree table by plot and species
surviving.trees <- surviving.trees[surviving.trees$DBH > 7.5,] # exclude small trees (it seems sometimes crews put smaller trees as saplings instead of surviving trees)
surviving.trees$ba <- (surviving.trees$DBH/2)^2*3.14
surviving.trees$count <- 1
surviving.trees.ag <- aggregate(surviving.trees[,c("count","ba")],by=list(species=surviving.trees$Species,Regen_Plot=surviving.trees$Regen_Plot),FUN=sum)
names(surviving.trees.ag) <- c("species","Regen_Plot","surviving.trees.count","surviving.trees.ba")

### merge regen and surviving tree DFs
regen.ag <- merge(regen.ag,surviving.trees.ag,by=c("Regen_Plot","species"),all=TRUE)


regen.ag$surviving.trees.count[is.na(regen.ag$surviving.trees.count)] <- 0
regen.ag$surviving.trees.ba[is.na(regen.ag$surviving.trees.ba)] <- 0

write.csv(regen.ag,"data_intermediate_processing_local/tree_summarized_sp_Davis.csv",row.names=FALSE)



#### 5. Integrate all plot and species data ####

### Read in raw data files (direct from DB export)
fire.years <- read.csv("data_fire/fire_years.csv",header=TRUE,stringsAsFactors=FALSE)
seed.tree <- read.csv("data_survey/Compiled/seed_tree_Davis.csv",header=TRUE,stringsAsFactors=FALSE)

### Read in the summarized data files
#plot.climate <- read.csv("../data_intermediate_processing_local/plot_climate_water_year.csv",stringsAsFactors=FALSE)
#plot.fire.data <- read.csv("data_intermediate_processing_local/plot_fire_data.csv",stringsAsFactors=FALSE)
plot.tree.sp <- read.csv("data_intermediate_processing_local/tree_summarized_sp_Davis.csv",stringsAsFactors=FALSE)

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

# 
# ## extract elevation for each plot
# dem <- raster("C:/Users/DYoung/Documents/UC Davis/GIS/CA abiotic layers/DEM/CA/new ncal/CAmerged12_albers2_firesmask_ncal.tif")
# plots$elev.m <- extract(dem,plots,method="bilinear")
# 
# ## extract march solar rad for each plot
# rad.march <- raster("C:/Users/DYoung/Documents/UC Davis/GIS/CA abiotic layers/solar rad/regen_srad_clear_sky/srad_clear_mar_regen_ncal.tif")
# plots$rad.march <- extract(rad.march,plots,method="bilinear")
# 
# #compile DF of geospatial data that was extracted for each plot
# plots.extracted <- plots[,c("Regen_Plot","elev.m","rad.march")]
# 

## merge the plot-level data into a single DF
plot.1 <- plot #merge(plot,plot.fire.data,by="Regen_Plot",all.x=TRUE)
plot.2 <- merge(plot.1,fire.years,by="Fire",all.x=TRUE)
plot.2$survey.years.post <- plot.2$Year - plot.2$fire.year

plot.3 <- plot.2 #merge(plot.2,plots.extracted,all.x=TRUE)

# ### get summarized climate data for each plot
# ##!! If want to look at weather beyond 4 years post-fire (for those fires that had more than 4 years), will need to make this relative to number of years post-fire
# plot.3.clim <- summarize.clim(plot.3,plot.climate,years.clim=1:3) #first three years after fire
# plot.3.clim2 <- summarize.clim(plot.3,plot.climate,years.clim=1:2) # for more than 3 years out, will need 2016 weather data
# names(plot.3.clim2) <- c(names(plot.3.clim2[1]),paste(names(plot.3.clim2)[2:ncol(plot.3.clim2)],".early",sep=""))
# plot.3.clim.both <- merge(plot.3.clim,plot.3.clim2)


### Find plots with planted trees and exclude them
seedl <- read.csv("data_survey/Compiled/tree_regen_Davis.csv")
seedl.p <- seedl[seedl$seed_veg_plant == "P",]
plots.planted <- unique(seedl.p$Regen_Plot)
plot.3 <- plot.3[!(plot.3$Regen_Plot %in% plots.planted),]


### get summarized regen data for each plot (also summarizes adults)
##!! NOTE that as written here, the code includes un-ageable species when tallying regen for ALL ages, but not for specific age classes
#plot.3.regen.old <- summarize.regen.ind(plot.3,plot.tree.sp,sp=c("ABCO","PSME","PIPO"),regen.ages="old",all.sp=TRUE,peryear=TRUE,incl.unk.age=TRUE)
#plot.3.regen.young <- summarize.regen.ind(plot.3,plot.tree.sp,sp=c("ABCO","PSME","PIPO"),regen.ages="young",all.sp=TRUE,peryear=TRUE)[,1:3] #only take the regen data (because funct also outupts adults data but we get that from the first call, the previous line)
plot.3.regen.all <- summarize.regen.ind(plot.3,plot.tree.sp,sp=c("ABCO","PSME","PIPO"),regen.ages="all",all.sp=TRUE,incl.unk.age=TRUE,peryear=FALSE)
# names(plot.3.regen.old)[3] <- "regen.count.old"
# names(plot.3.regen.young)[3] <- "regen.count.young"
# names(plot.3.regen.all)[3] <- "regen.count.all"

seedl$subsampled = ifelse(toupper(seedl$quadrants) %in% c("ALL","4","") | is.na(seedl$quadrants), FALSE, TRUE)
plot.3.regen.all = left_join(plot.3.regen.all,seedl %>% select(Regen_Plot, species = Species,subsampled))

plot.3.regen.all$subsampled = ifelse(is.na(plot.3.regen.all$subsampled),FALSE,TRUE)


plot.3.regen <- plot.3.regen.all


### species-level seed tree distance
## if a species within a plot has multiple seed trees listed, get the shortest distance
seed.tree.sp <- aggregate(seed.tree$Dist_m,by=list(seed.tree$Regen_Plot,seed.tree$Species),FUN=min)
names(seed.tree.sp) <- c("Regen_Plot","species","seed.tree.sp")
plot.3.regen <- merge(plot.3.regen,seed.tree.sp,all.x=TRUE)


### add the tallest seedling of each species
# 
# ## append seedling and resprout heights
# seedl.ht <- seedl[,c("Regen_Plot","Species","tallest_ht_cm")]
# respr.ht <- resprout[,c("Regen_Plot","Species","tallest_ht_cm")]
# regen.ht <- rbind(seedl.ht,respr.ht)
# 
# # in case there are multiple seedling entries per species-plot combination, aggregate
# regen.ht.agg <- aggregate(regen.ht$tallest_ht_cm,by=list(regen.ht$Regen_Plot,regen.ht$Species),FUN=max)
# names(regen.ht.agg) <- c("Regen_Plot","species","tallest_ht_cm")
# 
# ## Now do it for species groups
# pinus <- c("PIPO","PIJE","PILA","PIAT","PICO","PINUS","PISA2","PIMO3")
# regen.ht.pinus <- regen.ht[regen.ht$Species %in% pinus,]
# regen.ht.pinus.agg <- aggregate(regen.ht.pinus$tallest_ht_cm,by=list(regen.ht.pinus$Regen_Plot),FUN=max)
# names(regen.ht.pinus.agg) <- c("Regen_Plot","tallest_ht_cm")
# regen.ht.pinus.agg$species <- "PINUS.ALLSP"
# 
# shade <- c("ABCO","CADE27","ABIES","ABMA")
# shade <- c("ABCO","CADE27")
# regen.ht.shade <- regen.ht[regen.ht$Species %in% shade,]
# regen.ht.shade.agg <- aggregate(regen.ht.shade$tallest_ht_cm,by=list(regen.ht.shade$Regen_Plot),FUN=max)
# names(regen.ht.shade.agg) <- c("Regen_Plot","tallest_ht_cm")
# regen.ht.shade.agg$species <- "SHADE.ALLSP"
# 
# yellow <- c("PIPO","PIJE")
# regen.ht.yellow <- regen.ht[regen.ht$Species %in% yellow,]
# regen.ht.yellow.agg <- aggregate(regen.ht.yellow$tallest_ht_cm,by=list(regen.ht.yellow$Regen_Plot),FUN=max)
# names(regen.ht.yellow.agg) <- c("Regen_Plot","tallest_ht_cm")
# regen.ht.yellow.agg$species <- "PIPJ"
# 
# hdwd <- c("ALNUS","QUKE","QUCH2","ARME","LIDE3","CHCH","QUGA4","ACMA","ACME","CEMO2","CONU4","POTR5","QUBE5","QUJO3","QUWI","UMCA","ALRH")
# regen.ht.hdwd <- regen.ht[regen.ht$Species %in% hdwd,]
# regen.ht.hdwd.agg <- aggregate(regen.ht.hdwd$tallest_ht_cm,by=list(regen.ht.hdwd$Regen_Plot),FUN=max)
# names(regen.ht.hdwd.agg) <- c("Regen_Plot","tallest_ht_cm")
# regen.ht.hdwd.agg$species <- "HDWD.ALLSP"
# 
# 
# ## Append
# regen.ht.agg <- rbind.fill(regen.ht.agg,regen.ht.pinus.agg,regen.ht.shade.agg,regen.ht.yellow.agg,regen.ht.hdwd.agg)
# 
# 
# 
# # merge it in
# plot.3.regen <- merge(plot.3.regen,regen.ht.agg,all.x=TRUE)
# # species table ready for export


### add plot-level seed tree distance (shortest distance among all seed trees recorded for the plot)
seed.tree.any <- aggregate(seed.tree$Dist_m,by=list(seed.tree$Regen_Plot),FUN=min)
names(seed.tree.any) <- c("Regen_Plot","seed.tree.any")

### merge the plot, climate, and plot-level seed tree data
plot.clim <- plot.3
plot.clim.seedtree <- merge(plot.clim,seed.tree.any,all.x=TRUE)

## Correct misspelled FORBE
plot.clim.seedtree$FORB <- plot.clim.seedtree$FORBE
plot.clim.seedtree <- remove.vars(plot.clim.seedtree,"FORBE")


### Bring in dominant vegetation for Welch
dom.veg.welch <- read.csv("data_survey/Welch/DominantVegetation.txt",stringsAsFactors=FALSE)

dom.veg.df <- data.frame()

for(plt.name in unique(dom.veg.welch$Regen_Plot)) {
  
  dom.veg.plt <- dom.veg.welch[dom.veg.welch$Regen_Plot == plt.name,]
  dom.veg.plt.string <- paste(dom.veg.plt$Species,collapse=" ")
  dom.veg.plt.df <- data.frame(Regen_Plot = plt.name,dom.veg = dom.veg.plt.string)
  dom.veg.df <- rbind(dom.veg.df,dom.veg.plt.df)
  
}

plot.clim.seedtree.2 <- merge(plot.clim.seedtree,dom.veg.df,all.x=TRUE,by="Regen_Plot")


## Bring in dominant vegetation for other plots
plot.clim.seedtree.2$dom.veg.2016 <- paste(plot.clim.seedtree.2$dominant_tree_1,plot.clim.seedtree.2$dominant_tree_2,plot.clim.seedtree.2$dominant_tree_3,sep=" ")

plot.clim.seedtree.2$dom.veg.all <- paste(plot.clim.seedtree.2$dom.veg.2016,plot.clim.seedtree.2$dom.veg,sep=" ")






### write plot-level and species-level output files
write.csv(plot.clim.seedtree.2,"data_intermediate/plot_level_Davis.csv",row.names=FALSE)
write.csv(plot.3.regen,"data_intermediate/speciesXplot_level_Davis.csv",row.names=FALSE)




#### Outlier check ####

# 
# #test to see if CADE27 is really essentially only present in low sev plots
# 
# test <- merge(plot.tree.sp,plot.clim.seedtree[,c("Regen_Plot","FIRE_SEV")])
# 
# 
# 
# 
# 
# plot.3.regen <- read.csv("data_intermediate/speciesXplot_level.csv",header=TRUE,stringsAsFactors=FALSE)
# plot.clim.seedtree <- read.csv("data_intermediate/plot_level.csv",header=TRUE,stringsAsFactors=FALSE)
# 
# 
# clim.vars <- c("SHRUB","rad.march","ppt.normal","ppt.post","ppt.post.min","diff.norm.ppt.z","diff.norm.ppt.min.z", "seed_tree_distance_general","def.post","diff.norm.def.z","def.normal")
# 
# pairs(plot.clim.seedtree[,clim.vars])
# 
# plot.3.regen.conifer <- plot.3.regen[plot.3.regen$species=="CONIF.ALLSP",]
# 
# check.df <- merge(plot.clim.seedtree,plot.3.regen.conifer,by="Regen_Plot",all.x=TRUE)
# 
# all.vars <- c("SHRUB", "rad.march", "ppt.normal", "ppt.post", "ppt.post.min", "diff.norm.ppt.z", "diff.norm.ppt.min.z", "seed_tree_distance_general","regen.count.all","def.post","diff.norm.def.z","def.normal")
# pairs(check.df[,all.vars])
# 
# check.df.abbrev <- check.df[,c("Regen_Plot","seed_tree_distance_general","regen.count.all","SHRUB","diff.norm.ppt.z")]
# 
# check.df.abbrev[check.df.abbrev$seed_tree_distance_general>150 & check.df.abbrev$regen.count.all > 80,]
# check.df.abbrev[check.df.abbrev$regen.count.all > 250,]
# 
# pairs(check.df.abbrev[,-1])
# 
# check.df.abbrev[check.df.abbrev$seed_tree_distance_general > 300 & check.df.abbrev$regen.count.all > 5,]
