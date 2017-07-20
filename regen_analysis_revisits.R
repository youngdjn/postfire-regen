setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(party)
library(ggplot2)
library(brms)
library(pROC)
library(betareg)
library(car)
library(plyr)
library(data.table)
library(sf)

source("regen_analysis_functions.R")


#### 1. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

d.plot$Fire <- sapply(d.plot$Fire,FUN= function(x) simpleCap(tolower(x)))

d.plot[d.plot$Fire == "Btu Lightening","Fire"] <- "BTU Lightning"



## Look for plots with incomplete data specified in the comments
plots.exceptions <- grepl("#.*(INCOMPLETE|INCORRECT)",d.plot$NOTES)
d.plot <- d.plot[!plots.exceptions,]

# only keep the necessary columns
d.plot <- d.plot[,c("Regen_Plot","Fire","Year.of.Fire","Easting","Northing","aspect","slope","SHRUB","FORB","GRASS","HARDWOOD","CONIFER","FIRE_SEV","BA.Live1","Year","firesev","dist.to.low","fire.abbr","X5yr","fire.year","survey.years.post","elev.m","rad.march","tmean.post","ppt.post","ppt.post.min","tmean.normal","ppt.normal","seed.tree.any","diff.norm.ppt.z","diff.norm.ppt.min.z","seed_tree_distance_general","seed_tree_distance_conifer","seed_tree_distance_hardwood","diff.norm.ppt.z","diff.norm.ppt.min.z","tmean.post","ppt.post","ppt.post.min","perc.norm.ppt","perc.norm.ppt.min","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z","def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","diff.norm.def.max.z","diff.norm.aet.min.z","def.post","aet.post","def.post.max","aet.post.min","snow.post.min","snow.normal","snow.post","diff.norm.snow.z","diff.norm.snow.min.z","dominant_shrub_ht_cm","dom.veg.all")]

## Remove managed plots, plots in nonforested locations (e.g. exposed bedrock), etc.
plots.exclude <- read.csv("data_intermediate/plots_exclude.csv",header=T,stringsAsFactors=FALSE)
plot.ids.exclude <- plots.exclude[plots.exclude$Exclude != "",]$Regen_Plot
d.plot <- d.plot[!(d.plot$Regen_Plot %in% plot.ids.exclude),]

# # Remove any plots > 50m from seed source
# d.plot <- d.plot[which(d.plot$seed_tree_distance_general < 75),]

d.sp$regen.count.nonyoung <- d.sp$regen.count.all - d.sp$regen.count.young


# variable for regen presence/absence
d.sp$regen.presab.young <- ifelse(d.sp$regen.count.young > 0,TRUE,FALSE)
d.sp$regen.presab.old <- ifelse(d.sp$regen.count.old > 0,TRUE,FALSE)
d.sp$regen.presab.all <- ifelse(d.sp$regen.count.all > 0,TRUE,FALSE)
d.sp$regen.presab.nonyoung <- ifelse(d.sp$regen.count.nonyoung > 0,TRUE,FALSE)

# severity categories
high.sev <- c(3,4,5) # which field-assessed severity values are considered high severity
control <- c(0,1) # which field-assessed severity values are considered controls

# keep only high-sev
d.plot <- d.plot[d.plot$FIRE_SEV %in% high.sev,]

# list n plots in each survey year for each fire
d.plot.td <- as.data.table(d.plot)

fire.survey.yrs <- d.plot.td[,list(Fire=Fire,
                                    surv.post=survey.years.post,
                                    nplots=.N),
                              by=list(Fire,survey.years.post)]


fire.multsurveys <- fire.survey.yrs[,list(Fire=Fire,
                                          n.surv.yrs <- .N,
                                          surv.post= paste(surv.post,collapse=", ")
                                          ),
                                    by=list(Fire)]

fires.focal <- c("BTU Lightning","Straylor","Cub","Freds") #fires with enough revisit plots

d.plot <- d.plot[d.plot$Fire %in% fires.focal,]
d.plot <- d.plot[!(d.plot$Fire == "BTU Lightning" & d.plot$survey.years.post == 6),] ## don't want BTU 6-year plots

# d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Easting"] <- d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Easting"] - 100
# d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Northing"] <- d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Northing"] + 200




#export data to display in GIS
#write.csv(d.plot,"../geospatial output/plot.csv",row.names=FALSE)

coords <- as.matrix(d.plot[,c("Easting","Northing")])
plot.spatial <- st_multipoint(x=coords,dim="XY")
d.plot$plot.spatial <- plot.spatial

d.plot <- st_as_sf(d.plot,coords=c("Easting","Northing"),crs=26910)


## for each plot, get the closest nearby plot (of a diff survey year) and its distance

d.plot$revisit.plotname <- NA
d.plot$revisit.distance <- NA


if(sum(duplicated(d.plot$Regen_Plot)) > 0) {
  stop("Duplicated plot names")  
}


for (i in 1:length(unique(d.plot$Fire))) {
  
  fire <- unique(d.plot$Fire)[i]
  d.plot.fire <- d.plot[d.plot$Fire ==fire,]
  orig.year <- min(unique(d.plot.fire$survey.years.post))
  d.plot.fire.orig <- d.plot[d.plot$Fire ==fire & d.plot$survey.years.post==orig.year,]
  later.years <- unique(d.plot.fire$survey.years.post)
  later.years <- later.years[later.years > (orig.year + 1)]
  d.plot.fire.revisit <- d.plot[d.plot$Fire == fire & d.plot$survey.years.post %in% later.years,]
  
  for (plotname in unique(d.plot.fire.orig$Regen_Plot)) {
    
    d.plot.foc <- d.plot.fire.orig[d.plot.fire.orig$Regen_Plot == plotname,]
    
    # get the closest revisit plot on the same vire
    dists <- st_distance(d.plot.foc,d.plot.fire.revisit)
    
    mindist <- min(dists)
    closest.row <- which(dists == mindist)[1]
    
    closest.plot.name <- d.plot.fire.revisit[closest.row,"Regen_Plot"]
    st_geometry(closest.plot.name) <- NULL
    
    d.plot[(d.plot$Regen_Plot == plotname),c("revisit.plotname","revisit.distance")] <- cbind(closest.plot.name,mindist)

    
  }
}


d.plot$revisit.distance <- round(d.plot$revisit.distance)

d.plot.orig <- d.plot[!is.na(d.plot$revisit.distance),]
d.plot.orig <- d.plot.orig[d.plot.orig$revisit.distance < 100,]









