setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(party)
library(ggplot2)
library(brms)
library(pROC)
library(betareg)
library(car)
library(plyr)
library(data.table)
options(datatable.WhenJisSymbolThenCallingScope=TRUE)
library(sf)
library(viridis)

source("regen_analysis_functions.R")


#### 1. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp.ages <- read.csv("../data_intermediate_processing_local/tree_summarized_sp.csv",header=TRUE,stringsAsFactors=FALSE)

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
high.sev <- c(4,5) # which field-assessed severity values are considered high severity
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

fires.focal <- c("BTU Lightning","Straylor","Cub","Freds","Power") #fires with enough revisit plots

d.plot <- d.plot[d.plot$Fire %in% fires.focal,]
d.plot <- d.plot[!(d.plot$Fire == "BTU Lightning" & d.plot$survey.years.post == 6),] ## don't want BTU 6-year plots

# d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Easting"] <- d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Easting"] - 100
# d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Northing"] <- d.plot[d.plot$Fire == "Freds" & d.plot$survey.years.post == 8,"Northing"] + 200




#export data to display in GIS
# write.csv(d.plot,"../geospatial output/plot.csv",row.names=FALSE)

coords <- as.matrix(d.plot[,c("Easting","Northing")])
plot.spatial <- st_multipoint(x=coords,dim="XY")
d.plot$plot.spatial <- plot.spatial

d.plot <- st_as_sf(d.plot,coords=c("Easting","Northing"),crs=26910)




#### for each plot, get the closest nearby plot (of a diff survey year) and its distance ####

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
st_geometry(d.plot) <- NULL


#### Set max distance for considertation of nearby revisit ####
d.plot.orig <- d.plot[!is.na(d.plot$revisit.distance),]
d.plot.orig <- d.plot.orig[d.plot.orig$revisit.distance < 50,]





#### Get regen data from each plot (original and revisit), by species ####
d.plot.lookup <- d.plot.orig[,c("Regen_Plot","Fire","SHRUB","FIRE_SEV","revisit.plotname","revisit.distance")]

sp.opts <- c("PINUS.ALLSP","SHADE.ALLSP","PIPJ","ABCO","HDWD.ALLSP","QUKE")


plots.sps <- data.frame()


for(sp in sp.opts) {
  
  
  d.sp.focalsp <- d.sp[d.sp$species == sp,]
  
  
  for(i in 1:nrow(d.plot.lookup)) {
    
    orig.plt <- d.plot.lookup[i,"Regen_Plot"]
    revisit.plt <- d.plot.lookup[i,"revisit.plotname"]
    orig.sp <- d.sp.focalsp[d.sp.focalsp$Regen_Plot == orig.plt,]
    revisit.sp <- d.sp.focalsp[d.sp.focalsp$Regen_Plot == revisit.plt,]
    
    names(orig.sp) <- paste0(names(orig.sp),".orig")
    names(revisit.sp) <- paste0(names(revisit.sp),".revisit")
    
    plot.sp <- cbind(orig.sp,revisit.sp)
    plot.sp$species <- sp
    plot.sp$Fire <- d.plot.lookup[i,"Fire"]
    
    plots.sps <- rbind(plots.sps,plot.sp)
  
  }
}

# 
# 
# plot(regen.count.nonyoung.orig~regen.count.nonyoung.revisit,data=d.plot.orig)
# 
# table(d.plot.orig$regen.presab.all.orig,d.plot.orig$regen.presab.nonyoung.revisit)
# 
# 
# 
# 
# ggplot(d.plot.orig,aes(x=regen.presab.all.orig,y=regen.presab.nonyoung.revisit)) +
#   geom_bin2d() +
#   theme_bw()
# 
# 
# ggplot(d.plot.orig,aes(x=regen.presab.nonyoung.orig,fill=regen.presab.nonyoung.revisit)) +
#   geom_bar()
# 
# 




#### Aggregate regen data into original regen success categories (poor and good) ####

## prep

orig.type <- "all"
revisit.type <- "all"

orig.presab.var <- paste("regen.presab",orig.type,"orig",sep=".")
orig.count.var <- paste("regen.count",orig.type,"orig",sep=".")
revisit.presab.var <- paste("regen.presab",revisit.type,"revisit",sep=".")
revisit.count.var <- paste("regen.count",revisit.type,"revisit",sep=".")


plots.sps.keepvars <- plots.sps[,c(orig.presab.var,orig.count.var,revisit.presab.var,revisit.count.var,"species")]
names(plots.sps.keepvars) <- c("orig.presab","orig.count","revisit.presab","revisit.count","species")
plots.sps.keepvars <- data.table(plots.sps.keepvars)


## aggregate

revisit.agg <- plots.sps.keepvars[,list(orig.count = mean(orig.count),
                                        orig.count.low = quantile(orig.count,probs=.25),
                                        orig.count.high = quantile(orig.count,probs=.75),
                                        revisit.count = mean(revisit.count),
                                        revisit.count.low = quantile(revisit.count,probs=.25),
                                        revisit.count.high = quantile(revisit.count,probs=.75),
                                        revisit.prop = mean(revisit.presab),
                                        nplots = .N),
                                  by=list(species,orig.presab)]

revisit.agg.orig <- revisit.agg[,list(species,orig.presab,count = orig.count,nplots)]
revisit.agg.orig$survey <- "Initial"

revisit.agg.revisit <- revisit.agg[,list(species,orig.presab,count = revisit.count,prop=revisit.prop,nplots)]
revisit.agg.revisit$survey <- "Revisit"

d <- rbind.fill(revisit.agg.orig,revisit.agg.revisit)
d$labeltext <- paste0("n=",d$nplots)
d$labeltext[d$labeltext == "n=NA"] <- ""

#### Plot seedling counts ####
ggplot(d,aes(x=orig.presab,y=count,fill=survey)) +
  geom_bar(stat="identity",position="dodge",width=0.5) +
  facet_grid(~species) +
  theme_bw(14) +
  geom_text(aes(label=labeltext,y=-0.75)) +
  labs(x="Regeneration present at initial survey",y="Seedlings per plot")


#### Plot proportion plots with presence ####

## rename the species (only some of them--the others get dropped)



d.plot <- d[d$survey == "Revisit",]
ggplot(d.plot,aes(x=orig.presab,y=prop)) +
  geom_bar(stat="identity",position="dodge",width=0.5) +
  facet_grid(~species) +
  theme_bw(14) +
  geom_text(aes(label=labeltext,y=-0.02)) +
  labs(x="Regeneration present at initial survey",y="Proportion of plots with regeneration at revisit")

## patterns of older mimicking younger are stronger when considering all seedlings as opposed to "nonyoung" seedlings (in either or both survey times)




### for initial visit plot seedling counts by age for each fire and sp

d.plot.orig

d.sp.ages.int <- d.sp.ages[d.sp.ages$species %in% c("ABCO","PIPO"),]


#aggregate by sp and fire
ages.initial <- as.data.table(d.sp.ages.int)

ages.initial <- merge(ages.initial,d.plot.orig,by="Regen_Plot")

ages.init.agg <- ages.initial[,lapply(.SD,mean),by=list(species,Fire),.SDcols=3:29]

ages.init.melt <- melt(ages.init.agg,id.vars=c("species","Fire"),measure.vars=3:8)


ggplot(ages.init.melt,aes(x=variable,y=value)) +
  geom_bar(stat="identity") +
  facet_grid(species~Fire)


## now get the followup distrib for regen in plots that had and didn't in first survey

d.plot.orig

d.sp.ages.int <- d.sp.ages[d.sp.ages$species %in% c("ABCO","PIPO"),]


#aggregate by sp and fire

## need to load in whether it had regen of the given species (ABCO and PIPO first.)

ages.initial <- as.data.table(d.sp.ages.int)

ages.initial <- merge(ages.initial,d.plot.orig,by.x="Regen_Plot",by.y="revisit.plotname")

ages.init.agg <- ages.initial[,lapply(.SD,mean),by=list(species,Fire),.SDcols=3:29]

ages.init.melt <- melt(ages.init.agg,id.vars=c("species","Fire"),measure.vars=3:16)


orig.presab


ggplot(ages.init.melt,aes(x=variable,y=value)) +
  geom_bar(stat="identity") +
  facet_grid(species~Fire)





