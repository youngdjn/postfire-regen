setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(ggplot2)

source("regen_analysis_functions.R")

#### 0. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

# only keep the necessary columns
d.plot <- d.plot[,c("Regen_Plot","Fire","Year.of.Fire","Easting","Northing","aspect","slope","SHRUB","FORBE","GRASS","HARDWOOD","CONIFER","FIRE_SEV","BA.Live1","Year","firesev","dist.to.low","fire.abbr","X5yr","fire.year","survey.years.post","elev.m","rad.march","tmean.post","ppt.post","ppt.post.min","tmean.normal","ppt.normal","seed.tree.any")]

# thin to 5-year post fire plots #! need to allow 4-year fires once we add them and fix climate and regen summarization
d.plot <- d.plot[d.plot$survey.years.post == 5,]

# only Sierra Nevada fires #! need to add new fires to this list when they're in the dataset
sierra.fires <- c("STRAYLOR","CUB","RICH","DEEP","MOONLIGHT","ANTELOPE","BTU LIGHTENING","HARDING","BASSETTS","PENDOLA","AMERICAN RIVER","RALSTON","FREDS","SHOWERS","POWER")
d.plot <- d.plot[d.plot$Fire %in% sierra.fires,]

# fix incorrectly-named variable
d.plot$FORB <- d.plot$FORBE

# if no data on seed tree distance (or it was recorded as being at/beyond the limit of the laser) use remote-sensing-based value
d.plot$seed.tree.any.comb <- ifelse(is.na(d.plot$seed.tree.any) | (d.plot$seed.tree.any >= 150),d.plot$dist.to.low,d.plot$seed.tree.any)

# quadratic climate terms
d.plot$ppt.normal.sq <- d.plot$ppt.normal^2
d.plot$tmean.normal.sq <- d.plot$tmean.normal^2

# variable for regen presence/absence
d.sp$regen.presab.young <- ifelse(d.sp$regen.count.young > 0,TRUE,FALSE)
d.sp$regen.presab.old <- ifelse(d.sp$regen.count.old > 0,TRUE,FALSE)
d.sp$regen.presab.tot <- ifelse(d.sp$regen.count.tot > 0,TRUE,FALSE)

#! TEMPORARY: if there was no radiation data, set it equal to 0
d.plot$rad.march <- ifelse(is.na(d.plot$rad.march),0,d.plot$rad.march)

high.sev <- c(4,5) # which field-assessed severity values are considered high severity
control <- c(0,1) # which field-assessed severity values are considered controls




#### 1. Assign each plot a topoclimatic category ####

#! NOTE that when breaking down plots by factorial combinations topoclimatic variables,
# it might be important in the future to consider the interaction of variables.
# E.g., not all levels of a given variable may be available at all levels of another given variables
# E.g. at high precipitation, maybe there is only north aspects available. Doesn't seem to be the case here, but could potentially be with other variables.

fires <- unique(d.plot$Fire)
d.plot$precip.category
d.plot$rad.category #radiation

for(fire in fires) {
  
  ## Precipitation categories
  # determine what the precipitation breakpoints should be (here just using median) -- based on high severity plots only
  breaks <- quantile(d.plot[(d.plot$Fire == fire) & (d.plot$FIRE_SEV %in% high.sev),]$ppt.normal,probs=c(0.5),na.rm=TRUE)
  # categorize plots based on where they fall between the breakpoints  
  categories <- categorize(d.plot[d.plot$Fire==fire,]$ppt.normal,breaks,name="P")
  # store it into the plot data.frame
  d.plot[d.plot$Fire==fire,"precip.category"] <- categories
  
  ## Radiation categories #! note straylor cub do not have radiation yet
  # determine what the breakpoints should be (here just using median) -- based on high severity plots only
  breaks <- quantile(d.plot[(d.plot$Fire == fire) & (d.plot$FIRE_SEV %in% high.sev),]$rad.march,probs=c(0.5),na.rm=TRUE)
  # categorize plots based on where they fall between the breakpoints  
  categories <- categorize(d.plot[d.plot$Fire==fire,]$rad.march,breaks,name="R")
  # store it into the plot data.frame
  d.plot[d.plot$Fire==fire,"rad.category"] <- categories
  
  #! To-do: break down the categories further for the fires that have enough plots (may have different number of categories per fire). Goal of ~10 plots per category? Seed tree distance should probably be a category (or else simply exclude plots that are far from seed source)

}
  

## Create one variable reflecting the all-way factorial combination of topoclimatic categories
d.plot$topoclim.cat <- paste(d.plot$precip.category,d.plot$rad.category,sep="_")
  

### Plot relevant "topoclimate space" for each fire and see how the categories broke them down
## note that Cub and straylor do not have radiation (yet)

# look at high sev only
ggplot(d.plot[d.plot$FIRE_SEV %in% high.sev,],aes(x=ppt.normal,y=rad.march,col=topoclim.cat)) +
  geom_point() +
  facet_wrap(~Fire,scales="free")

# look at control only
ggplot(d.plot[d.plot$FIRE_SEV %in% control,],aes(x=ppt.normal,y=rad.march,col=topoclim.cat)) +
  geom_point() +
  facet_wrap(~Fire,scales="free")

#! obvious that we need more control plots for some of these topoclim categories (many are coming when we add this summer's data)





#### 2. Summarize (compute average) regen values (high-sev plots only) and adults (control plots only) by species across all plots in each topoclimatic category in each fire ####

## assign the trees by species their topoclimatic category and fire name. This also ensures that we only are looking at seedlings for whose plots we are interested (because with this merge operation, seedlings from plots not in d.plot will be dropped)
d.sp.cat <- merge(d.sp,d.plot[,c("Regen_Plot","topoclim.cat","Fire","FIRE_SEV")])

## preparing to aggregate tree data: get highsev and control plots only, each with only the columns relevant to it
d.sp.cat.highsev <- d.sp.cat[d.sp.cat$FIRE_SEV %in% high.sev,c("Fire","species","topoclim.cat","seed.tree.sp","regen.count.young","regen.count.old","regen.count.tot","regen.presab.young","regen.presab.old","regen.presab.tot")]
d.sp.cat.control <- d.sp.cat[d.sp.cat$FIRE_SEV %in% high.sev,c("Fire","species","topoclim.cat","adult.count","adult.ba")]

## aggregate tree data by species and topo category #! might want to also calculate SD or SE in a separate aggregate call? (would append the variable names with ".sd")
d.sp.agg.highsev <- aggregate(d.sp.cat.highsev[,c(-1,-2,-3)],by=list(d.sp.cat.highsev$species,d.sp.cat.highsev$topoclim.cat,d.sp.cat.highsev$Fire),FUN=mean)
names(d.sp.agg.highsev)[1:3] <- c("species","topoclim.cat","Fire")
d.sp.agg.control <- aggregate(d.sp.cat.control[,c(-1,-2,-3)],by=list(d.sp.cat.control$species,d.sp.cat.control$topoclim.cat,d.sp.cat.control$Fire),FUN=mean)
names(d.sp.agg.control)[1:3] <- c("species","topoclim.cat","Fire")

## merge the control (adults only) and highsev (seedlings only) tree data
d.sp.agg <- merge(d.sp.agg.highsev,d.sp.agg.control)

##preparing to aggregate plot data: label plots as highsev or control
# first remove the variables that are not useful
d.plot.c <- remove.vars(d.plot,c("Year.of.Fire","Easting","Northing","aspect","Year","precip.category","rad.category"))
# label plots as control or high sev
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% control,"control",NA)
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% high.sev,"highsev",d.plot.c$type)
# get rid of plots that are neither control nor high sev
d.plot.c <- d.plot.c[!is.na(d.plot.c$type),] 

## aggregate plots by Fire, topoclim category, and type (control or high sev)
d.plot.agg.mean <- aggregate(remove.vars(d.plot.c,c("Regen_Plot","topoclim.cat","type","fire.abbr","X5yr","Fire")),by=list(d.plot.c$Fire,d.plot.c$topoclim.cat,d.plot.c$type),FUN=mean)
names(d.plot.agg.mean)[1:3] <- c("Fire","topoclim.cat","type")
d.plot.agg.tot <- aggregate(d.plot.c["Regen_Plot"],by=list(d.plot.c$Fire,d.plot.c$topoclim.cat,d.plot.c$type),FUN=length)
names(d.plot.agg.tot) <- c("Fire","topoclim.cat","type","count")

## merge plot averages and counts
d.plot.agg <- merge(d.plot.agg.mean,d.plot.agg.tot)
# separate the control and high.sev plots and change the column names (appending ".control" and ".highsev")
d.plot.agg.control <- d.plot.agg[d.plot.agg$type == "control",]
d.plot.agg.highsev <- d.plot.agg[d.plot.agg$type == "highsev",]
names(d.plot.agg.control)[c(-1,-2,-3)] <- paste(names(d.plot.agg.control)[c(-1,-2,-3)],"control",sep=".")
names(d.plot.agg.highsev)[c(-1,-2,-3)] <- paste(names(d.plot.agg.highsev)[c(-1,-2,-3)],"highsev",sep=".")
# merge them back together
d.plot.agg.merged <- merge(remove.vars(d.plot.agg.highsev,"type"),d.plot.agg.control,all=TRUE)

## rename the objects to shorter names
d.plot.2 <- d.plot.agg.merged
d.sp.2 <- d.sp.agg

### The two objects above contain all the relevant data summarized by topoclimatic category within each fire
## NOTE that since we aggregated by computing the average, the presence-absence columns represent "% of plots with presence"
## also adult data only come from the control plots, and seedling data only from the highsev plots