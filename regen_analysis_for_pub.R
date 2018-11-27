# Analysis script for "Post-fire forest regeneration shows limited climate tracking and potential for drought-induced type conversion"
# By Derek J. N. Young, Chhaya M. Werner, Kevin R. Welch, Truman P. Young, Hugh D. Safford, Andrew M. Latimer
# This script performs data processing and statistical analyses to reproduce the analyses described in the article.
# The data files are described in the metadata document included in the Figshare repository.
# Correspondence to Derek Young: djyoung@ucdavis.edu

## Load required packages
library(ggplot2)
library(betareg)
library(plyr)
library(data.table)
library(gsubfn)
library(magrittr)
library(grid)
library(gridExtra)
library(viridis)
library(raster)
library(sf)
library(vegan)
library(ggvegan)

## Custom convenience functions used in script
source("regen_analysis_functions.R")

## Read in data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

## Create dataset of individual (non-aggregated) plots
d.plot.ind <- d.plot # object 

## label plots as control or high sev
# severity categories
high.sev <- c(4,5) # which field-assessed severity values are considered high severity
control <- c(0,1) # which field-assessed severity values are considered controls

d.plot.ind$type <- ifelse(d.plot.ind$FIRE_SEV %in% control,"control",NA)
d.plot.ind$type <- ifelse(d.plot.ind$FIRE_SEV %in% high.sev,"highsev",d.plot.ind$type)
# get rid of plots that are neither control nor high sev
d.plot.ind <- d.plot.ind[!is.na(d.plot.ind$type),]


#### 1. Assign each plot a topoclimatic category based on precipitation and radiation ####

d.plot.precat <- d.plot

# Remove climatic regions that don't have comparable control plots
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "American River") & (d.plot.precat$ppt.normal < 1700)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "BTU Lightning") & (d.plot.precat$ppt.normal > 2000)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "Bagley") & (d.plot.precat$rad.march < 4000)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "Cub") & (d.plot.precat$ppt.normal < 1300)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "Freds") & (d.plot.precat$ppt.normal > 1230)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "Power") & (d.plot.precat$rad.march < 6000)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "Ralston") & (d.plot.precat$ppt.normal < 1175)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "Rich") & (d.plot.precat$ppt.normal > 1080) & (d.plot.precat$rad.march < 5000)),]

fires <- unique(d.plot$Fire)

for(fire in fires) {
  
  ## Precipitation categories
  # determine what the precipitation breakpoints should be (here just using median) -- based on control plots only
  breaks <- quantile(d.plot.precat[(d.plot.precat$Fire == fire) & (d.plot.precat$FIRE_SEV %in% control),]$ppt.normal,probs=c(0.5),na.rm=TRUE)
  
  
  # for some fires with a small range of precip, override the precip breaks, so we just have one category per fire
  fires.small.precip.range <- c("American River","Antelope","Bagley","Bassetts","BTU Lightning","Harding","Straylor","Rich")
  if(fire %in% fires.small.precip.range) {
    breaks <- 9999
  }
  
  # categorize plots based on where they fall between the breakpoints  
  categories <- categorize(d.plot.precat[d.plot.precat$Fire==fire,]$ppt.normal,breaks,name="P")
  # store it into the plot data.frame
  d.plot.precat[d.plot.precat$Fire==fire,"precip.category"] <- categories
  
  ## Radiation categories
  # determine what the breakpoints should be (here just using median) -- based on high severity plots only
  breaks <- quantile(d.plot.precat[(d.plot.precat$Fire == fire) & (d.plot.precat$FIRE_SEV %in% high.sev),]$rad.march,probs=c(0.5),na.rm=TRUE)
  # categorize plots based on where they fall between the breakpoints
  
  #override the per-fire breaks
  breaks <- 6000
  
  categories <- categorize(d.plot.precat[d.plot.precat$Fire==fire,]$rad.march,breaks,name="R")
  # store it into the plot data.frame
  d.plot.precat[d.plot.precat$Fire==fire,"rad.category"] <- categories
  
}
  
## Create one variable reflecting the all-way factorial combination of topoclimatic categories
d.plot.precat$topoclim.cat <- paste(d.plot.precat$precip.category,d.plot.precat$rad.category,sep="_")

d.plot.cat <- d.plot.precat



#### 2. Summarize (compute average) regen values (high-sev plots only) and adults (control plots only) by species across all plots in each topoclimatic category in each fire ####

## assign the trees by species their topoclimatic category and fire name. This also ensures that we only are looking at seedlings for whose plots we are interested (because with this merge operation, seedlings from plots not in d.plot will be dropped)
d.sp.cat <- merge(d.sp,d.plot.cat[,c("Regen_Plot","topoclim.cat","Fire","FIRE_SEV","survey.years.post")])

## preparing to aggregate tree data: get highsev and control plots only, each with only the columns relevant to it
d.sp.cat.highsev <- d.sp.cat[d.sp.cat$FIRE_SEV %in% high.sev,c("Fire","species","topoclim.cat","regen.count.old","regen.count.all","regen.presab.old","regen.presab.all","survey.years.post")]
d.sp.cat.control <- d.sp.cat[d.sp.cat$FIRE_SEV %in% control,c("Fire","species","topoclim.cat","adult.count","adult.ba","survey.years.post")]

## for highsev plots (the ones where we're interested in regen), only consider plots surveyed 4-5 years post-fire
d.sp.cat.highsev <- d.sp.cat.highsev[d.sp.cat.highsev$survey.years.post %in% c(4,5),]

## aggregate tree data by species and topo category 
d.sp.agg.highsev <- aggregate(d.sp.cat.highsev[,c(-1,-2,-3)],by=list(d.sp.cat.highsev$species,d.sp.cat.highsev$topoclim.cat,d.sp.cat.highsev$Fire),FUN=mean,na.rm=TRUE)
names(d.sp.agg.highsev)[1:3] <- c("species","topoclim.cat","Fire")
d.sp.agg.control <- aggregate(d.sp.cat.control[,c(-1,-2,-3)],by=list(d.sp.cat.control$species,d.sp.cat.control$topoclim.cat,d.sp.cat.control$Fire),FUN=mean,na.rm=TRUE)
names(d.sp.agg.control)[1:3] <- c("species","topoclim.cat","Fire")

## merge the control (adults only) and highsev (seedlings only) tree data
d.sp.agg <- merge(d.sp.agg.highsev,d.sp.agg.control,all.x=TRUE,by=c("species","topoclim.cat","Fire"))

##preparing to aggregate plot (e.g. climate) data: label plots as highsev or control
# first remove the variables that are not useful
d.plot.c <- remove.vars(d.plot.cat,c("Year.of.Fire","Year","precip.category","rad.category"))
# label plots as control or high sev
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% control,"control",NA)
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% high.sev,"highsev",d.plot.c$type)
# get rid of plots that are neither control nor high sev
d.plot.c <- d.plot.c[!is.na(d.plot.c$type),]

## aggregate plots by Fire, topoclim category, and type (control or high sev)
d.plot.agg.mean <- aggregate(remove.vars(d.plot.c,c("Regen_Plot","topoclim.cat","type","fire.abbr","Fire","FIRE_SEV.cat")),by=list(d.plot.c$Fire,d.plot.c$topoclim.cat,d.plot.c$type),FUN=mean,na.rm=TRUE)
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

d.plot.agg.merged$count.control[is.na(d.plot.agg.merged$count.control)] <- 0
d.plot.agg.merged$count.highsev[is.na(d.plot.agg.merged$count.highsev)] <- 0

## rename the objects to shorter names
d.plot.2 <- d.plot.agg.merged
d.sp.2 <- d.sp.agg

### The two objects above contain all the relevant data summarized by topoclimatic category within each fire
## NOTE that since we aggregated by computing the average, the presence-absence columns represent "% of plots with presence"
## also adult data only come from the control plots, and seedling data only from the highsev plots





#### 3. Remove unneeded data rows and prep data frames for analysis ####

# Remove the topoclimatic categories with too few plots in either burned or control
d.plot.3 <- d.plot.2[which((d.plot.2$count.control > 4) & (d.plot.2$count.highsev > 0)),]

# only want to analyze high-severity plots burned 4-5 years post-fire
d.plot.c <- d.plot.c[(d.plot.c$survey.years.post %in% c(4,5)) & (d.plot.c$FIRE_SEV %in% c(4,5)),] 

# only want to analyze high-severity plots burned 4-5 years post-fire
d.plot.ind <- d.plot.ind[(d.plot.ind$survey.years.post %in% c(4,5)) & (d.plot.ind$FIRE_SEV %in% c(4,5)),] 







#### 4. Plot climate space of high-sev vs. reference plots ####

d.plot.keep <- merge(d.plot.cat,d.plot.3[,c("Fire","topoclim.cat","count.control","count.highsev")])
d.plot.keep <- d.plot.keep[d.plot.keep$count.control > 4 & d.plot.keep$count.highsev > 4,]

d.plot.keep$topoclim.cat <- gsub(".","",d.plot.keep$topoclim.cat,fixed=TRUE)

d.plot.keep$FIRE_SEV.cat <- gsub("control","Reference",d.plot.keep$FIRE_SEV.cat,fixed=TRUE)
d.plot.keep$FIRE_SEV.cat <- gsub("high.sev","High severity",d.plot.keep$FIRE_SEV.cat,fixed=TRUE)

p <- ggplot(d.plot.keep,aes(x=ppt.normal,y=rad.march,col=topoclim.cat,shape=FIRE_SEV.cat)) +
  geom_point(size=3) +
  facet_wrap(~Fire,scales="free") +
  theme_bw(16) +
  scale_shape_manual(values=c(1,3)) +
  labs(x="Normal annual precipitation (mm)",y="Solar exposure (Wh m-2 day-1)",shape="Plot type",color="Topoclimate category")

#Plot it  
p

#### 5. Plot climate space (Precip. normal vs. precip. anomalies) with fire labels (Fig. 2) ####

d.plot.ind$Fire <- as.factor(d.plot.ind$Fire)

### Plot of monitoring plots in climate space: minimum precip ### 

# For each fire, make a point at mean normal precip and mean precip anom 
fire.centers <- data.frame() 
fires <- unique(d.plot.ind$Fire) 
for(fire in fires) { 
  d.fire <- d.plot.ind[d.plot.ind$Fire == fire,] 
  norm.mean <- mean(d.fire$ppt.normal) 
  anom.mean <- mean(d.fire$diff.norm.ppt.min.z) 
  year <- d.fire$fire.year[1] 
  
  fire.center <- data.frame(Fire=fire,year=year,ppt.normal=norm.mean,diff.norm.ppt.min.z=anom.mean) 
  fire.centers <- rbind(fire.centers,fire.center)   
} 
fire.centers$fire.year <- paste0(fire.centers$Fire,", ",fire.centers$year) 


lines <- data.frame(rbind(
               #c(1900,-.27,2000,-0.42), # BTU
               #c(1650,-0.25,1500,-0.36), # CUB
               c(700,-0.8,750,-0.6) # moonlight
               ))
names(lines) <- c("x","y","xend","yend")

##shifting the fire name annotations so they look nice
fire.centers[fire.centers$Fire == "Antelope",c("ppt.normal","diff.norm.ppt.min.z")] <- c(500,-0.95)
fire.centers[fire.centers$Fire == "Freds",c("ppt.normal","diff.norm.ppt.min.z")] <- c(1170,-0.78)
fire.centers[fire.centers$Fire == "Power",c("ppt.normal","diff.norm.ppt.min.z")] <- c(1600,-0.85)
fire.centers[fire.centers$Fire == "Bassetts",c("ppt.normal","diff.norm.ppt.min.z")] <- c(2000,-0.95)
fire.centers[fire.centers$Fire == "Moonlight",c("ppt.normal","diff.norm.ppt.min.z")] <- c(750,-0.6)
fire.centers[fire.centers$Fire == "Straylor",c("ppt.normal","diff.norm.ppt.min.z")] <- c(400,-0.74)
fire.centers[fire.centers$Fire == "American River",c("ppt.normal","diff.norm.ppt.min.z")] <- c(1700,-0.14)
fire.centers[fire.centers$Fire == "Cub",c("ppt.normal","diff.norm.ppt.min.z")] <- c(1400,-0.40)
fire.centers[fire.centers$Fire == "BTU Lightning",c("ppt.normal","diff.norm.ppt.min.z")] <- c(2100,-0.37)
fire.centers[fire.centers$Fire == "Bagley",c("diff.norm.ppt.min.z")] <- -1.48
fire.centers[fire.centers$Fire == "Ralston",c("ppt.normal","diff.norm.ppt.min.z")] <- c(1400,-1.1)
fire.centers[fire.centers$Fire == "Chips",c("ppt.normal","diff.norm.ppt.min.z")] <- c(1250,-1.2)

p1 <-  ggplot(d.plot.ind,aes(y=diff.norm.ppt.min.z,x=ppt.normal,color=Fire)) + 
  geom_point(size=3) + 
  #geom_point(data=fire.centers,size=5,pch=1) + 
  geom_text(data=fire.centers,aes(label=fire.year),nudge_y=.04,size=6,color="darkgray") + 
  guides(color=FALSE) + 
  theme_bw(17) + 
  labs(x="Normal precipitation (mm)",y="Anomaly: post-fire minimum precipitation (SD)") + 
  scale_x_continuous(limits=c(180,2580))+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=lines,color="darkgray",size=1.0) +
  scale_color_viridis(discrete=TRUE)


### Plot of monitoring plots in climate space: mean precip ### 

# For each fire, make a point at mean normal precip and mean precip anom 
fire.centers <- data.frame() 
fires <- unique(d.plot.ind$Fire) 
for(fire in fires) { 
  d.fire <- d.plot.ind[d.plot.ind$Fire == fire,] 
  norm.mean <- mean(d.fire$ppt.normal) 
  anom.mean <- mean(d.fire$diff.norm.ppt.z) 
  year <- d.fire$fire.year[1] 
  
  fire.center <- data.frame(Fire=fire,year=year,ppt.normal=norm.mean,diff.norm.ppt.z=anom.mean) 
  fire.centers <- rbind(fire.centers,fire.center)   
} 
fire.centers$fire.year <- paste0(fire.centers$Fire,", ",fire.centers$year) 

fire.centers[fire.centers$Fire == "Antelope",c("diff.norm.ppt.z")] <- -0.49
fire.centers[fire.centers$Fire == "Moonlight",c("diff.norm.ppt.z")] <- -0.62
fire.centers[fire.centers$Fire == "Chips",c("ppt.normal","diff.norm.ppt.z")] <- c(1175,-0.90)
fire.centers[fire.centers$Fire == "Cub",c("diff.norm.ppt.z")] <- 0.22
fire.centers[fire.centers$Fire == "American River",c("diff.norm.ppt.z")] <- 0.44
fire.centers[fire.centers$Fire == "BTU Lightning",c("diff.norm.ppt.z")] <- 0.12
fire.centers[fire.centers$Fire == "Freds",c("diff.norm.ppt.z")] <- 0.35
fire.centers[fire.centers$Fire == "Power",c("diff.norm.ppt.z")] <- 0.61
fire.centers[fire.centers$Fire == "Bagley",c("diff.norm.ppt.z")] <- -.8
fire.centers[fire.centers$Fire == "Ralston",c("ppt.normal")] <- 1100


lines <- data.frame(rbind(c(1920,0.3,1920,0.2))) # BTU
names(lines) <- c("x","y","xend","yend")

p2 <- ggplot(d.plot.ind,aes(x=ppt.normal,y=diff.norm.ppt.z,color=Fire)) + 
  geom_point(size=3) + 
  #geom_point(data=fire.centers,size=5,pch=1) + 
  geom_text(data=fire.centers,aes(label=fire.year),nudge_y=.04,size=6,color="darkgray") + 
  guides(color=FALSE) + 
  theme_bw(17) + 
  labs(x="Normal precipitation (mm)",y="Anomaly: post-fire mean precipitation (SD)") + 
  scale_x_continuous(limits=c(180,2580)) +
  scale_color_viridis(discrete=TRUE)
  #geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=lines,color="darkgray",size=1.0)


blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob

#Plot it  
grid.arrange(p1,blank,p2,ncol=3,widths=c(0.49,0.02,0.49))




#### 6. Repeat cliamte space plots for AET (Fig. S3) ####

d.plot.ind$Fire <- as.factor(d.plot.ind$Fire)

### Plot of monitoring plots in climate space: minimum AET ### 

# For each fire, make a point at mean normal precip and mean precip anom 
fire.centers <- data.frame() 
fires <- unique(d.plot.ind$Fire) 
for(fire in fires) { 
  d.fire <- d.plot.ind[d.plot.ind$Fire == fire,] 
  norm.mean <- mean(d.fire$aet.normal) 
  anom.mean <- mean(d.fire$diff.norm.aet.min.z) 
  year <- d.fire$fire.year[1] 
  
  fire.center <- data.frame(Fire=fire,year=year,aet.normal=norm.mean,diff.norm.aet.min.z=anom.mean) 
  fire.centers <- rbind(fire.centers,fire.center)   
} 
fire.centers$fire.year <- paste0(fire.centers$Fire,", ",fire.centers$year) 


lines <- data.frame(rbind(
  #c(1900,-.27,2000,-0.42), # BTU
  #c(1650,-0.25,1500,-0.36), # CUB
  c(700,-0.8,750,-0.6) # moonlight
))
names(lines) <- c("x","y","xend","yend")

p1 <-  ggplot(d.plot.ind,aes(x=aet.normal,y=diff.norm.aet.min.z,color=Fire)) + 
  geom_point(size=3) + 
  #geom_point(data=fire.centers,size=5,pch=1) + 
  geom_text(data=fire.centers,aes(label=fire.year),nudge_y=.04,size=6,color="darkgray") + 
  guides(color=FALSE) + 
  theme_bw(17) + 
  labs(x="Normal AET (mm)",y="Anomaly: post-fire minimum AET (SD)") + 
  #scale_x_continuous(limits=c(180,2580))+
  #geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=lines,color="darkgray",size=1.0) +
  scale_color_viridis(discrete=TRUE)


### Plot of monitoring plots in climate space: mean AET ### 

# For each fire, make a point at mean normal precip and mean precip anom 
fire.centers <- data.frame() 
fires <- unique(d.plot.ind$Fire) 
for(fire in fires) { 
  d.fire <- d.plot.ind[d.plot.ind$Fire == fire,] 
  norm.mean <- mean(d.fire$aet.normal) 
  anom.mean <- mean(d.fire$diff.norm.aet.z) 
  year <- d.fire$fire.year[1] 
  
  fire.center <- data.frame(Fire=fire,year=year,aet.normal=norm.mean,diff.norm.aet.z=anom.mean) 
  fire.centers <- rbind(fire.centers,fire.center)   
} 
fire.centers$fire.year <- paste0(fire.centers$Fire,", ",fire.centers$year) 


lines <- data.frame(rbind(c(1920,0.3,1920,0.2))) # BTU
names(lines) <- c("x","y","xend","yend")


p2 <- ggplot(d.plot.ind,aes(x=aet.normal,y=diff.norm.aet.z,color=Fire)) + 
  geom_point(size=3) + 
  #geom_point(data=fire.centers,size=5,pch=1) + 
  geom_text(data=fire.centers,aes(label=fire.year),nudge_y=.04,size=6,color="darkgray") + 
  guides(color=FALSE) + 
  theme_bw(17) + 
  labs(x="Normal AET (mm)",y="Anomaly: post-fire mean AET (SD)") + 
  #scale_x_continuous(limits=c(180,2580)) +
  scale_color_viridis(discrete=TRUE)
#geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=lines,color="darkgray",size=1.0)


blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob

#Plot it  
grid.arrange(p1,blank,p2,ncol=3,widths=c(0.49,0.02,0.49))



#### 7. Plot-level regerssion analysis with GLM ####


### Center data frame ###

d <- d.plot.ind

# variables not to center
vars.leave <- c("SHRUB","GRASS","FIRE_SEV","Year","fire.year","survey.years.post","regen.count.old","regen.count.all","regen.presab.old","regen.presab.all","dominant_shrub_ht_cm","tallest_ht_cm")

# variables we use in models
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z", "def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z") # removed snow, adult.ba.agg
d <- d[complete.cases(d[,vars.focal]),]

d.c <- center.df(d,vars.leave)[["centered.df"]]

# data frame for making predictions from fitted models
d.center.dat <- center.df(d,vars.leave)[["center.data"]]

## Transformed variables for model fitting

d.c$SHRUB_c <- (d.c$SHRUB - mean(d.c$SHRUB))  / sd(d.c$SHRUB)

d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2
d.c$tmean.normal_c.sq <- d.c$tmean.normal_c^2
d.c$snow.normal_c.sq <- d.c$snow.normal_c^2

d.c$def.normal_c.sq <- d.c$def.normal_c^2
d.c$aet.normal_c.sq <- d.c$aet.normal_c^2

d.c$diff.norm.ppt.z_c.sq <- d.c$diff.norm.ppt.z_c^2
d.c$diff.norm.tmean.z_c.sq <- d.c$diff.norm.tmean.z_c^2
d.c$diff.norm.snow.z_c.sq <- d.c$diff.norm.snow.z_c^2

d.c$diff.norm.def.z_c.sq <- d.c$diff.norm.def.z_c^2
d.c$diff.norm.aet.z_c.sq <- d.c$diff.norm.aet.z_c^2

d.c$diff.norm.ppt.min.z_c.sq <- d.c$diff.norm.ppt.min.z_c^2
d.c$diff.norm.tmean.max.z_c.sq <- d.c$diff.norm.tmean.max.z_c^2

d.c$diff.norm.def.max.z_c.sq <- d.c$diff.norm.def.max.z_c^2
d.c$diff.norm.aet.min.z_c.sq <- d.c$diff.norm.aet.min.z_c^2

vars.focal.c <- paste0(vars.focal[-6],"_c")

## transform cover so it does not include 0 or 1 (for Beta distrib)
d.c$SHRUB.p <- d.c$SHRUB/100
d.c$SHRUB.pt <- (d.c$SHRUB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$GRASS.p <- d.c$GRASS/100
d.c$GRASS.pt <- (d.c$GRASS.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$Fire <- as.factor(d.c$Fire)

d.c.modfit <- d.c # the model fitting uses its own d.c that it modifies, so preserve original


## Set response options and loop through them

# species presence-absence
sp.opts <- c("HDWD.ALLSP","PIPJ","ABCO") # hardwood, yellow pine, white fir

# species cover
cover.opts <- c("COV.SHRUB","COV.GRASS")

# species height relative to shrubs ("height dominance")
ht.opts <- c("HT.PIPJ","HT.ABCO","HT.HDWD.ALLSP")

# all responses to fit models for
resp.opts <- c(sp.opts,cover.opts,ht.opts)

## variables to store outputs that are generated in the model-fitting loop
m.p <- list()
loos <- list()
dat.preds <- data.frame()
pred.obs <- data.frame()
d.loos.all <- data.frame()
d.loo.comps <- data.frame()
mods.best <- data.frame()
aucs <- data.frame()
pred.dat <- data.frame()
fit.dat <- data.frame()
d.maes.anoms <- data.frame()
rf.importance <- data.frame()
pred.rf <- data.frame()
fit.mods <- list()





sink("../run_output_V2F_temp.txt") # this file will store model fits etc
for(sp in resp.opts) {
  
  cat("\n\n#####")
  cat("Running model for: ",sp,"")
  cat("#####\n\n")

  ## select the correct plot-level and topoclimate category-level data for the given response variable (e.g., focal species)
  if(sp %in% cover.opts) {
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="PIPJ",] #pick any species; for cover it doesn't matter; just need to thin to one row per plot, because shrub cover is in the plot-level data (not the species-level data)
    d.sp.curr.plt <- d.sp[d.sp$species=="PIPJ",] 
  } else if(sp %in% ht.opts)  {
    sp.name <- substr(sp,4,100)
    d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp.name,]
    d.sp.curr.plt <- d.sp[d.sp$species==sp.name,]
  } else { # presence/absence
    d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp,]
    d.sp.curr.plt <- d.sp[d.sp$species==sp,]
  }
  
  #add the data for the current species response to the plot data
  d.c <- merge(d.c.modfit,d.sp.curr.plt,by=c("Regen_Plot"))
  
  # convert T/F to 0/1
  d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
  d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)
  
  ## Test whether seedling taller than shrub
  d.c$seedl.taller <- d.c$tallest_ht_cm > d.c$dominant_shrub_ht_cm
  d.c[d.c$regen.presab.old == FALSE,"seedl.taller"] <- 0 # where there is no seedling, it is not taller than shrub
  
  
  if(sp %in% cover.opts) {
    
    sp.cov <- substr(sp,5,100)
    sp.cov <- paste0(sp.cov,".pt")
    
    d.c$response.var <- d.c[,sp.cov]
    
  } else if(sp %in% ht.opts){
    
    d.c <- d.c[d.c$regen.presab.old == TRUE,] # this is where we select whether we want all plots or just plots where seedlings were present in the first place
    d.c$response.var <- d.c$seedl.taller
    
  } else  { # presence/absence
    
    d.c$response.var <- d.c$regen.presab.old.01
    
  }
  
  
    ### define candidate model formulas
    if(TRUE) { # to be able to collapse all these formulas
      formulas <- list()
      
      
      ### No seed tree ###
      
      formulas[["n0.a0"]] <- formula(response.var ~ 1)
      
      ## PPT
      
      formulas[["n0.aP"]] <- formula(response.var ~ diff.norm.ppt.z_c)
      formulas[["n0.aP2"]] <- formula(response.var ~ diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      
      
      formulas[["nP.a0"]] <- formula(response.var ~ ppt.normal_c)
      formulas[["nP.aP"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c)
      formulas[["nP.aP2"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2.a0"]] <- formula(response.var ~ ppt.normal_c + ppt.normal_c.sq)
      formulas[["nP2.aP"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2.aPdi"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq*diff.norm.ppt.z_c)
      formulas[["nP2.aP2"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      formulas[["nP.aPni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c)
      formulas[["nP.aP2ni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2.aPni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2.aP2ni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      


      ##DEF

      formulas[["n0.aD"]] <- formula(response.var ~ diff.norm.def.z_c)
      formulas[["n0.aD2"]] <- formula(response.var ~ diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)



      formulas[["nD.a0"]] <- formula(response.var ~ def.normal_c)
      formulas[["nD.aD"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c)
      formulas[["nD.aD2"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2.a0"]] <- formula(response.var ~ def.normal_c + def.normal_c.sq)
      formulas[["nD2.aD"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2.aDdi"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c.sq*diff.norm.def.z_c)
      formulas[["nD2.aD2"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
      formulas[["nD.aDni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c)
      formulas[["nD.aD2ni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2.aDni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2.aD2ni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



      ##AET


      formulas[["n0.aA"]] <- formula(response.var ~ diff.norm.aet.z_c)
      formulas[["n0.aA2"]] <- formula(response.var ~ diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)


      formulas[["nA.a0"]] <- formula(response.var ~     aet.normal_c)
      formulas[["nA.aA"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.z_c)
      formulas[["nA.aA2"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2.a0"]] <- formula(response.var ~     aet.normal_c +     aet.normal_c.sq)
      formulas[["nA2.aA"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq)
      formulas[["nA2.aAdi"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq*diff.norm.aet.z_c)
      formulas[["nA2.aA2"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)
      formulas[["nA.aAni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.z_c)
      formulas[["nA.aA2ni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2.aAni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c.sq)
      formulas[["nA2.aA2ni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)

      
      ##tmean
      
      formulas[["n0.aT"]] <- formula(response.var ~ diff.norm.tmean.z_c)
      formulas[["n0.aT2"]] <- formula(response.var ~ diff.norm.tmean.z_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
      
      
      
      formulas[["nT.a0"]] <- formula(response.var ~ tmean.normal_c)
      formulas[["nT.aT"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.z_c)
      formulas[["nT.aT2"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
      formulas[["nT2.a0"]] <- formula(response.var ~ tmean.normal_c + tmean.normal_c.sq)
      formulas[["nT2.aT"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq)
      formulas[["nT2.aTdi"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq*diff.norm.tmean.z_c)
      formulas[["nT2.aT2"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
      formulas[["nT.aTni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.z_c)
      formulas[["nT.aT2ni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
      formulas[["nT2.aTni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c.sq)
      formulas[["nT2.aT2ni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
      
      
      
      




      ## PPT
      formulas[["n0.aPmin"]] <- formula(response.var ~ diff.norm.ppt.min.z_c)
      formulas[["n0.aP2min"]] <- formula(response.var ~ diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)




      formulas[["nP.aPmin"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c)
      formulas[["nP.aP2min"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2.aPmin"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2.aPmindi"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq*diff.norm.ppt.min.z_c)
      formulas[["nP2.aP2min"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      formulas[["nP.aPminni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c)
      formulas[["nP.aP2minni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2.aPminni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2.aP2minni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)

       

      ##DEF

      formulas[["n0.aDmax"]] <- formula(response.var ~ diff.norm.def.max.z_c)
      formulas[["n0.aD2max"]] <- formula(response.var ~ diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

      formulas[["nD.aDmax"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c)
      formulas[["nD.aD2max"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2.aDmax"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2.aDmaxdi"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq*diff.norm.def.max.z_c)
      formulas[["nD2.aD2max"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
      formulas[["nD.aDmaxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c)
      formulas[["nD.aD2maxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2.aDmaxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2.aD2maxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



      ##AET

      formulas[["n0.aAmin"]] <- formula(response.var ~ diff.norm.aet.min.z_c)
      formulas[["n0.aA2min"]] <- formula(response.var ~ diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)


      formulas[["nA.aAmin"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.min.z_c)
      formulas[["nA.aA2min"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2.aAmin"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq)
      formulas[["nA2.aAmindi"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq*diff.norm.aet.min.z_c)
      formulas[["nA2.aA2min"]] <- formula(response.var ~     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)
      formulas[["nA.aAminni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.min.z_c)
      formulas[["nA.aA2minni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2.aAminni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c.sq)
      formulas[["nA2.aA2minni"]] <- formula(response.var ~     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)

      
      
      
      ##tmean
      
      formulas[["n0.aTmax"]] <- formula(response.var ~ diff.norm.tmean.max.z_c)
      formulas[["n0.aT2max"]] <- formula(response.var ~ diff.norm.tmean.max.z_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      
      formulas[["nT.aTmax"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.max.z_c)
      formulas[["nT.aT2max"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      formulas[["nT2.aTmax"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq)
      formulas[["nT2.aTmaxdi"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq*diff.norm.tmean.max.z_c)
      formulas[["nT2.aT2max"]] <- formula(response.var ~ tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)
      formulas[["nT.aTmaxni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.max.z_c)
      formulas[["nT.aT2maxni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      formulas[["nT2.aTmaxni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c.sq)
      formulas[["nT2.aT2maxni"]] <- formula(response.var ~ tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)
      
      
      
      
      
      

      ### Add seed tree ###      
      if(sp %in% c(sp.opts)) {
      

        formulas[["n0s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + 1)

        ## PPT
        
        formulas[["n0s.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.z_c)
        formulas[["n0s.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
        

        formulas[["nPs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c)
        formulas[["nPs.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c)
        formulas[["nPs.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
        formulas[["nP2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + ppt.normal_c.sq)
        formulas[["nP2s.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
        formulas[["nP2s.aPdi"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq*diff.norm.ppt.z_c)
        formulas[["nP2s.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
        formulas[["nPs.aPni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c)
        formulas[["nPs.aP2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
        formulas[["nP2s.aPni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
        formulas[["nP2s.aP2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
        
        
        
  
        ##DEF

        formulas[["n0s.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.z_c)
        formulas[["n0s.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)



        formulas[["nDs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c)
        formulas[["nDs.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c)
        formulas[["nDs.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
        formulas[["nD2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + def.normal_c.sq)
        formulas[["nD2s.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
        formulas[["nD2s.aDdi"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq*diff.norm.def.z_c)
        formulas[["nD2s.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
        formulas[["nDs.aDni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c)
        formulas[["nDs.aD2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
        formulas[["nD2s.aDni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
        formulas[["nD2s.aD2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



        ##AET


        formulas[["n0s.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.z_c)
        formulas[["n0s.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)


        formulas[["nAs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c)
        formulas[["nAs.aA"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c)
        formulas[["nAs.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
        formulas[["nA2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c +     aet.normal_c.sq)
        formulas[["nA2s.aA"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq)
        formulas[["nA2s.aAdi"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq*diff.norm.aet.z_c)
        formulas[["nA2s.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)
        formulas[["nAs.aAni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c)
        formulas[["nAs.aA2ni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
        formulas[["nA2s.aAni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c.sq)
        formulas[["nA2s.aA2ni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)

        
        ##tmean
        
        formulas[["n0s.aT"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.tmean.z_c)
        formulas[["n0s.aT2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.tmean.z_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
        
        
        
        formulas[["nTs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c)
        formulas[["nTs.aT"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c)
        formulas[["nTs.aT2"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
        formulas[["nT2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + tmean.normal_c.sq)
        formulas[["nT2s.aT"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq)
        formulas[["nT2s.aTdi"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq*diff.norm.tmean.z_c)
        formulas[["nT2s.aT2"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
        formulas[["nTs.aTni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c)
        formulas[["nTs.aT2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
        formulas[["nT2s.aTni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c.sq)
        formulas[["nT2s.aT2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
        
        
        
        
        

        
        
        
        ## PPT
        
        formulas[["n0s.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.min.z_c)
        formulas[["n0s.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        
        
        
        formulas[["nPs.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c)
        formulas[["nPs.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        formulas[["nP2s.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
        formulas[["nP2s.aPmindi"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq*diff.norm.ppt.min.z_c)
        formulas[["nP2s.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
        formulas[["nPs.aPminni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c)
        formulas[["nPs.aP2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        formulas[["nP2s.aPminni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
        formulas[["nP2s.aP2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
        
        
        
        ##DEF

        formulas[["n0s.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.max.z_c)
        formulas[["n0s.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

        formulas[["nDs.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c)
        formulas[["nDs.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
        formulas[["nD2s.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
        formulas[["nD2s.aDmaxdi"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq*diff.norm.def.max.z_c)
        formulas[["nD2s.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
        formulas[["nDs.aDmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c)
        formulas[["nDs.aD2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
        formulas[["nD2s.aDmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
        formulas[["nD2s.aD2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



        ##AET

        formulas[["n0s.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.min.z_c)
        formulas[["n0s.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)

        formulas[["nAs.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c)
        formulas[["nAs.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
        formulas[["nA2s.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq)
        formulas[["nA2s.aAmindi"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq*diff.norm.aet.min.z_c)
        formulas[["nA2s.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)
        formulas[["nAs.aAminni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c)
        formulas[["nAs.aA2minni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
        formulas[["nA2s.aAminni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c.sq)
        formulas[["nA2s.aA2minni"]] <- formula(response.var ~ seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)
      
        
        ##TMEAN
        
        formulas[["n0s.aTmax"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.tmean.max.z_c)
        formulas[["n0s.aT2max"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.tmean.max.z_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
        
        formulas[["nTs.aTmax"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c)
        formulas[["nTs.aT2max"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
        formulas[["nT2s.aTmax"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq)
        formulas[["nT2s.aTmaxdi"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq*diff.norm.tmean.max.z_c)
        formulas[["nT2s.aT2max"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)
        formulas[["nTs.aTmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c)
        formulas[["nTs.aT2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
        formulas[["nT2s.aTmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c.sq)
        formulas[["nT2s.aT2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)
        
        
        
        
        
      }
      
      
      ### With solar rad ###

      formulas[["n0r.a0"]] <- formula(response.var~ rad.march_c +  1)

      ## PPT

      formulas[["n0r.aP"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.z_c)
      formulas[["n0r.aP2"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)


      formulas[["nPr.a0"]] <- formula(response.var~ rad.march_c +  ppt.normal_c)
      formulas[["nPr.aP"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c)
      formulas[["nPr.aP2"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2r.a0"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + ppt.normal_c.sq)
      formulas[["nP2r.aP"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2r.aPdi"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq*diff.norm.ppt.z_c)
      formulas[["nP2r.aP2"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      formulas[["nPr.aPni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.z_c)
      formulas[["nPr.aP2ni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2r.aPni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2r.aP2ni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)



      ##DEF

      formulas[["n0r.aD"]] <- formula(response.var~ rad.march_c +  diff.norm.def.z_c)
      formulas[["n0r.aD2"]] <- formula(response.var~ rad.march_c +  diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)



      formulas[["nDr.a0"]] <- formula(response.var~ rad.march_c +  def.normal_c)
      formulas[["nDr.aD"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.z_c)
      formulas[["nDr.aD2"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2r.a0"]] <- formula(response.var~ rad.march_c +  def.normal_c + def.normal_c.sq)
      formulas[["nD2r.aD"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2r.aDdi"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.z_c + def.normal_c.sq*diff.norm.def.z_c)
      formulas[["nD2r.aD2"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
      formulas[["nDr.aDni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c)
      formulas[["nDr.aD2ni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2r.aDni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2r.aD2ni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



      ##AET


      formulas[["n0r.aA"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.z_c)
      formulas[["n0r.aA2"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)


      formulas[["nAr.a0"]] <- formula(response.var~ rad.march_c +      aet.normal_c)
      formulas[["nAr.aA"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.z_c)
      formulas[["nAr.aA2"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2r.a0"]] <- formula(response.var~ rad.march_c +      aet.normal_c +     aet.normal_c.sq)
      formulas[["nA2r.aA"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq)
      formulas[["nA2r.aAdi"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq*diff.norm.aet.z_c)
      formulas[["nA2r.aA2"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)
      formulas[["nAr.aAni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.z_c)
      formulas[["nAr.aA2ni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2r.aAni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.z_c +     aet.normal_c.sq)
      formulas[["nA2r.aA2ni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)

      
      ##TMEAN
      
      formulas[["n0r.aT"]] <- formula(response.var~ rad.march_c +  diff.norm.tmean.z_c)
      formulas[["n0r.aT2"]] <- formula(response.var~ rad.march_c +  diff.norm.tmean.z_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
      
      
      
      formulas[["nTr.a0"]] <- formula(response.var~ rad.march_c +  tmean.normal_c)
      formulas[["nTr.aT"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.z_c)
      formulas[["nTr.aT2"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
      formulas[["nT2r.a0"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + tmean.normal_c.sq)
      formulas[["nT2r.aT"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq)
      formulas[["nT2r.aTdi"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq*diff.norm.tmean.z_c)
      formulas[["nT2r.aT2"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
      formulas[["nTr.aTni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.z_c)
      formulas[["nTr.aT2ni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
      formulas[["nT2r.aTni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c.sq)
      formulas[["nT2r.aT2ni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
      
      
      
      

      
      

      ## PPT
      formulas[["n0r.aPmin"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.min.z_c)
      formulas[["n0r.aP2min"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)




      formulas[["nPr.aPmin"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c)
      formulas[["nPr.aP2min"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2r.aPmin"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2r.aPmindi"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq*diff.norm.ppt.min.z_c)
      formulas[["nP2r.aP2min"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      formulas[["nPr.aPminni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c)
      formulas[["nPr.aP2minni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2r.aPminni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2r.aP2minni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)



      ##DEF

      formulas[["n0r.aDmax"]] <- formula(response.var~ rad.march_c +  diff.norm.def.max.z_c)
      formulas[["n0r.aD2max"]] <- formula(response.var~ rad.march_c +  diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

      formulas[["nDr.aDmax"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c)
      formulas[["nDr.aD2max"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2r.aDmax"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2r.aDmaxdi"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq*diff.norm.def.max.z_c)
      formulas[["nD2r.aD2max"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
      formulas[["nDr.aDmaxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c)
      formulas[["nDr.aD2maxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2r.aDmaxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2r.aD2maxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



      ##AET

      formulas[["n0r.aAmin"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.min.z_c)
      formulas[["n0r.aA2min"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)


      formulas[["nAr.aAmin"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.min.z_c)
      formulas[["nAr.aA2min"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2r.aAmin"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq)
      formulas[["nA2r.aAmindi"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq*diff.norm.aet.min.z_c)
      formulas[["nA2r.aA2min"]] <- formula(response.var~ rad.march_c +      aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)
      formulas[["nAr.aAminni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.min.z_c)
      formulas[["nAr.aA2minni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2r.aAminni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c.sq)
      formulas[["nA2r.aA2minni"]] <- formula(response.var~ rad.march_c +      aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)

      
      
      ##TMEAN
      
      formulas[["n0r.aTmax"]] <- formula(response.var~ rad.march_c +  diff.norm.tmean.max.z_c)
      formulas[["n0r.aT2max"]] <- formula(response.var~ rad.march_c +  diff.norm.tmean.max.z_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      
      formulas[["nTr.aTmax"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.max.z_c)
      formulas[["nTr.aT2max"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      formulas[["nT2r.aTmax"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq)
      formulas[["nT2r.aTmaxdi"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq*diff.norm.tmean.max.z_c)
      formulas[["nT2r.aT2max"]] <- formula(response.var~ rad.march_c +  tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)
      formulas[["nTr.aTmaxni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.max.z_c)
      formulas[["nTr.aT2maxni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      formulas[["nT2r.aTmaxni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c.sq)
      formulas[["nT2r.aT2maxni"]] <- formula(response.var~ rad.march_c +  tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)
      
      
      
      

      ### Add seed tree ###
      if(sp %in% c(sp.opts)) {

        formulas[["n0sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + 1)

        ## PPT

        formulas[["n0sr.aP"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.z_c)
        formulas[["n0sr.aP2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)

        formulas[["nPsr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c)
        formulas[["nPsr.aP"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c)
        formulas[["nPsr.aP2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
        formulas[["nP2sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + ppt.normal_c.sq)
        formulas[["nP2sr.aP"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
        formulas[["nP2sr.aPdi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq*diff.norm.ppt.z_c)
        formulas[["nP2sr.aP2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
        formulas[["nPsr.aPni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c)
        formulas[["nPsr.aP2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
        formulas[["nP2sr.aPni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
        formulas[["nP2sr.aP2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)


        ##DEF

        formulas[["n0sr.aD"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.def.z_c)
        formulas[["n0sr.aD2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)

        formulas[["nDsr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c)
        formulas[["nDsr.aD"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c)
        formulas[["nDsr.aD2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
        formulas[["nD2sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + def.normal_c.sq)
        formulas[["nD2sr.aD"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
        formulas[["nD2sr.aDdi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq*diff.norm.def.z_c)
        formulas[["nD2sr.aD2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
        formulas[["nDsr.aDni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c)
        formulas[["nDsr.aD2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
        formulas[["nD2sr.aDni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
        formulas[["nD2sr.aD2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



        ##AET

        formulas[["n0sr.aA"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.z_c)
        formulas[["n0sr.aA2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)

        formulas[["nAsr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c)
        formulas[["nAsr.aA"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c)
        formulas[["nAsr.aA2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
        formulas[["nA2sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c +     aet.normal_c.sq)
        formulas[["nA2sr.aA"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq)
        formulas[["nA2sr.aAdi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c.sq*diff.norm.aet.z_c)
        formulas[["nA2sr.aA2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.z_c +     aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)
        formulas[["nAsr.aAni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c)
        formulas[["nAsr.aA2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
        formulas[["nA2sr.aAni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c.sq)
        formulas[["nA2sr.aA2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.z_c +     aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq +     aet.normal_c.sq)

        
        ##TMEAN
        
        formulas[["n0sr.aT"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.tmean.z_c)
        formulas[["n0sr.aT2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.tmean.z_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
        
        formulas[["nTsr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c)
        formulas[["nTsr.aT"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c)
        formulas[["nTsr.aT2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
        formulas[["nT2sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + tmean.normal_c.sq)
        formulas[["nT2sr.aT"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq)
        formulas[["nT2sr.aTdi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq*diff.norm.tmean.z_c)
        formulas[["nT2sr.aT2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
        formulas[["nTsr.aTni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c)
        formulas[["nTsr.aT2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq)
        formulas[["nT2sr.aTni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c.sq)
        formulas[["nT2sr.aT2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.z_c + tmean.normal_c + diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq)
        
        
        
        

        
        ## PPT

        formulas[["n0sr.aPmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.min.z_c)
        formulas[["n0sr.aP2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)


        formulas[["nPsr.aPmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c)
        formulas[["nPsr.aP2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        formulas[["nP2sr.aPmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
        formulas[["nP2sr.aPmindi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq*diff.norm.ppt.min.z_c)
        formulas[["nP2sr.aP2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
        formulas[["nPsr.aPminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c)
        formulas[["nPsr.aP2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        formulas[["nP2sr.aPminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
        formulas[["nP2sr.aP2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)



        ##DEF

        formulas[["n0sr.aDmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.def.max.z_c)
        formulas[["n0sr.aD2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

        formulas[["nDsr.aDmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c)
        formulas[["nDsr.aD2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
        formulas[["nD2sr.aDmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
        formulas[["nD2sr.aDmaxdi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq*diff.norm.def.max.z_c)
        formulas[["nD2sr.aD2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
        formulas[["nDsr.aDmaxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c)
        formulas[["nDsr.aD2maxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
        formulas[["nD2sr.aDmaxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
        formulas[["nD2sr.aD2maxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



        ##AET

        formulas[["n0sr.aAmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.min.z_c)
        formulas[["n0sr.aA2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)

        formulas[["nAsr.aAmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c)
        formulas[["nAsr.aA2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
        formulas[["nA2sr.aAmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq)
        formulas[["nA2sr.aAmindi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c.sq*diff.norm.aet.min.z_c)
        formulas[["nA2sr.aA2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c*diff.norm.aet.min.z_c +     aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)
        formulas[["nAsr.aAminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c)
        formulas[["nAsr.aA2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
        formulas[["nA2sr.aAminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c.sq)
        formulas[["nA2sr.aA2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c +     aet.normal_c + diff.norm.aet.min.z_c +     aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq +     aet.normal_c.sq)
      }
      
      ##TMEAN
      
      formulas[["n0sr.aTmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.tmean.max.z_c)
      formulas[["n0sr.aT2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.tmean.max.z_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      
      formulas[["nTsr.aTmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c)
      formulas[["nTsr.aT2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      formulas[["nT2sr.aTmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq)
      formulas[["nT2sr.aTmaxdi"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c.sq*diff.norm.tmean.max.z_c)
      formulas[["nT2sr.aT2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c*diff.norm.tmean.max.z_c + tmean.normal_c*diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)
      formulas[["nTsr.aTmaxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c)
      formulas[["nTsr.aT2maxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq)
      formulas[["nT2sr.aTmaxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c.sq)
      formulas[["nT2sr.aT2maxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + tmean.normal_c + diff.norm.tmean.max.z_c + tmean.normal_c + diff.norm.tmean.max.z_c.sq + diff.norm.tmean.max.z_c.sq + tmean.normal_c.sq)

    }
    
    
    
    ## for each formula (candidate model), fit the mdoel and evaluate fit using cross-validation
    cv.results <- data.frame()
    
    cat("\n")
    
    for(i in 1:length(formulas)) {
      
      #cat("\rEvaluating model ",i," of ",length(formulas))
      
      formula.name <- names(formulas)[i]
      cv.result <- cvfun.fire(formula=formulas[[i]],data=d.c) # function to fit model and evaluate fit using cross-validation
      cv.result.named <- data.frame(model=names(formulas)[i],cv.result,stringsAsFactors=FALSE)
      cv.results <- rbind(cv.results,cv.result.named)
      
    }
    
    cv.results$model <- as.character(cv.results$model)
    

    ### finding the best normal climate model and then the corresponding best post-fire weather (anomaly) model and best post-fire absolute weather models
    
    ## normal climate models
    search <- ".a0" # what text to search for in model names (a0 means there's no anomaly term, so it's a normal climate-only model)
    normal.model.names <- cv.results[grep(search,cv.results$model,fixed=TRUE),]$mod
    # normal.model.names <- normal.model.names[normal.model.names != "n0.a0"] # exclude null model

    norm.search.opts <- c(Pmean = "n(P|0)2?s?r?c?\\.a0",
                            Dmean = "n(D|0)2?s?r?c?\\.a0",
                            Amean = "n(A|0)2?s?r?c?\\.a0",
                            Tmean = "n(T|0)2?s?r?c?\\.a0",
                            Pmin = "n(P|0)2?s?r?c?\\.a0",
                            Dmax = "n(D|0)2?s?r?c?\\.a0",
                            Amin = "n(A|0)2?s?r?c?\\.a0",
                            Tmax = "n(T|0)2?s?r?c?\\.a0"
                            )


    anom.search.opts <- c(Pmean = "aP2?(ni)?(di)?$",
                          Dmean = "aD2?(ni)?(di)?$",
                          Amean = "aA2?(ni)?(di)?$",
                          Tmean = "aT2?(ni)?(di)?$",
                          Pmin = "aP2?min(ni)?(di)?$",
                          Dmax = "aD2?max(ni)?(di)?$",
                          Amin = "aA2?min(ni)?(di)?$",
                          Tmax = "aT2?max(ni)?(di)?$"
                          )
    

    
    
    d.maes.anoms.sp <- data.frame()
    
    ## for each of these categories of models, find the best ones
    for(i in 1:length(anom.search.opts)) {
      
      anom.name <- names(norm.search.opts)[i]
      norm.search <- norm.search.opts[i]
      anom.search <- anom.search.opts[i]
      null.names <- c("n0.a0","n0s.a0","n0r.a0","n0sr.a0","n0c.a0","n0sc.a0","n0rc.a0","n0src.a0") # models with no climare varis, but potentially some other env. vars like seed tree distance
      normal.model.names <- cv.results[grep(norm.search,cv.results$model,fixed=FALSE),]$model
      
      d.cv.normal <- cv.results[cv.results$model %in% normal.model.names,]
      best.normal.mod <- d.cv.normal[which(d.cv.normal$mae == min(d.cv.normal$mae,na.rm=TRUE)),]$mod[1]
      
      #what is the best anomaly model corresponding to that normal model, sticking with whichever normal we have
      best.normal.normal.part <- strsplit(as.character(best.normal.mod),".",fixed=TRUE)[[1]][1]
      best.normal.anom.search <- paste0(best.normal.normal.part,"\\.",anom.search)
      anom.model.names <- cv.results[grep(best.normal.anom.search,cv.results$model,fixed=FALSE),]$model
      d.cv.anom <- cv.results[cv.results$model %in% anom.model.names,]
      best.anom.mod <- d.cv.anom[which(d.cv.anom$mae == min(d.cv.anom$mae,na.rm=TRUE)),]$mod[1]
      
      #what is the best normal, excluding null models
      normal.model.names.nonull <- normal.model.names[!(normal.model.names %in% null.names)]
      d.cv.normal.nonull <- cv.results[cv.results$model %in% normal.model.names.nonull,]
      best.normal.nonull.mod <- d.cv.normal.nonull[which(d.cv.normal.nonull$mae == min(d.cv.normal.nonull$mae,na.rm=TRUE)),]$mod[1]
      
      #what is the best anomaly model all on its own?
      anom.search <- anom.search.opts[i]
      anom.own.model.names <- cv.results[grep(anom.search,cv.results$model,fixed=FALSE),]$model
      d.cv.anom.own <- cv.results[cv.results$model %in% anom.own.model.names,]
      best.anom.own.mod <- d.cv.anom.own[which(d.cv.anom.own$mae == min(d.cv.anom.own$mae,na.rm=TRUE)),]$mod[1]
      
 
      ## get maes (mean absolute errors) of best normal, best anomal, best normal.nonull, best.post
      best.anom.mae <- cv.results[cv.results$model == best.anom.mod,"mae"]
      best.anom.own.mae <- cv.results[cv.results$model == best.anom.own.mod,"mae"]
      best.anom.normal.mae <- cv.results[cv.results$model == best.normal.mod,"mae"]
      best.normal.nonull.mae <- cv.results[cv.results$model == best.normal.nonull.mod,"mae"]

      
      ### Store MAEs of the best models
      best.anomaly.mod <- best.anom.mod
      best.anomaly.own.mod <- best.anom.own.mod
      best.anom.normal <- best.normal.mod
      best.normal.nonull <- best.normal.nonull.mod

      
      # data frame row of the best anomaly models corresponding to the best normal model
      d.maes.anoms.sp.anom <- data.frame(best.anomaly.mod,best.anom.normal,best.normal.nonull,best.anomaly.own.mod,
                                         best.anom.mae,best.anom.normal.mae,best.normal.nonull.mae,best.anom.own.mae,
                                         sp,anom.name,stringsAsFactors=FALSE)
      
      d.maes.anoms.sp <- rbind(d.maes.anoms.sp,d.maes.anoms.sp.anom)
      
    }
    
    d.maes.anoms <- rbind(d.maes.anoms,d.maes.anoms.sp) # this stores the best models of each type (e.g., normal climate model and post-fire anomaly model) and their associated MAEs (mean absolute errors)

    
    ### Now for each anom model , make predictions
    
    # focal variables
    vars <- c("ppt.normal_c","ppt.normal_c.sq","diff.norm.ppt.z_c","diff.norm.ppt.z_c.sq","tmean.normal_c","tmean.normal_c.sq",
              "diff.norm.tmean.z_c","diff.norm.tmean.z_c.sq",
              "aet.normal_c","aet.normal_c.sq","diff.norm.aet.z_c","diff.norm.aet.z_c.sq","def.normal_c","def.normal_c.sq",
              "diff.norm.def.z_c","diff.norm.def.z_c.sq",
              "seed_tree_distance_general_c","rad.march_c","SHRUB_c"
    )
  
    ## middle, low, and high hypothetical values for making predictions, based on observed ranges in data
    mid.val <- sapply(vars,mid.val.fun,USE.NAMES=TRUE)
    
    low.val <- sapply(vars,low.val.fun,USE.NAMES=TRUE)
    names(low.val) <- vars
    
    high.val <- sapply(vars,high.val.fun,USE.NAMES=TRUE)
    names(high.val) <- vars
    
    
    ## create hypothetical values for prediction
    diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100) # for the anomaly vars, range from -1.5 SD to 1.5 SD
    diff.norm.seq.rev <- rev(diff.norm.seq)
    
    newdat <- data.frame(
      ppt.normal_c = c(rep(low.val["ppt.normal_c"],100),rep(mid.val["ppt.normal_c"],100),rep(high.val["ppt.normal_c"],100)),
      ppt.normal_c.sq = c(rep(low.val["ppt.normal_c"]^2,100),rep(mid.val["ppt.normal_c"]^2,100),rep(high.val["ppt.normal_c"]^2,100)),
      tmean.normal_c = c(rep(low.val["tmean.normal_c"],100),rep(mid.val["tmean.normal_c"],100),rep(high.val["tmean.normal_c"],100)),
      tmean.normal_c.sq = c(rep(low.val["tmean.normal_c"]^2,100),rep(mid.val["tmean.normal_c"]^2,100),rep(high.val["tmean.normal_c"]^2,100)),
      def.normal_c = c(rep(low.val["def.normal_c"],100),rep(mid.val["def.normal_c"],100),rep(high.val["def.normal_c"],100)),
      def.normal_c.sq = c(rep(low.val["def.normal_c"]^2,100),rep(mid.val["def.normal_c"]^2,100),rep(high.val["def.normal_c"]^2,100)),
      aet.normal_c = c(rep(low.val["aet.normal_c"],100),rep(mid.val["aet.normal_c"],100),rep(high.val["aet.normal_c"],100)),
      aet.normal_c.sq = c(rep(low.val["aet.normal_c"]^2,100),rep(mid.val["aet.normal_c"]^2,100),rep(high.val["aet.normal_c"]^2,100)),
      
      norm.level = c(rep("low",100),rep("mid",100),rep("high",100)),
      
      diff.norm.ppt.z_c = rep(diff.norm.seq,3),
      diff.norm.ppt.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.tmean.z_c = rep(diff.norm.seq,3),
      diff.norm.tmean.z_c.sq = rep(diff.norm.seq,3)^2,
      diff.norm.aet.z_c = rep(diff.norm.seq,3),
      diff.norm.aet.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.def.z_c = rep(diff.norm.seq,3),
      diff.norm.def.z_c.sq = rep(diff.norm.seq,3)^2,
      diff.norm.snow.z_c = rep(diff.norm.seq,3),
      diff.norm.snow.z_c.sq = rep(diff.norm.seq,3)^2,
      
      diff.norm.ppt.min.z_c = rep(diff.norm.seq,3),
      diff.norm.ppt.min.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.tmean.max.z_c = rep(diff.norm.seq,3),
      diff.norm.tmean.max.z_c.sq = rep(diff.norm.seq,3)^2,
      diff.norm.aet.min.z_c = rep(diff.norm.seq,3),
      diff.norm.aet.min.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.def.max.z_c = rep(diff.norm.seq,3),
      diff.norm.def.max.z_c.sq = rep(diff.norm.seq,3)^2,
      diff.norm.snow.min.z_c = rep(diff.norm.seq,3),
      diff.norm.snow.min.z_c.sq = rep(diff.norm.seq,3)^2,
      
      seed_tree_distance_general_c = mid.val["seed_tree_distance_general_c"],
      rad.march_c = mid.val["rad.march_c"],
      SHRUB_c = mid.val["SHRUB_c"]
      
    )
    
    # for each of the best models for each anomaly type, make predictions
    for(j in 1:nrow(d.maes.anoms.sp)) {
      
      d.maes.anoms.sp.row <- d.maes.anoms.sp[j,]
      best.anom.mod <- d.maes.anoms.sp.row$best.anomaly.mod
      best.anom.normal.mod <- d.maes.anoms.sp.row$best.anom.normal
      
  
      ##predict to new data
      
      if(sp %in% c(cover.opts)) {
      
        
        d.c.complete <- d.c[!is.na(d.c$response.var),]
        
        nboot <- 100
        nobs <- nrow(d.c.complete)
        npred <- nrow(newdat)
        preds.anom.boot <- matrix(nrow=npred,ncol=nboot)
        preds.norm.boot <- matrix(nrow=npred,ncol=nboot)
        
        for(k in 1:nboot) {
        
          ## use bootstrapping to fit betareg models to get CIs
          d.c.boot <- d.c.complete[sample(nobs,size=nobs,replace=TRUE),]
          
          #fit the anom and normal models
          mod.anom <- try(betareg(formulas[[best.anom.mod]],data=d.c.boot))
          if(class(mod.anom) == "try-error") {
            message <- paste0("Model error in bootstrapping for ",sp," ",formula.name,"\n")
            print(message)
            next()
          }
          
          mod.norm <- try(betareg(formulas[[best.anom.normal.mod]],data=d.c.boot))
          if(class(mod.norm) == "try-error") {
            message <- paste0("Model error in bootstrapping for ",sp," ",formula.name,"\n")
            print(message)
            next()
          }
          
          #predict the median
          pred.anom <-predict(mod.anom,newdat,type="response")
          pred.norm <- predict(mod.norm,newdat,type="response")
          
          preds.anom.boot[,k] <- pred.anom
          preds.norm.boot[,k] <- pred.norm
        
        }
        
        pred.anom.mid <- apply(preds.anom.boot,1,quantile,probs=0.5,na.rm=TRUE)
        pred.anom.lwr <- apply(preds.anom.boot,1,quantile,probs=0.025,na.rm=TRUE)
        pred.anom.upr <- apply(preds.anom.boot,1,quantile,probs=0.975,na.rm=TRUE)
        
        pred.anom <- data.frame(pred.low=pred.anom.lwr,pred.mid=pred.anom.mid,pred.high=pred.anom.upr)
        
        pred.norm.mid <- apply(preds.norm.boot,1,quantile,probs=0.5,na.rm=TRUE)
        pred.norm.lwr <- apply(preds.norm.boot,1,quantile,probs=0.025,na.rm=TRUE)
        pred.norm.upr <- apply(preds.norm.boot,1,quantile,probs=0.975,na.rm=TRUE)
        
        pred.norm <- data.frame(pred.low=pred.norm.lwr,pred.mid=pred.norm.mid,pred.high=pred.norm.upr)
        
  
      } else { # presence-absence
        
        d.c.complete <- d.c[complete.cases(d.c$response.var),]
        
        #fit the anom and normal models
        mod.anom <- glm(formulas[[best.anom.mod]],data=d.c,family="binomial")
        mod.norm <- glm(formulas[[best.anom.normal.mod]],data=d.c,family="binomial")
        
        p <- predict(mod.anom,newdat,type="link",se.fit=TRUE)
        low <- p$fit - 1.96*p$se
        high <- p$fit + 1.96*p$se
        high[high>400] <- 400
        pred.anom <- data.frame(pred.low=inv.logit(low),pred.mid=inv.logit(p$fit),pred.high=inv.logit(high))
        
        p <- predict(mod.norm,newdat,type="link",se.fit=TRUE)
        low <- p$fit - 1.96*p$se
        high <- p$fit + 1.96*p$se
        high[high>400] <- 400
        pred.norm <- data.frame(pred.low=inv.logit(low),pred.mid=inv.logit(p$fit),pred.high=inv.logit(high))
        
      }
      

      print("Anom mod")
      print(coef(mod.anom))
      print("Norm mod")
      print(coef(mod.norm))
      
      names(pred.anom) <- c("pred.low","pred.mid","pred.high")
      names(pred.norm) <- c("pred.low","pred.mid","pred.high")
      
      pred.anom.dat <- cbind(pred.anom,newdat)
      pred.anom.dat$type <- "anom"
      pred.anom.dat$mod <- d.maes.anoms.sp.row$best.anomaly.mod
      
      pred.norm.dat <- cbind(pred.norm,newdat)
      pred.norm.dat$type <- "norm"
      pred.norm.dat$mod <- d.maes.anoms.sp.row$best.anom.normal
      
      pred.dat.sp <- rbind(pred.norm.dat,pred.anom.dat)
      #pred.dat.sp <- pred.anom.dat
      
      pred.dat.sp$anom <- d.maes.anoms.sp.row$anom.name
      pred.dat.sp$sp <- d.maes.anoms.sp.row$sp
      pred.dat.sp$rad.level <- d.maes.anoms.sp.row$rad.level

      
      pred.dat <- rbind(pred.dat,pred.dat.sp)  

      ## get fitted and observed values (for figure, etc.)
      d.c.complete <- d.c[!is.na(d.c$response.var),]
      
      if(sp %in% c(sp.opts,ht.opts)) {
        
        mod.anom <- glm(formulas[[best.anom.mod]],data=d.c.complete,family="binomial")
        mod.norm <- glm(formulas[[best.anom.normal.mod]],data=d.c.complete,family="binomial")
        

      } else  { # cover
        
        mod.anom <- betareg(formulas[[best.anom.mod]],data=d.c.complete)
        mod.norm <- betareg(formulas[[best.anom.normal.mod]],data=d.c.complete)
        
      } 
      
      # store the fitted model to use for later operations
      fit.mods[[paste0(sp,"_",best.anom.mod)]] <- mod.anom
      fit.mods[[paste0(sp,"_",best.anom.normal.mod)]] <- mod.norm
      
      # compute and store fitted values
      fit.anom <- as.data.frame(predict(mod.anom))
      fit.norm <- as.data.frame(predict(mod.norm))
      names(fit.anom) <- names(fit.norm) <- c("fitted")
      fit.anom$type <- "anom"
      fit.norm$type <- "norm"
      
      d.c.complete <- d.c[complete.cases(d.c$response.var),]
      
      # store fitted values along with predictors
      fit.anom.dat <- cbind(fit.anom,d.c.complete)
      fit.norm.dat <- cbind(fit.norm,d.c.complete)
      
      fit.dat.sp <- rbind(fit.anom.dat,fit.norm.dat)
      fit.dat.sp$anom <- d.maes.anoms.sp.row$anom.name
      fit.dat.sp$sp <- sp
      
      if(sp %in% sp.opts) {
        
        fit.dat.sp$fitted <- inv.logit(fit.dat.sp$fitted)
        
      }
      
      fit.dat <- rbind.fill(fit.dat,fit.dat.sp)
      

    }
   
}
sink(file=NULL)


save.image("../model_workspace_V2F3.Rwkspc")
load("../model_workspace_V2F3.Rwkspc")




#### 8. Make overall exploration figure of model fits for each response (Fig. S1) ####

### Prep  ###

#for each sp-rad combo, find which was the best anomaly term, which anomalies were better than their corresponding normals, and plot symbols indicating
d.maes.anoms$anom.better <- d.maes.anoms$best.anom.mae < d.maes.anoms$best.anom.normal.mae
d.maes.anoms$anom.improvement <-  d.maes.anoms$best.anom.normal.mae - d.maes.anoms$best.anom.mae

d.maes.anoms$best.of.species <- ""
d.maes.anoms$most.improved <- ""

d.maes.anoms.short <- d.maes.anoms[,c("sp","anom.name","anom.better","anom.improvement","best.of.species","most.improved")]

pred.dat.comb <- merge(pred.dat,d.maes.anoms.short,all.x=TRUE,by.x=c("sp","anom"),by.y=c("sp","anom.name"))

pred.dat.comb$anom.improvement <- round(pred.dat.comb$anom.improvement,4)
#pred.dat.comb[pred.dat.comb$anom.improvement < -0.01 , "anom.improvement"] <- NA


#Compute total MAE across all species, to see which is the best #
d.mae.agg <- aggregate(d.maes.anoms[,c("best.anom.mae","best.anom.normal.mae","anom.improvement")],by=list(d.maes.anoms$anom.name),FUN=mean,na.rm=TRUE)



### Plot counterfactuals (predictions for hypothetical values of the predictor variables) ###

pred.dat.plotting <- pred.dat.comb



pred.dat.plotting <- as.data.frame(pred.dat.plotting)

#for every column in predictions, if it ends in _c, uncenter it

for(col.name in names(pred.dat.plotting)) {
  if(grepl("_c$",col.name)) { #it's a centered col
    col.name.uncentered <- substr(col.name,1,nchar(col.name)-2)
    
    col.mean <- d.center.dat[d.center.dat$var == col.name,"var.mean"]
    col.sd <- d.center.dat[d.center.dat$var == col.name,"var.sd"]
    
    pred.dat.plotting[,as.character(col.name.uncentered)] <- pred.dat.plotting[,col.name] * col.sd + col.mean
    
    
  }
}

pred.dat.plotting <- data.table(pred.dat.plotting)

pred.dat.plotting$pred.val <- NA

# Turn predictor variable names into intelligible names for plotting
center.rename <- list("Pmin"="diff.norm.ppt.min.z_c","Pmean"="diff.norm.ppt.z_c","Amin"="diff.norm.aet.min.z_c","Amean"="diff.norm.aet.z_c","Dmax"="diff.norm.def.max.z_c","Dmean"="diff.norm.def.z_c",   "Tmax"="diff.norm.tmean.max.z_c","Tmean"="diff.norm.tmean.z_c")
uncenter.rename <- list("Pmin"="diff.norm.ppt.min.z","Pmean"="diff.norm.ppt.z","Amin"="diff.norm.aet.min.z","Amean"="diff.norm.aet.z","Dmax"="diff.norm.def.max.z","Dmean"="diff.norm.def.z",   "Tmax"="diff.norm.tmean.max.z","Tmean"="diff.norm.tmean.z")

pred.dat.plotting$pred.center.name <- gsubfn("\\S+",center.rename,as.character(pred.dat.plotting$anom))
pred.dat.plotting$pred.uncenter.name <- gsubfn("\\S+",uncenter.rename,as.character(pred.dat.plotting$anom))

##drop columns that weren't translated (deficit)
pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$pred.center.name %in% center.rename,]
pred.dat.plotting <- as.data.table(pred.dat.plotting)

#### if the anom model was better, plot the anom; otherwise, plot the baseline
pred.dat.plotting <- pred.dat.plotting[(pred.dat.plotting$anom.better & pred.dat.plotting$type == "anom") |
                                         (!pred.dat.plotting$anom.better & pred.dat.plotting$type == "norm"),]


pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$sp %in% c("ABCO","COV.GRASS","COV.SHRUB","COV.FORB","HDWD.ALLSP","CONIF.ALLSP","HT.CONIF.ALLSP","HT.ABCO","HT.HDWD.ALLSP","HT.PIPJ","PINUS.ALLSP","PIPJ","SHADE.ALLSP"),]
#pred.dat.plotting <- pred.dat.plotting[!(pred.dat.plotting$anom %in% c("Dmax","Dmean")),]

## make species names intelligible for plotting
spsub <- list("PIPJ" = "Yellow pine\npresence\n(% of plots)","ABCO" = "White fir\npresence\n(% of plots)","CONIF.ALLSP" = "Conifer\npresence\n(% of plots)","PSME" = "Douglas-fir\npresence\n(% of plots)","QUKE"="Black oak\npresence\n(% of plots)","PINUS.ALLSP" = "All pines\npresence\n(% of plots)","SHADE.ALLSP" = "Shade-tolerant conifers\npresence\n(% of plots)","HDWD.ALLSP" = "Broadleaved trees\npresence\n(% of plots)","COV.SHRUB" = "Shrubs\n(% cover)","COV.GRASS" = "Graminoids\n(% cover)","COV.FORB" = "Forbs\n(% cover)","HT.PIPJ"="Yellow pine\nheight dominance\n(% of plots)","HT.ABCO"="White fir\nheight dominance\n(% of plots)","HT.PINUS.ALLSP"="Pines\nheight dominance\n(% of plots)","HT.SHADE.ALLSP"="Shde-tolerant conifers\nheight dominance\n(% of plots)","HT.HDWD.ALLSP"="Broadleaved trees\nheight dominance\n(% of plots)","HT.CONIF.ALLSP"="Conifer\nheight dominance\n(% of plots)")
pred.dat.plotting$sp <- gsubfn("\\S+",spsub,as.character(pred.dat.plotting$sp))
pred.dat.plotting$sp <- factor(pred.dat.plotting$sp,levels=spsub)

## make clim var names intelligible for plotting
anomsub <- list("Pmin" = "Min precip","Pmean" = "Mean precip","Amin" = "Min AET","Amean" = "Mean AET","Dmax" = "Max CWD","Dmean" = "Mean CWD",     "Tmax" = "Max temp","Tmean" = "Mean temp")
pred.dat.plotting$anom <- gsubfn("\\S+",anomsub,as.character(pred.dat.plotting$anom))
pred.dat.plotting$anom <- factor(pred.dat.plotting$anom,levels=anomsub)

## if the anom model was worse than the null model, make the lines dashed and remove the confidence bands
pred.dat.plotting[which(pred.dat.plotting$anom.improvement < 0),c("pred.low","pred.high")] <- NA
pred.dat.plotting$linestyle = 1 # solid, by default
pred.dat.plotting[which(pred.dat.plotting$anom.improvement < 0),"linestyle"] <- 2 # dashed if the anomaly model was worse (i.e., the prediction plot does not include the anomaly)
pred.dat.plotting <- pred.dat.plotting[!is.na(pred.dat.plotting$sp),]



### for those normal models that are not null (i.e., include normal climate predictors), include low and high hypothetical climate normal levels; for those that are, include only mid
pred.dat.plotting$norm.beg <- substr(pred.dat.plotting$mod,1,2)

pred.dat.plotting <- as.data.table(pred.dat.plotting)

pred.dat.plotting <- pred.dat.plotting[(norm.beg != "n0" & norm.level %in% c("low","high")) | (norm.beg == "n0" & norm.level == "mid")]


pred.dat.plotting[,c("pred.mid","pred.low","pred.high")] <- 100 * pred.dat.plotting[,c("pred.mid","pred.low","pred.high")] # convert to percentages

levels(pred.dat.plotting$norm.level)
levels(pred.dat.plotting$norm.level) <- c("High","Low","High and low")

## do not plot CWD or temp predictors
pred.dat.plotting = pred.dat.plotting[pred.dat.plotting$anom %in% c("Min precip","Mean precip","Min AET","Mean AET"),]



p <- ggplot(pred.dat.plotting,aes(x=diff.norm.ppt.z_c,y=pred.mid,color=norm.level,fill=norm.level)) +
  geom_line(aes(linetype=as.character(linestyle),size=linestyle)) +
  geom_ribbon(aes(ymin=pred.low,ymax=pred.high),alpha=0.3,color=NA) +
  facet_grid(anom~sp) +
  #geom_text(aes(0,0.8,label=anom.improvement),size=4,color="black") +
  theme_bw(12) +
  labs(x="Anomaly value",y="Model prediction (%)",color="Normal climate\n(precip. or AET)",fill="Normal climate\n(precip. or AET)") +
  scale_color_manual(values=c("High"="turquoise4","Low"="darkorange1","High and low" = "gray26")) +
  scale_fill_manual(values=c("High"="turquoise4","Low"="darkorange1","High and low" = "gray26")) +
  theme(panel.grid.minor = element_blank(),strip.background = element_blank(), panel.border = element_rect(colour = "black",size=0.6), strip.text = element_text(size = 9)) +
  theme(plot.margin = unit(c(-.1,0.5,0,0.5), "cm"),legend.key = element_rect(size = 0.5),legend.key.size = unit(1.2, 'lines')) +
  geom_vline(xintercept=0,linetype="longdash") +
  scale_linetype_manual(values=c(1,2),guide=FALSE) +
  scale_size(range=c(1.5,1),guide=FALSE)
  
#Plot it  
p



#### 9. Plot prediction fits for the 5 main responses (Fig. 3) ####

anom.var <- "Pmin"
pred.dat.comb <- merge(pred.dat,d.maes.anoms.short,all.x=TRUE,by.x=c("sp","anom"),by.y=c("sp","anom.name"))


## uncenter the predictor vars in the predictions data frame
pred.dat.comb <- as.data.frame(pred.dat.comb)
pred.dat.comb <- pred.dat.comb[(pred.dat.comb$anom.better & pred.dat.comb$type == "anom") |
                                         (!pred.dat.comb$anom.better & pred.dat.comb$type == "norm"),]



#for every column in predictions, if it ends in _c, uncenter it
for(col.name in names(pred.dat.comb)) {
  if(grepl("_c$",col.name)) { #it's a centered col
    col.name.uncentered <- substr(col.name,1,nchar(col.name)-2)
    
    col.mean <- d.center.dat[d.center.dat$var == col.name,"var.mean"]
    col.sd <- d.center.dat[d.center.dat$var == col.name,"var.sd"]
    
    pred.dat.comb[,as.character(col.name.uncentered)] <- pred.dat.comb[,col.name] * col.sd + col.mean
    
    
  }
}


pred.dat.comb <- data.table(pred.dat.comb)
pred.dat.comb <- pred.dat.comb[anom==anom.var]

pred.dat.comb <- pred.dat.comb[(pred.dat.comb$anom.better & pred.dat.comb$type == "anom") |
                                 (!pred.dat.comb$anom.better & pred.dat.comb$type == "norm"),]


pred.dat.comb[which(pred.dat.comb$anom.improvement < 0),c("pred.low","pred.high")] <- NA
pred.dat.comb$linestyle = 1 # solid, by default
pred.dat.comb[which(pred.dat.comb$anom.improvement < 0),"linestyle"] <- 2 # dashed if the anomaly model was worse (i.e., the prediction plot does not include the anomaly)
pred.dat.comb <- pred.dat.comb[!is.na(pred.dat.comb$sp),]


if(anom.var == "Pmin") {
  anom.var.c <- "diff.norm.ppt.min.z_c"
  anom.var.nc <- "diff.norm.ppt.min.z"
}

if(anom.var == "Pmean") {
  anom.var.c <- "diff.norm.ppt.z_c"
  anom.var.nc <- "diff.norm.ppt.z"
}

if(anom.var == "Amin") {
  anom.var.c <- "diff.norm.aet.min.z_c"
  anom.var.nc <- "diff.norm.aet.min.z"
}

if(anom.var == "Amean") {
  anom.var.c <- "diff.norm.aet.z_c"
  anom.var.nc <- "diff.norm.aet.z"
}


#for plotting a vertical line at the average anomaly value across all surveyed plots
anom.mid <- d.center.dat[d.center.dat$var == anom.var.c,"var.mean"]


### for those normal models that are not null, include low and high norm levels; for those that are, include only mid
pred.dat.comb$norm.beg <- substr(pred.dat.comb$mod,1,2)

pred.dat.comb <- pred.dat.comb[(norm.beg != "n0" & norm.level %in% c("low","high")) | (norm.beg == "n0" & norm.level == "mid")]

pred.dat.comb[,c("pred.mid","pred.low","pred.high")] <- 100 * pred.dat.comb[,c("pred.mid","pred.low","pred.high")] # convert to percentages

pred.dat.comb$anom.var <- pred.dat.comb[,anom.var.nc,with=FALSE]


levels(pred.dat.comb$norm.level)
levels(pred.dat.comb$norm.level) <- c("High","Low","High and low")


plot.cats <- c("presab",
               "cov")
plot.sps <- list(c("ABCO","PIPJ","HDWD.ALLSP"),
                 c("COV.SHRUB","COV.GRASS"))

p <- list()

for(i in 1:length(plot.cats)) {
  
  plot.sp <- plot.sps[[i]]

  # former pred.dat.plotting <- pred.dat.comb[type == "anom" & sp %in% plot.sp,]
  pred.dat.plotting <- pred.dat.comb[sp %in% plot.sp,]
  
  pred.dat.plotting$sp <- factor(pred.dat.plotting$sp,plot.sp)
  
  if(plot.cats[[i]] == "presab") {
    ylab <- "Percentage of plots\nwith regeneration"
    levels(pred.dat.plotting$sp) <- c("White fir","Yellow pine","Broadleaved trees")
  } else if(plot.cats[[i]] == "ht") {
    ylab <- "Percentage of plots where\ndominant in height"
  } else if(plot.cats[[i]] == "cov") {
    ylab <- "Percent cover"
    levels(pred.dat.plotting$sp) <- c("Shrubs","Graminoids","Forbs")
    

    
  }
  

  p[[i]] <- ggplot(pred.dat.plotting,aes(x=diff.norm.ppt.min.z,y=pred.mid,color=norm.level,fill=norm.level)) +
    geom_line(aes(linetype=as.character(linestyle),size=linestyle)) +
    geom_ribbon(aes(ymin=pred.low,ymax=pred.high),alpha=0.3,color=NA) +
    facet_wrap(~sp) +
    theme_bw(16) +
    labs(x="Anomaly: post-fire minimum precipitation (SD)",y=ylab,color="Normal\nprecipitation",fill="Normal\nprecipitation") +
    scale_color_manual(values=c("High"="turquoise4","Low"="darkorange1","High and low" = "gray26"),drop=FALSE) +
    scale_fill_manual(values=c("High"="turquoise4","Low"="darkorange1","High and low" = "gray26"),drop=FALSE) +
    theme(panel.grid.minor = element_blank(),strip.background = element_blank(), panel.border = element_rect(colour = "black",size=0.6), strip.text = element_text(size = 16,vjust=0)) +
    theme(plot.margin = unit(c(-.1,0.5,0,0.5), "cm"),legend.key.size = unit(1.5, 'lines')) +
    geom_vline(xintercept=anom.mid,linetype="longdash") +
    scale_linetype_manual(values=c(1,2),guide=FALSE) +
    scale_size(range=c(1.5,1.1), guide=FALSE)
  
  
  if(i == 1) {
    p[[i]] <- p[[i]] +
      #theme(plot.margin = unit(c(-.1,0.5,1,0.5), "cm")) +
      guides(fill=FALSE,color=FALSE)
  } else {
    p[[i]] <- p[[i]] +
    theme(plot.margin = unit(c(0.5,4,1,0), "cm")) +
    theme(legend.position=c(1.26,0.5)) +
    guides(color = guide_legend(override.aes = list(size=1.5)))
  }

  
}


a <- ggplot_gtable(ggplot_build(p[[1]]))
b <-  ggplot_gtable(ggplot_build(p[[2]]))

lay = rbind(c(1,1,1),
            c(NA,2,NA))

#Plot it  
grid.arrange(a,b,ncol=3,layout_matrix = lay,heights=c(1,1.1),widths=c(0.055,1,0.205))


## Compute what the normal precip hypothetical values are (for reporting in pub) ##

ppt.low <- low.val.fun("ppt.normal_c")
ppt.high <- high.val.fun("ppt.normal_c")

ppt.low <- (ppt.low * d.center.dat[d.center.dat$var=="ppt.normal_c","var.sd"]) + d.center.dat[d.center.dat$var=="ppt.normal_c","var.mean"]
ppt.high <- (ppt.high * d.center.dat[d.center.dat$var=="ppt.normal_c","var.sd"]) + d.center.dat[d.center.dat$var=="ppt.normal_c","var.mean"]




#### 10. Table of cross-validation errors (Table S1) ####
d.maes.anoms <- as.data.table(d.maes.anoms)
cv.errors <- d.maes.anoms[,.(Species=sp,Variable=anom.name,Baseline=best.anom.normal.mae,Anomaly=best.anom.mae)]
cv.err.melt <- melt(cv.errors,id.vars=c("Species","Variable"),measure.vars=c("Baseline","Anomaly"),variable.name="type")
cv.err.cast <- dcast(cv.err.melt,Species~Variable+type,value.var="value",fun=mean)
cv.err.cast[,2:13] <- round(cv.err.cast[,2:13],3)

response.rename <- list("PIPJ" = "Yellow pine presence",
                        "ABCO" = "White fir presence",

                        "HDWD.ALLSP" = "Broadleaved species presence",
                        
                        "COV.SHRUB" = "Shrub cover",
                        "COV.GRASS" = "Graminoid cover",
                        
                        "HT.PIPJ" = "Yellow pine height dominance",
                        "HT.ABCO" = "White fir height dominance",
                        "HT.HDWD.ALLSP" = "Broadleaved species height dominance"
                        )

cv.err.cast$Species <- gsubfn(pattern="\\S+",replacement=response.rename,x=cv.err.cast$Species)

cv.err.cast$Species <- factor(cv.err.cast$Species,levels=response.rename)
cv.err.cast <- cv.err.cast[!is.na(cv.err.cast$Species)]

cv.err.cast <- cv.err.cast[order(Species),]


colnames <- names(cv.err.cast)

colnames <- gsub("mean_Baseline","_Baseline",colnames)
colnames <- gsub("min_Baseline","_Baseline",colnames)
colnames <- gsub("max_Baseline","_Baseline",colnames)

colnames <- gsub("Amean","Mean AET",colnames)
colnames <- gsub("Amin","Minimum AET",colnames)
colnames <- gsub("Dmean","Mean CWD",colnames)
colnames <- gsub("Dmax","Maximum CWD",colnames)
colnames <- gsub("Pmean","Mean precpitation",colnames)
colnames <- gsub("Pmin","Minimum precipitation",colnames)

colnames <- gsub("P_Baseline","Precipitation baseline",colnames)
colnames <- gsub("A_Baseline","AET baseline",colnames)
colnames <- gsub("D_Baseline","CWD baseline",colnames)

colnames <- gsub("Anomaly","anomaly",colnames)
colnames <- gsub("_"," ",colnames)
colnames <- gsub("Species","Regeneration response variable",colnames)

names(cv.err.cast) <- colnames
cv.err.cast <- cv.err.cast[,which(!duplicated(names(cv.err.cast))),with=FALSE]


write.csv(cv.err.cast,"../tables/cv_err_V2F.csv")



#### 11. Table of model coefs (Tables S3-S5) ####

## For each anomaly type and model type (anom or norm), extract the coefs

coefs <- data.frame()

for(i in 1:nrow(d.maes.anoms)) {
  
  cv.row <- d.maes.anoms[i,]
  sp <- cv.row$sp
  norm.mod <- cv.row$best.anom.normal
  anom.mod <- cv.row$best.anomaly.mod
  
  m.norm <- fit.mods[[paste0(sp,"_",norm.mod)]]
  m.anom <- fit.mods[[paste0(sp,"_",anom.mod)]]
  
  if(sp %in% c(sp.opts,ht.opts)) { # glm regression (not betareg)
    
    anom.coefs <- as.data.frame(summary(m.anom)$coefficients[,1:2],optional=TRUE)
    anom.coefs$var.name <- row.names(anom.coefs)
    anom.coefs$type <- "Anomaly"
    anom.coefs$variable <- cv.row$anom.name
    
    norm.coefs <- as.data.frame(summary(m.norm)$coefficients[,1:2,drop=FALSE])
    norm.coefs$var.name <- row.names(norm.coefs)
    norm.coefs$type <- "Baseline"
    norm.coefs$variable <- cv.row$anom.name
    
    norm.coefs.prec <- NULL
    anom.coefs.prec <- NULL
    
  } else { #it's a cover resonse (beta regression)
    
    anom.coefs <- as.data.frame(summary(m.anom)$coefficients$mean[,1:2],optional=TRUE)
    anom.coefs$var.name <- row.names(anom.coefs)
    anom.coefs$type <- "Anomaly"
    anom.coefs$variable <- cv.row$anom.name
    
    anom.coefs.prec <- as.data.frame(summary(m.anom)$coefficients$precision[,1:2,drop=FALSE],optional=TRUE)
    anom.coefs.prec$var.name <- row.names(anom.coefs.prec)
    anom.coefs.prec$type <- "Anomaly"
    anom.coefs.prec$variable <- cv.row$anom.name
    
    
    norm.coefs <- as.data.frame(summary(m.norm)$coefficients$mean[,1:2,drop=FALSE])
    norm.coefs$var.name <- row.names(norm.coefs)
    norm.coefs$type <- "Baseline"
    norm.coefs$variable <- cv.row$anom.name
    
    norm.coefs.prec <- as.data.frame(summary(m.norm)$coefficients$precision[,1:2,drop=FALSE],optional=TRUE)
    norm.coefs.prec$var.name <- row.names(anom.coefs.prec)
    norm.coefs.prec$type <- "Baseline"
    norm.coefs.prec$variable <- cv.row$anom.name
    

  }
  
  coefs.row <- rbind(norm.coefs,norm.coefs.prec,anom.coefs,anom.coefs.prec)
  coefs.row$sp <- sp
  coefs <- rbind(coefs,coefs.row)
  
  
}

coefs <- as.data.table(coefs)

coefs[,c("Estimate","Std. Error")] <- coefs[,lapply(.SD,round,digits=3),.SDcols=c("Estimate","Std. Error")]

coefs$var.name <- gsub(pattern="aet\\.|ppt\\.|def\\.",replacement="",x=coefs$var.name)
coefs$var.name <- gsub(pattern="min\\.|max\\.",replacement="",x=coefs$var.name)
coefs$est.se <- paste0(coefs$Estimate," (",coefs$`Std. Error`,")")

coefs.cast <- dcast(coefs,sp+variable+type~var.name,value.var="est.se")
coefs.cast <- as.data.table(coefs.cast)

dummy.df <- data.table("sp"="a","variable"="a","type"="a","(Intercept)"=1,"normal_c"=1,"normal_c.sq"=1,"seed_tree_distance_general_c"=1,"rad.march_c"=1,"diff.norm.z_c"=1,"diff.norm.z_c.sq"=1,"normal_c:diff.norm.z_c"=1,"normal_c:diff.norm.z_c.sq"=1,"diff.norm.z_c:normal_c.sq"=1,"(phi)"=1)
coefs.cast <- rbind.fill(coefs.cast,dummy.df)
coefs.cast <- as.data.table(coefs.cast)

coefs.cast <- coefs.cast[,c("sp","variable","type",Intercept="(Intercept)","normal_c","normal_c.sq",seed.tree="seed_tree_distance_general_c","rad.march_c","diff.norm.z_c","diff.norm.z_c.sq","normal_c:diff.norm.z_c","normal_c:diff.norm.z_c.sq","diff.norm.z_c:normal_c.sq","(phi)"),with=FALSE]


responses.keep <- c("ABCO","PIPJ","HDWD.ALLSP","COV.SHRUB","COV.GRASS","HT.PIPJ","HT.ABCO","HT.HDWD.ALLSP")
anoms.keep <- c("Amin","Amean","Pmin","Pmean")

coefs.cast <- coefs.cast[coefs.cast$sp %in% responses.keep & coefs.cast$variable %in% anoms.keep,]


anomsub <- list("Pmean"="Mean Precipitation","Pmin" = "Minimum Precipitation","Amean" = "Mean AET","Amin"="Minimum AET")
coefs.cast$variable <- gsubfn("\\S+",anomsub,coefs.cast$variable)

colsub <- list("sp"="Response","variable" = "Anomaly","type"="Model","(Intercept)" = "Intercept","normal_c"="Normal climate","normal_c.sq"="Normal climate^2","seed_tree_distance_general_c"="Seed tree dist.","rad.march_c"="Solar exposure","diff.norm.z_c"="Anomaly","diff.norm.z_c.sq"="Anomaly^2","normal_c:diff.norm.z_c" = "Normal climate * Anomaly","normal_c:diff.norm.z_c.sq" = "Normal climate * Anomaly^2","diff.norm.z_c:normal_c.sq" = "Normal climate^2 * Anomaly","(phi)" = "Phi")
names(coefs.cast) <- gsubfn("\\S+",colsub,names(coefs.cast))

coefs.cast[is.na(coefs.cast)] <- ""

#change the first word of a two-word string to "Mean/Min"
fixfront <- function(x) {
  second <- strsplit(x," ",fixed=TRUE)[[1]][2]
  return(second)
}

coefs.cast[Model=="Baseline","Anomaly"] <- apply(coefs.cast[Model=="Baseline","Anomaly"],FUN=fixfront,MARGIN=1)
coefs.cast$Model <- paste(coefs.cast$Anomaly,coefs.cast$Model,sep=" ")


modelsub <- c("Precipitation Baseline" = "Precipitation baseline","Mean Precipitation Anomaly"="Mean precipitation anomaly","Minimum Precipitation Anomaly"="Minimum precipitation anomaly","AET Baseline"="AET baseline","Mean AET Anomaly"="Mean AET anomaly","Minimum AET Anomaly"="Minimum AET anomaly")
coefs.cast$Model <- mapvalues(x=coefs.cast$Model,from=names(modelsub),to=modelsub)

coefs.cast$Model <- as.factor(coefs.cast$Model)
coefs.cast$Model <- factor(coefs.cast$Model,levels=modelsub)

coefs.cast <- coefs.cast[!duplicated(coefs.cast),]
coefs.cast <- coefs.cast[,-2] #remove "anomaly" column


### Presence/absence
keep.vars <- c("PIPJ","ABCO","HDWD.ALLSP")
coefs.presab <- coefs.cast[coefs.cast$Response %in% keep.vars,]
spsub <- list("PIPJ" = "Yellow pine","ABCO" = "White fir","HDWD.ALLSP" = "Broadleaved trees")
coefs.presab$Response <- gsubfn("\\S+",spsub,coefs.presab$Response)
coefs.presab$Response <- factor(coefs.presab$Response,levels=spsub)
coefs.presab <- coefs.presab[,-"Phi"]
coefs.presab <- coefs.presab[order(Response,Model)]

coefs.presab[coefs.presab == ""] <- "--"

write.csv(coefs.presab,"../tables/coefs_presab_V2F.csv")


### Height dominaince
keep.vars <- c("HT.PIPJ","HT.ABCO","HT.HDWD.ALLSP")
coefs.dom <- coefs.cast[coefs.cast$Response %in% keep.vars,]
spsub <- list("HT.PIPJ" = "Yellow pine","HT.ABCO" = "White fir","HT.PINUS.ALLSP" = "Pines","HT.SHADE.ALLSP" = "Shade-tolerant conifers","HT.HDWD.ALLSP" = "Broadleaved trees")
coefs.dom$Response <- gsubfn("\\S+",spsub,coefs.dom$Response)
coefs.dom$Response <- factor(coefs.dom$Response,levels=spsub)
coefs.dom <- coefs.dom[,-c("Phi","Seed tree dist.")]
coefs.dom <- coefs.dom[order(Response,Model)]

coefs.dom[coefs.dom == ""] <- "--"

write.csv(coefs.dom,"../tables/coefs_dom_V2F.csv")

### Cover
keep.vars <- c("COV.SHRUB","COV.GRASS")
coefs.cov <- coefs.cast[coefs.cast$Response %in% keep.vars,]
spsub <- list("COV.SHRUB"="Shrubs","COV.GRASS"="Graminoids")
coefs.cov$Response <- gsubfn("\\S+",spsub,coefs.cov$Response)
coefs.cov$Response <- factor(coefs.cov$Response,levels=spsub)
coefs.cov <- coefs.cov[,-c("Seed tree dist.")]
coefs.cov <- coefs.cov[order(Response,Model)]

coefs.cov[coefs.cov == ""] <- "--"

write.csv(coefs.cov,"../tables/coefs_cov_V2F.csv")




#### 12. For each species and rad group, plot fitted vs. observed, for normal and anom side by side (Figs. S4-S5) ####

fit.dat.plotting <- fit.dat

fit.dat.plotting$resid <- fit.dat.plotting$response.var-fit.dat.plotting$fitted

## PPT min only
fit.dat.ppt <- as.data.table(fit.dat.plotting[fit.dat$anom == "Pmin",])

## summarize by fire, sp

fit.dat.ppt.fire <- fit.dat.ppt[,list(fitted=mean(fitted),observed=mean(response.var),anom.var=mean(diff.norm.ppt.min.z_c),resid=mean(resid)),by=.(Fire,type,sp)]

fit.dat.ppt.fire[,c("observed","fitted")] <- 100 * fit.dat.ppt.fire[,c("observed","fitted")]


resps.keep <- c("PIPJ","ABCO","HDWD.ALLSP","COV.SHRUB","COV.GRASS")

fit.dat.ppt.fire <- fit.dat.ppt.fire[fit.dat.ppt.fire$sp %in% resps.keep,]

fit.dat.ppt.fire$sp <- as.factor(fit.dat.ppt.fire$sp)
fit.dat.ppt.fire$sp <- factor(fit.dat.ppt.fire$sp,resps.keep)

levels(fit.dat.ppt.fire$sp)
levels(fit.dat.ppt.fire$sp) <- c("Yellow pine regeneration\n(% of plots)","White fir regeneration\n(% of plots)","Broadleaf tree\nregeneration\n(% of plots)","Shrubs\n(% cover)","Graminoids\n(% cover)")

fit.dat.ppt.fire$type <- as.factor(fit.dat.ppt.fire$type)
levels(fit.dat.ppt.fire$type)
fit.dat.ppt.fire$type <- factor(fit.dat.ppt.fire$type,c("norm","anom"))
levels(fit.dat.ppt.fire$type)
levels(fit.dat.ppt.fire$type) <- c("Baseline","Post-fire\nanomaly")

dummy.df <- fit.dat.ppt.fire[,list(observed=c(0,max(observed,fitted)),fitted=c(0,max(observed,fitted))),by=list(sp,type,sp)]


p <- ggplot(fit.dat.ppt.fire,aes(x=observed,y=fitted,color=type)) +
  geom_point(size=2) +
  geom_abline(slope=1,intercept=0,size=1,color="gray") +
  geom_point(data=dummy.df,alpha=0) +
  facet_wrap(~sp,scales="free") +
  labs(x="Observed value",y="Fitted value",color="Model type") +
  theme_bw(18) +
  theme(plot.title=element_text(hjust=0.5)) +
  scale_color_manual(values=c("turquoise4","darkorange1")) +
  theme(panel.grid.minor = element_blank(),strip.background = element_blank(), panel.border = element_rect(colour = "black",size=0.6), strip.text = element_text(size = 16,vjust=0)) +
  theme(plot.margin = unit(c(-.1,0.5,0.1,0.5), "cm")) +
  theme(legend.position=c(.83,.25),legend.key.size = unit(2, 'lines'))

p


###
### Repeat for min AET
###

fit.dat.plotting <- fit.dat

fit.dat.plotting$resid <- fit.dat.plotting$response.var-fit.dat.plotting$fitted


## min anom only
fit.dat.ppt <- as.data.table(fit.dat.plotting[fit.dat$anom == "Amin",])

## summarize by fire, sp

fit.dat.ppt.fire <- fit.dat.ppt[,list(fitted=mean(fitted),observed=mean(response.var),anom.var=mean(diff.norm.ppt.min.z_c),resid=mean(resid)),by=.(Fire,type,sp)]
# fit.dat.ppt.fire <- fit.dat.ppt[,list(fitted=fitted,observed=response.var,anom.var=diff.norm.ppt.min.z_c,Fire,type,sp)]

fit.dat.ppt.fire[,c("observed","fitted")] <- 100 * fit.dat.ppt.fire[,c("observed","fitted")]


resps.keep <- c("PIPJ","ABCO","HDWD.ALLSP","COV.SHRUB","COV.GRASS")

fit.dat.ppt.fire <- fit.dat.ppt.fire[fit.dat.ppt.fire$sp %in% resps.keep,]

fit.dat.ppt.fire$sp <- as.factor(fit.dat.ppt.fire$sp)
fit.dat.ppt.fire$sp <- factor(fit.dat.ppt.fire$sp,resps.keep)

levels(fit.dat.ppt.fire$sp)
levels(fit.dat.ppt.fire$sp) <- c("Yellow pine regeneration\n(% of plots)","White fir regeneration\n(% of plots)","Broadleaf tree\nregeneration\n(% of plots)","Shrubs\n(% cover)","Graminoids\n(% cover)")

fit.dat.ppt.fire$type <- as.factor(fit.dat.ppt.fire$type)
levels(fit.dat.ppt.fire$type)
fit.dat.ppt.fire$type <- factor(fit.dat.ppt.fire$type,c("norm","anom"))
levels(fit.dat.ppt.fire$type)
levels(fit.dat.ppt.fire$type) <- c("Baseline","Post-fire\nanomaly")

dummy.df <- fit.dat.ppt.fire[,list(observed=c(0,max(observed,fitted)),fitted=c(0,max(observed,fitted))),by=list(sp,type,sp)]


p <- ggplot(fit.dat.ppt.fire,aes(x=observed,y=fitted,color=type)) +
  geom_point(size=2) +
  geom_abline(slope=1,intercept=0,size=1,color="gray") +
  geom_point(data=dummy.df,alpha=0) +
  facet_wrap(~sp,scales="free") +
  labs(x="Observed value",y="Fitted value",color="Model type") +
  theme_bw(18) +
  theme(plot.title=element_text(hjust=0.5)) +
  scale_color_manual(values=c("turquoise4","darkorange1")) +
  theme(panel.grid.minor = element_blank(),strip.background = element_blank(), panel.border = element_rect(colour = "black",size=0.6), strip.text = element_text(size = 16,vjust=0)) +
  theme(plot.margin = unit(c(-.1,0.5,0.1,0.5), "cm")) +
  theme(legend.position=c(.83,.25),legend.key.size = unit(2, 'lines'))


#Plot it  
p





#### 13. Topoclimate category-level multivariate analysis (CCA) ####

focal.sp <- c("PIPJ","ABCO","PILA","PSME")

focal.cols <- c("Fire","topoclim.cat","species","regen.count.old","regen.count.all","adult.ba")
d.sp.simp <- d.sp.2[d.sp.2$species %in% focal.sp,focal.cols]
names(d.sp.simp)[4:6] <- c("r.old","r.all","a.ba")
d.sp.simp <- as.data.table(d.sp.simp)
d.sp.cast <- dcast(d.sp.simp,topoclim.cat + Fire~species,value.var=c("r.old","r.all","a.ba"))

regen.var <- "r.old"

d.all <- merge(d.plot.3,d.sp.cast,by=c("Fire","topoclim.cat"))
d.all <- d.all[(d.all$count.highsev > 4) & (d.all$count.control > 4),]

sp.cols <- grep(regen.var,names(d.all))
d.all.sp <- d.all[,sp.cols]
d.all.sp.tot <- rowSums(d.all.sp)
keep.rows <- d.all.sp.tot > 0
d.all.sp <- d.all.sp[keep.rows,]
d.all <- d.all[keep.rows,]


### remove the prefixes from the regen and adult species names

adult.cols <- grep("^a[.]ba",names(d.all))
adult.colnames.simp <- substr(names(d.all)[adult.cols],6,100)
names(d.all)[adult.cols] <- adult.colnames.simp

regen.cols <- grep("^r[.]",names(d.all.sp))
regen.colnames.simp <- substr(names(d.all.sp)[regen.cols],7,100)
names(d.all.sp)[regen.cols] <- regen.colnames.simp

### change species names to common name abbreviations

sublist <- c("PILA" = "SP",
             "ABCO" = "WF",
             "PSME" = "DF",
             "PIPJ" = "YP",
             "rad.march.highsev" = "solar.rad",
             "tmean.normal.highsev" = "temp.norm",
             "ppt.normal.highsev" = "ppt.norm",
             "diff.norm.ppt.min.z.highsev" = "ppt.anom",
             "diff.norm.tmean.max.z.highsev" = "temp.anom")

d.names <- names(d.all)
d.names.new <- match(d.names,names(sublist))
d.names.match <- !is.na(d.names.new)
d.names[d.names.match] <- sublist[d.names.new][d.names.match]
names(d.all) <- d.names

d.names <- names(d.all.sp)
d.names.new <- match(d.names,names(sublist))
d.names.match <- !is.na(d.names.new)
d.names[d.names.match] <- sublist[d.names.new][d.names.match]
names(d.all.sp) <- d.names



# Better predictor variable names
d.all$`Solar Radiation` <- d.all$solar.rad
d.all$`Normal Precip` <- d.all$ppt.norm
d.all$`Normal Temp` <- d.all$temp.norm

d.all$`Precip Anomaly` <- d.all$ppt.anom
d.all$TempAnom <- d.all$temp.anom


## Test different predictor variable sets

cc.full <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Solar Radiation` + `Precip Anomaly` + WF+DF+SP+YP,data=d.all)
cc.np <- cca(d.all.sp ~ `Normal Temp` + `Solar Radiation` + `Precip Anomaly` + WF+DF+SP+YP,data=d.all)
cc.nt <- cca(d.all.sp ~ `Normal Precip` +`Solar Radiation` + `Precip Anomaly` + WF+DF+SP+YP,data=d.all)
cc.sr <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Precip Anomaly` + WF+DF+SP+YP,data=d.all)
cc.pa <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Solar Radiation` + WF+DF+SP+YP,data=d.all)
cc.wf <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Solar Radiation` + `Precip Anomaly` + DF+SP+YP,data=d.all)
cc.df <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Solar Radiation` + `Precip Anomaly` + WF+SP+YP,data=d.all)
cc.sp <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Solar Radiation` + `Precip Anomaly` + WF+DF+YP,data=d.all)
cc.yp <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Solar Radiation` + `Precip Anomaly` + WF+DF+SP,data=d.all)
cc.nosp <- cca(d.all.sp ~ `Normal Precip` + `Normal Temp` + `Solar Radiation` + `Precip Anomaly`,data=d.all)

## Check proportion of inertia explained
cc.full
cc.np
cc.nt
cc.sr
cc.pa
cc.wf
cc.df
cc.sp
cc.yp
cc.nosp



#### 14. Plot the CCA (Fig. 4) ####
ccadata <- fortify(cc.full)

ccadata$Label <- gsub("`","",ccadata$Label)

ccsites <- ccadata[ccadata$Score %in% "sites",]
ccsp <- ccadata[ccadata$Score %in% "species",]
ccvects <- ccadata[ccadata$Score %in% "biplot",]


ccvects[,3:4] <- 3.8 * ccvects[,3:4]

ccvects.1 <- ccvects
ccvects.1[,3:4] <- 1.2 * ccvects[,3:4]

p <- ggplot(ccsites,aes(x=CCA1,y=CCA2)) +
  geom_point(color="darkorange1") +
  geom_hline(yintercept=0,color="grey") +
  geom_vline(xintercept=0,color="grey") +
  theme_bw(14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_segment(data=ccvects,aes(x=0,xend=CCA1,y=0,yend=CCA2), arrow = arrow(length = unit(0.4, "cm")),size=0.7,colour="grey") +
  geom_text(data=ccsp,aes(label=Label),color="turquoise4",size=5.5,fontface=2) +
  #lims(x=c(-4,4),y=c(-4.1,4.1)) +
  geom_text(data=ccvects.1,aes(x=CCA1,y=CCA2,label=Label),size=5) +
  labs(x="Axis 1",y="Axis 2") +
  scale_x_continuous(limits=c(-1.7,3.2))

#Plot it  
p
