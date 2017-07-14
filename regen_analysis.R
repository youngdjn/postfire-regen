setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(party)
library(ggplot2)
library(brms)
library(pROC)
library(betareg)


source("regen_analysis_functions.R")

#### 0. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)



## Look for plots with incomplete data specified in the comments
plots.exceptions <- grepl("#.*(INCOMPLETE|INCORRECT)",d.plot$NOTES)
d.plot <- d.plot[!plots.exceptions,]


# only keep the necessary columns
d.plot <- d.plot[,c("Regen_Plot","Fire","Year.of.Fire","Easting","Northing","aspect","slope","SHRUB","FORB","GRASS","HARDWOOD","CONIFER","FIRE_SEV","BA.Live1","Year","firesev","dist.to.low","fire.abbr","X5yr","fire.year","survey.years.post","elev.m","rad.march","tmean.post","ppt.post","ppt.post.min","tmean.normal","ppt.normal","seed.tree.any","diff.norm.ppt.z","diff.norm.ppt.min.z","seed_tree_distance_general","seed_tree_distance_conifer","seed_tree_distance_hardwood","diff.norm.ppt.z","diff.norm.ppt.min.z","tmean.post","ppt.post","ppt.post.min","perc.norm.ppt","perc.norm.ppt.min","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z","def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","diff.norm.def.max.z","diff.norm.aet.min.z","def.post","aet.post","def.post.max","aet.post.min","snow.post.min","snow.normal","snow.post","diff.norm.snow.z","diff.norm.snow.min.z","dominant_shrub_ht_cm","dom.veg.all")]

# only Sierra Nevada fires #! removed DEEP
sierra.fires <- c("STRAYLOR","CUB","RICH","MOONLIGHT","ANTELOPE","BTU LIGHTENING","HARDING","BASSETTS","PENDOLA","AMERICAN RIVER","RALSTON","FREDS","SHOWERS","POWER","BAGLEY","PEAK","CHIPS")
d.plot <- d.plot[d.plot$Fire %in% sierra.fires,]

d.plot <- d.plot[d.plot$Regen_Plot != "SHR0900015",]

## Remove managed plots, plots in nonforested locations (e.g. exposed bedrock), etc.
plots.exclude <- read.csv("data_intermediate/plots_exclude.csv",header=T,stringsAsFactors=FALSE)
plot.ids.exclude <- plots.exclude[plots.exclude$Exclude != "",]$Regen_Plot
d.plot <- d.plot[!(d.plot$Regen_Plot %in% plot.ids.exclude),]




# ## Removed this code: because most of the time NA means seed source not visible, and excluding all plots >50 m from seed source
# # if no data on seed tree distance (or it was recorded as being at/beyond the limit of the laser) use remote-sensing-based value
# d.plot$seed.tree.any.comb <- ifelse(is.na(d.plot$seed_tree_distance_general) | (d.plot$seed_tree_distance_general >= 150),d.plot$dist.to.low,d.plot$seed_tree_distance_general)

## Replaced with this:
# Remove any plots > 50m from seed source
d.plot <- d.plot[which(d.plot$seed_tree_distance_general < 75),]


# quadratic climate terms
d.plot$ppt.normal.sq <- d.plot$ppt.normal^2
d.plot$tmean.normal.sq <- d.plot$tmean.normal^2
d.plot$snow.normal.sq <- d.plot$snow.normal^2

# variable for regen presence/absence
d.sp$regen.presab.young <- ifelse(d.sp$regen.count.young > 0,TRUE,FALSE)
d.sp$regen.presab.old <- ifelse(d.sp$regen.count.old > 0,TRUE,FALSE)
d.sp$regen.presab.all <- ifelse(d.sp$regen.count.all > 0,TRUE,FALSE)

# #! TEMPORARY: if there was no radiation data, set it equal to 0
# d.plot$rad.march <- ifelse(is.na(d.plot$rad.march),0,d.plot$rad.march)

high.sev <- c(4,5) # which field-assessed severity values are considered high severity
control <- c(0,1) # which field-assessed severity values are considered controls




#### 1. Assign each plot a topoclimatic category ####

#! NOTE that when breaking down plots by factorial combinations topoclimatic variables,
# it might be important in the future to consider the interaction of variables.
# E.g., not all levels of a given variable may be available at all levels of another given variables
# E.g. at high precipitation, maybe there is only north aspects available. Doesn't seem to be the case here, but could potentially be with other variables.


## get what average radiation is across all the fires
avg.rad <- mean(d.plot$rad.march)


fires <- unique(d.plot$Fire)
d.plot$precip.category
d.plot$rad.category #radiation

for(fire in fires) {
  
  ## Precipitation categories
  # determine what the precipitation breakpoints should be (here just using median) -- based on high severity plots only
  breaks <- quantile(d.plot[(d.plot$Fire == fire) & (d.plot$FIRE_SEV %in% high.sev),]$ppt.normal,probs=c(0.5),na.rm=TRUE)
  
  
  # for some fires with a small range of precip, override the precip breaks, so we just have one category per fire
  fires.small.precip.range <- c("AMERICAN RIVER","ANTELOPE","BAGLEY","BASSETTS","BTU LIGHTENING","HARDING","STRAYLOR")
  if(fire %in% fires.small.precip.range) {
    breaks <- 0
  }
  
  # categorize plots based on where they fall between the breakpoints  
  categories <- categorize(d.plot[d.plot$Fire==fire,]$ppt.normal,breaks,name="P")
  # store it into the plot data.frame
  d.plot[d.plot$Fire==fire,"precip.category"] <- categories
  
  ## Radiation categories
  # determine what the breakpoints should be (here just using median) -- based on high severity plots only
  breaks <- quantile(d.plot[(d.plot$Fire == fire) & (d.plot$FIRE_SEV %in% high.sev),]$rad.march,probs=c(0.5),na.rm=TRUE)
  # categorize plots based on where they fall between the breakpoints
  
  #override the per-fire breaks
  breaks <- avg.rad
  
  categories <- categorize(d.plot[d.plot$Fire==fire,]$rad.march,breaks,name="R")
  # store it into the plot data.frame
  d.plot[d.plot$Fire==fire,"rad.category"] <- categories
  
  #! To-do: break down the categories further for the fires that have enough plots (may have different number of categories per fire). Goal of ~10 plots per category? Seed tree distance should probably be a category (or else simply exclude plots that are far from seed source)

}
  

## Create one variable reflecting the all-way factorial combination of topoclimatic categories
d.plot$topoclim.cat <- paste(d.plot$precip.category,d.plot$rad.category,sep="_")
# 
# Removed now that only have one precip category
# ## Make an exception for Harding: only two categories
# d.plot[d.plot$Fire == "HARDING","topoclim.cat"] <- ifelse(d.plot[d.plot$Fire=="HARDING","rad.march"] > 6250,"P.1_R.2","P.1_R.1")




  

### Plot relevant "topoclimate space" for each fire and see how the categories broke them down
## note that Cub and straylor do not have radiation (yet)

## remove PEAK and PENDOLA
d.plot <- d.plot[!(d.plot$Fire %in% c("PEAK","PENDOLA","SHOWERS")),]

library(car)
d.plot$FIRE_SEV.cat <- recode(d.plot$FIRE_SEV,"control='control';high.sev='high.sev';else=NA")

d.plot.cat <- d.plot

# plot high sev and control
ggplot(d.plot.precat[!is.na(d.plot.precat$FIRE_SEV.cat),],aes(x=ppt.normal,y=rad.march,col=FIRE_SEV.cat)) +
  geom_point(size=3) +
  facet_wrap(~Fire,scales="fixed") +
  geom_hline(yintercept=avg.rad) +
  theme_bw(16)

# look at control only
ggplot(d.plot.precat[d.plot.precat$FIRE_SEV %in% control,],aes(x=ppt.normal,y=rad.march,col=topoclim.cat)) +
  geom_point() +
  facet_wrap(~Fire,scales="free") +
  theme_bw(16)

# look at high sev only
ggplot(d.plot.precat[d.plot.precat$FIRE_SEV %in% high.sev,],aes(x=ppt.normal,y=rad.march,col=topoclim.cat)) +
  geom_point() +
  facet_wrap(~Fire,scales="free") +
  theme_bw(16)



## Remove climatic regions that do not have comparable controls and high sev
## NOTE: this is only necessary and recommended when doing an analysis that involves comparing control adults with high sev regen

d.plot.precat <- d.plot
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "BTU LIGHTENING") & (d.plot.precat$ppt.normal > 2000)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "BAGLEY") & (d.plot.precat$rad.march < 4000)),]
d.plot.precat <- d.plot.precat[!((d.plot.precat$Fire == "CUB") & (d.plot.precat$ppt.normal < 1300)),]

d.plot.domveg <- d.plot.precat # save this data frame at this stage for extracting dominant vegetation data (doesn't matter what severity or years post-fire)

# remove plots that are high severity but not surveyed in years 4-5 post-fire
d.plot.precat <- d.plot.precat[!(d.plot.precat$FIRE_SEV.cat == "high.sev" & !(d.plot.precat$survey.years.post %in% c(4,5))),]



write.csv(d.plot.precat,"../../all_plots.csv",row.names=FALSE)



#### 2. Summarize (compute average) regen values (high-sev plots only) and adults (control plots only) by species across all plots in each topoclimatic category in each fire ####

## assign the trees by species their topoclimatic category and fire name. This also ensures that we only are looking at seedlings for whose plots we are interested (because with this merge operation, seedlings from plots not in d.plot will be dropped)
d.sp.cat <- merge(d.sp,d.plot.precat[,c("Regen_Plot","topoclim.cat","Fire","FIRE_SEV","survey.years.post")])

## preparing to aggregate tree data: get highsev and control plots only, each with only the columns relevant to it
d.sp.cat.highsev <- d.sp.cat[d.sp.cat$FIRE_SEV %in% high.sev,c("Fire","species","topoclim.cat","seed.tree.sp","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all","survey.years.post")]
d.sp.cat.control <- d.sp.cat[d.sp.cat$FIRE_SEV %in% control,c("Fire","species","topoclim.cat","adult.count","adult.ba","survey.years.post")]

## for highsev plots (the ones where we're interested in regen), only consider plots surveyed 4-5 years post-fire
d.sp.cat.highsev <- d.sp.cat.highsev[d.sp.cat.highsev$survey.years.post %in% c(4,5),]

## aggregate tree data by species and topo category #! might want to also calculate SD or SE in a separate aggregate call? (would append the variable names with ".sd")
d.sp.agg.highsev <- aggregate(d.sp.cat.highsev[,c(-1,-2,-3)],by=list(d.sp.cat.highsev$species,d.sp.cat.highsev$topoclim.cat,d.sp.cat.highsev$Fire),FUN=mean,na.rm=TRUE)
names(d.sp.agg.highsev)[1:3] <- c("species","topoclim.cat","Fire")
d.sp.agg.control <- aggregate(d.sp.cat.control[,c(-1,-2,-3)],by=list(d.sp.cat.control$species,d.sp.cat.control$topoclim.cat,d.sp.cat.control$Fire),FUN=mean,na.rm=TRUE)
names(d.sp.agg.control)[1:3] <- c("species","topoclim.cat","Fire")

## merge the control (adults only) and highsev (seedlings only) tree data
d.sp.agg <- merge(d.sp.agg.highsev,d.sp.agg.control,all.x=TRUE,by=c("species","topoclim.cat","Fire"))

##preparing to aggregate plot (e.g. climate) data: label plots as highsev or control
# first remove the variables that are not useful
d.plot.c <- remove.vars(d.plot.precat,c("Year.of.Fire","Easting","Northing","aspect","Year","precip.category","rad.category"))
# label plots as control or high sev
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% control,"control",NA)
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% high.sev,"highsev",d.plot.c$type)
# get rid of plots that are neither control nor high sev
d.plot.c <- d.plot.c[!is.na(d.plot.c$type),]

## aggregate plots by Fire, topoclim category, and type (control or high sev)
d.plot.agg.mean <- aggregate(remove.vars(d.plot.c,c("Regen_Plot","topoclim.cat","type","fire.abbr","X5yr","Fire","FIRE_SEV.cat")),by=list(d.plot.c$Fire,d.plot.c$topoclim.cat,d.plot.c$type),FUN=mean,na.rm=TRUE)
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


#### 3. Steps required prior to any analysis ####

# Remove the topoclimatic categories with too few plots in either burned or control
d.plot.3 <- d.plot.2[which((d.plot.2$count.control > 4) & (d.plot.2$count.highsev > 0)),]

# Compute additional variables
d.sp.2$proportion.young <- d.sp.2$regen.count.young / d.sp.2$regen.count.all



### Compute dominant vegetation for each remaining topoclimate category, based on stated observed nearby dominants. use all plots regardless of severity

non.tree.dom.veg <- c("ARPA6","CEIN3","CECO","CEPR") # species to exclude from list


for(fire in unique(d.plot.3$Fire)) {
  
  d.plot.fire <- d.plot.3[d.plot.3$Fire == fire,]
  
  for(topoclim.cat in unique(d.plot.fire$topoclim.cat)) {
    
    
    d.plot.fire.cat <- d.plot.domveg[(d.plot.domveg$Fire == fire) & (d.plot.domveg$topoclim.cat == topoclim.cat),]
    
    if (nrow(d.plot.fire.cat) < 5) {
      cat("Less than 5 plots in ",fire," ",topoclim.cat,". Skipping.\n")
      next()
    }
    
    
    fire.cat.sp <- paste(d.plot.fire.cat$dom.veg.all,collapse = " ")
    
    # use white space to split string into vector
    sps <- scan(text=fire.cat.sp,what="")
    sps <- sps[!is.na(sps)]
    sps <- sps[!sps %in% non.tree.dom.veg]
    sps[sps %in% c("PIJE","PIPO")] <- "PIPJ"
    sps.tab <- table(sps)
    sps.tab <- sort(sps.tab,decreasing=TRUE)
    sps.tab <- sps.tab/sum(sps.tab)
    
    #print(sps.tab)
    
    # keep all species that got > 10% of the mentions
    sps.tab <- sps.tab[sps.tab > .2]
    
    # keep the top 4 species
    sps.tab <- sps.tab[1:4]
    
    # remove NAs (generated if list of species was < 4 prior to last line)
    sps.tab <- sps.tab[!is.na(sps.tab)]
    
    # get names
    dom.sp <- names(sps.tab)
    
    dom.sp <- paste(dom.sp,collapse=", ")
    
    
    d.plot.3[(d.plot.3$Fire == fire) & (d.plot.3$topoclim.cat == topoclim.cat),"dom.tree.sp.obs"] <- dom.sp
    
  }
  
}








### Plot histogram of seed tree distance ####

d.plot.highsev <- d.plot[d.plot$FIRE_SEV %in% c(4,5),]

ggplot(d.plot.highsev,aes(seed_tree_distance_general)) +
  geom_histogram() +
  facet_wrap(~Fire)

library(plyr)

a <- ddply(d.plot.highsev,~Fire,summarise,mean=mean(seed_tree_distance_general),sd=sd(seed_tree_distance_general))

ggplot(a,aes(x=Fire,y=mean)) +
  geom_point() +
  geom_errorbar(ymax=a$mean+a$sd,ymin=a$mean-a$sd) +
  scale_y_continuous(limits=c(0,50))
  



#### 4. Determine dominant adult tree species (from control plots) within each category ####
library(reshape)

d.sp.pre <- d.sp.2[,c("species","topoclim.cat","Fire","adult.ba")]

d.fire.sp <- cast(d.sp.pre,Fire + topoclim.cat ~ species,value="adult.ba")

d.fire.sp$PIPJ <- d.fire.sp$PIPO + d.fire.sp$PIJE

sp.names <- names(d.fire.sp)
names.drop <- c("ALL","ALNUS","CONIF.ALLSP","CONIFER","HDWD.ALLSP","JUNIPERUS","PINUS.ALLSP","SHADE.ALLSP","ABIES","PIPO","PIJE")
sp.names <- sp.names[!(sp.names %in% names.drop)]
d.fire.sp <- d.fire.sp[,sp.names]
d.fire.sp[,4:length(names(d.fire.sp))] <- d.fire.sp[,4:length(names(d.fire.sp))] / rowSums(d.fire.sp[,4:length(names(d.fire.sp))])

d.fire.sp <- d.fire.sp[complete.cases(d.fire.sp),]

## for each one, get a list of all species with > 20% BA

d.fire.sp$dom.tree.sp.ba <- ""

for (i in 1:nrow(d.fire.sp)) {
  
  row.ba <- d.fire.sp[i,4:ncol(d.fire.sp)]
  
  row.ba <- sort(row.ba,decreasing=TRUE)
  
  row.ba <- unlist(row.ba)
  names.store <- names(row.ba)
  row.ba <- as.numeric(row.ba)
  names(row.ba) <- names.store
  
  row.ba <- row.ba[row.ba > .2]
  
  row.ba <- row.ba[1:4]
  row.ba <- row.ba[!is.na(row.ba)]
  
  dom.sp.row <- paste(names(row.ba),collapse=", ")
  
  d.fire.sp[i,"dom.tree.sp.ba"] <- dom.sp.row
}

d.fire.sp <- d.fire.sp[,c("Fire","topoclim.cat","dom.tree.sp.ba")]

d.plot.3 <- merge(d.plot.3,d.fire.sp,all.x=TRUE)


d.view <- d.plot.3[,c("Fire","topoclim.cat","dom.tree.sp.ba","dom.tree.sp.obs","count.control")]






#### 6.5 Plot-level analysis with GLM; also run randomForest to get importance scores ####

d.plot.3 <- d.plot.3[d.plot.3$count.control > 4 & d.plot.3$count.highsev > 0,]


## Next, fit models ## 

library(brms)
library(loo)

d.plot <- d.plot.c

d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV %in% c(4,5)),]
# 
# ## remove an outlier plot with 20 abco that is preventing model conversion
# d.plot <- d.plot[!(d.plot$Regen_Plot == "CHI1248"),]

## remove another potential outlier plot: extremely high normal precip and high numbers of ABCO way above other plots with similar precip
d.plot <- d.plot[!(d.plot$Regen_Plot == "BTU1300185"),]

# ## remove a plot that has no grass cover data
# d.plot <- d.plot[!(d.plot$Regen_Plot == "AMR1300445"),]


## remove an NA fire value
d.plot <- d.plot[!is.na(d.plot$Fire),]

## other potential plots to remove: BTU1300050, CUB1300305

sp.opts <- c("PIPO","ABCO","PILA","SHADE.ALLSP","QUKE","PINUS.ALLSP","HDWD.ALLSP","PSME","CADE27","PIPJ")
cover.opts <- c("COV.SHRUB","COV.GRASS","COV.FORB")

# All
ht.opts <- c("HT.PIPJ","HT.HDWD.ALLSP","HT.PINUS.ALLSP","HT.SHADE.ALLSP")
htabs.opts <- c("HTABS.PIPJ","HTABS.PIPO","HTABS.SHRUB","HTABS.ABCO","HTABS.PSME","HTABS.PILA","HTABS.QUKE","HTABS.CADE27","HTABS.HDWD.ALLSP","HTABS.PINUS.ALLSP","HTABS.SHADE.ALLSP")
prop.opts <- c("PROP.CONIF","PROP.PINUS","PROP.SHADE")
prop.opts <- NULL

# # sp grps
# ht.opts <- c("HT.HDWD.ALLSP","HT.PINUS.ALLSP","HT.SHADE.ALLSP")
# htabs.opts <- c("HTABS.HDWD.ALLSP","HTABS.PINUS.ALLSP","HTABS.SHRUB","HTABS.SHADE.ALLSP")
# prop.opts <- c("PROP.CONIF","PROP.PINUS")
# 



resp.opts <- c(prop.opts,htabs.opts,sp.opts,cover.opts,ht.opts)


do.regression <-TRUE

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



# #Remove Dry fires
# dry.fires <- c("STRAYLOR","HARDING","ANTELOPE","MOONLIGHT")
# d.plot <- d.plot[!(d.plot$Fire %in% dry.fires),]

# Remove Bassetts
# d.plot <- d.plot[d.plot$Fire != "BASSETTS",]

# d.plot <- d.plot[!is.na(d.plot$diff.norm.snow.z),] #! need to enable this line if running code for snow models

for(sp in resp.opts) {
  
  cat("\n\n#####")
  cat("Running model for: ",sp,"")
  cat("#####\n\n")

  if(sp %in% cover.opts) {
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="PIPO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
    d.sp.curr.plt <- d.sp[d.sp$species=="PIPO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
  } else if(sp %in% ht.opts)  {
    sp.name <- substr(sp,4,100)
    d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp.name,]
    d.sp.curr.plt <- d.sp[d.sp$species==sp.name,]
  } else if(sp == "HTABS.SHRUB") {
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="PIPO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
    d.sp.curr.plt <- d.sp[d.sp$species=="PIPO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
  } else if(sp %in% htabs.opts) {
    sp.name <- substr(sp,7,100)
    d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp.name,]
    d.sp.curr.plt <- d.sp[d.sp$species==sp.name,]
  }  else if(sp == "PROP.PINUS") {
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="PINUS.ALLSP",]
    d.sp.curr.plt <- d.sp[d.sp$species=="PINUS.ALLSP",]
    d.sp.conif.plt <- d.sp[d.sp$species=="CONIF.ALLSP",]
    
    
    regen.pinus.old <- d.sp.curr.plt$regen.count.old / d.sp.conif.plt$regen.count.old
    regen.pinus.all <- d.sp.curr.plt$regen.count.all / d.sp.conif.plt$regen.count.all
    
    regen.pinus.old <- ifelse(is.nan(regen.pinus.old),NA,regen.pinus.old)
    regen.pinus.all <- ifelse(is.nan(regen.pinus.all),NA,regen.pinus.all)
    
    d.sp.curr.plt$prop.regen.pinus.old <- regen.pinus.old
    d.sp.curr.plt$prop.regen.pinus.all <- regen.pinus.all
    
  } else if(sp == "PROP.SHADE") {
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="SHADE.ALLSP",]
    d.sp.curr.plt <- d.sp[d.sp$species=="SHADE.ALLSP",]
    d.sp.conif.plt <- d.sp[d.sp$species=="CONIF.ALLSP",]
    
    
    regen.shade.old <- d.sp.curr.plt$regen.count.old / d.sp.conif.plt$regen.count.old
    regen.shade.all <- d.sp.curr.plt$regen.count.all / d.sp.conif.plt$regen.count.all
    
    regen.shade.old <- ifelse(is.nan(regen.shade.old),NA,regen.shade.old)
    regen.shade.all <- ifelse(is.nan(regen.shade.all),NA,regen.shade.all)
    
    d.sp.curr.plt$prop.regen.shade.old <- regen.shade.old
    d.sp.curr.plt$prop.regen.shade.all <- regen.shade.all
    
  } else if(sp == "PROP.CONIF") {
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="CONIF.ALLSP",]
    d.sp.curr.plt <- d.sp[d.sp$species=="CONIF.ALLSP",]
    d.sp.conif.plt <- d.sp[d.sp$species=="ALL",]
    
    
    regen.conif.old <- d.sp.curr.plt$regen.count.old / d.sp.conif.plt$regen.count.old
    regen.conif.all <- d.sp.curr.plt$regen.count.all / d.sp.conif.plt$regen.count.all
    
    regen.conif.old <- ifelse(is.nan(regen.conif.old),NA,regen.conif.old)
    regen.conif.all <- ifelse(is.nan(regen.conif.all),NA,regen.conif.all)
    
    d.sp.curr.plt$prop.regen.conif.old <- regen.conif.old
    d.sp.curr.plt$prop.regen.conif.all <- regen.conif.all
    
  }  else  {
    
    d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp,]
    d.sp.curr.plt <- d.sp[d.sp$species==sp,]
  }
  
  d <- merge(d.plot,d.plot.3,by=c("Fire","topoclim.cat")) # this effectively thins to plots that belong to a topoclimate category that has enough plots in it
  
  d <- merge(d,d.sp.curr.plt,by=c("Regen_Plot"))
  
  names(d.sp.curr.agg) <- paste0(names(d.sp.curr.agg),".agg")
  d <- merge(d,d.sp.curr.agg,by.x=c("Fire","topoclim.cat"),by.y=c("Fire.agg","topoclim.cat.agg")) # add the aggregated species data (from this we just want adult BA from the control plot)
  
  
  vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all","dominant_shrub_ht_cm","tallest_ht_cm","prop.regen.pinus.old","prop.regen.pinus.all","prop.regen.shade.old","prop.regen.hdwd.old","prop.regen.hdwd.old","prop.regen.conif.old")
  vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z", "def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","def.post","aet.post","adult.ba.agg","snow.post")
  d <- d[complete.cases(d[,vars.focal]),]
  d.c <- center.df(d,vars.leave)
  
  d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2
  d.c$tmean.normal_c.sq <- d.c$tmean.normal_c^2
  d.c$snow.normal_c.sq <- d.c$snow.normal_c^2
  
  d.c$def.normal_c.sq <- d.c$def.normal_c^2
  d.c$aet.normal_c.sq <- d.c$aet.normal_c^2
  
  # ####!!!! trick model: make diff.norm into diff.norm.min
  # d.c$diff.norm.ppt.z_c <- d.c$diff.norm.ppt.min.z_c
  # d.c$diff.norm.tmean.z_c <- d.c$diff.norm.tmean.max.z_c
  # d.c$diff.norm.aet.z_c <- d.c$diff.norm.aet.min.z_c
  # d.c$diff.norm.def.z_c <- d.c$diff.norm.def.max.z_c
  # #### end trick model
  
  
  d.c$diff.norm.ppt.z_c.sq <- d.c$diff.norm.ppt.z_c^2
  d.c$diff.norm.tmean.z_c.sq <- d.c$diff.norm.tmean.z_c^2
  d.c$diff.norm.snow.z_c.sq <- d.c$diff.norm.snow.z_c^2
  
  d.c$diff.norm.def.z_c.sq <- d.c$diff.norm.def.z_c^2
  d.c$diff.norm.aet.z_c.sq <- d.c$diff.norm.aet.z_c^2
  
  
  
  d.c$diff.norm.ppt.min.z_c.sq <- d.c$diff.norm.ppt.min.z_c^2
  d.c$diff.norm.tmean.max.z_c.sq <- d.c$diff.norm.tmean.max.z_c^2
  d.c$diff.norm.snow.min.z_c.sq <- d.c$diff.norm.snow.min.z_c^2
  
  
  d.c$diff.norm.def.max.z_c.sq <- d.c$diff.norm.def.max.z_c^2
  d.c$diff.norm.aet.min.z_c.sq <- d.c$diff.norm.aet.min.z_c^2
  
  
  
  
  d.c$ppt.post_c.sq <- d.c$ppt.post_c^2
  d.c$tmean.post_c.sq <- d.c$tmean.post_c^2
  d.c$snow.post_c.sq <- d.c$snow.post_c^2
  
  d.c$def.post_c.sq <- d.c$def.post_c^2
  d.c$aet.post_c.sq <- d.c$aet.post_c^2
  
  d.c$ppt.post.min_c.sq <- d.c$ppt.post.min_c^2
  #d.c$tmean.post.max_c.sq <- d.c$tmean.post.max_c^2
  d.c$snow.post.min_c.sq <- d.c$snow.post.min_c^2
  
  d.c$def.post.max_c.sq <- d.c$def.post.max_c^2
  d.c$aet.post.min_c.sq <- d.c$aet.post.min_c^2
  
  d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
  d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)
  
  
  vars.focal.c <- paste0(vars.focal[-6],"_c")
  
  ## transform cover so it does not include 0 or 1 (for Beta distrib)
  d.c$SHRUB.p <- d.c$SHRUB/100
  d.c$SHRUB.pt <- (d.c$SHRUB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
  
  d.c$GRASS.p <- d.c$GRASS/100
  d.c$GRASS.pt <- (d.c$GRASS.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
  
  d.c$HARDWOOD.p <- d.c$HARDWOOD/100
  d.c$HARDWOOD.pt <- (d.c$HARDWOOD.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
  
  d.c$FORB.p <- d.c$FORB/100
  d.c$FORB.pt <- (d.c$FORB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
  
  d.c$CONIFER.p <- d.c$CONIFER/100
  d.c$CONIFER.pt <- (d.c$CONIFER.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
  
  
  d.c$Fire <- as.factor(d.c$Fire)
  
  #d.c <- d.c[!(d.c$Fire == "RICH"),]
  
  d.c$regen.count.all.int <- ceiling(d.c$regen.count.all)
  d.c$regen.count.old.int <- ceiling(d.c$regen.count.old)
  
  
  ## Test whether seedling taller than shrub
  d.c$seedl.taller <- d.c$tallest_ht_cm > d.c$dominant_shrub_ht_cm
  
  
  
  
  
  
  if(sp %in% cover.opts) {
    
    sp.cov <- substr(sp,5,100)
    sp.cov <- paste0(sp.cov,".pt")
    
    d.c$response.var <- d.c[,sp.cov]
    
  } else if(sp %in% ht.opts){
    
    d.c <- d.c[d.c$regen.presab.old == TRUE,] # this is where we select whether we want all plots where the species was present or just old seedlings
    d.c$response.var <- d.c$seedl.taller
    
  } else if(sp == "HTABS.SHRUB") {
    d.c$response.var <- d.c$dominant_shrub_ht_cm
    d.c <- d.c[!is.na(d.c$response.var),]
  } else if(sp %in% htabs.opts) {
    d.c <- d.c[d.c$regen.presab.old == TRUE, ] # this is where we select whether we want all plots where the species was present or just old seedlings
    d.c$response.var <- d.c$tallest_ht_cm
  } else if(sp == "PROP.PINUS") {
    
    prop.regen.pinus.old.pt <- (d.c$prop.regen.pinus.old*(nrow(d.c)-1) + 0.5) / nrow(d.c)
    d.c$response.var <- prop.regen.pinus.old.pt
    
  } else if(sp %in% "PROP.SHADE") {
    
    prop.regen.shade.old.pt <- (d.c$prop.regen.shade.old*(nrow(d.c)-1) + 0.5) / nrow(d.c)
    d.c$response.var <- prop.regen.shade.old.pt
    
  } else if(sp %in% "PROP.CONIF") {
    
    prop.regen.conif.old.pt <- (d.c$prop.regen.conif.old*(nrow(d.c)-1) + 0.5) / nrow(d.c)
    d.c$response.var <- prop.regen.conif.old.pt
    
  } else if(sp %in% "PROP.HDWD") {
    
    prop.regen.hdwd.old.pt <- (d.c$prop.regen.hdwd.old*(nrow(d.c)-1) + 0.5) / nrow(d.c)
    d.c$response.var <- prop.regen.hdwd.old.pt
    
  }  else  {
    
    d.c$response.var <- d.c$regen.presab.old.01
    
  }
  
  if(do.regression) { 
    
    if(TRUE) {
      
      
      formulas <- list()
      
      
      #### No seed tree ####
      
      formulas[["n0.a0"]] <- formula(response.var ~ 1)
      formulas[["n0a.a0"]] <- formula(response.var ~ adult.ba.agg_c)
      
      ## PPT
      
      formulas[["n0a.aP"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.ppt.z_c)
      formulas[["n0.aP"]] <- formula(response.var ~ diff.norm.ppt.z_c)
      formulas[["n0a.aP2"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["n0.aP2"]] <- formula(response.var ~ diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      
      formulas[["nPa.a0"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c)
      formulas[["nPa.aP"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c)
      formulas[["nPa.aP2"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["pPa"]] <- formula(response.var ~ adult.ba.agg_c + ppt.post_c)
      formulas[["nP2a.a0"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + ppt.normal_c.sq)
      formulas[["nP2a.aP"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2a.aP2"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2a"]] <- formula(response.var ~ adult.ba.agg_c + ppt.post_c + ppt.post_c.sq)
      formulas[["nPa.aPni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c)
      formulas[["nPa.aP2ni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2a.aPni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2a.aP2ni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      
      formulas[["nP.a0"]] <- formula(response.var ~ ppt.normal_c)
      formulas[["nP.aP"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c)
      formulas[["nP.aP2"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["pP"]] <- formula(response.var ~ ppt.post_c)
      formulas[["nP2.a0"]] <- formula(response.var ~ ppt.normal_c + ppt.normal_c.sq)
      formulas[["nP2.aP"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2.aP2"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2"]] <- formula(response.var ~ ppt.post_c + ppt.post_c.sq)
      formulas[["nP.aPni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c)
      formulas[["nP.aP2ni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2.aPni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2.aP2ni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      

      ## snow

# 
#       formulas[["n0a.aS"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.snow.z_c)
#       formulas[["n0.aS"]] <- formula(response.var ~ diff.norm.snow.z_c)
#       formulas[["n0a.aS2"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.snow.z_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["n0.aS2"]] <- formula(response.var ~ diff.norm.snow.z_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
# 
#       formulas[["nSa.a0"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c)
#       formulas[["nSa.aS"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c)
#       formulas[["nSa.aS2"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["pSa"]] <- formula(response.var ~ adult.ba.agg_c + snow.post_c)
#       formulas[["nS2a.a0"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + snow.normal_c.sq)
#       formulas[["nS2a.aS"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2a.aS2"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2a"]] <- formula(response.var ~ adult.ba.agg_c + snow.post_c + snow.post_c.sq)
#       formulas[["nSa.aSni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c)
#       formulas[["nSa.aS2ni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["nS2a.aSni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2a.aS2ni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)
# 
#       formulas[["nS.a0"]] <- formula(response.var ~ snow.normal_c)
#       formulas[["nS.aS"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.z_c)
#       formulas[["nS.aS2"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["pS"]] <- formula(response.var ~ snow.post_c)
#       formulas[["nS2.a0"]] <- formula(response.var ~ snow.normal_c + snow.normal_c.sq)
#       formulas[["nS2.aS"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2.aS2"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2"]] <- formula(response.var ~ snow.post_c + snow.post_c.sq)
#       formulas[["nS.aSni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.z_c)
#       formulas[["nS.aS2ni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["nS2.aSni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2.aS2ni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)
# 






      ##DEF

      formulas[["n0a.aD"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.def.z_c)
      formulas[["n0a.aD2"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["n0.aD"]] <- formula(response.var ~ diff.norm.def.z_c)
      formulas[["n0.aD2"]] <- formula(response.var ~ diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)


      formulas[["nDa.a0"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c)
      formulas[["nDa.aD"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.z_c)
      formulas[["nDa.aD2"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["pDa"]] <- formula(response.var ~ adult.ba.agg_c + def.post_c)
      formulas[["nD2a.a0"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + def.normal_c.sq)
      formulas[["nD2a.aD"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2a.aD2"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
      formulas[["pD2a"]] <- formula(response.var ~ adult.ba.agg_c + def.post_c + def.post_c.sq)
      formulas[["nDa.aDni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.z_c)
      formulas[["nDa.aD2ni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2a.aDni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2a.aD2ni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)

      formulas[["nD.a0"]] <- formula(response.var ~ def.normal_c)
      formulas[["nD.aD"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c)
      formulas[["nD.aD2"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["pD"]] <- formula(response.var ~ def.post_c)
      formulas[["nD2.a0"]] <- formula(response.var ~ def.normal_c + def.normal_c.sq)
      formulas[["nD2.aD"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2.aD2"]] <- formula(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
      formulas[["pD2"]] <- formula(response.var ~ def.post_c + def.post_c.sq)
      formulas[["nD.aDni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c)
      formulas[["nD.aD2ni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2.aDni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2.aD2ni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



      ##AET

      formulas[["n0a.aA"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.aet.z_c)
      formulas[["n0a.aA2"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)

      formulas[["n0.aA"]] <- formula(response.var ~ diff.norm.aet.z_c)
      formulas[["n0.aA2"]] <- formula(response.var ~ diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)

      formulas[["nAa.a0"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c)
      formulas[["nAa.aA"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c)
      formulas[["nAa.aA2"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["pAa"]] <- formula(response.var ~ adult.ba.agg_c + aet.post_c)
      formulas[["nA2a.a0"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + aet.normal_c.sq)
      formulas[["nA2a.aA"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2a.aA2"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)
      formulas[["pA2a"]] <- formula(response.var ~ adult.ba.agg_c + aet.post_c + aet.post_c.sq)
      formulas[["nAa.aAni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c)
      formulas[["nAa.aA2ni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2a.aAni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2a.aA2ni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)

      formulas[["nA.a0"]] <- formula(response.var ~ aet.normal_c)
      formulas[["nA.aA"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.z_c)
      formulas[["nA.aA2"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["pA"]] <- formula(response.var ~ aet.post_c)
      formulas[["nA2.a0"]] <- formula(response.var ~ aet.normal_c + aet.normal_c.sq)
      formulas[["nA2.aA"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2.aA2"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)
      formulas[["pA2"]] <- formula(response.var ~ aet.post_c + aet.post_c.sq)
      formulas[["nA.aAni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.z_c)
      formulas[["nA.aA2ni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2.aAni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2.aA2ni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)




      
      
      ## PPT
      
      formulas[["n0a.aPmin"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.ppt.min.z_c)
      formulas[["n0a.aP2min"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["n0.aPmin"]] <- formula(response.var ~ diff.norm.ppt.min.z_c)
      formulas[["n0.aP2min"]] <- formula(response.var ~ diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      
      
      
      formulas[["nPa.aPmin"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c)
      formulas[["nPa.aP2min"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["pPmina"]] <- formula(response.var ~ adult.ba.agg_c + ppt.post.min_c)
      formulas[["nP2a.aPmin"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2a.aP2min"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2mina"]] <- formula(response.var ~ adult.ba.agg_c + ppt.post.min_c + ppt.post.min_c.sq)
      formulas[["nPa.aPminni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c)
      formulas[["nPa.aP2minni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2a.aPminni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2a.aP2minni"]] <- formula(response.var ~ adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      
      formulas[["nP.aPmin"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c)
      formulas[["nP.aP2min"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["pPmin"]] <- formula(response.var ~ ppt.post.min_c)
      formulas[["nP2.aPmin"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2.aP2min"]] <- formula(response.var ~ ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2min"]] <- formula(response.var ~ ppt.post.min_c + ppt.post.min_c.sq)
      formulas[["nP.aPminni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c)
      formulas[["nP.aP2minni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2.aPminni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2.aP2minni"]] <- formula(response.var ~ ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      
      
      ## Snow
# 
#       formulas[["n0a.aSmin"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.snow.min.z_c)
#       formulas[["n0a.aS2min"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.snow.min.z_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["n0.aSmin"]] <- formula(response.var ~ diff.norm.snow.min.z_c)
#       formulas[["n0.aS2min"]] <- formula(response.var ~ diff.norm.snow.min.z_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
# 
# 
# 
#       formulas[["nSa.aSmin"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c)
#       formulas[["nSa.aS2min"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["pSmina"]] <- formula(response.var ~ adult.ba.agg_c + snow.post.min_c)
#       formulas[["nS2a.aSmin"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2a.aS2min"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2mina"]] <- formula(response.var ~ adult.ba.agg_c + snow.post.min_c + snow.post.min_c.sq)
#       formulas[["nSa.aSminni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c)
#       formulas[["nSa.aS2minni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["nS2a.aSminni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2a.aS2minni"]] <- formula(response.var ~ adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
# 
#       formulas[["nS.aSmin"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.min.z_c)
#       formulas[["nS.aS2min"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["pSmin"]] <- formula(response.var ~ snow.post.min_c)
#       formulas[["nS2.aSmin"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2.aS2min"]] <- formula(response.var ~ snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2min"]] <- formula(response.var ~ snow.post.min_c + snow.post.min_c.sq)
#       formulas[["nS.aSminni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.min.z_c)
#       formulas[["nS.aS2minni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["nS2.aSminni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2.aS2minni"]] <- formula(response.var ~ snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
# 
# 

      ##DEF

      formulas[["n0a.aDmax"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.def.max.z_c)
      formulas[["n0a.aD2max"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

      formulas[["n0.aDmax"]] <- formula(response.var ~ diff.norm.def.max.z_c)
      formulas[["n0.aD2max"]] <- formula(response.var ~ diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

      formulas[["nDa.aDmax"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c)
      formulas[["nDa.aD2max"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["pDmaxa"]] <- formula(response.var ~ adult.ba.agg_c + def.post.max_c)
      formulas[["nD2a.aDmax"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2a.aD2max"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
      formulas[["pD2maxa"]] <- formula(response.var ~ adult.ba.agg_c + def.post.max_c + def.post.max_c.sq)
      formulas[["nDa.aDmaxni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c)
      formulas[["nDa.aD2maxni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2a.aDmaxni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2a.aD2maxni"]] <- formula(response.var ~ adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)

      formulas[["nD.aDmax"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c)
      formulas[["nD.aD2max"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["pDmax"]] <- formula(response.var ~ def.post.max_c)
      formulas[["nD2.aDmax"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2.aD2max"]] <- formula(response.var ~ def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
      formulas[["pD2max"]] <- formula(response.var ~ def.post.max_c + def.post.max_c.sq)
      formulas[["nD.aDmaxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c)
      formulas[["nD.aD2maxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2.aDmaxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2.aD2maxni"]] <- formula(response.var ~ def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



      ##AET

      formulas[["n0a.aAmin"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.aet.min.z_c)
      formulas[["n0a.aA2min"]] <- formula(response.var ~ adult.ba.agg_c + diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)

      formulas[["n0.aAmin"]] <- formula(response.var ~ diff.norm.aet.min.z_c)
      formulas[["n0.aA2min"]] <- formula(response.var ~ diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)


      formulas[["nAa.aAmin"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c)
      formulas[["nAa.aA2min"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["pAmina"]] <- formula(response.var ~ adult.ba.agg_c + aet.post.min_c)
      formulas[["nA2a.aAmin"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2a.aA2min"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)
      formulas[["pA2mina"]] <- formula(response.var ~ adult.ba.agg_c + aet.post.min_c + aet.post.min_c.sq)
      formulas[["nAa.aAminni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c)
      formulas[["nAa.aA2minni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2a.aAminni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2a.aA2minni"]] <- formula(response.var ~ adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)

      formulas[["nA.aAmin"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.min.z_c)
      formulas[["nA.aA2min"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["pAmin"]] <- formula(response.var ~ aet.post.min_c)
      formulas[["nA2.aAmin"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2.aA2min"]] <- formula(response.var ~ aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)
      formulas[["pA2min"]] <- formula(response.var ~ aet.post.min_c + aet.post.min_c.sq)
      formulas[["nA.aAminni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.min.z_c)
      formulas[["nA.aA2minni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2.aAminni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2.aA2minni"]] <- formula(response.var ~ aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)

      
      #### With seed tree ####
      
      formulas[["n0s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + 1)
      formulas[["n0as.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c)
      
      ## PPT
      
      formulas[["n0as.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.ppt.z_c)
      formulas[["n0s.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.z_c)
      formulas[["n0as.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["n0s.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      
      formulas[["nPas.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c)
      formulas[["nPas.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c)
      formulas[["nPas.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["pPas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.post_c)
      formulas[["nP2as.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + ppt.normal_c.sq)
      formulas[["nP2as.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2as.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2as"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.post_c + ppt.post_c.sq)
      formulas[["nPas.aPni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c)
      formulas[["nPas.aP2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2as.aPni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2as.aP2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      
      formulas[["nPs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c)
      formulas[["nPs.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c)
      formulas[["nPs.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["pPs"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.post_c)
      formulas[["nP2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + ppt.normal_c.sq)
      formulas[["nP2s.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2s.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2s"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.post_c + ppt.post_c.sq)
      formulas[["nPs.aPni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c)
      formulas[["nPs.aP2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["nP2s.aPni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2s.aP2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.z_c + ppt.normal_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      
      
      


      ## Snow
# 
#       formulas[["n0as.aS"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.snow.z_c)
#       formulas[["n0s.aS"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.snow.z_c)
#       formulas[["n0as.aS2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.snow.z_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["n0s.aS2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.snow.z_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
# 
#       formulas[["nSas.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c)
#       formulas[["nSas.aS"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c)
#       formulas[["nSas.aS2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["pSas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.post_c)
#       formulas[["nS2as.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + snow.normal_c.sq)
#       formulas[["nS2as.aS"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2as.aS2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2as"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.post_c + snow.post_c.sq)
#       formulas[["nSas.aSni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c)
#       formulas[["nSas.aS2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["nS2as.aSni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2as.aS2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)
# 
#       formulas[["nSs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c)
#       formulas[["nSs.aS"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.z_c)
#       formulas[["nSs.aS2"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["pSs"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.post_c)
#       formulas[["nS2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + snow.normal_c.sq)
#       formulas[["nS2s.aS"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2s.aS2"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.z_c + snow.normal_c*diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2s"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.post_c + snow.post_c.sq)
#       formulas[["nSs.aSni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.z_c)
#       formulas[["nSs.aS2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq)
#       formulas[["nS2s.aSni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c.sq)
#       formulas[["nS2s.aS2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.z_c + snow.normal_c + diff.norm.snow.z_c.sq + diff.norm.snow.z_c.sq + snow.normal_c.sq)





      ##DEF

      formulas[["n0as.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.def.z_c)
      formulas[["n0as.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["n0s.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.z_c)
      formulas[["n0s.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)


      formulas[["nDas.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c)
      formulas[["nDas.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.z_c)
      formulas[["nDas.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["pDas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.post_c)
      formulas[["nD2as.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + def.normal_c.sq)
      formulas[["nD2as.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2as.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
      formulas[["pD2as"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.post_c + def.post_c.sq)
      formulas[["nDas.aDni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.z_c)
      formulas[["nDas.aD2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2as.aDni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2as.aD2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)

      formulas[["nDs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c)
      formulas[["nDs.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c)
      formulas[["nDs.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["pDs"]] <- formula(response.var ~ seed_tree_distance_general_c + def.post_c)
      formulas[["nD2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + def.normal_c.sq)
      formulas[["nD2s.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2s.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
      formulas[["pD2s"]] <- formula(response.var ~ seed_tree_distance_general_c + def.post_c + def.post_c.sq)
      formulas[["nDs.aDni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c)
      formulas[["nDs.aD2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2s.aDni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2s.aD2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



      ##AET

      formulas[["n0as.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.aet.z_c)
      formulas[["n0as.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)

      formulas[["n0s.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.z_c)
      formulas[["n0s.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)

      formulas[["nAas.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c)
      formulas[["nAas.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c)
      formulas[["nAas.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["pAas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.post_c)
      formulas[["nA2as.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + aet.normal_c.sq)
      formulas[["nA2as.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2as.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)
      formulas[["pA2as"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.post_c + aet.post_c.sq)
      formulas[["nAas.aAni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c)
      formulas[["nAas.aA2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2as.aAni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2as.aA2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)

      formulas[["nAs.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c)
      formulas[["nAs.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c)
      formulas[["nAs.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["pAs"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.post_c)
      formulas[["nA2s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + aet.normal_c.sq)
      formulas[["nA2s.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2s.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)
      formulas[["pA2s"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.post_c + aet.post_c.sq)
      formulas[["nAs.aAni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c)
      formulas[["nAs.aA2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2s.aAni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2s.aA2ni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)



      
      
      
      ## PPT
      
      formulas[["n0as.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.ppt.min.z_c)
      formulas[["n0as.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["n0s.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.min.z_c)
      formulas[["n0s.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      
      
      
      formulas[["nPas.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c)
      formulas[["nPas.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["pPminas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.post.min_c)
      formulas[["nP2as.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2as.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2minas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.post.min_c + ppt.post.min_c.sq)
      formulas[["nPas.aPminni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c)
      formulas[["nPas.aP2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2as.aPminni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2as.aP2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      
      formulas[["nPs.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c)
      formulas[["nPs.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["pPmins"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.post.min_c)
      formulas[["nP2s.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2s.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2mins"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.post.min_c + ppt.post.min_c.sq)
      formulas[["nPs.aPminni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c)
      formulas[["nPs.aP2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2s.aPminni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2s.aP2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      
      
      


      ## Snow
# 
#       formulas[["n0as.aSmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.snow.min.z_c)
#       formulas[["n0as.aS2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.snow.min.z_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["n0s.aSmin"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.snow.min.z_c)
#       formulas[["n0s.aS2min"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.snow.min.z_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
# 
# 
# 
#       formulas[["nSas.aSmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c)
#       formulas[["nSas.aS2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["pSminas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.post.min_c)
#       formulas[["nS2as.aSmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2as.aS2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2minas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.post.min_c + snow.post.min_c.sq)
#       formulas[["nSas.aSminni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c)
#       formulas[["nSas.aS2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["nS2as.aSminni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2as.aS2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
# 
#       formulas[["nSs.aSmin"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.min.z_c)
#       formulas[["nSs.aS2min"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["pSmins"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.post.min_c)
#       formulas[["nS2s.aSmin"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2s.aS2min"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c*diff.norm.snow.min.z_c + snow.normal_c*diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
#       formulas[["pS2mins"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.post.min_c + snow.post.min_c.sq)
#       formulas[["nSs.aSminni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.min.z_c)
#       formulas[["nSs.aS2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq)
#       formulas[["nS2s.aSminni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c.sq)
#       formulas[["nS2s.aS2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + snow.normal_c + diff.norm.snow.min.z_c + snow.normal_c + diff.norm.snow.min.z_c.sq + diff.norm.snow.min.z_c.sq + snow.normal_c.sq)
# 



      ##DEF

      formulas[["n0as.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.def.max.z_c)
      formulas[["n0as.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

      formulas[["n0s.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.max.z_c)
      formulas[["n0s.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

      formulas[["nDas.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c)
      formulas[["nDas.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["pDmaxas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.post.max_c)
      formulas[["nD2as.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2as.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
      formulas[["pDmax2as"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.post.max_c + def.post.max_c.sq)
      formulas[["nDas.aDmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c)
      formulas[["nDas.aD2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2as.aDmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2as.aD2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)

      formulas[["nDs.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c)
      formulas[["nDs.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["pDmaxs"]] <- formula(response.var ~ seed_tree_distance_general_c + def.post.max_c)
      formulas[["nD2s.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2s.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
      formulas[["pD2maxs"]] <- formula(response.var ~ seed_tree_distance_general_c + def.post.max_c + def.post.max_c.sq)
      formulas[["nDs.aDmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c)
      formulas[["nDs.aD2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2s.aDmaxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2s.aD2maxni"]] <- formula(response.var ~ seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



      ##AET

      formulas[["n0as.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.aet.min.z_c)
      formulas[["n0as.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)

      formulas[["n0s.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.min.z_c)
      formulas[["n0s.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)


      formulas[["nAas.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c)
      formulas[["nAas.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["pAminas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.post.min_c)
      formulas[["nA2as.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2as.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)
      formulas[["pA2minas"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.post.min_c + aet.post.min_c.sq)
      formulas[["nAas.aAminni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c)
      formulas[["nAas.aA2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2as.aAminni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2as.aA2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)

      formulas[["nAs.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c)
      formulas[["nAs.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["pAmins"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.post.min_c)
      formulas[["nA2s.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2s.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)
      formulas[["pA2mins"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.post.min_c + aet.post.min_c.sq)
      formulas[["nAs.aAminni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c)
      formulas[["nAs.aA2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2s.aAminni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2s.aA2minni"]] <- formula(response.var ~ seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)

    }
    
    
    
    
    
    
    cv.results <- data.frame()
    
    cat("\n")
    
    for(i in 1:length(formulas)) {
      
      #cat("\rEvaluating model ",i," of ",length(formulas))
      
      formula.name <- names(formulas)[i]
      cv.result <- cvfun.fire(formulas[[i]],data=d.c)
      cv.result.named <- data.frame(model=names(formulas)[i],cv.result,stringsAsFactors=FALSE)
      cv.results <- rbind(cv.results,cv.result.named)
      
    }
    
    cv.results$model <- as.character(cv.results$model)
    
    
    
    
    
    
    
    ## finding the best normal and then the corresponding best anomaly and post
    
    ## normal models
    search <- ".a0"
    normal.model.names <- cv.results[grep(search,cv.results$model,fixed=TRUE),]$mod
    # normal.model.names <- normal.model.names[normal.model.names != "n0.a0"] # exclude null model

    if(sp %in% cover.opts) { # don't allow a normal with adults in it if this is a cover model
      norm.search.opts <- c(Pmean = "n(P|0)2?\\.a0",
                            Dmean = "n(D|0)2?\\.a0",
                            Amean = "n(A|0)2?\\.a0",
                            # Smean = "n(S|0)2?\\.a0",
                            Pmin = "n(P|0)2?\\.a0",
                            Dmax = "n(D|0)2?\\.a0",
                            Amin = "n(A|0)2?\\.a0"
                            
                            # ,
                            # Smin = "n(S|0)2?\\.a0"
                            )
    } else {
      norm.search.opts <- c(Pmean = "n(P|0)2?a?s?\\.a0",
                            Dmean = "n(D|0)2?a?s?\\.a0",
                            Amean = "n(A|0)2?a?s?\\.a0",
                            # Smean = "n(S|0)2?a?s?\\.a0",
                            Pmin = "n(P|0)2?a?s?\\.a0",
                            Dmax = "n(D|0)2?a?s?\\.a0",
                            Amin = "n(A|0)2?a?s?\\.a0"
                            # ,
                            # Smin = "n(S|0)2?a?s?\\.a0"
                            )
    }

    anom.search.opts <- c(Pmean = "aP2?(ni)?$",
                          Dmean = "aD2?(ni)?$",
                          Amean = "aA2?(ni)?$",
                          # Smean = "aS2?(ni)?$",
                          Pmin = "aP2?min(ni)?$",
                          Dmax = "aD2?max(ni)?$",
                          Amin = "aA2?min(ni)?$"
                          # ,
                          # Smin = "aS2?min(ni)?$"
    )


    post.search.opts <- c(Pmean = "pP2?a?s?$",
                          Dmean = "pD2?a?s?$",
                          Amean = "pA2?a?s?$",
                          # Smean = "pS2?a?s?$",
                          Pmin = "pP2?mina?s?$",
                          Dmax = "pD2?maxa?s?$",
                          Amin = "pA2?mina?s?$"
                          # ,
                          # Smin = "pS2?mina?s?$"
    )


    
    
    # ## For ppt only
    # 
    # if(sp %in% cover.opts) { # don't allow a normal with adults in it if this is a cover model
    #   norm.search.opts <- c(Pmean = "n(P|0)2?\\.a0",
    #                         Pmin = "n(P|0)2?\\.a0")
    # } else {
    #   norm.search.opts <- c(Pmean = "n(P|0)2?a?s?\\.a0",
    #                         Pmin = "n(P|0)2?a?s?\\.a0")
    # }
    # 
    # anom.search.opts <- c(Pmean = "aP2?(ni)?$",
    #                       Pmin = "aP2?min(ni)?$")
    # 
    # 
    # post.search.opts <- c(Pmean = "pP2?a?s?$",
    #                       Pmin = "pP2?mina?s?$")
    # 





    
    
    d.maes.anoms.sp <- data.frame()
    
    
    for(i in 1:length(anom.search.opts)) {
      
      anom.name <- names(norm.search.opts)[i]
      norm.search <- norm.search.opts[i]
      anom.search <- anom.search.opts[i]
      null.names <- c("n0.a0","n0a.a0","n0s.a0","n0as.a0")
      normal.model.names <- cv.results[grep(norm.search,cv.results$model,fixed=FALSE),]$model
      
      d.cv.normal <- cv.results[cv.results$model %in% normal.model.names,]
      best.normal.mod <- d.cv.normal[which(d.cv.normal$mae == min(d.cv.normal$mae,na.rm=TRUE)),]$mod[1]
      
      #what is the best anomaly model for that normal, sticking with whichever anom we have
      best.normal.normal.part <- strsplit(as.character(best.normal.mod),".",fixed=TRUE)[[1]][1]
      best.normal.anom.search <- paste0(best.normal.normal.part,"\\.",anom.search)
      anom.model.names <- cv.results[grep(best.normal.anom.search,cv.results$model,fixed=FALSE),]$model
      d.cv.anom <- cv.results[cv.results$model %in% anom.model.names,]
      best.anom.mod <- d.cv.anom[which(d.cv.anom$mae == min(d.cv.anom$mae,na.rm=TRUE)),]$mod[1]
      
      #what is the best normal, excluding null models
      normal.model.names.nonull <- normal.model.names[!(normal.model.names %in% null.names)]
      d.cv.normal.nonull <- cv.results[cv.results$model %in% normal.model.names.nonull,]
      best.normal.nonull.mod <- d.cv.normal.nonull[which(d.cv.normal.nonull$mae == min(d.cv.normal.nonull$mae,na.rm=TRUE)),]$mod[1]
      
      
  
      
      
      #what is the best post model for that normal, sticking with whichever anom we have
      post.mod.name.1 <- gsub("^n","p",best.normal.nonull.mod)
      post.mod.name.2 <- substr(post.mod.name.1,1,nchar(post.mod.name.1)-3)
      extreme.text <- "(min)?(max)?"
      post.mod.name.p1 <- substr(post.mod.name.2,1,2)
      post.mod.name.p2 <- substr(post.mod.name.2,3,3)
      post.mod.name.p3 <- substr(post.mod.name.2,4,100)
      post.mod.search <- paste0(post.mod.name.p1,extreme.text,post.mod.name.p2,extreme.text,post.mod.name.p3)
      post.model.names <- cv.results[grep(post.mod.search,cv.results$model,fixed=FALSE),]$model
      d.cv.post <- cv.results[cv.results$model %in% post.model.names,]
      best.post.mod <- d.cv.post[which(d.cv.post$mae == min(d.cv.post$mae,na.rm=TRUE)),]$mod[1]
      
      
      ## get maes of best normal, best anomal, best normal.nonull, best.post
      best.anom.mae <- cv.results[cv.results$model == best.anom.mod,"mae"]
      best.anom.normal.mae <- cv.results[cv.results$model == best.normal.mod,"mae"]
      best.normal.nonull.mae <- cv.results[cv.results$model == best.normal.nonull.mod,"mae"]
      best.post.mae <- cv.results[cv.results$model == best.post.mod,"mae"]
      
      
      ### Store MAEs of the best models
      best.anomaly.mod <- best.anom.mod
      best.anom.normal <- best.normal.mod
      best.normal.nonull <- best.normal.nonull.mod
      best.post <- best.post.mod
      
      
      d.maes.anoms.sp.anom <- data.frame(best.anomaly.mod,best.anom.normal,best.normal.nonull,best.post,
                                         best.anom.mae,best.anom.normal.mae,best.normal.nonull.mae,best.post.mae,
                                         sp,anom.name,stringsAsFactors=FALSE)
      
      d.maes.anoms.sp <- rbind(d.maes.anoms.sp,d.maes.anoms.sp.anom)
      
    }
    
    d.maes.anoms <- rbind(d.maes.anoms,d.maes.anoms.sp)
    
    
    
    
    
    
    
    
    
    
    
    ## Now for each anom model , make predictions
    
    
    
    
    ### Function to find the center of a "fixed" predictor variable within the species distribution
    mid.val.fun <- function(var) {
      mid.val <- mean(d.c[,var],na.rm=TRUE)
      return(mid.val)
    }
    
    low.val.fun <- function(var) {
      a <- quantile(d.c[,var],0.25,na.rm=TRUE)
      #a <- min(d.c[,var])
      return(a)
    }
    
    high.val.fun <- function(var) {
      a <- quantile(d.c[,var],0.75,na.rm=TRUE)
      #a <- max(d.c[,var])
      return(a)
    }
    
    vars <- c("ppt.normal_c","ppt.normal_c.sq","diff.norm.ppt.z_c","diff.norm.ppt.z_c.sq","tmean.normal_c","tmean.normal_c.sq",
              "diff.norm.tmean.z_c","diff.norm.tmean.z_c.sq",
              "snow.normal_c","diff.norm.snow.z_c","diff.norm.snow.z_c.sq","diff.norm.snow.min.z_c","diff.norm.snow.min.z_c.sq",
              "aet.normal_c","aet.normal_c.sq","diff.norm.aet.z_c","diff.norm.aet.z_c.sq","def.normal_c","def.normal_c.sq",
              "diff.norm.def.z_c","diff.norm.def.z_c.sq",
              "seed_tree_distance_general_c","rad.march_c","adult.ba.agg_c"
    )
    
    
    
    mid.val <- sapply(vars,mid.val.fun,USE.NAMES=TRUE)
    
    low.val <- sapply(vars,low.val.fun,USE.NAMES=TRUE)
    names(low.val) <- vars
    
    high.val <- sapply(vars,high.val.fun,USE.NAMES=TRUE)
    names(high.val) <- vars
    
    
    
    diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100)
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
      snow.normal_c = c(rep(low.val["snow.normal_c"],100),rep(mid.val["snow.normal_c"],100),rep(high.val["snow.normal_c"],100)),
      snow.normal_c.sq = c(rep(low.val["snow.normal_c"]^2,100),rep(mid.val["snow.normal_c"]^2,100),rep(high.val["snow.normal_c"]^2,100)),
      
      norm.level = c(rep("low",100),rep("mid",100),rep("high",100)),
      
      diff.norm.ppt.z_c = rep(diff.norm.seq,3),
      diff.norm.ppt.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.tmean.z_c = rep(diff.norm.seq.rev,3),
      diff.norm.tmean.z_c.sq = rep(diff.norm.seq.rev,3)^2,
      diff.norm.aet.z_c = rep(diff.norm.seq,3),
      diff.norm.aet.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.def.z_c = rep(diff.norm.seq.rev,3),
      diff.norm.def.z_c.sq = rep(diff.norm.seq.rev,3)^2,
      diff.norm.snow.z_c = rep(diff.norm.seq.rev,3),
      diff.norm.snow.z_c.sq = rep(diff.norm.seq.rev,3)^2,
      
      diff.norm.ppt.min.z_c = rep(diff.norm.seq,3),
      diff.norm.ppt.min.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.tmean.max.z_c = rep(diff.norm.seq.rev,3),
      diff.norm.tmean.max.z_c.sq = rep(diff.norm.seq.rev,3)^2,
      diff.norm.aet.min.z_c = rep(diff.norm.seq,3),
      diff.norm.aet.min.z_c.sq = rep(diff.norm.seq^2,3),
      diff.norm.def.max.z_c = rep(diff.norm.seq.rev,3),
      diff.norm.def.max.z_c.sq = rep(diff.norm.seq.rev,3)^2,
      diff.norm.snow.min.z_c = rep(diff.norm.seq.rev,3),
      diff.norm.snow.min.z_c.sq = rep(diff.norm.seq.rev,3)^2,
      
      seed_tree_distance_general_c = mid.val["seed_tree_distance_general_c"],
      rad.march_c = mid.val["rad.march_c"],
      adult.ba.agg_c = mid.val["adult.ba.agg_c"]
      
    )
    
    for(j in 1:nrow(d.maes.anoms.sp)) {
      
      d.maes.anoms.sp.row <- d.maes.anoms.sp[j,]
      best.anom.mod <- d.maes.anoms.sp.row$best.anomaly.mod
      best.anom.normal.mod <- d.maes.anoms.sp.row$best.anom.normal
      
  
      ##predict to new data
      
      if(sp %in% c(cover.opts,prop.opts)) {
      
        
        d.c.complete <- d.c[!is.na(d.c$response.var),]
        
        nboot <- 30
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
        
  
      } else if (sp %in% htabs.opts) {
        
        
        d.c.complete <- d.c[!is.na(d.c$response.var),]
        
        #fit the anom and normal models
        mod.anom <- glm(formulas[[best.anom.mod]],data=d.c,family="gaussian")
        mod.norm <- glm(formulas[[best.anom.normal.mod]],data=d.c,family="gaussian")
        
        p <- predict(mod.anom,newdat,type="link",se.fit=TRUE)
        low <- p$fit - 1.96*p$se
        high <- p$fit + 1.96*p$se
        high[high>400] <- 400
        pred.anom <- data.frame(pred.low=(low),pred.mid=(p$fit),pred.high=(high))
        
        p <- predict(mod.norm,newdat,type="link",se.fit=TRUE)
        low <- p$fit - 1.96*p$se
        high <- p$fit + 1.96*p$se
        high[high>400] <- 400
        pred.norm <- data.frame(pred.low=(low),pred.mid=(p$fit),pred.high=(high))
        
      } else {
        
        d.c.complete <- d.c[!is.na(d.c$response.var),]
        
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
      
      
      
      ## get fitted and observed
      fit.anom <- as.data.frame(predict(mod.anom))
      fit.norm <- as.data.frame(predict(mod.norm))
      names(fit.anom) <- names(fit.norm) <- c("fitted")
      fit.anom$type <- "anom"
      fit.norm$type <- "norm"
      
      d.c.complete <- d.c[!is.na(d.c$response.var),]
      
      fit.anom.dat <- cbind(fit.anom,d.c.complete)
      fit.norm.dat <- cbind(fit.norm,d.c.complete)
      
      fit.dat.sp <- rbind(fit.anom.dat,fit.norm.dat)
      fit.dat.sp$anom <- d.maes.anoms.sp.row$anom.name
      fit.dat.sp$sp <- sp
      
      fit.dat <- rbind.fill(fit.dat,fit.dat.sp)
      
      if(! (sp %in% c(prop.opts,cover.opts,htabs.opts))) {
        
        fit.dat$fitted <- inv.logit(fit.dat$fitted)
        
      }
      
      
      
      
      
    }
  }
  
  # 
  # ### run random forest
  # 
  # 
  # pred.names <- c("seed_tree_distance_general_c", "adult.ba.agg_c", "ppt.normal_c", "diff.norm.ppt.z_c", "diff.norm.ppt.min.z_c")
  # 
  # data <- d.c[,c(pred.names,"response.var")]
  # data <- data[complete.cases(data),]
  # 
  # if(!(sp %in% c(htabs.opts,prop.opts,cover.opts))) { # if it's not a continuous variable response
  #   data$response.var <- as.factor(data$response.var)
  # }
  # 
  # m2.mean <- cforest(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.z_c,data=data)
  # 
  # rf.importance.sp <- as.data.frame(t(varimp(m2.mean,conditional=TRUE)))
  # rf.importance.sp$sp <- sp
  # rf.importance.sp$anom <- "mean"
  # rf.importance <- rbind.fill(rf.importance,rf.importance.sp)
  # 
  # m2.min <- cforest(response.var ~ seed_tree_distance_general_c + adult.ba.agg_c + ppt.normal_c + diff.norm.ppt.min.z_c,data=data)
  # 
  # rf.importance.sp <- as.data.frame(t(varimp(m2.min,conditional=TRUE)))
  # rf.importance.sp$sp <- sp
  # rf.importance.sp$anom <- "min"
  # 
  # 
  # rf.importance <- rbind.fill(rf.importance,rf.importance.sp)
  # 
  # 
  # ### predict for min
  # 
  # preds <- predict(m2.min,newdata=newdat,type="prob")
  # preds.df <- data.frame()
  # 
  # for(k in 1:length(preds)) {
  #   preds.df <- rbind(preds.df,preds[[k]])  
  # }
  # 
  # if(!(sp %in% c(htabs.opts,prop.opts,cover.opts,ht.opts))) { # if it's not a continuous variable response
  #   pred <- preds.df$response.var.1
  # } else if (sp %in% ht.opts) {
  #   pred <- preds.df$response.var.TRUE
  # } else {
  #   pred <- preds.df$response.var
  # } 
  # 
  # pred.newdat <- cbind(pred,newdat)
  # 
  # pred.newdat$anom <- "min"
  # pred.newdat$species <- sp
  # 
  # pred.rf <- rbind(pred.rf,pred.newdat)
  # 
  # 
  # 
  # ### predict for mean
  # 
  # preds <- predict(m2.mean,newdata=newdat,type="prob")
  # preds.df <- data.frame()
  # 
  # for(k in 1:length(preds)) {
  #   preds.df <- rbind(preds.df,preds[[k]])  
  # }
  # 
  # if(!(sp %in% c(htabs.opts,prop.opts,cover.opts,ht.opts))) { # if it's not a continuous variable response
  #   pred <- preds.df$response.var.1
  # } else if (sp %in% ht.opts) {
  #   pred <- preds.df$response.var.TRUE
  # }  else {
  #   pred <- preds.df$response.var
  # }
  # 
  # pred.newdat <- cbind(pred,newdat)
  # 
  # pred.newdat$anom <- "mean"
  # pred.newdat$species <- sp
  # 
  # pred.rf <- rbind(pred.rf,pred.newdat)
  # 
}


#write.csv(rf.importance,"rf_importance.csv")

#### Prep for plots etc ####

## make a single factor that combines species and rad level
pred.dat$sp.rad <- paste(pred.dat$sp,pred.dat$rad.level,sep=".")


#for each sp-rad combo, find which was the best anomaly, which anomalies were better than their corresponding normals, and plot symbols indicating
d.maes.anoms$anom.better <- d.maes.anoms$best.anom.mae < d.maes.anoms$best.anom.normal.mae
d.maes.anoms$anom.improvement <-  d.maes.anoms$best.anom.normal.mae - d.maes.anoms$best.anom.mae

d.maes.anoms$best.of.species <- ""
d.maes.anoms$most.improved <- ""

for(sp in unique(d.maes.anoms$sp)) {
  d.maes.anoms.sp <- d.maes.anoms[d.maes.anoms$sp == sp,]
  for(rad in d.maes.anoms.sp$rad.level) {
    d.maes.anoms.sp.rad <- d.maes.anoms.sp[d.maes.anoms.sp$rad.level == rad,]
    
    best.anom <- d.maes.anoms.sp.rad[d.maes.anoms.sp.rad$best.anom.mae == min(d.maes.anoms.sp.rad$best.anom.mae)[1],"anom.name"]
    d.maes.anoms[(d.maes.anoms$sp==sp) & (d.maes.anoms$rad.level == rad) & (d.maes.anoms$anom.name == best.anom),"best.of.species"] <- ""# "best anom"
    
    most.improved <- d.maes.anoms.sp.rad[d.maes.anoms.sp.rad$anom.improvement == max(d.maes.anoms.sp.rad$anom.improvement)[1],"anom.name"]
    d.maes.anoms[(d.maes.anoms$sp==sp) & (d.maes.anoms$rad.level == rad) & (d.maes.anoms$anom.name == most.improved),"most.improved"] <- ""# "largest improvement"
  }
}

d.maes.anoms.short <- d.maes.anoms[,c("sp","anom.name","anom.better","anom.improvement","best.of.species","most.improved")]

pred.dat.comb <- merge(pred.dat,d.maes.anoms.short,all.x=TRUE,by.x=c("sp","anom"),by.y=c("sp","anom.name"))
pred.dat.comb$anom.improvement <- round(pred.dat.comb$anom.improvement,2)
pred.dat.comb[pred.dat.comb$anom.improvement < 0 , "anom.improvement"] <- NA


#### Compute total MAE across all species, to see which is the best ####
d.mae.agg <- aggregate(d.maes.anoms[,c("best.anom.mae","best.anom.normal.mae","anom.improvement")],by=list(d.maes.anoms$anom.name),FUN=mean,na.rm=TRUE)



#### Plot counterfactuals ####

## get rid of mid and get rid of normal predictions
pred.dat.plotting <- pred.dat.comb[pred.dat.comb$norm.level %in% c("low","high"),]
#pred.dat.plotting <- pred.dat.comb[pred.dat.comb$norm.level %in% c("mid"),]
pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$type == "anom",]

#pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$sp=="PILA",]
#pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$anom=="Amin",]

ggplot(pred.dat.plotting,aes(x=diff.norm.ppt.z_c,y=pred.mid,color=norm.level,fill=norm.level)) +
  geom_point() +
  geom_ribbon(aes(ymin=pred.low,ymax=pred.high),alpha=0.3,color=NA) +
  facet_wrap(anom~sp.rad,nrow=8,scales="free_y") +
  geom_text(aes(0,1,label=best.of.species),size=3,color="black") +
  geom_text(aes(0,0.9,label=most.improved),size=3,color="black") +
  geom_text(aes(0,0.8,label=anom.improvement),size=3,color="black")




#### testinttgggg


pred.rf.plot <- pred.rf[pred.rf$norm.level %in% c("low","high"),]

ggplot(pred.rf.plot,aes(x=diff.norm.ppt.z_c,y=pred,color=norm.level,fill=norm.level)) +
  geom_line() +
  facet_wrap(anom~species,scales="free_y",nrow=2)




#### Display the normal vs. post analysis ####
library(reshape)

##calc how much better the post is than the normal
d.maes.anoms$post.v.norm <- d.maes.anoms$best.normal.nonull.mae - d.maes.anoms$best.post.mae
d.maes.anoms$post.v.norm[d.maes.anoms$post.v.norm < 0] <- NA

##make a table: species x anom, where value is the normal-post improvement

# positive means post has lower error
post.table <- cast(d.maes.anoms,anom.name~sp,value='post.v.norm')
post.table







#### For each species and rad group, plot fitted vs. observed, for normal and anom side by side ####

fit.dat$resid <- fit.dat$response.var-fit.dat$fitted




## AET min only
fit.dat.ppt <- fit.dat[fit.dat$anom == "Pmin",]



fit.obs.plots <- list()

for(sp in unique(fit.dat$sp)) {
  #for(rad.level in unique(fit.dat$rad.level)) {
    
    #fit.dat.sp <- fit.dat.aet[fit.dat.aet$sp == sp & fit.dat.aet$rad.level == rad.level,]
  fit.dat.sp <- fit.dat.ppt[fit.dat.ppt$sp == sp,]
    
    fit.dat.sp$type <- factor(fit.dat.sp$type,levels=c("norm","anom"))
    
    # dummy points
    max.xy <- max(c(fit.dat.sp$response.var,fit.dat.sp$fitted))
    dummy <- data.frame(response.var=c(0,max.xy),fitted=c(0,max.xy))
    
    plot.name <- paste0(sp,".",rad.level)
    
    fit.obs.plots[[plot.name]] <- ggplot(fit.dat.sp,aes(x=response.var,y=fitted)) +
      geom_point() +
      facet_grid(.~type) +
      geom_abline(slope=1,intercept=0) +
      geom_point(data=dummy,alpha=0) +
      labs(x="Observed",y="Fitted",title=plot.name) +
      theme_bw() +
      theme(plot.title=element_text(hjust=0.5))
    
  #}
}

library(gridExtra)
n <- length(fit.obs.plots)
nCol <- 4
do.call("grid.arrange",c(fit.obs.plots,ncol=nCol))




#### Plot residuals by fire ####
fit.dat.plot <- fit.dat[fit.dat$type=="norm",]
fit.dat.plot <- fit.dat.plot[fit.dat.plot$anom == "Pmin",]

# fit.dat.plot$rad.level <- as.factor(fit.dat.plot$rad.level)


plot(fit.dat.plot$diff.norm.aet.min.z.highsev_c~fit.dat.plot$Fire)


fit.dat.plot$Fire <- factor(fit.dat.plot$Fire,c("BAGLEY","CHIPS","RALSTON","BASSETTS","MOONLIGHT","ANTELOPE","HARDING","FREDS","POWER","STRAYLOR","RICH","BTU LIGHTENING","CUB","AMERICAN RIVER"))


ggplot(fit.dat.plot,aes(x=Fire,y=resid)) +
  geom_boxplot(position="dodge",width=0.5) +
  facet_grid(sp~.) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=1))

ggplot(fit.dat.plot,aes(y=diff.norm.ppt.min.z.highsev_c,x=fit.dat.plot$Fire)) +
  geom_point()










#### 6.8 Plot-level analaysis decision trees (PARTY)) ####



library(party)

### prep the data frame






sp <- "SHADE.ALLSP"



d.plot <- d.plot.c

d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV %in% c(4,5)),]

## remove an outlier plot with 20 abco that is preventing model conversion
d.plot <- d.plot[!(d.plot$Regen_Plot == "CHI1248"),]

## remove another potential outlier plot: extremely high normal precip and high numbers of ABCO way above other plots with similar precip
d.plot <- d.plot[!(d.plot$Regen_Plot == "BTU1300185"),]

## remove a plot that has no grass cover data
d.plot <- d.plot[!(d.plot$Regen_Plot == "AMR1300445"),]


## remove an NA fire value
d.plot <- d.plot[!is.na(d.plot$Fire),]
# 
# 
# #Dry fires
# dry.fires <- c("STRAYLOR","HARDING","ANTELOPE","MOONLIGHT")
# d.plot <- d.plot[!(d.plot$Fire %in% dry.fires),]
# 



if(sp %in% cover.opts) {
  d.sp.curr.agg <- d.sp.2[d.sp.2$species=="ABCO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
  d.sp.curr.plt <- d.sp[d.sp$species=="ABCO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
} else {
  d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp,]
  d.sp.curr.plt <- d.sp[d.sp$species==sp,]
}

d <- merge(d.plot,d.plot.3,by=c("Fire","topoclim.cat")) # this effectively thins to plots that belong to a topoclimate category that has enough plots in it

d <- merge(d,d.sp.curr.plt,by=c("Regen_Plot"))

names(d.sp.curr.agg) <- paste0(names(d.sp.curr.agg),".agg")
d <- merge(d,d.sp.curr.agg,by.x=c("Fire","topoclim.cat"),by.y=c("Fire.agg","topoclim.cat.agg")) # add the aggregated species data (from this we just want adult BA from the control plot)




vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z", "def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","def.post","aet.post","adult.ba.agg")
d <- d[complete.cases(d[,vars.focal]),]
d.c <- center.df(d,vars.leave)

d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2
d.c$tmean.normal_c.sq <- d.c$tmean.normal_c^2
d.c$snow.normal_c.sq <- d.c$snow.normal_c^2

d.c$def.normal_c.sq <- d.c$def.normal_c^2
d.c$aet.normal_c.sq <- d.c$aet.normal_c^2

# ####!!!! trick model: make diff.norm into diff.norm.min
# d.c$diff.norm.ppt.z_c <- d.c$diff.norm.ppt.min.z_c
# d.c$diff.norm.tmean.z_c <- d.c$diff.norm.tmean.max.z_c
# d.c$diff.norm.aet.z_c <- d.c$diff.norm.aet.min.z_c
# d.c$diff.norm.def.z_c <- d.c$diff.norm.def.max.z_c
# #### end trick model


d.c$diff.norm.ppt.z_c.sq <- d.c$diff.norm.ppt.z_c^2
d.c$diff.norm.tmean.z_c.sq <- d.c$diff.norm.tmean.z_c^2
d.c$diff.norm.snow.z_c.sq <- d.c$diff.norm.snow.z_c^2

d.c$diff.norm.def.z_c.sq <- d.c$diff.norm.def.z_c^2
d.c$diff.norm.aet.z_c.sq <- d.c$diff.norm.aet.z_c^2



d.c$diff.norm.ppt.min.z_c.sq <- d.c$diff.norm.ppt.min.z_c^2
d.c$diff.norm.tmean.max.z_c.sq <- d.c$diff.norm.tmean.max.z_c^2
d.c$diff.norm.snow.min.z_c.sq <- d.c$diff.norm.snow.min.z_c^2


d.c$diff.norm.def.max.z_c.sq <- d.c$diff.norm.def.max.z_c^2
d.c$diff.norm.aet.min.z_c.sq <- d.c$diff.norm.aet.min.z_c^2




d.c$ppt.post_c.sq <- d.c$ppt.post_c^2
d.c$tmean.post_c.sq <- d.c$tmean.post_c^2

d.c$def.post_c.sq <- d.c$def.post_c^2
d.c$aet.post_c.sq <- d.c$aet.post_c^2

d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)


vars.focal.c <- paste0(vars.focal[-6],"_c")

## transform cover so it does not include 0 or 1 (for Beta distrib)
d.c$SHRUB.p <- d.c$SHRUB/100
d.c$SHRUB.pt <- (d.c$SHRUB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$GRASS.p <- d.c$GRASS/100
d.c$GRASS.pt <- (d.c$GRASS.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$HARDWOOD.p <- d.c$HARDWOOD/100
d.c$HARDWOOD.pt <- (d.c$HARDWOOD.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$FORB.p <- d.c$FORB/100
d.c$FORB.pt <- (d.c$FORB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$CONIFER.p <- d.c$CONIFER/100
d.c$CONIFER.pt <- (d.c$CONIFER.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)


d.c$Fire <- as.factor(d.c$Fire)

#d.c <- d.c[!(d.c$Fire == "RICH"),]

d.c$regen.count.all.int <- ceiling(d.c$regen.count.all)
d.c$regen.count.old.int <- ceiling(d.c$regen.count.old)


if(sp %in% cover.opts) {
  
  sp.cov <- substr(sp,5,100)
  sp.cov <- paste0(sp.cov,".pt")
  
  d.c$response.var <- d.c[,sp.cov]
  
} else {
  
  d.c$response.var <- d.c$regen.presab.old.01
  
}




d.c$response.var <- ifelse(d.c$response.var==1,"present","absent")
d.c$response.var <- as.factor(d.c$response.var)

a <- ctree(response.var ~ ppt.normal_c + tmean.normal_c + snow.normal_c + aet.normal_c + def.normal_c +
             diff.norm.ppt.z_c + diff.norm.tmean.z_c + diff.norm.snow.z_c + diff.norm.aet.z_c + diff.norm.def.z_c +
             seed_tree_distance_general_c + adult.ba.agg_c, data=d.c)


a <- ctree(response.var ~ def.normal_c +
             diff.norm.def.z_c +
             seed_tree_distance_general_c + adult.ba.agg_c + rad.march_c, data=d.c)


a <- ctree(response.var ~ ppt.normal_c + 
           diff.norm.ppt.min.z_c +
           seed_tree_distance_general_c + adult.ba.agg_c, data=d.c)



#control=ctree_control(testtype=c("Univariate"),mincriterion=0.80)

plot(a)























#### 8. Cluster-level multivariate analysis ####




#d.plot.3
#d.sp.2
library(data.table)



focal.sp <- c("PIPO","ABCO","PILA","QUKE") #   "ABMA","QUCH2","ARME") #ABMA, PIPO, QUCH2, 
focal.sp <- c("PIPO","ABCO","PILA","QUKE","QUCH2","LIDE3","ABMA")
#focal.sp <- c("PIPO","ABCO","PILA","PSME","QUKE","LIDE3","QUCH2","CADE27")
#focal.sp <- c("PINUS.ALLSP","SHADE.ALLSP","HDWD.ALLSP")
focal.cols <- c("Fire","topoclim.cat","species","regen.presab.old","regen.presab.all","adult.ba")
#focal.cols <- c("Fire","topoclim.cat","species","regen.count.old","regen.count.all","adult.ba")


d.sp.simp <- d.sp.2[d.sp.2$species %in% focal.sp,focal.cols]

names(d.sp.simp)[4:6] <- c("r.old","r.all","a.ba")

d.sp.simp <- as.data.table(d.sp.simp)

d.sp.cast <- dcast(d.sp.simp,topoclim.cat + Fire~species,value.var=c("r.old","r.all","a.ba"))


regen.var <- "r.old"

d.all <- merge(d.plot.3,d.sp.cast,by=c("Fire","topoclim.cat"))

# #Dry fires
# dry.fires <- c("STRAYLOR","HARDING","ANTELOPE","MOONLIGHT")
# d.all <- d.all[!(d.all$Fire %in% dry.fires),]

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

sublist <- c("LIDE3" = "Tanoak",
             "PILA" = "SP",
             "ABCO" = "WF",
             "ABMA" = "RF",
             "QUKE" = "BlkOak",
             "QUCH2" = "CynOak",
             "PIPO" = "PP",
             "rad.march.highsev" = "solar.rad",
             "tmean.normal.highsev" = "temp.norm",
             "ppt.normal.highsev" = "ppt.norm",
             "diff.norm.ppt.z.highsev" = "ppt.anom",
             "diff.norm.tmean.z.highsev" = "temp.anom")

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




library(vegan)


### try adding diff predictors
# 
# c0 <- cca(d.all.sp ~ 1,data=d.all)
# c1 <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev + seed_tree_distance_general.highsev + BA.Live1.control,data=d.all)
# c1.fire <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev + seed_tree_distance_general.highsev + BA.Live1.control + Fire,data=d.all)
# c2  <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev + seed_tree_distance_general.highsev + BA.Live1.control + diff.norm.ppt.z.highsev + diff.norm.tmean.z.highsev,data=d.all)
# 
# c3  <- cca(d.all.sp ~ ppt.normal.highsev + rad.march.highsev + diff.norm.ppt.z.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
# c3.nosp  <- cca(d.all.sp ~ ppt.normal.highsev + rad.march.highsev + diff.norm.ppt.z.highsev,data=d.all)
# c4 <- cca(d.all.sp ~ ppt.normal.highsev + rad.march.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
# 
# c3.def  <- cca(d.all.sp ~ def.normal.highsev + rad.march.highsev + diff.norm.def.z.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
# c4.def <- cca(d.all.sp ~ def.normal.highsev + rad.march.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
# 
# c.sp.def  <- cca(d.all.sp ~ def.normal.highsev + aet.normal.highsev + rad.march.highsev + diff.norm.def.z.highsev + diff.norm.aet.z.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
# c.nosp.def<- cca(d.all.sp ~ def.normal.highsev + aet.normal.highsev + rad.march.highsev + diff.norm.def.z.highsev + diff.norm.aet.z.highsev,data=d.all)
# 
# c.sp  <- cca(d.all.sp ~ tmean.normal.highsev + ppt.normal.highsev + rad.march.highsev + diff.norm.tmean.z.highsev + diff.norm.ppt.z.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PILA + a.ba_QUKE,data=d.all)
# c.nosp<- cca(d.all.sp ~ tmean.normal.highsev + ppt.normal.highsev + rad.march.highsev + diff.norm.tmean.z.highsev + diff.norm.ppt.z.highsev,data=d.all)
# 



d.all$SolarRadiation <- d.all$solar.rad
d.all$NormalPrecip <- d.all$ppt.norm
d.all$NormalTemp <- d.all$temp.norm

d.all$PrecipAnom <- d.all$ppt.anom
d.all$TempAnom <- d.all$temp.anom





# all non-anomaly, excluding species
cc1 <- cca(d.all.sp ~ NormalPrecip + NormalTemp + SolarRadiation,data=d.all)

# all non-anomaly, including species
cc2 <- cca(d.all.sp ~ NormalPrecip + NormalTemp + SolarRadiation + WF + PP + SP + BlkOak + CynOak + Tanoak + RF,data=d.all)

# non-anamoly and anomaly, excluding species
#cc2b <- cca(d.all.sp ~ NormalPrecipal.highsev + tmean.normal.highsev + diff.norm.ppt.z.highsev + diff.norm.tmean.z.highsev + rad.march.highsev,data=d.all)

# all including anomaly, including species
cc3 <- cca(d.all.sp ~ NormalPrecip + NormalTemp + PrecipAnom + TempAnom + SolarRadiation + WF + PP + SP + BlkOak + CynOak + Tanoak + RF,data=d.all)


cc1
cc2
#cc2b
cc3

plot(cc1)
plot(cc2,choices=c(1,2))
#plot(cc2b,choices=c(1,2))
plot(cc3,choices=c(1,2))




vars.focal.nmds <- c("NormalPrecip","NormalTemp","PrecipAnom","TempAnom","SolarRadiation","WF","PP","SP","BlkOak","CynOak","Tanoak","RF")
d.focal.nmds <- d.all[,vars.focal.nmds]


a <- metaMDS(d.all.sp)
b <- envfit(a,d.focal.nmds)


plot(a,type="t")
plot(b)














##
library(devtools)
install_github('fawda123/ggord')
library(ggord)

png("plot1.png", 800, 300);


ggord(cc1, ptslab=TRUE,size=2,addsize=5,txt=5,margins=0) +
  theme(legend.position = 'top') +
  theme(plot.margin=unit(c(0,0,0,0),"cm"))

dev.off();



### Redo CCA for wet plots only
median.ppt <- median(d.all$ppt.norm)
keep <- which(d.all$ppt.norm > median.ppt)
d.all.wet <- d.all[keep,]
d.all.sp.wet <- d.all.sp[keep,]

cc3.wet <- cca(d.all.sp.wet ~ ppt.anom + temp.anom,data=d.all.wet)
plot(cc3.wet)


### Redo CCA for dry plots only
median.ppt <- median(d.all$ppt.norm)
keep <- which((d.all$ppt.norm < median.ppt) & !(d.all$Fire == "MOONLIGHT" & d.all$topoclim.cat == "P.2_R.1")) #second part is to remove an outlier category that had red fir even though these are supposed to be dry plots
d.all.wet <- d.all[keep,]
d.all.sp.wet <- d.all.sp[keep,]


cc3.wet <- cca(d.all.sp.wet ~ ppt.anom + temp.anom + solar.rad,data=d.all.wet)
plot(cc3.wet)





# 
# #### repeat, but by cover types
# 
# 
# 
# # all non-anomaly, excluding species
# cc1 <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev,data=d.all)
# 
# # all non-anomaly, including species
# cc2 <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev + a.ba_PINUS.ALLSP + a.ba_SHADE.ALLSP + a.ba_HDWD.ALLSP,data=d.all)
# 
# # non-anamoly and anomaly, excluding species
# cc2b <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + diff.norm.ppt.z.highsev + diff.norm.tmean.z.highsev + rad.march.highsev,data=d.all)
# 
# # all including anomaly, including species
# cc3 <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + diff.norm.ppt.z.highsev + diff.norm.tmean.z.highsev + rad.march.highsev + a.ba_PINUS.ALLSP + a.ba_SHADE.ALLSP + a.ba_HDWD.ALLSP,data=d.all)
# 
# cc1
# cc2
# cc3



Cairo(file=paste0("../Figures/Fig5_CCA1_",Sys.Date(),".png"),width=1300,height=1200,ppi=200,res=200,dpi=200)
plot(cc1,choices=c(1,2))
dev.off()

Cairo(file=paste0("../Figures/Fig5_CCA2_",Sys.Date(),".png"),width=1300,height=1200,ppi=200,res=200,dpi=200)
plot(cc2,choices=c(1,2))
dev.off()

Cairo(file=paste0("../Figures/Fig5_CCA3_",Sys.Date(),".png"),width=1300,height=1200,ppi=200,res=200,dpi=200)
plot(cc3,choices=c(1,2))
dev.off()

Cairo(file=paste0("../Figures/Fig5_CCA4wet_",Sys.Date(),".png"),width=1200,height=1200,ppi=200,res=200,dpi=200)
plot(cc3.wet,choices=c(1,2))
dev.off()








#### 9. Mantel test to explain regen species comp with control species comp ####

regen.var <- "r.old"
control.var <- "a.ba"
d.all <- merge(d.plot.3,d.sp.cast,by=c("Fire","topoclim.cat"))
sp.cols <- grep(regen.var,names(d.all)) #regen
a.cols <- grep(control.var,names(d.all)) #adult (control)
d.all.sp <- d.all[,sp.cols]
d.all.a <- d.all[,a.cols]
d.all.sp.tot <- rowSums(d.all.sp)
d.all.a.tot <- rowSums(d.all.a)
keep.rows <- ((d.all.sp.tot > 0) & (d.all.a.tot > 0))
d.all.sp <- d.all.sp[keep.rows,]
d.all.a <- d.all.a[keep.rows,]
d.all <- d.all[keep.rows,]

#d.all.sp, d.all.a

regen.dist <- vegdist(d.all.sp,method="jaccard")
adult.dist <- vegdist(d.all.a,method="jaccard")

mantel(regen.dist,adult.dist)


## significant relationship between the two, relatively high r








