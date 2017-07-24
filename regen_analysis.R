setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(party)
library(ggplot2)
library(brms)
library(pROC)
library(betareg)
library(car)
library(plyr)
library(data.table)

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

# only northern Sierra Nevada fires
sierra.fires <- c("Straylor","Cub","Rich","Moonlight","Antelope","BTU Lightning","Harding","Bassetts","American River","Ralston","Freds","Power","Bagley","Chips")
d.plot <- d.plot[d.plot$Fire %in% sierra.fires,]

## Remove managed plots, plots in nonforested locations (e.g. exposed bedrock), etc.
plots.exclude <- read.csv("data_intermediate/plots_exclude.csv",header=T,stringsAsFactors=FALSE)
plot.ids.exclude <- plots.exclude[plots.exclude$Exclude != "",]$Regen_Plot
d.plot <- d.plot[!(d.plot$Regen_Plot %in% plot.ids.exclude),]


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

# severity categories
high.sev <- c(4,5) # which field-assessed severity values are considered high severity
control <- c(0,1) # which field-assessed severity values are considered controls

# categorize severities
d.plot$FIRE_SEV.cat <- recode(d.plot$FIRE_SEV,"control='control';high.sev='high.sev';else=NA")
d.plot <- d.plot[!is.na(d.plot$FIRE_SEV.cat),]

# remove plots that are high severity but not surveyed in years 4-5 post-fire
d.plot <- d.plot[which(!(d.plot$FIRE_SEV.cat == "high.sev" & !(d.plot$survey.years.post %in% c(4,5)))),]


# must have shrub cover
d.plot <- d.plot[!is.na(d.plot$SHRUB),]




# 
# 
# 
# d.plot <- d.plot[!((d.plot$Fire == "American River") & (d.plot$ppt.normal < 1700)),]
# d.plot <- d.plot[!((d.plot$Fire == "BTU Lightning") & (d.plot$ppt.normal > 2000)),]
# d.plot <- d.plot[!((d.plot$Fire == "Bagley") & (d.plot$rad.march < 4000)),]
# d.plot <- d.plot[!((d.plot$Fire == "Cub") & (d.plot$ppt.normal < 1300)),]
# d.plot <- d.plot[!((d.plot$Fire == "Freds") & (d.plot$ppt.normal > 1230)),]
# d.plot <- d.plot[!((d.plot$Fire == "Power") & (d.plot$rad.march < 6000)),]
# d.plot <- d.plot[!((d.plot$Fire == "Ralston") & (d.plot$ppt.normal < 1175)),]
# d.plot <- d.plot[!((d.plot$Fire == "Rich") & (d.plot$ppt.normal > 1080) & (d.plot$rad.march < 5000)),]
# 







#### 1. Assign each plot a topoclimatic category ####
fires <- unique(d.plot$Fire)

for(fire in fires) {
  
  ## Precipitation categories
  # determine what the precipitation breakpoints should be (here just using median) -- based on control plots only
  breaks <- quantile(d.plot[(d.plot$Fire == fire) & (d.plot$FIRE_SEV %in% control),]$ppt.normal,probs=c(0.5),na.rm=TRUE)
  
  
  # for some fires with a small range of precip, override the precip breaks, so we just have one category per fire
  fires.small.precip.range <- c("American River","Antelope","Bagley","Bassetts","BTU Lightning","Harding","Straylor","Rich")
  if(fire %in% fires.small.precip.range) {
    breaks <- 9999
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
  breaks <- 6000
  
  categories <- categorize(d.plot[d.plot$Fire==fire,]$rad.march,breaks,name="R")
  # store it into the plot data.frame
  d.plot[d.plot$Fire==fire,"rad.category"] <- categories
  
}
  
## Create one variable reflecting the all-way factorial combination of topoclimatic categories
d.plot$topoclim.cat <- paste(d.plot$precip.category,d.plot$rad.category,sep="_")



### Plot relevant "topoclimate space" for each fire and see how the categories broke them down

d.plot.precat <- d.plot


# look at both together
ggplot(d.plot.precat,aes(x=ppt.normal,y=rad.march,col=topoclim.cat,shape=FIRE_SEV.cat)) +
  geom_point() +
  facet_wrap(~Fire,scales="free") +
  theme_bw(16) +
  scale_shape_manual(values=c(1,16))



## Remove climatic regions that do not have comparable controls and high sev
## NOTE: this is only necessary and recommended when doing an analysis that involves comparing control adults with high sev regen




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




#### 3. Remove irrelevant data rows ####

# Remove the topoclimatic categories with too few plots in either burned or control
d.plot.3 <- d.plot.2[which((d.plot.2$count.control > 4) & (d.plot.2$count.highsev > 0)),]

# only want to analyze high-severity plots burned 4-5 years post-fire
d.plot.c <- d.plot.c[(d.plot.c$survey.years.post %in% c(4,5)) & (d.plot.c$FIRE_SEV %in% c(4,5)),] 












#### 4. Compute dominant vegetation for each remaining topoclimate category, based on stated observed nearby dominants. use all plots regardless of severity ####

non.tree.dom.veg <- c("ARPA6","CEIN3","CECO","CEPR") # species to exclude from list

for(fire in unique(d.plot.3$Fire)) {
  
  d.plot.fire <- d.plot.3[d.plot.3$Fire == fire,]
  
  #for(topoclim.cat in unique(d.plot.fire$topoclim.cat)) {
    
    
    # d.plot.fire.cat <- d.plot.c[(d.plot.c$Fire == fire) & (d.plot.c$topoclim.cat == topoclim.cat),]
  d.plot.fire.cat <- d.plot.c[(d.plot.c$Fire == fire),] 
    
    if (nrow(d.plot.fire.cat) < 5) {
      cat("Less than 5 plots in ",fire," ",topoclim.cat,". Skipping.\n")
      next()
    }
  
    if (nrow(d.plot.fire.cat) < 5) {
      cat("Less than 5 plots in ",fire," ",". Skipping.\n")
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
    
    
    #d.plot.3[(d.plot.3$Fire == fire) & (d.plot.3$topoclim.cat == topoclim.cat),"dom.tree.sp.obs"] <- dom.sp
    d.plot.3[(d.plot.3$Fire == fire),"dom.tree.sp.obs"] <- dom.sp
    
  #}
  
}








#### 5. Determine dominant adult tree species (from control plots) within each category ####
library(reshape)

d.sp.pre <- d.sp.2[,c("species","topoclim.cat","Fire","adult.ba")]

#d.fire.sp <- cast(d.sp.pre,Fire + topoclim.cat ~ species,value="adult.ba")
d.fire.sp <- cast(d.sp.pre,Fire ~ species,value="adult.ba",fun.aggregate=sum)

d.fire.sp$PIPJ <- d.fire.sp$PIPO + d.fire.sp$PIJE

sp.names <- names(d.fire.sp)
names.drop <- c("ALL","ALNUS","CONIF.ALLSP","CONIFER","HDWD.ALLSP","JUNIPERUS","PINUS.ALLSP","SHADE.ALLSP","ABIES","PIPO","PIJE")
sp.names <- sp.names[!(sp.names %in% names.drop)]
d.fire.sp <- d.fire.sp[,sp.names]
d.fire.sp[,3:length(names(d.fire.sp))] <- d.fire.sp[,3:length(names(d.fire.sp))] / rowSums(d.fire.sp[,3:length(names(d.fire.sp))])

d.fire.sp <- d.fire.sp[complete.cases(d.fire.sp),]

## for each one, get a list of all species with > 20% BA

d.fire.sp$dom.tree.sp.ba <- ""

for (i in 1:nrow(d.fire.sp)) {
  
  row.ba <- d.fire.sp[i,3:ncol(d.fire.sp)]
  
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

#d.fire.sp <- d.fire.sp[,c("Fire","topoclim.cat","dom.tree.sp.ba")]
d.fire.sp <- d.fire.sp[,c("Fire","dom.tree.sp.ba")]

d.plot.3 <- merge(d.plot.3,d.fire.sp,all.x=TRUE)


d.domsp <- d.plot.3[,c("Fire","topoclim.cat","dom.tree.sp.ba","dom.tree.sp.obs","count.control")]


### Combine dom sp lists: take ba-based and append any sp from obs-based that are not in list.

d.domsp$domsp.comb <- NULL

for(i in 1:nrow(d.domsp)) {
  
  row <- d.domsp[i,]
  
  ba <- row$dom.tree.sp.ba
  obs <- row$dom.tree.sp.obs
  
  ba.list <- strsplit(ba,split=", ")[[1]]
  obs.list <- strsplit(obs,split=", ")[[1]]
  
  diff <- setdiff(obs.list,ba.list)
  
  full <- paste(c(ba.list,diff),collapse=", ")
  
  d.domsp[i,"domsp.comb"] <- full
  
  
}

# add to plot-level data

d.domsp <- d.domsp[!duplicated(d.domsp$Fire),]
d.domsp <- remove.vars(d.domsp,c("topoclim.cat","count.control"))

d.plot.c <- merge(d.plot.c,d.domsp,all.x=TRUE,by=c("Fire"))






#### 6. Compile summary statistics and tables ####



### fire-level ###

d <- as.data.table(d.plot.c)
#d <- merge(d,d.plot.3[,c("Fire","topoclim.cat")]) # keep only the plot data that has corresponding topoclim cats that are being kept (i.e. correct severity, etc.)


d.fire <- d[,list(fire.year=mean(fire.year),
                    years.post=paste(sort(unique(survey.years.post)),collapse=", "),
                    nplots=.N,
                    ppt.norm=summary.string(ppt.normal,decimals=0),
                    ppt.anom=summary.string(diff.norm.ppt.min.z,decimals=2),
                    dom.sp=first(domsp.comb)
                  ),
                  by=Fire]

d.fire$Fire <- sapply(d.fire$Fire,FUN= function(x) simpleCap(tolower(x)))

d.fire[d.fire$Fire == "Btu Lightening","Fire"] <- "BTU Lightning"




### topoclim-cat level ###

#! keep only the species data that has corresponding plots that are being kept (i.e. correct severity, etc.)
d <- d.sp.2

keep.sp <- c("ABCO","PIPJ","PILA","PSME","PINUS.ALLSP","SHADE.ALLSP","HDWD.ALLSP")
keep.cols <- c("species","regen.presab.old","adult.ba","topoclim.cat","Fire")
d <- d[d$species %in% keep.sp,]

d <- as.data.table(d)

d.sp.cast <- dcast(d,Fire + topoclim.cat~species,value.var=c("regen.presab.old","adult.ba"))

# add the other necessary columns
keep.cols <- c("Fire","topoclim.cat","SHRUB.highsev","GRASS.highsev","FORB.highsev","count.highsev","count.control")
d2 <- d.plot.3[,keep.cols]

#merge all necessary data together
d.merge <- merge(d2,d.sp.cast,by=c("Fire","topoclim.cat"))

d.merge <- as.data.table(d.merge)

d.cat <- d.merge[,list(Fire=Fire,
                       topo.cat=topoclim.cat,
                       
                       n.highsev = count.highsev,
                       n.ref= count.control,
                       
                       regen.ABCO = 100*round(regen.presab.old_ABCO,digits=2),
                       regen.PILA = 100*round(regen.presab.old_PILA,digits=2),
                       regen.PIPJ = 100*round(regen.presab.old_PIPJ,digits=2),
                       regen.PSME = 100*round(regen.presab.old_PSME,digits=2),
                       regen.Pinus = 100*round(regen.presab.old_PINUS.ALLSP,digits=2),
                       regen.Shade = 100*round(regen.presab.old_SHADE.ALLSP,digits=2),
                       regen.Hdwd = 100*round(regen.presab.old_HDWD.ALLSP,digits=2),
                       
                       ref.ABCO = round(adult.ba_ABCO,digits=0),
                       ref.PILA = round(adult.ba_PILA,digits=0),
                       ref.PIPJ = round(adult.ba_PIPJ,digits=0),
                       ref.PSME = round(adult.ba_PSME,digits=0),
                       ref.Pinus = round(adult.ba_PINUS.ALLSP,digits=0),
                       ref.Shade = round(adult.ba_SHADE.ALLSP,digits=0),
                       ref.Hdwd = round(adult.ba_HDWD.ALLSP,digits=0),
                       
                       cov.shrub=round(SHRUB.highsev,digits=0),
                       cov.forb=round(FORB.highsev,digits=0),
                       cov.grass=round(GRASS.highsev,digits=0)
                       )
                  ]


#### 7. Plot highsev vs. reference ####


d.plot.keep <- merge(d.plot,d.plot.3[,c("Fire","topoclim.cat")])


d.plot.keep$topoclim.cat <- gsub(".","",d.plot.keep$topoclim.cat,fixed=TRUE)

d.plot.keep$FIRE_SEV.cat <- gsub("control","Reference",d.plot.keep$FIRE_SEV.cat,fixed=TRUE)
d.plot.keep$FIRE_SEV.cat <- gsub("high.sev","High severity",d.plot.keep$FIRE_SEV.cat,fixed=TRUE)



p <- ggplot(d.plot.keep,aes(x=ppt.normal,y=rad.march,col=topoclim.cat,shape=FIRE_SEV.cat)) +
  geom_point(size=3) +
  facet_wrap(~Fire,scales="free") +
  theme_bw(16) +
  scale_shape_manual(values=c(1,3)) +
  labs(x="Normal annual precipitation (mm)",y="Solar exposure (Wh m-2 day-1)",shape="Plot type",color="Topoclimate category")

Cairo(file=paste0("../Figures/FigS1_ref_vs_highsev",Sys.Date(),".png"),width=2800,height=2000,ppi=200,res=200,dpi=200)
p
dev.off()


#### 8. Plot climate space with fire labels ####


d.plot.c$Fire <- as.factor(d.plot.c$Fire)

### Plot of monitoring plots in climate space: minimum precip ### 



# For each fire, make a point at mean normal precip and mean precip anom 
fire.centers <- data.frame() 
fires <- unique(d.plot.c$Fire) 
for(fire in fires) { 
  d.fire <- d.plot.c[d.plot.c$Fire == fire,] 
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

##shifts

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


p1 <-  ggplot(d.plot.c,aes(x=ppt.normal,y=diff.norm.ppt.min.z,color=Fire)) + 
  geom_point(size=3) + 
  #geom_point(data=fire.centers,size=5,pch=1) + 
  geom_text(data=fire.centers,aes(label=fire.year),nudge_y=.04,size=6,color="darkgray") + 
  guides(color=FALSE) + 
  theme_bw(17) + 
  labs(x="Normal precipitation (mm)",y="Postfire minimum precipitation anomaly (SD)") + 
  scale_x_continuous(limits=c(180,2580))+
  geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=lines,color="darkgray",size=1.0) +
  scale_color_viridis(discrete=TRUE)



### Plot of monitoring plots in climate space: mean precip ### 

# For each fire, make a point at mean normal precip and mean precip anom 
fire.centers <- data.frame() 
fires <- unique(d.plot.c$Fire) 
for(fire in fires) { 
  d.fire <- d.plot.c[d.plot.c$Fire == fire,] 
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




p2 <- ggplot(d.plot.c,aes(x=ppt.normal,y=diff.norm.ppt.z,color=Fire)) + 
  geom_point(size=3) + 
  #geom_point(data=fire.centers,size=5,pch=1) + 
  geom_text(data=fire.centers,aes(label=fire.year),nudge_y=.04,size=6,color="darkgray") + 
  guides(color=FALSE) + 
  theme_bw(17) + 
  labs(x="Normal precipitation (mm)",y="Postfire mean precipitation anomaly (SD)") + 
  scale_x_continuous(limits=c(180,2580)) +
  scale_color_viridis(discrete=TRUE)
  #geom_segment(aes(x=x,y=y,xend=xend,yend=yend),data=lines,color="darkgray",size=1.0)


library(gridExtra)
blank<-rectGrob(gp=gpar(col="white")) # make a white spacer grob

Cairo(file=paste0("../Figures/Fig1_climateSpace_2pptanom_",Sys.Date(),".png"),width=3100,height=1500,ppi=200,res=200,dpi=200) 
grid.arrange(p1,blank,p2,ncol=3,widths=c(0.49,0.02,0.49))
dev.off()




#### 9. Graphical summary of abiotic environment comparing reference vs. highsev ####

keep.vars <- c("Fire","topoclim.cat","ppt.normal","rad.march","diff.norm.ppt.min.z","diff.norm.ppt.z","FIRE_SEV.cat") # removed seed_tree_distance_general

d <- as.data.table(merge(d.plot[,keep.vars],d.plot.3[,c("Fire","topoclim.cat")]))

d.med <- d[,lapply(.SD,median),by=list(Fire,topoclim.cat,FIRE_SEV.cat)]
d.med$level <- "mid"

d.low <- d[,lapply(.SD,quantile,probs=0.25),by=list(Fire,topoclim.cat,FIRE_SEV.cat)]
d.low$level <- "low"

d.high <- d[,lapply(.SD,quantile,probs=0.75),by=list(Fire,topoclim.cat,FIRE_SEV.cat)]
d.high$level <- "high"

d.agg <- rbind(d.med,d.low,d.high)

d.melt <- melt(d.agg,id.vars=c("Fire","topoclim.cat","level","FIRE_SEV.cat"))
d.cast <- dcast(d.melt,Fire + topoclim.cat + FIRE_SEV.cat + variable ~ level)

d.cast$variable <- gsub("ppt.normal","Normal precipitation (mm)",d.cast$variable)
d.cast$variable <- gsub("rad.march","Solar exposure (Wh m-2 day-1)",d.cast$variable)
d.cast$variable <- gsub("diff.norm.ppt.min.z","Minimum precipitation anomaly (SD)",d.cast$variable)
d.cast$variable <- gsub("diff.norm.ppt.z","Mean precipitation anomaly (SD)",d.cast$variable)
d.cast$variable <- gsub("seed_tree_distance_general","Seed tree distance (m)",d.cast$variable)

d.cast$FIRE_SEV.cat <- gsub("control","Reference",d.cast$FIRE_SEV.cat)
d.cast$FIRE_SEV.cat <- gsub("high.sev","High severity",d.cast$FIRE_SEV.cat)

d.cast$fire.cat <- paste(d.cast$Fire,d.cast$topoclim.cat,sep=": ")

plots <- list()

for(i in 1:length(unique(d.cast$variable))) {
  
  var <- unique(d.cast$variable)[i]
  
  d.var <- d.cast[d.cast$variable == var,]
  
  plots[[var]] <-ggplot(d.var,aes(x=fire.cat,y=mid,color=FIRE_SEV.cat)) +
    geom_point(position=position_dodge(width=.5),size=2) +
    geom_errorbar(aes(ymin=low,ymax=high),position=position_dodge(width=.5),width=0,size=1) +
    theme_bw(12) +
    theme(axis.title.y = element_text(size=8))
  

  
  if(i == length(unique(d.cast$variable))) {
    
    plots[[var]] <- plots[[var]] +
      theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
      labs(x="Fire and topoclimatic category",y=var,color="Plot type")
    
  } else {
    
    plots[[var]] <- plots[[var]] +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      labs(y=var,color="Plot type")
  }
  


}

a <- ggplot_gtable(ggplot_build(plots[[1]]))
b <-  ggplot_gtable(ggplot_build(plots[[2]]))
c <-  ggplot_gtable(ggplot_build(plots[[3]]))
d <-  ggplot_gtable(ggplot_build(plots[[4]]))

maxWidth = unit.pmax(a$widths[2:3], b$widths[2:3], c$widths[2:3], d$widths[2:3])

a$widths[2:3] <- maxWidth
b$widths[2:3] <- maxWidth
c$widths[2:3] <- maxWidth
d$widths[2:3] <- maxWidth




Cairo(file=paste0("../Figures/Fig1_anom_normal_abiotic_",Sys.Date(),".png"),width=1500,height=2000,ppi=200,res=200,dpi=200) 
grid.arrange(a,b,c,d,ncol=1,heights=c(1,1,1,1.7))
dev.off()




### 9. Plot histogram of seed tree distance ####

d.plot.highsev <- d.plot[d.plot.c$FIRE_SEV %in% c(4,5),]

ggplot(d.plot.highsev,aes(seed_tree_distance_general)) +
  geom_histogram() +
  facet_wrap(~Fire)

library(plyr)

a <- ddply(d.plot.highsev,~Fire,summarise,mean=mean(seed_tree_distance_general),sd=sd(seed_tree_distance_general))

ggplot(a,aes(x=Fire,y=mean)) +
  geom_point() +
  geom_errorbar(ymax=a$mean+a$sd,ymin=a$mean-a$sd) +
  scale_y_continuous(limits=c(0,50))
















#### 10. Plot-level analysis with GLM; also run randomForest to get importance scores ####




### Center data frame ###

d <- d.plot.c


vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all","dominant_shrub_ht_cm","tallest_ht_cm","prop.regen.pinus.old","prop.regen.pinus.all","prop.regen.shade.old","prop.regen.hdwd.old","prop.regen.hdwd.old","prop.regen.conif.old","regen.count.broader.old")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z", "def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","def.post","aet.post") # removed snow, adult.ba.agg
d <- d[complete.cases(d[,vars.focal]),]

d.c <- center.df(d,vars.leave)[["centered.df"]]
d.center.dat <- center.df(d,vars.leave)[["center.data"]]

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

d.c$ppt.post_c.sq <- d.c$ppt.post_c^2
d.c$tmean.post_c.sq <- d.c$tmean.post_c^2

d.c$def.post_c.sq <- d.c$def.post_c^2
d.c$aet.post_c.sq <- d.c$aet.post_c^2

d.c$ppt.post.min_c.sq <- d.c$ppt.post.min_c^2

d.c$def.post.max_c.sq <- d.c$def.post.max_c^2
d.c$aet.post.min_c.sq <- d.c$aet.post.min_c^2

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

d.c.modfit <- d.c # because the model fitting uses its own d.c




## Set response options and loop through them

sp.opts <- c("ABCO","SHADE.ALLSP","PINUS.ALLSP","HDWD.ALLSP","PIPJ")
sp.opts <- c("SHADE.ALLSP","PINUS.ALLSP","HDWD.ALLSP")
cover.opts <- c("COV.SHRUB","COV.GRASS","COV.FORB")

# All
ht.opts <- c("HT.HDWD.ALLSP","HT.PINUS.ALLSP","HT.SHADE.ALLSP")
htabs.opts <- c("HTABS.PIPJ","HTABS.SHRUB","HTABS.ABCO","HTABS.PSME","HTABS.PILA","HTABS.QUKE","HTABS.CADE27","HTABS.HDWD.ALLSP","HTABS.PINUS.ALLSP","HTABS.SHADE.ALLSP")
htabs.opts <- c("HTABS.PIPJ","HTABS.SHRUB","HTABS.ABCO","HTABS.HDWD.ALLSP","HTABS.PINUS.ALLSP","HTABS.SHADE.ALLSP")

prop.opts <- c("PROP.CONIF","PROP.PINUS","PROP.SHADE")
prop.opts <- NULL
htabs.opts <- NULL 
count.opts <- c("COUNT.PINUS.ALLSP","COUNT.SHADE.ALLSP","COUNT.HDWD.ALLSP")

# # sp grps
# ht.opts <- c("HT.HDWD.ALLSP","HT.PINUS.ALLSP","HT.SHADE.ALLSP")
# htabs.opts <- c("HTABS.HDWD.ALLSP","HTABS.PINUS.ALLSP","HTABS.SHRUB","HTABS.SHADE.ALLSP")
# prop.opts <- c("PROP.CONIF","PROP.PINUS")
# 



resp.opts <- c(count.opts,prop.opts,htabs.opts,sp.opts,cover.opts,ht.opts)
resp.opts <- c(sp.opts,cover.opts)


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

fit.mods <- list()





sink("run_output.txt")
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
  }  else if(sp %in% count.opts) {
    sp.name <- substr(sp,7,100)
    d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp.name,]
    d.sp.curr.plt <- d.sp[d.sp$species==sp.name,]
  } else if(sp == "PROP.PINUS") {
    d.sp.curr.plt <- d.sp[d.sp$species=="PINUS.ALLSP",]
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="PINUS.ALLSP",]
    d.sp.all.plt <- d.sp[d.sp$species=="CONIF.ALLSP",]
    d.sp.curr.plt$regen.count.broader.old <- d.sp.all.plt$regen.count.old
    
  } else if(sp == "PROP.SHADE") {
    d.sp.curr.plt <- d.sp[d.sp$species=="SHADE.ALLSP",]
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="SHADE.ALLSP",]
    d.sp.all.plt <- d.sp[d.sp$species=="CONIF.ALLSP",]
    d.sp.curr.plt$regen.count.broader.old <- d.sp.all.plt$regen.count.old
    
  } else if(sp == "PROP.CONIF") {
    d.sp.curr.plt <- d.sp[d.sp$species=="CONIF.ALLSP",]
    d.sp.curr.agg <- d.sp.2[d.sp.2$species=="CONIF.ALLSP",]
    d.sp.all.plt <- d.sp[d.sp$species=="ALL",]
    d.sp.curr.plt$regen.count.broader.old <- d.sp.all.plt$regen.count.old
    
  }  else  {
    
    d.sp.curr.agg <- d.sp.2[d.sp.2$species==sp,]
    d.sp.curr.plt <- d.sp[d.sp$species==sp,]
  }
  
  #d <- merge(d.plot.c,d.plot.3,by=c("Fire","topoclim.cat")) # this effectively thins to plots that belong to a topoclimate category that has enough plots in it

  #add the data for the current species response to the plot data
  d.c <- merge(d.c.modfit,d.sp.curr.plt,by=c("Regen_Plot"))
  
  # names(d.sp.curr.agg) <- paste0(names(d.sp.curr.agg),".agg")
  # d <- merge(d,d.sp.curr.agg,by.x=c("Fire","topoclim.cat"),by.y=c("Fire.agg","topoclim.cat.agg")) # add the aggregated species data (from this we just want adult BA from the control plot)
  

  d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
  d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)
  
  d.c$regen.count.all.int <- ceiling(d.c$regen.count.all)
  d.c$regen.count.old.int <- ceiling(d.c$regen.count.old)
  
  ## Test whether seedling taller than shrub
  d.c$seedl.taller <- d.c$tallest_ht_cm > d.c$dominant_shrub_ht_cm
  d.c[d.c$regen.presab.old == FALSE,"seedl.taller"] <- 0 # where there is no seedling, it is not taller than shrub
  
  
  
  
  if(sp %in% cover.opts) {
    
    sp.cov <- substr(sp,5,100)
    sp.cov <- paste0(sp.cov,".pt")
    
    d.c$response.var <- d.c[,sp.cov]
    
  } else if(sp %in% ht.opts){
    
    # no longer looking only at plots where sp existed in the first place d.c <- d.c[d.c$regen.presab.old == TRUE,] # this is where we select whether we want all plots where the species was present or just old seedlings
    d.c$response.var <- d.c$seedl.taller
    
  } else if(sp == "HTABS.SHRUB") {
    d.c$response.var <- d.c$dominant_shrub_ht_cm
    d.c <- d.c[!is.na(d.c$response.var),]
  } else if(sp %in% htabs.opts) {
    d.c <- d.c[d.c$regen.presab.old == TRUE, ] # this is where we select whether we want all plots where the species was present or just old seedlings
    d.c$response.var <- d.c$tallest_ht_cm
  } else if(sp %in% count.opts) {
    d.c <- d.c[d.c$regen.presab.old == TRUE, ] # this is where we select whether we want all plots where the species was present or just old seedlings
    d.c$response.var <- log(d.c$regen.count.old)
  } else if(sp %in% prop.opts) {
    
    response.var <- cbind(d.c$regen.count.old,d.c$regen.count.broader.old-d.c$regen.count.old)
    response.var <- response.var * 3
    
    d.c$response.var <- response.var
    
  } else  {
    
    d.c$response.var <- d.c$regen.presab.old.01
    
  }
  
  
    
    if(TRUE) {
      
      
      formulas <- list()
      
      
      ### No seed tree ###
      
      formulas[["n0.a0"]] <- formula(response.var ~ 1)
      
      ## PPT
      
      formulas[["n0.aP"]] <- formula(response.var ~ diff.norm.ppt.z_c)
      formulas[["n0.aP2"]] <- formula(response.var ~ diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      
      
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
      


      ##DEF

      formulas[["n0.aD"]] <- formula(response.var ~ diff.norm.def.z_c)
      formulas[["n0.aD2"]] <- formula(response.var ~ diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)



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


      formulas[["n0.aA"]] <- formula(response.var ~ diff.norm.aet.z_c)
      formulas[["n0.aA2"]] <- formula(response.var ~ diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)


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
      formulas[["n0.aPmin"]] <- formula(response.var ~ diff.norm.ppt.min.z_c)
      formulas[["n0.aP2min"]] <- formula(response.var ~ diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)




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

       

      ##DEF

      formulas[["n0.aDmax"]] <- formula(response.var ~ diff.norm.def.max.z_c)
      formulas[["n0.aD2max"]] <- formula(response.var ~ diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

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

      formulas[["n0.aAmin"]] <- formula(response.var ~ diff.norm.aet.min.z_c)
      formulas[["n0.aA2min"]] <- formula(response.var ~ diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)


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


      ### Add seed tree ###      
      if(sp %in% c(sp.opts,count.opts)) {
      

        formulas[["n0s.a0"]] <- formula(response.var ~ seed_tree_distance_general_c + 1)

        ## PPT
        
        formulas[["n0s.aP"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.z_c)
        formulas[["n0s.aP2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
        

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
        
        
        
  
        ##DEF

        formulas[["n0s.aD"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.z_c)
        formulas[["n0s.aD2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.z_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)



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


        formulas[["n0s.aA"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.z_c)
        formulas[["n0s.aA2"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)


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
        
        formulas[["n0s.aPmin"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.min.z_c)
        formulas[["n0s.aP2min"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        
        
        
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
        
        
        
        ##DEF

        formulas[["n0s.aDmax"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.max.z_c)
        formulas[["n0s.aD2max"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

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

        formulas[["n0s.aAmin"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.min.z_c)
        formulas[["n0s.aA2min"]] <- formula(response.var ~ seed_tree_distance_general_c + diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)

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
      
      
      ### With solar rad ###

      formulas[["n0r.a0"]] <- formula(response.var~ rad.march_c +  1)

      ## PPT

      formulas[["n0r.aP"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.z_c)
      formulas[["n0r.aP2"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)


      formulas[["nPr.a0"]] <- formula(response.var~ rad.march_c +  ppt.normal_c)
      formulas[["nPr.aP"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c)
      formulas[["nPr.aP2"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
      formulas[["pPr"]] <- formula(response.var~ rad.march_c +  ppt.post_c)
      formulas[["nP2r.a0"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + ppt.normal_c.sq)
      formulas[["nP2r.aP"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
      formulas[["nP2r.aP2"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2r"]] <- formula(response.var~ rad.march_c +  ppt.post_c + ppt.post_c.sq)
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
      formulas[["pDr"]] <- formula(response.var~ rad.march_c +  def.post_c)
      formulas[["nD2r.a0"]] <- formula(response.var~ rad.march_c +  def.normal_c + def.normal_c.sq)
      formulas[["nD2r.aD"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2r.aD2"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
      formulas[["pD2r"]] <- formula(response.var~ rad.march_c +  def.post_c + def.post_c.sq)
      formulas[["nDr.aDni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c)
      formulas[["nDr.aD2ni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
      formulas[["nD2r.aDni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
      formulas[["nD2r.aD2ni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



      ##AET


      formulas[["n0r.aA"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.z_c)
      formulas[["n0r.aA2"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)


      formulas[["nAr.a0"]] <- formula(response.var~ rad.march_c +  aet.normal_c)
      formulas[["nAr.aA"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.z_c)
      formulas[["nAr.aA2"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["pAr"]] <- formula(response.var~ rad.march_c +  aet.post_c)
      formulas[["nA2r.a0"]] <- formula(response.var~ rad.march_c +  aet.normal_c + aet.normal_c.sq)
      formulas[["nA2r.aA"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2r.aA2"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)
      formulas[["pA2r"]] <- formula(response.var~ rad.march_c +  aet.post_c + aet.post_c.sq)
      formulas[["nAr.aAni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.z_c)
      formulas[["nAr.aA2ni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
      formulas[["nA2r.aAni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.z_c + aet.normal_c.sq)
      formulas[["nA2r.aA2ni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)



      ## PPT
      formulas[["n0r.aPmin"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.min.z_c)
      formulas[["n0r.aP2min"]] <- formula(response.var~ rad.march_c +  diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)




      formulas[["nPr.aPmin"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c)
      formulas[["nPr.aP2min"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["pPminr"]] <- formula(response.var~ rad.march_c +  ppt.post.min_c)
      formulas[["nP2r.aPmin"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2r.aP2min"]] <- formula(response.var~ rad.march_c +  ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
      formulas[["pP2minr"]] <- formula(response.var~ rad.march_c +  ppt.post.min_c + ppt.post.min_c.sq)
      formulas[["nPr.aPminni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c)
      formulas[["nPr.aP2minni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
      formulas[["nP2r.aPminni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
      formulas[["nP2r.aP2minni"]] <- formula(response.var~ rad.march_c +  ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)



      ##DEF

      formulas[["n0r.aDmax"]] <- formula(response.var~ rad.march_c +  diff.norm.def.max.z_c)
      formulas[["n0r.aD2max"]] <- formula(response.var~ rad.march_c +  diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

      formulas[["nDr.aDmax"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c)
      formulas[["nDr.aD2max"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["pDmaxr"]] <- formula(response.var~ rad.march_c +  def.post.max_c)
      formulas[["nD2r.aDmax"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2r.aD2max"]] <- formula(response.var~ rad.march_c +  def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
      formulas[["pD2maxr"]] <- formula(response.var~ rad.march_c +  def.post.max_c + def.post.max_c.sq)
      formulas[["nDr.aDmaxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c)
      formulas[["nDr.aD2maxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
      formulas[["nD2r.aDmaxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
      formulas[["nD2r.aD2maxni"]] <- formula(response.var~ rad.march_c +  def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



      ##AET

      formulas[["n0r.aAmin"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.min.z_c)
      formulas[["n0r.aA2min"]] <- formula(response.var~ rad.march_c +  diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)


      formulas[["nAr.aAmin"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.min.z_c)
      formulas[["nAr.aA2min"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["pAminr"]] <- formula(response.var~ rad.march_c +  aet.post.min_c)
      formulas[["nA2r.aAmin"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2r.aA2min"]] <- formula(response.var~ rad.march_c +  aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)
      formulas[["pA2minr"]] <- formula(response.var~ rad.march_c +  aet.post.min_c + aet.post.min_c.sq)
      formulas[["nAr.aAminni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.min.z_c)
      formulas[["nAr.aA2minni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
      formulas[["nA2r.aAminni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c.sq)
      formulas[["nA2r.aA2minni"]] <- formula(response.var~ rad.march_c +  aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)


      ### Add seed tree ###
      if(sp %in% c(sp.opts,count.opts)) {

        formulas[["n0sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + 1)

        ## PPT

        formulas[["n0sr.aP"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.z_c)
        formulas[["n0sr.aP2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)

        formulas[["nPsr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c)
        formulas[["nPsr.aP"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c)
        formulas[["nPsr.aP2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq)
        formulas[["pPsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.post_c)
        formulas[["nP2sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + ppt.normal_c.sq)
        formulas[["nP2sr.aP"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq)
        formulas[["nP2sr.aP2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq)
        formulas[["pP2sr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.post_c + ppt.post_c.sq)
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
        formulas[["pDsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.post_c)
        formulas[["nD2sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + def.normal_c.sq)
        formulas[["nD2sr.aD"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c.sq)
        formulas[["nD2sr.aD2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)
        formulas[["pD2sr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.post_c + def.post_c.sq)
        formulas[["nDsr.aDni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c)
        formulas[["nDsr.aD2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq)
        formulas[["nD2sr.aDni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c.sq)
        formulas[["nD2sr.aD2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.z_c + def.normal_c + diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq)



        ##AET

        formulas[["n0sr.aA"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.z_c)
        formulas[["n0sr.aA2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.z_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)

        formulas[["nAsr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c)
        formulas[["nAsr.aA"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c)
        formulas[["nAsr.aA2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
        formulas[["pAsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.post_c)
        formulas[["nA2sr.a0"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + aet.normal_c.sq)
        formulas[["nA2sr.aA"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq)
        formulas[["nA2sr.aA2"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)
        formulas[["pA2sr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.post_c + aet.post_c.sq)
        formulas[["nAsr.aAni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c)
        formulas[["nAsr.aA2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq)
        formulas[["nA2sr.aAni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c.sq)
        formulas[["nA2sr.aA2ni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.z_c + aet.normal_c + diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq)


        
        ## PPT

        formulas[["n0sr.aPmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.min.z_c)
        formulas[["n0sr.aP2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.ppt.min.z_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)


        formulas[["nPsr.aPmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c)
        formulas[["nPsr.aP2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        formulas[["pPminsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.post.min_c)
        formulas[["nP2sr.aPmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c.sq)
        formulas[["nP2sr.aP2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c*diff.norm.ppt.min.z_c + ppt.normal_c*diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)
        formulas[["pP2minsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.post.min_c + ppt.post.min_c.sq)
        formulas[["nPsr.aPminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c)
        formulas[["nPsr.aP2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq)
        formulas[["nP2sr.aPminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c.sq)
        formulas[["nP2sr.aP2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + ppt.normal_c + diff.norm.ppt.min.z_c + ppt.normal_c + diff.norm.ppt.min.z_c.sq + diff.norm.ppt.min.z_c.sq + ppt.normal_c.sq)



        ##DEF

        formulas[["n0sr.aDmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.def.max.z_c)
        formulas[["n0sr.aD2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.def.max.z_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)

        formulas[["nDsr.aDmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c)
        formulas[["nDsr.aD2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
        formulas[["pDmaxsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.post.max_c)
        formulas[["nD2sr.aDmax"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c.sq)
        formulas[["nD2sr.aD2max"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c*diff.norm.def.max.z_c + def.normal_c*diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)
        formulas[["pD2maxsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.post.max_c + def.post.max_c.sq)
        formulas[["nDsr.aDmaxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c)
        formulas[["nDsr.aD2maxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq)
        formulas[["nD2sr.aDmaxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c.sq)
        formulas[["nD2sr.aD2maxni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + def.normal_c + diff.norm.def.max.z_c + def.normal_c + diff.norm.def.max.z_c.sq + diff.norm.def.max.z_c.sq + def.normal_c.sq)



        ##AET

        formulas[["n0sr.aAmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.min.z_c)
        formulas[["n0sr.aA2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + diff.norm.aet.min.z_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)

        formulas[["nAsr.aAmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c)
        formulas[["nAsr.aA2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
        formulas[["pAminsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.post.min_c)
        formulas[["nA2sr.aAmin"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c.sq)
        formulas[["nA2sr.aA2min"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c*diff.norm.aet.min.z_c + aet.normal_c*diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)
        formulas[["pA2minsr"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.post.min_c + aet.post.min_c.sq)
        formulas[["nAsr.aAminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c)
        formulas[["nAsr.aA2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq)
        formulas[["nA2sr.aAminni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c.sq)
        formulas[["nA2sr.aA2minni"]] <- formula(response.var~ rad.march_c +  seed_tree_distance_general_c + aet.normal_c + diff.norm.aet.min.z_c + aet.normal_c + diff.norm.aet.min.z_c.sq + diff.norm.aet.min.z_c.sq + aet.normal_c.sq)
      }
      
      
    #   
    # ### Take all these formulas and add shrub cover to them (as an additional formula)
    #   
    #   if(sp %in% c(sp.opts,ht.opts)) {
    #   
    #     for(k in 1:length(formulas)) {
    #       
    #       form.name <- names(formulas)[k]
    #       name.parts <- strsplit(form.name,split=".",fixed=TRUE)[[1]]
    #       
    #       form <- formulas[[k]]
    #       form.char <- as.character(form)[3]
    #       form.char <- paste("response.var ~ ",form.char,"SHRUB_c",sep=" + ")
    #       form.new <- as.formula(form.char)
    #       
    #       form.name.new <- paste0(name.parts[1],"c")
    #       form.name.new <- paste(form.name.new,name.parts[2],sep=".")
    #       
    #       formulas[[form.name.new]] <- form.new
    #       
    #     }
    #   
    #   }
      
      
      
      
      
      
      
      
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

    norm.search.opts <- c(Pmean = "n(P|0)2?s?r?c?\\.a0",
                            Dmean = "n(D|0)2?s?r?c?\\.a0",
                            Amean = "n(A|0)2?s?r?c?\\.a0",
                            Pmin = "n(P|0)2?s?r?c?\\.a0",
                            Dmax = "n(D|0)2?s?r?c?\\.a0",
                            Amin = "n(A|0)2?s?r?c?\\.a0"
                            )


    anom.search.opts <- c(Pmean = "aP2?(ni)?$",
                          Dmean = "aD2?(ni)?$",
                          Amean = "aA2?(ni)?$",
                          Pmin = "aP2?min(ni)?$",
                          Dmax = "aD2?max(ni)?$",
                          Amin = "aA2?min(ni)?$"
                          )


    post.search.opts <- c(Pmean = "pP2?s?r?c?$",
                          Dmean = "pD2?s?r?c?$",
                          Amean = "pA2?s?r?c?$",
                          Pmin = "pP2?mins?r?c?$",
                          Dmax = "pD2?maxs?r?c?$",
                          Amin = "pA2?mins?r?c?$"
                          )


    
    
    d.maes.anoms.sp <- data.frame()
    
    
    for(i in 1:length(anom.search.opts)) {
      
      anom.name <- names(norm.search.opts)[i]
      norm.search <- norm.search.opts[i]
      anom.search <- anom.search.opts[i]
      null.names <- c("n0.a0","n0s.a0","n0r.a0","n0sr.a0","n0c.a0","n0sc.a0","n0rc.a0","n0src.a0")
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
      
      
  
      
      
      #what is the best post model for that normal, sticking with whichever anom we have from the normal
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
    
    
    
    
    vars <- c("ppt.normal_c","ppt.normal_c.sq","diff.norm.ppt.z_c","diff.norm.ppt.z_c.sq","tmean.normal_c","tmean.normal_c.sq",
              "diff.norm.tmean.z_c","diff.norm.tmean.z_c.sq",
              "aet.normal_c","aet.normal_c.sq","diff.norm.aet.z_c","diff.norm.aet.z_c.sq","def.normal_c","def.normal_c.sq",
              "diff.norm.def.z_c","diff.norm.def.z_c.sq",
              "seed_tree_distance_general_c","rad.march_c","SHRUB_c"
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
      SHRUB_c = mid.val["SHRUB_c"]
      
    )
    
    for(j in 1:nrow(d.maes.anoms.sp)) {
      
      d.maes.anoms.sp.row <- d.maes.anoms.sp[j,]
      best.anom.mod <- d.maes.anoms.sp.row$best.anomaly.mod
      best.anom.normal.mod <- d.maes.anoms.sp.row$best.anom.normal
      
  
      ##predict to new data
      
      if(sp %in% c(cover.opts,prop.opts)) {
      
        
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
        
  
      } else if (sp %in% c(htabs.opts,count.opts)) {
        
        
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
        
      } else { # presab, ht
        
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
      

      #cat("###")
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
      
      
      
      ## get fitted and observed
      
      d.c.complete <- d.c[!is.na(d.c$response.var),]
      
      if(sp %in% c(sp.opts,ht.opts)) {
        
        mod.anom <- glm(formulas[[best.anom.mod]],data=d.c.complete,family="binomial")
        mod.norm <- glm(formulas[[best.anom.normal.mod]],data=d.c.complete,family="binomial")
        

        
      } else if(sp %in% c(cover.opts,prop.opts)) {
        
        mod.anom <- betareg(formulas[[best.anom.mod]],data=d.c.complete)
        mod.norm <- betareg(formulas[[best.anom.normal.mod]],data=d.c.complete)
        
      } else { #htabs, count
        
        mod.anom <- glm(formulas[[best.anom.mod]],data=d.c.complete,family="gaussian")
        mod.norm <- glm(formulas[[best.anom.normal.mod]],data=d.c.complete,family="gaussian")
        
      }
      
      fit.mods[[paste0(sp,"_",best.anom.mod)]] <- mod.anom
      fit.mods[[paste0(sp,"_",best.anom.normal.mod)]] <- mod.norm
      

      fit.anom <- as.data.frame(predict(mod.anom))
      fit.norm <- as.data.frame(predict(mod.norm))
      names(fit.anom) <- names(fit.norm) <- c("fitted")
      fit.anom$type <- "anom"
      fit.norm$type <- "norm"
      
      d.c.complete <- d.c[complete.cases(d.c$response.var),]
      
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







#### 11. Plot overall exploration plot ####

### Prep  ###

#for each sp-rad combo, find which was the best anomaly, which anomalies were better than their corresponding normals, and plot symbols indicating
d.maes.anoms$anom.better <- d.maes.anoms$best.anom.mae < d.maes.anoms$best.anom.normal.mae
d.maes.anoms$anom.improvement <-  d.maes.anoms$best.anom.normal.mae - d.maes.anoms$best.anom.mae

d.maes.anoms$best.of.species <- ""
d.maes.anoms$most.improved <- ""

d.maes.anoms.short <- d.maes.anoms[,c("sp","anom.name","anom.better","anom.improvement","best.of.species","most.improved")]

pred.dat.comb <- merge(pred.dat,d.maes.anoms.short,all.x=TRUE,by.x=c("sp","anom"),by.y=c("sp","anom.name"))
pred.dat.comb$anom.improvement <- round(pred.dat.comb$anom.improvement,2)
pred.dat.comb[pred.dat.comb$anom.improvement < 0 , "anom.improvement"] <- NA


#Compute total MAE across all species, to see which is the best #
d.mae.agg <- aggregate(d.maes.anoms[,c("best.anom.mae","best.anom.normal.mae","anom.improvement")],by=list(d.maes.anoms$anom.name),FUN=mean,na.rm=TRUE)



### Plot counterfactuals ###

## get rid of mid and get rid of normal predictions
pred.dat.plotting <- pred.dat.comb[pred.dat.comb$norm.level %in% c("low","high"),]
#pred.dat.plotting <- pred.dat.comb[pred.dat.comb$norm.level %in% c("mid"),]
pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$type == "anom",]

pred.dat.plotting$anom.improvement <- round(pred.dat.plotting$anom.improvement,2)

#pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$sp=="PILA",]
#pred.dat.plotting <- pred.dat.plotting[pred.dat.plotting$anom=="Amin",]


ggplot(pred.dat.plotting,aes(x=diff.norm.ppt.z_c,y=pred.mid,color=norm.level,fill=norm.level)) +
  geom_point() +
  geom_ribbon(aes(ymin=pred.low,ymax=pred.high),alpha=0.3,color=NA) +
  facet_wrap(anom~sp,nrow=6,scales="free_y") +
  geom_text(aes(0,0.8,label=anom.improvement),size=3,color="black")




#### 12. Plot pub-quality counterfactual fits for the 9 main responses ####
pred.dat.comb <- merge(pred.dat,d.maes.anoms.short,all.x=TRUE,by.x=c("sp","anom"),by.y=c("sp","anom.name"))


## uncenter the predictor vars in the predictions data frame

pred.dat.comb <- as.data.frame(pred.dat.comb)

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
pred.dat.comb <- pred.dat.comb[anom=="Pmin"]

### for those normal models that are not null, include low and high norm levels; for those that are, include only mid
pred.dat.comb$norm.beg <- substr(pred.dat.comb$mod,1,2)

pred.dat.comb <- pred.dat.comb[(norm.beg != "n0" & norm.level %in% c("low","high")) | (norm.beg == "n0" & norm.level == "mid")]


pred.dat.comb[pred.dat.comb$sp %in% cover.opts,c("pred.mid","pred.low","pred.high")] <- 100 * pred.dat.comb[pred.dat.comb$sp %in% cover.opts,c("pred.mid","pred.low","pred.high")]



levels(pred.dat.comb$norm.level)
levels(pred.dat.comb$norm.level) <- c("High","Low","High and low")




plot.cats <- c("presab",
               #"ht",
               "cov")
plot.sps <- list(c("PINUS.ALLSP","SHADE.ALLSP"),
                 #c("HT.PINUS.ALLSP","HT.SHADE.ALLSP","HT.HDWD.ALLSP"),
                 c("COV.SHRUB","COV.GRASS"))

p <- list()

for(i in 1:length(plot.cats)) {
  
  plot.sp <- plot.sps[[i]]

  pred.dat.plotting <- pred.dat.comb[type == "anom" & sp %in% plot.sp,]
  
  pred.dat.plotting$sp <- factor(pred.dat.plotting$sp,plot.sp)
  
  if(plot.cats[[i]] == "presab") {
    ylab <- "Proportion of plots\nwith regeneration"
    levels(pred.dat.plotting$sp) <- c("Pines","Shade-tolerant\nconifers","Broadleaved trees")
  } else if(plot.cats[[i]] == "ht") {
    ylab <- "Proportion plots where\ndominant in height"
  } else if(plot.cats[[i]] == "cov") {
    ylab <- "Percent cover"
    levels(pred.dat.plotting$sp) <- c("Shrubs","Graminoids","Forbs")
    

    
  }
  

  p[[i]] <- ggplot(pred.dat.plotting,aes(x=diff.norm.ppt.min.z,y=pred.mid,color=norm.level,fill=norm.level)) +
    geom_line(size=1.5) +
    geom_ribbon(aes(ymin=pred.low,ymax=pred.high),alpha=0.3,color=NA) +
    facet_wrap(~sp) +
    theme_bw(16) +
    labs(x="Postfire minimum precipitation anomaly (SD)",y=ylab,color="Normal\nprecipitation",fill="Normal\nprecipitation") +
    scale_color_manual(values=c("High"="turquoise4","Low"="darkorange1","High and low" = "gray44")) +
    scale_fill_manual(values=c("High"="turquoise4","Low"="darkorange1","High and low" = "gray44")) +
    theme(panel.grid.minor = element_blank(),strip.background = element_blank(), panel.border = element_rect(colour = "black",size=0.6), strip.text = element_text(size = 16,vjust=0)) +
    theme(plot.margin = unit(c(-.1,0.5,0,0.5), "cm"))
  
  
  if(i == 2) {
    p[[i]] <- p[[i]] + theme(legend.position=c(0.98,0.95),legend.justification=c(1,1),legend.background = element_rect(fill="transparent")) +
      theme(plot.margin = unit(c(-.1,0.5,.50,0.5), "cm"))
  } else {
    p[[i]] <- p[[i]] + guides(fill=FALSE,color=FALSE) + labs(x="")
  }

  
}



a <- ggplot_gtable(ggplot_build(p[[1]]))
b <-  ggplot_gtable(ggplot_build(p[[2]]))

maxWidth = unit.pmax(a$widths[2:3], b$widths[2:3])

a$widths[2:3] <- maxWidth
b$widths[2:3] <- maxWidth


tiff(file=paste0("../Figures/FigX_prediction_plots_",Sys.Date(),".tiff"),width=1600,height=1600,res=200) 
grid.arrange(a,b,ncol=1,heights=c(1,1))
dev.off()


#### Compute what the normal precip hypothetical values are ####

ppt.low <- low.val.fun("ppt.normal_c")
ppt.high <- high.val.fun("ppt.normal_c")

ppt.low <- (ppt.low * d.center.dat[d.center.dat$var=="ppt.normal_c","var.sd"]) + d.center.dat[d.center.dat$var=="ppt.normal_c","var.mean"]
ppt.high <- (ppt.high * d.center.dat[d.center.dat$var=="ppt.normal_c","var.sd"]) + d.center.dat[d.center.dat$var=="ppt.normal_c","var.mean"]




#### 14. Table of cross-validation errors ####
d.maes.anoms <- as.data.table(d.maes.anoms)
cv.errors <- d.maes.anoms[,.(Species=sp,Variable=anom.name,Baseline=best.anom.normal.mae,Anomaly=best.anom.mae)]

cv.err.melt <- melt(cv.errors,id.vars=c("Species","Variable"),measure.vars=c("Baseline","Anomaly"),variable.name="type")

cv.err.cast <- dcast(cv.err.melt,Species~Variable+type,value.var="value",fun=mean)

cv.err.cast[,-1] <- round(cv.err.cast[,-1],2)





#### 15. Table of normal vs. post MAEs ####
library(reshape)

##calc how much better the post is than the normal
d.maes.anoms$post.v.norm <- d.maes.anoms$best.normal.nonull.mae - d.maes.anoms$best.post.mae
d.maes.anoms$post.v.norm[d.maes.anoms$post.v.norm < 0] <- NA

##make a table: species x anom, where value is the normal-post improvement

# positive means post has lower error
post.table <- cast(d.maes.anoms,anom.name~sp,value='post.v.norm')
post.table


#### 16. Table of model coefs ####

## For each anomaly type and model type (anom or norm), extract the coefs

coefs <- data.frame()

for(i in 1:nrow(d.maes.anoms)) {
  
  cv.row <- d.maes.anoms[i,]
  sp <- cv.row$sp
  norm.mod <- cv.row$best.anom.normal
  anom.mod <- cv.row$best.anomaly.mod
  
  m.norm <- fit.mods[[paste0(sp,"_",norm.mod)]]
  m.anom <- fit.mods[[paste0(sp,"_",anom.mod)]]
  
  if(sp %in% c(sp.opts,ht.opts)) {
    
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

dummy.df <- data.table("sp"="a","variable"="a","type"="a","(Intercept)"=1,"normal_c"=1,"normal_c.sq"=1,"seed_tree_distance_general_c"=1,"rad.march_c"=1,"diff.norm.z_c"=1,"diff.norm.z_c.sq"=1,"normal_c:diff.norm.z_c"=1,"normal_c:diff.norm.z_c.sq"=1,"(phi)"=1)
coefs.cast <- rbind.fill(coefs.cast,dummy.df)
coefs.cast <- as.data.table(coefs.cast)

coefs.cast <- coefs.cast[,c("sp","variable","type",Intercept="(Intercept)","normal_c","normal_c.sq",seed.tree="seed_tree_distance_general_c","rad.march_c","diff.norm.z_c","diff.norm.z_c.sq","normal_c:diff.norm.z_c","normal_c:diff.norm.z_c.sq","(phi)"),with=FALSE]





#### 13. For each species and rad group, plot fitted vs. observed, for normal and anom side by side ####

fit.dat$resid <- fit.dat$response.var-fit.dat$fitted


## PPT min only
fit.dat.ppt <- as.data.table(fit.dat[fit.dat$anom == "Pmin",])

## summarize by fire, sp

fit.dat.ppt.fire <- fit.dat.ppt[,list(fitted=mean(fitted),observed=mean(response.var),anom.var=mean(diff.norm.ppt.min.z_c),resid=mean(resid)),by=.(Fire,type,sp)]
# fit.dat.ppt.fire <- fit.dat.ppt[,list(fitted=fitted,observed=response.var,anom.var=diff.norm.ppt.min.z_c,Fire,type,sp)]

fit.dat.ppt.fire[sp %in% cover.opts,c("observed","fitted")] <- 100* fit.dat.ppt.fire[sp %in% cover.opts,c("observed","fitted")]

fit.dat.ppt.fire$sp <- as.factor(fit.dat.ppt.fire$sp)


fit.dat.ppt.fire$sp <- factor(fit.dat.ppt.fire$sp,c("PINUS.ALLSP","SHADE.ALLSP","HDWD.ALLSP","COV.SHRUB","COV.GRASS","COV.FORB"))

levels(fit.dat.ppt.fire$sp)
levels(fit.dat.ppt.fire$sp) <- c("Pine regeneration\n(% of plots)","Shade tolerant conifer\nspecies regeneration\n(% of plots)","Broadleaved species\nregeneration\n(% of plots)","Shrubs\n(% cover)","Graminoids\n(% cover)","Forbs\n(% cover)")

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
  theme(plot.margin = unit(c(-.1,0.5,0.1,0.5), "cm"))

tiff(file=paste0("../Figures/FigX_fitted_vs_observed_",Sys.Date(),".tiff"),width=2100,height=1600,res=200) 
p
dev.off()



#### 14. Plot residuals of regen fits (of baseline non-anomaly models) by fire: are we missing anything? ####


fit.dat.plot <- fit.dat.ppt.fire[fit.dat.ppt.fire$type=="Baseline",]

ggplot(fit.dat.plot,aes(y=anom.var,x=Fire)) +
  geom_point()

fit.dat.plot <- fit.dat.plot[!is.na(fit.dat.plot$Fire),]

fit.dat.plot$Fire <- factor(fit.dat.plot$Fire,c("Bagley","Chips","Harding","Ralston","Bassetts","Moonlight","Antelope","Freds","Power","Straylor","Rich","Btu Lightning","Cub","American River"))

ggplot(fit.dat.plot,aes(x=Fire,y=resid)) +
  geom_point(size=2) +
  facet_grid(sp~.) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))






#### 6.8 Plot-level analaysis decision trees (PARTY)) ####



library(party)

### prep the data frame

sp <- "PINUS.ALLSP"

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


a <- ctree(response.var ~ ppt.normal_c + 
           diff.norm.ppt.min.z_c +
           seed_tree_distance_general_c + adult.ba.agg_c + rad.march_c, data=d.c)



#control=ctree_control(testtype=c("Univariate"),mincriterion=0.80)

plot(a)













#### 8. Cluster-level multivariate analysis ####




#d.plot.3
#d.sp.2
library(data.table)



focal.sp <- c("PIPJ","ABCO","PILA","PSME")
# focal.sp <- c("PINUS.ALLSP","HDWD.ALLSP","SHADE.ALLSP")

focal.cols <- c("Fire","topoclim.cat","species","regen.presab.old","regen.presab.all","adult.ba")
d.sp.simp <- d.sp.2[d.sp.2$species %in% focal.sp,focal.cols]
names(d.sp.simp)[4:6] <- c("r.old","r.all","a.ba")
d.sp.simp <- as.data.table(d.sp.simp)
d.sp.cast <- dcast(d.sp.simp,topoclim.cat + Fire~species,value.var=c("r.old","r.all","a.ba"))

regen.var <- "r.old"

d.all <- merge(d.plot.3,d.sp.cast,by=c("Fire","topoclim.cat"))
d.all <- d.all[(d.all$count.highsev > 4) & (d.all$count.control > 4),]

# ## wet plots only
# ppt.cutoff <- mean(d.all$ppt.normal.highsev)
# d.all <- d.all[d.all$ppt.normal.highsev > ppt.cutoff,]

# ## dry plots only
# ppt.cutoff <- mean(d.all$ppt.normal.highsev)
# d.all <- d.all[d.all$ppt.normal.highsev < ppt.cutoff,]



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
             "PSME" = "DF",
             "PIPJ" = "YP",
             "CADE27" = "IC",
             "HDWD.ALLSP" = "HDWD",
             "SHADE.ALLSP" = "SHADE",
             "PINUS.ALLSP" = "PINUS",
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




library(vegan)


d.all$SolarRadiation <- d.all$solar.rad
d.all$NormalPrecip <- d.all$ppt.norm
d.all$NormalTemp <- d.all$temp.norm

d.all$PrecipAnom <- d.all$ppt.anom
d.all$TempAnom <- d.all$temp.anom





# all non-anomaly, excluding species
cc1 <- cca(d.all.sp ~ NormalPrecip + NormalTemp + SolarRadiation,data=d.all)

# all non-anomaly, including species
cc2 <- cca(d.all.sp ~ NormalPrecip + NormalTemp + SolarRadiation + WF+DF+SP+YP,data=d.all)

# all including anomaly, including species
cc3 <- cca(d.all.sp ~ NormalPrecip + NormalTemp + PrecipAnom + SolarRadiation + WF+DF+SP+YP,data=d.all)


#plot(cc1)
#plot(cc2,choices=c(1,2))
#plot(cc2b,choices=c(1,2))
plot(cc3,choices=c(1,2))




vars.focal.nmds <- c("NormalPrecip","NormalTemp","PrecipAnom","SolarRadiation","DF","YP","WF","SP")
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








