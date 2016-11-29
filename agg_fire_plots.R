#comment from derek
library(reshape2)
library(data.table)
source("v_dodge.R")
setwd("C:/Users/DYoung/Dropbox/Research projects/Post-fire regen shared")

deg2rad <- function(deg) {
  return(deg * (3.141593/180))
}

d <- read.csv("Data/regen_clim_full.csv",header=TRUE)
d.save <- d


##take the first row for each fire and save to a file for plotting general fire locations
first.fire <- !duplicated(d$Fire)
d.first.fire <- d[first.fire,]
write.csv(d.first.fire,"regen_one_point_per_fire.csv",row.names=FALSE)




##(optional) thin to Sierra
sierra.fires <- c("STRAYLOR","CUB","RICH","DEEP","MOONLIGHT","ANTELOPE","BTU LIGHTENING","HARDING","BASSETTS","PENDOLA","AMERICAN RIVER","RALSTON","FREDS","SHOWERS","POWER")
d$in.sierra <- ifelse(d$Fire %in% sierra.fires,TRUE,FALSE)
d <- d[d$in.sierra,]



d$FORB <- d$FORBE

d$firesev <- d$FIRE_SEV # use field-assessed fire sev

d <- d[!(d$Fire %in% c("ELK","CHINABACK","ZACAM","SHOWERS")),] #,"DEEP","STRAYLOR","HARDING","MOONLIGHT","ANTELOPE")),] # rich has some uncertainty in fire and sampling year could be added later; just removed because of uncertainty in fire and sampling year
# kevin says deep, straylor, harding, showers were anomalous (showers was already removed, not sure why)
# also removing moonlight, antelope because so dry
d$Fire <- as.factor(d$Fire)

# ## remove plots with over 2000 mm precip (some BTU plots)
# d <- d[d$ppt.normal < 2000,]

d$seedtr.comb <- ifelse((is.na(d$seedtr.dist) | (d$seedtr.dist >= 150)),d$dist.to.low,d$seedtr.dist)
d$seedtr.any.comb <- ifelse(is.na(d$seed.tree.any) | (d$seed.tree.any >= 150),d$dist.to.low,d$seed.tree.any)

## calculate northness and eastness
d$northness <- cos(deg2rad(d$aspect))
d$eastness <- sin(deg2rad(d$aspect))

d$ppt.normal.sq <- d$ppt.normal^2
d$tmean.normal.sq <- d$tmean.normal^2

d$regen.presab <- ifelse(d$regen.count > 0,TRUE,FALSE)


agg.fire <- function(d,focal.vars,species.count.vars,species.presab.vars) {
  
  species.vars <- unique(c(species.count.vars,species.presab.vars))
  all.focal.vars <- unique(c(focal.vars,species.count.vars,species.presab.vars))
  fire.mean <- aggregate(d[,all.focal.vars],by=list(d$Fire,d$species),FUN=mean,na.rm=TRUE)
  fire.high <- aggregate(d[,all.focal.vars],by=list(d$Fire,d$species),quantile,probs=.75,na.rm=TRUE)
  fire.low <- aggregate(d[,all.focal.vars],by=list(d$Fire,d$species),quantile,probs=.25,na.rm=TRUE)
  
  fire.mean$summary.type <- "mean"
  fire.high$summary.type <- "high"
  fire.low$summary.type <- "low"
  
  fire.agg <- rbind(fire.mean,fire.high,fire.low)
  names(fire.agg)[1:2] <- c("Fire","species")
  
  fire.agg.m <- melt(fire.agg,id.vars=c("Fire","species","summary.type"))
  d.fire.m2 <- dcast(fire.agg.m,Fire+species+variable~summary.type)
  
  ## remove high and low from presence/absence vars (only mean is meaningful as proportion of plots that had presence)
  d.fire.m2[d.fire.m2$variable %in% species.presab.vars,c("high","low")] <- NA
  
  
  ##do not replicate the non-species variables for each species
  #get first species
  sp1 <- unique(d.fire.m2$species)[1]
  d.fire.m2.nonsp <- d.fire.m2[(d.fire.m2$species==sp1) & !(d.fire.m2$variable %in% species.vars),]
  d.fire.m2.sp <- d.fire.m2[(d.fire.m2$variable %in% species.vars),]
  d.fire.m2.sp$variable <- paste(d.fire.m2.sp$species,d.fire.m2.sp$variable,sep=".")
  
  d.fire.m2.full <- rbind(d.fire.m2.nonsp,d.fire.m2.sp)
  
  #remove species column
  d.fire.m2.full <- d.fire.m2.full[,-grep("species",names(d.fire.m2.full))]
  
  # add nplots
  d.sp1 <- d[d$sp==sp1,]
  fire.count <- aggregate(d.sp1$Regen_Plot,by=list(d.sp1$Fire),FUN=length)
  #make fake data frame to append to the main data frame
  fire.count.df <- data.frame(Fire=fire.count$Group.1,variable="nplots",high=NA,low=NA,mean=fire.count$x)
  
  d.fire.m2.full <- rbind(d.fire.m2.full,fire.count.df)
  
  return(d.fire.m2.full)
  
}


focal.vars <- c("GRASS","SHRUB","diff.norm.ppt.min.z","diff.norm.ppt.z","northness","seedtr.any.comb","ppt.normal")
species.count.vars <- c("regen.count")
species.presab.vars <- c("regen.presab")


d.highsev <- d[(d$firesev %in% c(4,5)) & (d$survey.years.post == 5),] # only medium-high and high severity; use only 5-year postfire plots. this dataset will represent post-fire regen
d.lowsev <- d[d$firesev %in% c(0,1,2,3),] # this dataset will represent pre-fire adult vegetation



d.highsev.fire <- agg.fire(d.highsev,focal.vars,species.count.vars,species.presab.vars)
d.lowsev.fire <- agg.fire(d.lowsev,focal.vars,species.count.vars,species.presab.vars)

d.highsev.fire$sev <- "high"
d.lowsev.fire$sev <- "low"
d.fire <- rbind(d.highsev.fire,d.lowsev.fire)



#define interesting regen variables
interesting <- c(
  
  "SHRUB","high","shrub","regen",
  "SHRUB","low","shrub","control",
  
  "CONIFER.regen.count","high","conif","regen",
  "CONIFER.surviving.regen.count","low","conif","control",
  "CONIFER.regen.presab","high","conif.perc","regen",
  "CONIFER.surviving.regen.presab","low","conif.perc","control",
  
  "HARDWOOD.regen.count","high","hw","regen",
  "HARDWOOD.surviving.regen.count","low","hw","control",
  "HARDWOOD.regen.presab","high","hw.perc","regen",
  "HARDWOOD.surviving.regen.presab","low","hw.perc","control",
  
  "PIPO.regen.count","high","PIPO","regen",
  "PIPO.surviving.regen.count","low","PIPO","control",
  "PIPO.regen.presab","high","PIPO.perc","regen",
  "PIPO.surviving.regen.presab","low","PIPO.perc","control",
  
  "ABCO.regen.count","high","ABCO","regen",
  "ABCO.surviving.regen.count","low","ABCO","control",
  "ABCO.regen.presab","high","ABCO.perc","regen",
  "ABCO.surviving.regen.presab","low","ABCO.perc","control"
)


interesting.mat <- matrix(interesting,ncol=4,byrow=TRUE)
interesting.df <- data.frame(interesting.mat)
names(interesting.df) <- c("variable", "sev", "gen.name","treatment")

d.interesting <- merge(d.fire,interesting.df,by=c("variable","sev"))


d.interesting



#### Plot mean, low, high val for each fire (burned and unburned side-by-side)
# 
# focal.vars.plot <- c("diff.norm.ppt.z","diff.norm.ppt.min.z","ppt.normal","PIPO.regen.count","PIPO.regen.presab","CONIFER.regen.count","CONIFER.regen.presab","HARDWOOD.regen.count","HARDWOOD.regen.presab")
# d.fire.focal <- d.fire[d.fire$variable %in% focal.vars.plot,]

library(ggplot2)

ggplot(d.interesting,aes(x=mean,y=Fire,color=treatment)) +
  geom_point(size=5) + #,position=position_dodgev(height=0.5)
  facet_grid(.~gen.name,scales="free_x") +
  theme_bw(20) +
  geom_errorbarh(aes(xmin=low,xmax=high)) #position=position_dodgev(height=0.5)
  
  

### plot non-regen (i.e., climate, environment, sample size) vars for each plot
non.regen.vars <- c("diff.norm.ppt.min.z","diff.norm.ppt.z","northness","seedtr.any.comb","ppt.normal","nplots")
d.non.regen <- d.fire[d.fire$variable %in% non.regen.vars,]

ggplot(d.non.regen,aes(x=mean,y=Fire,color=sev)) +
  geom_point(size=5) + #,position=position_dodgev(height=0.5)
  facet_grid(.~variable,scales="free_x") +
  theme_bw(20) +
  geom_errorbarh(aes(xmin=low,xmax=high)) #,position=position_dodgev(height=0.5)










