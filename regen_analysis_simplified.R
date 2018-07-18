setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")
source("regen_analysis_functions.R")

library(tidyverse)
library(car)
library(lme4)
library(betareg)



#### 1. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

# make fire names capitalized as proper nouns
d.plot$Fire <- sapply(d.plot$Fire,FUN= function(x) simpleCap(tolower(x)))

# rename an incorrectly named fire
d.plot[d.plot$Fire == "Btu Lightening","Fire"] <- "BTU Lightning"

## Omit plots with incomplete data specified in the comments
plots.exceptions <- grepl("#.*(INCOMPLETE|INCORRECT)",d.plot$NOTES)
d.plot <- d.plot[!plots.exceptions,]

d.plot$tmean.post.max <- NA # not currently using this variable

# only keep the necessary columns
d.plot <- d.plot[,c("Regen_Plot","Fire","Year.of.Fire","Easting","Northing","aspect","slope","SHRUB","FORB","GRASS","HARDWOOD","CONIFER","FIRE_SEV","BA.Live1","Year","firesev","dist.to.low","fire.abbr","X5yr","fire.year","survey.years.post","elev.m","rad.march","tmean.post","ppt.post","ppt.post.min","tmean.normal","ppt.normal","seed.tree.any","diff.norm.ppt.z","diff.norm.ppt.min.z","seed_tree_distance_general","seed_tree_distance_conifer","seed_tree_distance_hardwood","diff.norm.ppt.z","diff.norm.ppt.min.z","tmean.post","tmean.post.max","ppt.post","ppt.post.min","perc.norm.ppt","perc.norm.ppt.min","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z","def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","diff.norm.def.max.z","diff.norm.aet.min.z","def.post","aet.post","def.post.max","aet.post.min","snow.post.min","snow.normal","snow.post","diff.norm.snow.z","diff.norm.snow.min.z","dominant_shrub_ht_cm","dom.veg.all")]

# only northern Sierra Nevada fires
sierra.fires <- c("Straylor","Cub","Rich","Moonlight","Antelope","BTU Lightning","Harding","Bassetts","American River","Ralston","Freds","Power","Bagley","Chips")
d.plot <- d.plot[d.plot$Fire %in% sierra.fires,]

## Remove managed plots, plots in nonforested locations (e.g. exposed bedrock), etc.
plots.exclude <- read.csv("data_intermediate/plots_exclude.csv",header=T,stringsAsFactors=FALSE)
plot.ids.exclude <- plots.exclude[plots.exclude$Exclude != "",]$Regen_Plot
d.plot <- d.plot[!(d.plot$Regen_Plot %in% plot.ids.exclude),]

# Remove any plots > 75m from seed source
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
d.plot$FIRE_SEV.cat <- car::recode(d.plot$FIRE_SEV,"control='control';high.sev='high.sev';else=NA")
d.plot <- d.plot[!is.na(d.plot$FIRE_SEV.cat),]

# remove plots that are high severity but not surveyed in years 4-5 post-fire (but if they are control plots they can be from any year post-fire)
d.plot <- d.plot[which(!(d.plot$FIRE_SEV.cat == "high.sev" & !(d.plot$survey.years.post %in% c(4,5)))),]

# must have shrub cover
d.plot <- d.plot[!is.na(d.plot$SHRUB),]

# Dataset of individual (non-aggregated) plots
d.plot.ind <- d.plot

# first remove the variables that are not useful

# label plots as control or high sev
d.plot.ind$type <- ifelse(d.plot.ind$FIRE_SEV %in% control,"control",NA)
d.plot.ind$type <- ifelse(d.plot.ind$FIRE_SEV %in% high.sev,"highsev",d.plot.ind$type)
# get rid of plots that are neither control nor high sev
d.plot.ind <- d.plot.ind[!is.na(d.plot.ind$type),]

# filter to high-severity plots only
d.plot.ind <- d.plot.ind %>%
  filter(FIRE_SEV.cat == "high.sev")



#### 2. Center and standardize predictors ####

d <- d.plot.ind

#predictors not to transform
vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all","dominant_shrub_ht_cm","tallest_ht_cm","prop.regen.pinus.old","prop.regen.pinus.all","prop.regen.shade.old","prop.regen.hdwd.old","prop.regen.hdwd.old","prop.regen.conif.old","regen.count.broader.old")

# candidate predictor variables
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z", "def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","def.post","aet.post") # removed snow, adult.ba.agg
d <- d[complete.cases(d[,vars.focal]),]

# center predictors
d.c <- center.df(d,vars.leave)[["centered.df"]]
d.center.dat <- center.df(d,vars.leave)[["center.data"]]

d.c$SHRUB_c <- (d.c$SHRUB - mean(d.c$SHRUB))  / sd(d.c$SHRUB)

# compute quadratic forms of the centered predictors
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


#### 3. Fit models (examples) ####


### Example: PIPO seedling presence-absence
focal.sp <- "PIPJ"

d.sp.focal <- d.sp %>%
  filter(species == focal.sp)

# Join the species data to the plot data
d.mod <- left_join(d.c.highsev,d.sp.focal)

## this is the best-fit mean precipitation anomaly model as determined by the cross-validation procedure (ref manuscript)
m <- glm(regen.presab.old ~ rad.march_c + diff.norm.ppt.z_c, data=d.mod, family="binomial")
summary(m)




### Example: ABCO seedling presence-absence
focal.sp <- "ABCO"

d.sp.focal <- d.sp %>%
  filter(species == focal.sp)

d.c.highsev <- d.c %>%
  filter(FIRE_SEV.cat == "high.sev")

# Join the species data to the plot data
d.mod <- left_join(d.c.highsev,d.sp.focal)

## this is the best-fit mean precipitation anomaly model as determined by the cross-validation procedure (ref manuscript)
m <- glm(regen.presab.old ~ ppt.normal_c + ppt.normal_c.sq + seed_tree_distance_general_c + rad.march_c + diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq, data=d.mod, family="binomial")
summary(m)




### Example: shrub cover mean precipitation anomaly model
d.mod <- d.c.highsev
m <- betareg(SHRUB.pt ~ diff.norm.ppt.z_c + diff.norm.ppt.z_c.sq, data=d.mod)
summary(m)

