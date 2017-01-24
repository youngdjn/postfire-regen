setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(ggplot2)

source("regen_analysis_functions.R")

#### 0. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

# only keep the necessary columns
d.plot <- d.plot[,c("Regen_Plot","Fire","Year.of.Fire","Easting","Northing","aspect","slope","SHRUB","FORB","GRASS","HARDWOOD","CONIFER","FIRE_SEV","BA.Live1","Year","firesev","dist.to.low","fire.abbr","X5yr","fire.year","survey.years.post","elev.m","rad.march","tmean.post","ppt.post","ppt.post.min","tmean.normal","ppt.normal","seed.tree.any","diff.norm.ppt.z","diff.norm.ppt.min.z","seed_tree_distance_general","seed_tree_distance_conifer","seed_tree_distance_hardwood","diff.norm.ppt.z","diff.norm.ppt.min.z","tmean.post","ppt.post","ppt.post.min")]

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
d.plot <- d.plot[which(d.plot$seed_tree_distance_general < 50),]


# quadratic climate terms
d.plot$ppt.normal.sq <- d.plot$ppt.normal^2
d.plot$tmean.normal.sq <- d.plot$tmean.normal^2

# variable for regen presence/absence
d.sp$regen.presab.young <- ifelse(d.sp$regen.count.young > 0,TRUE,FALSE)
d.sp$regen.presab.old <- ifelse(d.sp$regen.count.old > 0,TRUE,FALSE)
d.sp$regen.presab.all <- ifelse(d.sp$regen.count.all > 0,TRUE,FALSE)

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


## Make an exception for Harding: only two categories
d.plot[d.plot$Fire == "HARDING","topoclim.cat"] <- ifelse(d.plot[d.plot$Fire=="HARDING","rad.march"] > 6250,"P.1_R.2","P.1_R.1")




  

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




#### 2. Summarize (compute average) regen values (high-sev plots only) and adults (control plots only) by species across all plots in each topoclimatic category in each fire ####

## assign the trees by species their topoclimatic category and fire name. This also ensures that we only are looking at seedlings for whose plots we are interested (because with this merge operation, seedlings from plots not in d.plot will be dropped)
d.sp.cat <- merge(d.sp,d.plot[,c("Regen_Plot","topoclim.cat","Fire","FIRE_SEV","survey.years.post")])

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
d.plot.c <- remove.vars(d.plot,c("Year.of.Fire","Easting","Northing","aspect","Year","precip.category","rad.category"))
# label plots as control or high sev
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% control,"control",NA)
d.plot.c$type <- ifelse(d.plot.c$FIRE_SEV %in% high.sev,"highsev",d.plot.c$type)
# get rid of plots that are neither control nor high sev
d.plot.c <- d.plot.c[!is.na(d.plot.c$type),]
# get rid of regen (high sev) plots that are not surveyed 4-5 years post
d.plot.c <- d.plot.c[(d.plot.c$type == "control") | (d.plot.c$survey.years.post %in% c(4,5)),]

## aggregate plots by Fire, topoclim category, and type (control or high sev)
d.plot.agg.mean <- aggregate(remove.vars(d.plot.c,c("Regen_Plot","topoclim.cat","type","fire.abbr","X5yr","Fire")),by=list(d.plot.c$Fire,d.plot.c$topoclim.cat,d.plot.c$type),FUN=mean,na.rm=TRUE)
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
d.plot.3 <- d.plot.2[which((d.plot.2$count.control > 3) & (d.plot.2$count.highsev > 3)),] #! removed restriction on control count

# Compute additional variables
d.sp.2$proportion.young <- d.sp.2$regen.count.young / d.sp.2$regen.count.all



#### 4. Exploration of pairwise correlations between response variables and predictor variables ####

# Analysis for a single focal species (or species group, such as CONIFER)
focalsp <- "PIPO"
d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
d.mod[is.nan(d.mod)] <- NA

# List interesting variables and make pairs plots
response <- "regen.presab.old"
preds <- c("SHRUB.highsev", "rad.march.highsev", "ppt.normal.highsev", "ppt.post.highsev", "ppt.post.min.highsev", "diff.norm.ppt.z.highsev", "diff.norm.ppt.min.z.highsev", "seed_tree_distance_general.highsev", "adult.ba", "adult.count","proportion.young")
#pairs(d.mod[,c(response,preds)])

# Get the pearson correlation of each predictor with the response (individually)
corr <- NULL
for(i in 1:length(preds)) {
  cat(preds[i])
  corr[preds[i]] <- cor(d.mod[,response],d.mod[,preds[i]],use="complete.obs")
}
corr


### Make a heatmap of pairwise correlations for different predictors and responses

# Define the interesting set of responses
sps <- c("CONIF.ALLSP","PIPO","ABCO","HDWD.ALLSP","PSME","PILA","QUKE","SHADE.ALLSP","PINUS.ALLSP")
responses <- c("regen.presab.old","regen.presab.all","regen.count.all","proportion.young")

opts <- expand.grid(responses,sps,stringsAsFactors=FALSE)
names(opts) <- c("response.opt","sp.opt")

opts.names <- paste(opts$sp.opt,opts$response.opt,sep="-")

#interesting predictors
#early climate
preds <- c("SHRUB.highsev", "rad.march.highsev", "ppt.normal.highsev", "ppt.post.highsev", "ppt.post.min.highsev", "diff.norm.ppt.z.highsev", "diff.norm.ppt.min.z.highsev", "seed_tree_distance_general.highsev", "adult.ba", "adult.count")
# full climate
preds <- c("SHRUB.highsev", "rad.march.highsev", "ppt.normal.highsev", "ppt.post.highsev", "ppt.post.min.highsev", "diff.norm.ppt.z.highsev", "diff.norm.ppt.min.z.highsev", "seed_tree_distance_general.highsev", "adult.ba", "adult.count")

# data frame to store correlation values
corr.df <- data.frame(opt=NA,pred=NA,cor=NA,pval=NA)

for(i in 1:nrow(opts)) {
  
  opt <- opts[i,]
  
  focalsp <- opt$sp.opt
  d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
  d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above
  
  d.mod[d.mod==Inf] <- NA
  d.mod[is.nan(d.mod)] <- NA
  
  newfires <- c("BAGLEY","PEAK","CHIPS")
  #! test removing new fires
  #d.mod <- d.mod[!(d.mod$Fire %in% newfires),]
  
  
  # List interesting variables and make pairs plots
  response <- opt$response.opt

  

  for(j in 1:length(preds)) {
    corr <- cor(d.mod[,response],d.mod[,preds[j]],use="complete.obs")
    
    mod.df <- data.frame(y=d.mod[,response],pred=d.mod[,preds[j]])
    
    #if the response is all the same number (e.g. all plots had no hardwoods), can't fit model; skip
    if(length(unique(mod.df$pred)) < 2) {
      next()
    }
    
    mod.df[is.nan(mod.df)] <- NA
    mod.df[mod.df==Inf] <- NA

    m <- lm(y~pred,data=mod.df,na.action=na.omit)
    p <- summary(m)$coefficients["pred","Pr(>|t|)"]

    corr.df <- rbind(corr.df,data.frame(opt=opts.names[i],pred=preds[j],cor=corr,pval=p))
    
    
  }
}
# remove first row (blank)
corr.df <- corr.df[-1,]

## Plot heatmap
library(ggplot2)


corr.df$cor.sig <- ifelse(corr.df$pval < 0.05, corr.df$cor,NA) #sed non-sig correlations to NA
corr.df$sig <- ifelse(corr.df$pval < 0.05,"*","")

ggplot(corr.df,aes(x=opt,y=pred)) +
  geom_tile(aes(fill=cor)) +
  scale_fill_gradientn(colours=c("red","white","blue"),limits=c(1,-1)) +
  geom_text(aes(label=sig)) +
  labs(x="response",y="predictor") +
  theme_bw(15) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  



# Summary of correlation results: patterns we would expect, and pretty strong! Comparing regen count and regen pres/ab, all of the site factors are more correlated with regen presab, EXCEPT seed tree and adult ba/count
# Shrubs are not important for old seedlings, but they're important for all seedlings. So shrubs become more important later!
# When looking at all seedlings instead of just old seedlings, all predictors are even stronger (weird!) even post-fire weather anomaly. Is it because shrubs do poorly when seedlings initially do well? Hard to tease that apart with this dataset; maybe with a multiple regression Exception to this is adult BA and count better explain old regen than all regen.
# Patterns even stronger when looking at "CONIFER" (especially radiation--weird!), generally even stronger with using all seedlings instead of just old seedlings
# Post precip is a better predictor than normal precip


# Picking the more interesting variables for a multiple regression
# 

focalsp <- "CONIFER"
d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above

newfires <- c("BAGLEY","PEAK","CHIPS")
#! test removing new fires
#d.mod <- d.mod[!(d.mod$Fire %in% newfires),]


resp <- "SHRUB.highsev"
preds <- c("ppt.normal.highsev","rad.march.highsev","diff.norm.ppt.z.highsev","adult.ba")
preds <- c("adult.ba")
preds <- c("rad.march.highsev","diff.norm.ppt.z.highsev","adult.ba")
preds <- c("ppt.normal.highsev","diff.norm.ppt.z.highsev")
preds <- c("ppt.normal.highsev","rad.march.highsev","diff.norm.ppt.z.highsev","adult.ba","diff.norm.ppt.min.z.highsev")
# look for autocorrelation among predictors
pairs(d.mod[,preds])
# not bad

m <- lm(SHRUB.highsev~ppt.normal.highsev*diff.norm.ppt.z.highsev,data=d.mod[,c(resp,preds)])
summary(m)







