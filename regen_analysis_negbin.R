setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(ggplot2)

source("regen_analysis_functions.R")

#### 0. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

# only keep the necessary columns
d.plot <- d.plot[,c("Regen_Plot","Fire","Year.of.Fire","Easting","Northing","aspect","slope","SHRUB","FORB","GRASS","HARDWOOD","CONIFER","FIRE_SEV","BA.Live1","Year","firesev","dist.to.low","fire.abbr","X5yr","fire.year","survey.years.post","elev.m","rad.march","tmean.post","ppt.post","ppt.post.min","tmean.normal","ppt.normal","seed.tree.any","diff.norm.ppt.z","diff.norm.ppt.min.z","seed_tree_distance_general","seed_tree_distance_conifer","seed_tree_distance_hardwood","diff.norm.ppt.z","diff.norm.ppt.min.z","tmean.post","ppt.post","ppt.post.min","perc.norm.ppt","perc.norm.ppt.min","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z","def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","diff.norm.def.max.z","diff.norm.aet.min.z","def.post","aet.post")]

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

pairs(d.mod[,preds])


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
  





#### 4.1 Statistical models with cluster-level data ####

# Analysis for a single focal species (or species group, such as CONIFER)
focalsp <- "ABCO"
d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above



vars.leave <- c("","FORB.highsev","SHRUB.highsev","GRASS.highsev","CONIFER.highsev","HARDWOOD.highsev","FIRE_SEV.highsev","firesev.highsev","fire.year.highsev","survey.years.post.highsev","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal.highsev","diff.norm.ppt.z.highsev","ppt.normal.sq.highsev","rad.march.highsev","seed_tree_distance_general.highsev","SHRUB.highsev","adult.count","adult.ba")
d.mod <- d.mod[complete.cases(d.mod[,vars.focal]),]
d.c <- center.df(d.mod,vars.leave)
# 
# ## transform cover so it does not include 0 or 1 (for Beta distrib)
# d.c$SHRUB.p <- d.c$SHRUB/100
# d.c$SHRUB.pt <- (d.c$SHRUB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
# 
# d.c$GRASS.p <- d.c$GRASS/100
# d.c$GRASS.pt <- (d.c$GRASS.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
# 
# d.c$HARDWOOD.p <- d.c$HARDWOOD/100
# d.c$HARDWOOD.pt <- (d.c$HARDWOOD.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
# 
# d.c$FORB.p <- d.c$FORB/100
# d.c$FORB.pt <- (d.c$FORB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
# 
# d.c$CONIFER.p <- d.c$CONIFER/100
# d.c$CONIFER.pt <- (d.c$CONIFER.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)
# 

d.c$Fire <- as.factor(d.c$Fire)

d.c <- d.c[!(d.c$Fire == "RICH"),]

d.c$regen.count.all.int <- ceiling(d.c$regen.count.all)
d.c$regen.count.old.int <- ceiling(d.c$regen.count.old)


d.c$ppt.normal.highsev_c.sq <- d.c$ppt.normal.highsev_c^2
d.c$tmean.normal.highsev_c.sq <- d.c$tmean.normal.highsev_c^2


vars <- c("ppt.normal.highsev_c" , "ppt.normal.highsev_c.sq" , "tmean.normal.highsev_c" , "tmean.normal.highsev_c.sq" , "seed_tree_distance_general.highsev_c" , "rad.march.highsev_c")
d.foc <- d.c[,vars]
pairs(d.foc)

d.c$regen.count.all.nz <- ifelse(d.c$regen.count.all == 0,0.1,d.c$regen.count.all)

m.nofire <- brm(regen.count.all.nz ~ ppt.normal.highsev_c + ppt.normal.highsev_c.sq + seed_tree_distance_general.highsev_c + rad.march.highsev_c,family="Gamma",data=d.c,iter=6000,control = list(adapt_delta = 0.90),cores=3,chains=3)
m.fire <- brm(regen.count.all.nz ~ ppt.normal.highsev_c + ppt.normal.highsev_c.sq + seed_tree_distance_general.highsev_c + rad.march.highsev_c + (1|Fire),family="Gamma",data=d.c,iter=6000,control = list(adapt_delta = 0.90),cores=3,chains=3)
m.adult <- brm(regen.count.all.nz ~ ppt.normal.highsev_c + ppt.normal.highsev_c.sq + seed_tree_distance_general.highsev_c + rad.march.highsev_c + adult.count_c,family="Gamma",data=d.c,iter=6000,control = list(adapt_delta = 0.90),cores=3,chains=3)
m.adultba <- brm(regen.count.all.nz ~ ppt.normal.highsev_c + ppt.normal.highsev_c.sq + seed_tree_distance_general.highsev_c + rad.march.highsev_c + adult.ba_c,family="Gamma",data=d.c,iter=6000,control = list(adapt_delta = 0.90),cores=3,chains=3)
m.adultba.fire <- brm(regen.count.all.nz ~ ppt.normal.highsev_c + ppt.normal.highsev_c.sq + seed_tree_distance_general.highsev_c + rad.march.highsev_c + adult.ba_c + (1|Fire),family="Gamma",data=d.c,iter=6000,control = list(adapt_delta = 0.90),cores=3,chains=3)

loos <- loo(m.nofire,m.fire,m.adult,m.adultba,m.adultba.fire)
loos

summary(m)



###+ For cluster-level analysis, regardless of species, adding fire random effect dramatically improves model (but probably because very few data points per fire), regardless of species




# Summary of correlation results: patterns we would expect, and pretty strong! Comparing regen count and regen pres/ab, all of the site factors are more correlated with regen presab, EXCEPT seed tree and adult ba/count
# Shrubs are not important for old seedlings, but they're important for all seedlings. So shrubs become more important later!
# When looking at all seedlings instead of just old seedlings, all predictors are even stronger (weird!) even post-fire weather anomaly. Is it because shrubs do poorly when seedlings initially do well? Hard to tease that apart with this dataset; maybe with a multiple regression Exception to this is adult BA and count better explain old regen than all regen.
# Patterns even stronger when looking at "CONIFER" (especially radiation--weird!), generally even stronger with using all seedlings instead of just old seedlings
# Post precip is a better predictor than normal precip









#### 5. Plot-level exploration ####


d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV > 3),]

sp <- "PIPO"
d.sp.curr <- d.sp[d.sp$species==sp,]
d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","def.normal","diff.norm.def.z")
d <- d[complete.cases(d[,vars.focal]),]
d.c <- center.df(d,vars.leave)

vars.focal.c <- paste0(vars.focal[-6],"_c")
pairs(d.c[,vars.focal.c])

d.c$SHRUB.p <- d.c$SHRUB/100
d.c$SHRUB.pt <- (d.c$SHRUB.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$GRASS.p <- d.c$GRASS/100
d.c$GRASS.pt <- (d.c$GRASS.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)

d.c$HARDWOOD.p <- d.c$HARDWOOD/100
d.c$HARDWOOD.pt <- (d.c$HARDWOOD.p*(nrow(d.c)-1) + 0.5) / nrow(d.c)


d.c$Fire <- as.factor(d.c$Fire)

d.c <- d.c[!(d.c$Fire == "RICH"),]


ggplot(d.c,aes(x=diff.norm.ppt.z_c,y=ppt.normal_c,color=Fire,shape=regen.presab.all)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(1,19)) +
  #facet_wrap(~Fire,scales="free") +
  theme_bw(20)












#### 6. Plot-level analysis with BRMS ####



## First, data frame for making predictions

diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100)

newdat.ppt <- data.frame(
  ppt.normal_c = c(rep(-1,100),rep(1,100)),
  ppt.normal_c.sq = c(rep(1,100),rep(1,100)),
  norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.ppt.z_c = rep(diff.norm.seq,2),
  diff.norm.ppt.z_c.sq = rep(diff.norm.seq^2,2),
  tmean.normal_c = 0,
  tmean.normal_c.sq = 0,
  diff.norm.tmean.z_c = 0,
  diff.norm.tmean.z_c.sq = 0,
  scenario = "ppt",
  
  aet.normal_c = c(rep(-1,100),rep(1,100)),
  aet.normal_c.sq = c(rep(1,100),rep(1,100)),
  diff.norm.aet.z_c = rep(diff.norm.seq,2),
  diff.norm.aet.z_c.sq = rep(diff.norm.seq^2,2),
  def.normal_c = 0,
  def.normal_c.sq = 0,
  diff.norm.def.z_c = 0,
  diff.norm.def.z_c.sq = 0
)

newdat.tmean <- data.frame(
  tmean.normal_c = c(rep(-1,100),rep(1,100)),
  tmean.normal_c.sq = c(rep(1,100),rep(1,100)),
  norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.tmean.z_c = rep(diff.norm.seq,2),
  diff.norm.tmean.z_c.sq = rep(diff.norm.seq^2,2),
  ppt.normal_c = 0,
  ppt.normal_c.sq = 0,
  diff.norm.ppt.z_c = 0,
  diff.norm.ppt.z_c.sq = 0,
  scenario = "tmean",
  
  def.normal_c = c(rep(-1,100),rep(1,100)),
  def.normal_c.sq = c(rep(1,100),rep(1,100)),
  diff.norm.def.z_c = rep(diff.norm.seq,2),
  diff.norm.def.z_c.sq = rep(diff.norm.seq^2,2),
  aet.normal_c = 0,
  aet.normal_c.sq = 0,
  diff.norm.aet.z_c = 0,
  diff.norm.aet.z_c.sq = 0
)

newdat <- rbind(newdat.ppt,newdat.tmean)
interact.df <- data.frame("Fake"=rep(NA,nrow(newdat)))
for(i in 1:ncol(newdat)) { # for each col
  for(j in 1:ncol(newdat)) {
    name.i <- names(newdat)[i]
    name.j <- names(newdat)[j]
    name.inter <- paste(name.i,name.j,sep=":")
    val.inter <- newdat[,i] * newdat[,j]
    interact.df <- cbind(interact.df,val.inter)
    names(interact.df)[ncol(interact.df)] <- name.inter
  }
}

newdat <- cbind(newdat,interact.df)
newdat$"(Intercept)" <- 1 # this is the "predictor" value to multiple the intercept by




## Next, fit models ## 

library(brms)
library(loo)

d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV > 3),]


sp.opts <- c("PINUS.ALLSP","SHADE.ALLSP","HDWD.ALLSP","PIPO","ABCO","ABMA","CONIF.ALLSP","PSME","PILA","CADE27","PIJE","PIPJ")
sp.opts <- c("PINUS.ALLSP","SHADE.ALLSP","HDWD.ALLSP","PIPO","ABCO","PSME","PILA","PIPJ") # reduced


cover.opts <- c("COV.SHRUB","COV.GRASS","COV.FORB","COV.HARDWOOD","COV.CONIFER")
cover.opts <- c("COV.SHRUB","COV.GRASS") # reduced
#cover.opts <- NULL
sp.opts <- c(sp.opts,cover.opts)

m.p <- list()

loos <- list()

dat.preds <- data.frame()
pred.obs <- data.frame()

d.loos.all <- data.frame()
d.loo.comps <- data.frame()

for(sp in sp.opts) { # about 1 hr per species
  
      cat("\n\n#####")
      cat("Running model for: ",sp,"")
      cat("#####\n\n")
      
      
      if(sp %in% cover.opts) {
        d.sp.curr <- d.sp[d.sp$species=="PIPO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
      } else {
        d.sp.curr <- d.sp[d.sp$species==sp,]
      }
      
      d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
      vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
      vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z", "def.normal","aet.normal","diff.norm.def.z","diff.norm.aet.z","def.post","aet.post")
      d <- d[complete.cases(d[,vars.focal]),]
      d.c <- center.df(d,vars.leave)

      d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2
      d.c$tmean.normal_c.sq <- d.c$tmean.normal_c^2
      
      d.c$def.normal_c.sq <- d.c$def.normal_c^2
      d.c$aet.normal_c.sq <- d.c$aet.normal_c^2
      
      # ####!!!! trick model: make diff.norm into diff.norm.min
      # d.c$diff.norm.ppt.z_c <- d.c$diff.norm.ppt.min.z_c
      # d.c$diff.norm.tmean.z_c <- d.c$diff.norm.tmean.max.z_c
      # #### end trick model
      
      
      d.c$diff.norm.ppt.z_c.sq <- d.c$diff.norm.ppt.z_c^2
      d.c$diff.norm.tmean.z_c.sq <- d.c$diff.norm.tmean.z_c^2
      
      d.c$diff.norm.def.z_c.sq <- d.c$diff.norm.def.z_c^2
      d.c$diff.norm.aet.z_c.sq <- d.c$diff.norm.aet.z_c^2
      
      d.c$ppt.post_c.sq <- d.c$ppt.post_c^2
      d.c$tmean.post_c.sq <- d.c$tmean.post_c^2
      
      d.c$def.post_c.sq <- d.c$def.post_c^2
      d.c$aet.post_c.sq <- d.c$aet.post_c^2
      
      d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
      d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)
      
      
      vars.focal.c <- paste0(vars.focal[-6],"_c")
      pairs(d.c[,vars.focal.c])
      
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
      
      d.c <- d.c[!(d.c$Fire == "RICH"),]
      
      d.c$regen.count.all.int <- ceiling(d.c$regen.count.all)
      d.c$regen.count.old.int <- ceiling(d.c$regen.count.old)
      
      
      if(sp %in% cover.opts) {
        
        sp.cov <- substr(sp,5,100)
        sp.cov <- paste0(sp.cov,".pt")
        
        d.c$cov.response <- d.c[,sp.cov]
  
        d.c$response.var <- d.c$cov.response
        
        mod.family <- "Beta"
        
      } else {

        d.c$response.var <- d.c$regen.count.old.int
        
        mod.family <- "zero_inflated_negbinomial"
        
      }
      
      m <- list()
      
      m[["n0.a0"]] <- brm(response.var ~ 1 + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nPT.a0"]] <- brm(response.var ~ ppt.normal_c + ppt.normal_c.sq + tmean.normal_c + tmean.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nPT.aPT"]] <- brm(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nPT.aPT2"]] <- brm(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq + tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["pPT"]] <- brm(response.var ~ ppt.post_c + ppt.post_c.sq + tmean.post_c + tmean.post_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nP.a0"]] <- brm(response.var ~ ppt.normal_c + ppt.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nP.aP"]] <- brm(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nP.aP2"]] <- brm(response.var ~ ppt.normal_c*diff.norm.ppt.z_c + ppt.normal_c*diff.norm.ppt.z_c.sq + diff.norm.ppt.z_c.sq + ppt.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["pP"]] <- brm(response.var ~ ppt.post_c + ppt.post_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nT.p0"]] <- brm(response.var ~ tmean.normal_c + tmean.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nT.aT"]] <- brm(response.var ~ tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nT.aT2"]] <- brm(response.var ~ tmean.normal_c*diff.norm.tmean.z_c + tmean.normal_c*diff.norm.tmean.z_c.sq + diff.norm.tmean.z_c.sq + tmean.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["pT"]] <- brm(response.var ~ tmean.post_c + tmean.post_c.sq + tmean.post_c + tmean.post_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      
      m[["n0.a0"]] <- brm(response.var ~ 1 + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nAD.a0"]] <- brm(response.var ~ aet.normal_c + aet.normal_c.sq + def.normal_c + def.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nAD.aAD"]] <- brm(response.var ~ aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq + def.normal_c*diff.norm.def.z_c + def.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nAD.aAD2"]] <- brm(response.var ~ aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq + def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["pAD"]] <- brm(response.var ~ aet.post_c + aet.post_c.sq + def.post_c + def.post_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nA.a0"]] <- brm(response.var ~ aet.normal_c + aet.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nA.aA"]] <- brm(response.var ~ aet.normal_c*diff.norm.aet.z_c + aet.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nA.aA2"]] <- brm(response.var ~ aet.normal_c*diff.norm.aet.z_c + aet.normal_c*diff.norm.aet.z_c.sq + diff.norm.aet.z_c.sq + aet.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["pA"]] <- brm(response.var ~ aet.post_c + aet.post_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nD.p0"]] <- brm(response.var ~ def.normal_c + def.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nD.aD"]] <- brm(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["nD.aD2"]] <- brm(response.var ~ def.normal_c*diff.norm.def.z_c + def.normal_c*diff.norm.def.z_c.sq + diff.norm.def.z_c.sq + def.normal_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      m[["pD"]] <- brm(response.var ~ def.post_c + def.post_c.sq + def.post_c + def.post_c.sq + (1|Fire),family=mod.family,data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      

      
      #get individual loos (will get specific pairs later)
      
      d.loos <- data.frame()
      
      for(mod in names(m)) {
        cat(mod," ")
        loo.mod <- loo(m[[mod]])
        d.mod <- data.frame(sp=sp,mod=mod,LOOIC=loo.mod$looic,SE=loo.mod$se_looic)
        d.loos <- rbind(d.loos,d.mod)
        
      }
      
      d.loos$upr <- d.loos$LOOIC + 1.96*d.loos$SE
      
      d.loos.all <- rbind(d.loos.all,d.loos)
      
      ### Specify model categories
      
      ## normal models
      search <- ".a0"
      normal.model.names <- d.loos[grep(search,d.loos$mod,fixed=TRUE),]$mod
      
      ## anomaly models
      search <- ".a"
      anomaly.model.names <- d.loos[grep(search,d.loos$mod,fixed=TRUE),]$mod
      anomaly.model.names <- anomaly.model.names[!(anomaly.model.names %in% normal.model.names)]
      
      ### Find the best normal model
      
      d.loos.normal <- d.loos[d.loos$mod %in% normal.model.names,]
      best.normal.mod <- d.loos.normal[d.loos.normal$upr == min(d.loos.normal$upr),]$mod[1]
      
      ### Compare null to best normal
      comp.null.normal <- loo(m[["n0.a0"]],m[[best.normal.mod]])$ic_diffs__
      
      ### Get post corresponding to best normal
      norm.part <- strsplit(as.character(best.normal.mod),".",fixed=TRUE)[[1]][1]
      norm.post.mod <- sub("n","p",norm.part)
      
      ### Compare best normal to corresponding post
      comp.normal.post <- loo(m[[best.normal.mod]],m[[norm.post.mod]])$ic_diffs__
      
      
      ### Get anomaly models corresponding to best normal
      
      normal.part <- strsplit(as.character(best.normal.mod),".",fixed=TRUE)[[1]][1]
      search <- paste0(normal.part,".")
      normal.matches <- grepl(search,d.loos$mod,fixed=TRUE)
      search <- paste0(normal.part,".a0")
      non.anom.match <- grepl(search,d.loos$mod,fixed=TRUE)
      normal.anom <- normal.matches & (!non.anom.match)
      normal.anom.names <- d.loos[normal.anom,]$mod
      
      ### Get the best corresponding anomaly model
      
      d.loos.anom <- d.loos[d.loos$mod %in% normal.anom.names,]
      best.normal.anom <- d.loos.anom[d.loos.anom$upr == min(d.loos.anom$upr),]$mod[1]
      
      ### Compare best normal with corresponding best anomaly
      comp.normal.anom <- loo(m[[best.normal.mod]],m[[best.normal.anom]])$ic_diffs__
      
      
      ### Get the best anomaly, independent of best normal
      d.loos.anomaly <- d.loos[d.loos$mod %in% anomaly.model.names,]
      best.anomaly.mod <- d.loos.anomaly[d.loos.anomaly$upr == min(d.loos.anomaly$upr),]$mod[1]
      
      ### Get the normal corresponding to the best anomaly
      best.anom.normal.part <- strsplit(as.character(best.anomaly.mod),".",fixed=TRUE)[[1]][1]
      best.anom.normal <- paste0(best.anom.normal.part,".a0")
      
      ### Compare corresponding normal to best anomaly
      comp.anom.normal <- loo(m[[best.anom.normal]],m[[best.anomaly.mod]])$ic_diffs__
      
      
      
      ### Store loos comps
      d.loo.comps.sp <- rbind(comp.null.normal,comp.normal.post,comp.normal.anom,comp.anom.normal)
      d.loo.comps.sp <- as.data.frame(d.loo.comps.sp)
      d.loo.comps.sp$comp <-c("null.normal","normal.post","normal.anom","anom.normal")
      d.loo.comps.sp$sp <- sp
      
      d.loo.comps <- rbind(d.loo.comps,d.loo.comps.sp)
      
      
      
      
      
      
      
      ### For best anomaly model, store predictions for counterfactual plots
      m.sp <- m[[best.anomaly.mod]]
      m.p.sp <- posterior_samples(m.sp,pars="^b")
      
      names(m.p.sp) <- sapply(names(m.p.sp),substr,start=3,stop=1000)
      names(m.p.sp)[1] <- "(Intercept)"
      
      mod.vars <- colnames(m.p.sp)
      
      newdat.ordered <- newdat[,mod.vars]
      
      # for each row of newdat.ordered, multiply data by all coef values, sum across, and compute median and conf int
      preds <- matrix(nrow=nrow(newdat.ordered),ncol=3)
      for(i in 1:nrow(newdat.ordered)) {
        
        newdat.row <- as.matrix(newdat.ordered[i,])
        
        prod <- sweep(m.p.sp,MARGIN=2,newdat.row,`*`)
        sums <- rowSums(prod)
        
        fit <- median(sums)
        upr <- quantile(sums,probs=0.975)
        lwr <- quantile(sums,probs=0.025)
        
        preds[i,] <- c(fit,lwr,upr)
        
      }
      
      if(sp %in% cover.opts) {
        preds <- inv.logit(preds)
      } else {
        preds <- exp(preds)
      }
      
      
      colnames(preds) <- c("fit","lwr","upr")
      preds <- as.data.frame(preds)
      
      sp.grp <- sp
      dat.preds.sp <- cbind(newdat,preds,sp.grp,best.anomaly.mod)
      dat.preds <- rbind(dat.preds,dat.preds.sp) # make sure this is working right to store all species in th e same DF
      

      
      
      
      
      ### store a predicted vs. observed data frame here: one for the best anomaly model and one for the corresponding normal model ###
      predicted.normal <- predict(m[[best.anom.normal]])[,1]
      predicted.anomaly <- predict(m[[best.anomaly.mod]])[,1]
      
      observed <- d.c$response.var
      

      pred.obs.normal <- data.frame(pred=predicted.normal,obs=observed,mod="normal",sp=sp)
      pred.obs.anomaly <- data.frame(pred=predicted.anomaly,obs=observed,mod="anomaly",sp=sp)
      
      pred.obs.sp <- rbind(pred.obs.normal,pred.obs.anomaly)
      
      if(sp %in% cover.opts) {
        pred.obs.sp$pred <- inv.logit(pred.obs.sp$pred)
      } else {
        pred.obs.sp$pred <- exp(pred.obs.sp$pred)
      }
      
      pred.obs <- rbind(pred.obs,pred.obs.sp)
      
      
      
      
      
      
}




dat.preds <- read.csv("../data_model_summaries/summaries_wb_1/counterfactual_df.csv",header=TRUE)



## plot predictions over precip anomaly

dat.pred <- dat.preds[dat.preds$scenario=="ppt",]


ggplot(dat.pred,aes(x=diff.norm.ppt.z_c,y=fit,color=norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip")) +
  facet_wrap(~sp.grp,scales="free",ncol=5) +
  scale_y_continuous(limits=c(0,1))


## plot predictions over tmean anomaly

dat.pred <- dat.preds[dat.preds$scenario=="tmean",]

ggplot(dat.pred,aes(x=diff.norm.tmean.z_c,y=fit,color=norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal tmean"),color=guide_legend(title="Normal tmean")) +
  facet_wrap(~sp.grp,scales="free",ncol=5) +
  scale_y_continuous(limits=c(0,10))






## save counterfactual and loo dataframes, and predicted vs. observeds

write.csv(dat.preds,"counterfactual_df.csv",row.names=FALSE)
write.csv(m.loos,"loos_df.csv",row.names=FALSE)
write.csv(pred.obs,"pred_obs.csv",row.names=FALSE)





### plot looic diffs

## null vs: Pn, PTn
## Pn vs: PTn
## Pn vs: Pp
## PTn vs: PTp
## Pn vs: Pna, Pna2
## PTn vs: PTna, PTna2
## Pna vs: Pna2
## PTna vs: PTna2

m.loos <- read.csv("summaries_2/loos_df.csv",header=TRUE,stringsAsFactors=TRUE)

comps <- c("m.null - m.Pn",
           "m.null - m.PTn",
           "m.Pn - m.PTn",
           "m.Pn - m.Pp",
           "m.PTn - m.PTp",
           "m.Pn - m.Pna",
           "m.Pn - m.Pna2",
           "m.PTn - m.PTna",
           "m.PTn - m.PTna2",
           "m.Pna - m.Pna2",
           "m.PTna - m.PTna2")

m.loos.thin <- m.loos[m.loos$comp %in% comps,]

m.loos.thin$upr <- m.loos.thin$LOOdiff + m.loos.thin$SE*1.96
m.loos.thin$lwr <- m.loos.thin$LOOdiff - m.loos.thin$SE*1.96

m.loos.thin$comp <- factor(m.loos.thin$comp,rev(comps))

ggplot(m.loos.thin,aes(y=comp,x=LOOdiff)) +
  geom_point() +
  geom_errorbarh(aes(xmax=upr,xmin=lwr)) +
  facet_wrap(~sp,scales="free_x",ncol=2) +
  theme_bw(12) +
  geom_vline(xintercept=0,color="red")





#### 7. Test how much better a model is with FIRE ####


d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV > 3),]


sp <- "ABCO"
d.sp.curr <- d.sp[d.sp$species==sp,]
d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z")
d <- d[complete.cases(d[,vars.focal]),]
d.c <- center.df(d,vars.leave)

d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2
d.c$tmean.normal_c.sq <- d.c$tmean.normal_c^2

# ####!!!! trick model: make diff.norm into diff.norm.min
# d.c$diff.norm.ppt.z_c <- d.c$diff.norm.ppt.min.z_c
# d.c$diff.norm.tmean.z_c <- d.c$diff.norm.tmean.max.z_c
# #### end trick model


d.c$diff.norm.ppt.z_c.sq <- d.c$diff.norm.ppt.z_c^2
d.c$diff.norm.tmean.z_c.sq <- d.c$diff.norm.tmean.z_c^2

d.c$ppt.post_c.sq <- d.c$ppt.post_c^2
d.c$tmean.post_c.sq <- d.c$tmean.post_c^2

d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)


vars.focal.c <- paste0(vars.focal[-6],"_c")
pairs(d.c[,vars.focal.c])

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

d.c <- d.c[!(d.c$Fire == "RICH"),]

d.c$regen.count.all.int <- ceiling(d.c$regen.count.all)
d.c$regen.count.old.int <- ceiling(d.c$regen.count.old)


library(brms)

m.nofire <- brm(regen.count.all.int ~ ppt.normal_c + ppt.normal_c.sq + tmean.normal_c + tmean.normal_c.sq + seed_tree_distance_general_c + rad.march_c,family="zero_inflated_negbinomial",data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
m.fire <- brm(regen.count.all.int ~ ppt.normal_c + ppt.normal_c.sq + tmean.normal_c + tmean.normal_c.sq + seed_tree_distance_general_c + rad.march_c + (1|Fire),family="zero_inflated_negbinomial",data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)


loos <- loo(m.nofire,m.fire)
loos




###+ Summary of results: with PIPO, using all non-anomaly predictors without fire random effect is worse than including fire random effect, but with ABCO, adding fire random effect provides little improvement. Is this legit though?








#### 8. Cluster-level multivariate analysis ####




#d.plot.3
#d.sp.2
library(data.table)



focal.sp <- c("PIPJ","ABCO","PILA","PSME","QUKE","CADE27") #   "ABMA","QUCH2","ARME") #ABMA, PIPO, QUCH2, 
#focal.sp <- c("PIPO","ABCO","PILA","PSME","QUKE","LIDE3","QUCH2","CADE27")
#focal.sp <- c("PINUS.ALLSP","SHADE.ALLSP","HDWD.ALLSP")
focal.cols <- c("Fire","topoclim.cat","species","regen.count.old","regen.count.all","adult.ba")


d.sp.simp <- d.sp.2[d.sp.2$species %in% focal.sp,focal.cols]

names(d.sp.simp)[4:6] <- c("r.old","r.all","a.ba")

d.sp.simp <- as.data.table(d.sp.simp)
d.sp.cast <- dcast(d.sp.simp,topoclim.cat + Fire~species,value.var=c("r.old","r.all","a.ba"))


regen.var <- "r.all"

d.all <- merge(d.plot.3,d.sp.cast,by=c("Fire","topoclim.cat"))
sp.cols <- grep(regen.var,names(d.all))
d.all.sp <- d.all[,sp.cols]
d.all.sp.tot <- rowSums(d.all.sp)
keep.rows <- d.all.sp.tot > 0
d.all.sp <- d.all.sp[keep.rows,]
d.all <- d.all[keep.rows,]


library(vegan)


### try adding diff predictors
# 
# c0 <- cca(d.all.sp ~ 1,data=d.all)
# c1 <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev + seed_tree_distance_general.highsev + BA.Live1.control,data=d.all)
# c1.fire <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev + seed_tree_distance_general.highsev + BA.Live1.control + Fire,data=d.all)
# c2  <- cca(d.all.sp ~ ppt.normal.highsev + tmean.normal.highsev + rad.march.highsev + seed_tree_distance_general.highsev + BA.Live1.control + diff.norm.ppt.z.highsev + diff.norm.tmean.z.highsev,data=d.all)

c3  <- cca(d.all.sp ~ ppt.normal.highsev + rad.march.highsev + diff.norm.ppt.z.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
c4 <- cca(d.all.sp ~ ppt.normal.highsev + rad.march.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)

c3.def  <- cca(d.all.sp ~ def.normal.highsev + rad.march.highsev + diff.norm.def.z.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
c4.def <- cca(d.all.sp ~ def.normal.highsev + rad.march.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)
c5.def  <- cca(d.all.sp ~ def.normal.highsev + aet.normal.highsev + rad.march.highsev + diff.norm.def.z.highsev + diff.norm.aet.z.highsev + a.ba_ABCO + a.ba_PIPJ + a.ba_PSME + a.ba_PILA + a.ba_QUKE,data=d.all)


# cca1.plot <- plot(c1,choices=c(1,2))
# cca1.fire.plot <- plot(c1.fire,choices=c(1,2))
# cca2.plot <- plot(c2,choices=c(1,2))

cca3.plot <- plot(c3,choices=c(1,2))

cca3.def.plot <- plot(c3.def,choices=c(1,2))
cca5.def.plot <- plot(c5.def,choices=c(1,2))


extractAIC(c1)


summary(c2)


#### Mantel test to explain regen species comp with control species comp ####

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







d.sp.curr <- d.sp[d.sp$species==sp,]
d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z")
d <- d[complete.cases(d[,vars.focal]),]
d.c <- center.df(d,vars.leave)

d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2
d.c$tmean.normal_c.sq <- d.c$tmean.normal_c^2

# ####!!!! trick model: make diff.norm into diff.norm.min
# d.c$diff.norm.ppt.z_c <- d.c$diff.norm.ppt.min.z_c
# d.c$diff.norm.tmean.z_c <- d.c$diff.norm.tmean.max.z_c
# #### end trick model


d.c$diff.norm.ppt.z_c.sq <- d.c$diff.norm.ppt.z_c^2
d.c$diff.norm.tmean.z_c.sq <- d.c$diff.norm.tmean.z_c^2

d.c$ppt.post_c.sq <- d.c$ppt.post_c^2
d.c$tmean.post_c.sq <- d.c$tmean.post_c^2

d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)


vars.focal.c <- paste0(vars.focal[-6],"_c")
pairs(d.c[,vars.focal.c])

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

d.c <- d.c[!(d.c$Fire == "RICH"),]

d.c$regen.count.all.int <- ceiling(d.c$regen.count.all)
d.c$regen.count.old.int <- ceiling(d.c$regen.count.old)












