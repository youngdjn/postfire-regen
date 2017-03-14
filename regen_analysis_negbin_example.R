setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(ggplot2)

source("regen_analysis_functions.R")

#### 0. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

# only keep the necessary columns
d.plot <- d.plot[,c("Regen_Plot","Fire","Year.of.Fire","Easting","Northing","aspect","slope","SHRUB","FORB","GRASS","HARDWOOD","CONIFER","FIRE_SEV","BA.Live1","Year","firesev","dist.to.low","fire.abbr","X5yr","fire.year","survey.years.post","elev.m","rad.march","tmean.post","ppt.post","ppt.post.min","tmean.normal","ppt.normal","seed.tree.any","diff.norm.ppt.z","diff.norm.ppt.min.z","seed_tree_distance_general","seed_tree_distance_conifer","seed_tree_distance_hardwood","diff.norm.ppt.z","diff.norm.ppt.min.z","tmean.post","ppt.post","ppt.post.min","perc.norm.ppt","perc.norm.ppt.min","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z")]

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




#### Climate space of topocimatic categories ####

d.plot.3$fire.topoclim.cat <- paste(d.plot.3$Fire,d.plot.3$topoclim.cat,sep="-")

ggplot(d.plot.3,aes(x=ppt.normal.highsev,y=diff.norm.ppt.z.highsev,color=Fire)) +
  geom_point(size=3) +
  labs(x="Normal precip",y="Postfire precip anomaly (3-year average)") +
  theme_bw(15)




#### Model selection ####

# Define the interesting set of responses
sps <- c("CONIF.ALLSP","PIPO","ABCO","HDWD.ALLSP","PSME","PILA","QUKE","SHADE.ALLSP","PINUS.ALLSP")
responses <- c("regen.presab.old","regen.presab.all","regen.count.all","proportion.young")

opts <- expand.grid(responses,sps,stringsAsFactors=FALSE)
names(opts) <- c("response.opt","sp.opt")

opts.names <- paste(opts$sp.opt,opts$response.opt,sep="-")

preds <- list(
  c("nullmod"),
  c("adult.count"),
  c("seed_tree_distance_general.highsev"),
  c("rad.march.highsev"),
  c("ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("adult.count","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("adult.count","seed_tree_distance_general.highsev"),
  c("adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("adult.count","ppt.normal.highsev","ppt.normal.sq.highsev"),
  
  c("diff.norm.ppt.z.highsev"),
  c("diff.norm.ppt.z.highsev","adult.count"),
  c("diff.norm.ppt.z.highsev","seed_tree_distance_general.highsev"),
  c("diff.norm.ppt.z.highsev","rad.march.highsev"),
  c("diff.norm.ppt.z.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.z.highsev","adult.count","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.z.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.z.highsev","adult.count","seed_tree_distance_general.highsev"),
  c("diff.norm.ppt.z.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.z.highsev","adult.count","ppt.normal.highsev","ppt.normal.sq.highsev"),
  
  c("diff.norm.ppt.min.z.highsev"),
  c("diff.norm.ppt.min.z.highsev","adult.count"),
  c("seed_tree_distance_general.highsev"),
  c("diff.norm.ppt.min.z.highsev","rad.march.highsev"),
  c("diff.norm.ppt.min.z.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.min.z.highsev","adult.count","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.min.z.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.min.z.highsev","adult.count","seed_tree_distance_general.highsev"),
  c("diff.norm.ppt.min.z.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev"),
  c("diff.norm.ppt.min.z.highsev","adult.count","ppt.normal.highsev","ppt.normal.sq.highsev")
)


aic.m <- matrix(nrow=nrow(opts),ncol=length(preds))
for(j in 1:nrow(opts)) {

  for(i in 1:length(preds)) {
    
    opts.row <- opts[j,]
    
    focalsp <- opts.row$sp.opt
    d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
    d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above
    d.mod$nullmod <- seq(from=0,to=1,length.out=nrow(d.mod))
    d.mod$y <- d.mod[,opts.row$response.opt]
    
    d.mod.curr <- d.mod[,c("y",preds[[i]])]
    
    m <- lm(y~.,data=d.mod.curr)
    aic.m[j,i] <- AIC(m)
    
  }
  
  aic.min <- min(aic.m[j,])
  aic.m[j,] <- aic.m[j,] - aic.min

}

rownames(aic.m) <- opts.names
aic.m <- round(aic.m)

write.csv(aic.m,"data_analysis_output/model_selection.csv",row.names=TRUE)


#### GLMnet ####

library(glmnet)

response <- c("regen.abcount.all","regen.prescount.all")
predictors <- c("diff.norm.ppt.z.highsev","diff.norm.ppt.min.z.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev")

#alternative
#predictors <- c("diff.norm.ppt.min.z.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev")


focalsp <- "SHADE.ALLSP"
d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above

d.mod$regen.prescount.all <- d.mod$regen.presab.all * d.mod$count.highsev
d.mod$regen.abcount.all <- d.mod$count.highsev - d.mod$regen.prescount.all

#d.mod$interaction <- d.mod$diff.norm.ppt.min.z.highsev * d.mod$ppt.normal.highsev
#predictors <- c(predictors,"interaction")

d.mod.glmnet <- d.mod[,c(response,predictors)]
#d.mod.glmnet <- complete.cases(d.mod.glmnet)
d.mod.glmnet <- as.matrix(d.mod.glmnet)




lambdas = NULL
for (i in 1:100)
{
  fit <- cv.glmnet(d.mod.glmnet[,predictors],d.mod.glmnet[,response],nfolds=10,family="binomial")
  errors = data.frame(fit$lambda,fit$cvm)
  lambdas <- rbind(lambdas,errors)
}
# take mean cvm for each lambda
lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)

# select the best one
bestindex = which(lambdas[2]==min(lambdas[2]))
bestlambda = lambdas[bestindex,1]
bestlambda
# and now run glmnet once more with it
fit <- glmnet(d.mod.glmnet[,predictors],d.mod.glmnet[,response],family="binomial",lambda=bestlambda)
coef(fit)


n.bootstraps <- 1000
boot.coefs <- matrix(nrow=length(predictors)+1,ncol=n.bootstraps)
### bootstrap the fits. Sample the data 1000 times, with replacement
for(i in 1:n.bootstraps) {
  d.mod.glmnet.samp <- d.mod.glmnet[sample(nrow(d.mod.glmnet),size=nrow(d.mod.glmnet),replace=TRUE),]
  fit <- glmnet(d.mod.glmnet.samp[,predictors],d.mod.glmnet.samp[,response],family="binomial",lambda=bestlambda)
  coefs <- coef(fit)
  boot.coefs[,i] <- as.vector(coefs)
}

rownames(boot.coefs) <- rownames(coefs)

boot.coefs[boot.coefs==0] <- NA





m <- glm()




#### GLMNET with interactions ####


response <- c("regen.presab.all")
predictors <- c("diff.norm.ppt.z.highsev","diff.norm.ppt.min.z.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev")

#alternative
#predictors <- c("diff.norm.ppt.min.z.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev")


focalsp <- "PINUS.ALLSP"
d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above

#d.mod$interaction <- d.mod$diff.norm.ppt.min.z.highsev * d.mod$ppt.normal.highsev
#predictors <- c(predictors,"interaction")

d.mod.glmnet <- d.mod[,predictors]
d.mod.glmnet$y <- d.mod[,response]
#d.mod.glmnet <- complete.cases(d.mod.glmnet)



f <- as.formula(y ~ .*.)
y <- d.mod.glmnet$y
x <- model.matrix(f,d.mod.glmnet)[,-1]

cvfit <- cv.glmnet(x,y)
plot(cvfit)
coef(cvfit,s="lambda.min")




#### Model selection and averaging using MuMIn ####

library(snow)
library(MuMIn)


focalsp <- "CONIF.ALLSP"
d.sp.2.singlesp <- d.sp.2[d.sp.2$species==focalsp,]
d.mod <- merge(d.plot.3,d.sp.2.singlesp,all.x=TRUE) # data frame for modeling. Has regen-specific and plot-specific data for the species (or species group) specified above


### plotting climate vs. regen ###

ggplot(d.mod,aes(y=ppt.normal.highsev,x=diff.norm.ppt.z.highsev,color=regen.presab.all)) +
  geom_point(size=3) +
  labs(y="Normal precip",x="Postfire precip anomaly (3-year average)") +
  theme_bw(15)











predictors <- c("diff.norm.ppt.z.highsev","diff.norm.ppt.min.z.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev","count.highsev","SHRUB.highsev","GRASS.highsev")
#predictors <- c("perc.norm.ppt.highsev","perc.norm.ppt.min.highsev","adult.count","seed_tree_distance_general.highsev","rad.march.highsev","ppt.normal.highsev","ppt.normal.sq.highsev","count.highsev")

d.mod <- d.mod[,c("regen.presab.all","regen.presab.old","regen.count.all",predictors)]
d.mod <- d.mod[complete.cases(d.mod),]
 
d.mod$count.present.all <- d.mod$regen.presab.all * d.mod$count.highsev
d.mod$count.absent.all <- d.mod$count.highsev - d.mod$count.present.all

d.mod$count.present.old <- d.mod$regen.presab.old * d.mod$count.highsev
d.mod$count.absent.old <- d.mod$count.highsev - d.mod$count.present.old

d.mod$GRASS.highsev.prop <- d.mod$GRASS.highsev/100
d.mod$SHRUB.highsev.prop <- d.mod$SHRUB.highsev/100

d.mod$regen.count.all <- round(d.mod$regen.count.all)

#m.full <- lm(regen.presab.all ~ (diff.norm.ppt.z.highsev + diff.norm.ppt.min.z.highsev + ppt.normal.highsev)^2 + ppt.normal.sq.highsev + rad.march.highsev + adult.count,data=d.mod,na.action="na.fail")
#m.full <- glm(cbind(count.present.all,count.absent.all) ~ (perc.norm.ppt.highsev + perc.norm.ppt.min.highsev + ppt.normal.highsev)^2 + ppt.normal.sq.highsev + rad.march.highsev + adult.count,data=d.mod,na.action="na.fail",family="binomial")
#m.full <- glm(cbind(count.present.all,count.absent.all) ~ (diff.norm.ppt.z.highsev + diff.norm.ppt.min.z.highsev + ppt.normal.highsev)^2 + ppt.normal.sq.highsev + rad.march.highsev + adult.count,data=d.mod,na.action="na.fail",family="binomial")
#m.full <- betareg(GRASS.highsev.prop ~ (diff.norm.ppt.z.highsev + diff.norm.ppt.min.z.highsev + ppt.normal.highsev)^2 + ppt.normal.sq.highsev + rad.march.highsev + adult.count,data=d.mod,na.action="na.fail")
m.full <- glm(regen.count.all ~ (diff.norm.ppt.z.highsev + diff.norm.ppt.min.z.highsev + ppt.normal.highsev)^2 + ppt.normal.sq.highsev + rad.march.highsev + adult.count,data=d.mod,na.action="na.fail",family=poisson)



cl <- makeCluster(4, type = "SOCK")
clusterExport(cl,c("d.mod","betareg"))
a <- pdredge(m.full,cluster=cl, subset=(!(`ppt.normal.sq.highsev` & !`ppt.normal.highsev`)
                                            & !(`diff.norm.ppt.z.highsev` & `diff.norm.ppt.min.z.highsev`)
                                          & !`diff.norm.ppt.min.z.highsev`
                                        ))
stopCluster(cl)

imp <- data.frame(var=names(importance(a)),imp=importance(a))
best.model.num <- as.character(row.names(a[1,]))
best.model <- get.models(a,best.model.num)[[1]]
summary(best.model)
#best.model <- get.models(a,"29")[[1]]
























#### Make predictions ####

## Variation over normal precip, with two levels of post-fire precip

ppt.norm.seq <- seq(from=500,to=2000,length.out=100)

# newdat.pptnorm <- data.frame(
#                     ppt.normal.highsev = rep(ppt.norm.seq,2),
#                     ppt.normal.sq.highsev = rep(ppt.norm.seq^2,2),
#                     diff.norm.ppt.z.highsev = c(rep(-0.8,100),rep(0.55,100)),
#                     diff.norm.level = c(rep("low",100),rep("high",100)),
#                     adult.count = 1,
#                     rad.march.highsev = 6000)
                    
newdat.pptnorm <- data.frame(
  ppt.normal.highsev = rep(ppt.norm.seq,2),
  ppt.normal.sq.highsev = rep(ppt.norm.seq^2,2),
  diff.norm.ppt.min.z.highsev = c(rep(-1.5,100),rep(-0.2,100)),
  diff.norm.level = c(rep("low",100),rep("high",100)),
  adult.count = 1,
  rad.march.highsev = 6000)



ppt.norm.pred.obj <- predict(best.model,newdata=newdat.pptnorm,type="link",se.fit=TRUE)
ppt.norm.pred <- data.frame(fit=ppt.norm.pred.obj)
ppt.norm.pred$high <- ppt.norm.pred$fit + ppt.norm.pred$se.fit * 1.96
ppt.norm.pred$low <- ppt.norm.pred$fit - ppt.norm.pred$se.fit * 1.96


ppt.norm.pred$fit <- best.model$family$linkinv(ppt.norm.pred$fit)
ppt.norm.pred$low <- best.model$family$linkinv(ppt.norm.pred$low)
ppt.norm.pred$high <- best.model$family$linkinv(ppt.norm.pred$high)

pptnorm.dat.pred <- cbind(newdat.pptnorm,ppt.norm.pred)

ggplot(pptnorm.dat.pred,aes(x=ppt.normal.highsev,y=fit,color=diff.norm.level)) +
  geom_line(size=1) +
  #geom_ribbon(aes(ymin=low,ymax=high,fill=diff.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Postfire precip anomaly"),color=guide_legend(title="Postfire precip anomaly"))
  
  
## Variation over postfire anomaly, with two levels of normal precip

diff.norm.seq <- seq(from=-0.8,to=0.5,length.out=100)

newdat.diffnorm <- data.frame(
  ppt.normal.highsev = c(rep(750,100),rep(1750,100)),
  ppt.normal.sq.highsev = c(rep(750^2,100),rep(1750^2,100)),
  ppt.norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.ppt.z.highsev = rep(diff.norm.seq,2),
  adult.count = 1,
  rad.march.highsev = 6000)

ppt.norm.pred.obj <- predict(best.model,newdata=newdat.diffnorm,type="link",se.fit=TRUE)
#ppt.norm.pred <- data.frame(fit=ppt.norm.pred.obj)
ppt.norm.pred <- ppt.norm.pred.obj
ppt.norm.pred$high <- ppt.norm.pred$fit + ppt.norm.pred$se.fit * 1.96
ppt.norm.pred$low <- ppt.norm.pred$fit - ppt.norm.pred$se.fit * 1.96

ppt.norm.pred$fit <- best.model$family$linkinv(ppt.norm.pred$fit)
ppt.norm.pred$low <- best.model$family$linkinv(ppt.norm.pred$low)
ppt.norm.pred$high <- best.model$family$linkinv(ppt.norm.pred$high)

diffnorm.dat.pred <- cbind(newdat.diffnorm,ppt.norm.pred)

ggplot(diffnorm.dat.pred,aes(x=diff.norm.ppt.z.highsev,y=fit,color=ppt.norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=low,ymax=high,fill=ppt.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip"))












## predictions from GLMnet ## (not yet working)


diff.norm.seq <- seq(from=-1.5,to=0,length.out=100)

newdat.diffnorm <- data.frame(
    diff.norm.ppt.z.highsev = 0,
    diff.norm.ppt.min.z.highsev = rep(diff.norm.seq,2),
    ppt.normal.highsev = c(rep(750,100),rep(1750,100)))


diff.norm.pred <- predict(best.model,newdata=newdat.diffnorm,se.fit=TRUE) #a.avg, best.model, cvfit
pred.mat <- as.matrix(newdat.diffnorm)
diff.norm.pred <- predict(cvfit,newx=pred.mat,s="lambda.min") # cvfit
diff.norm.pred$low <- diff.norm.pred$fit + diff.norm.pred$se.fit * 1.96
diff.norm.pred$high <- diff.norm.pred$fit - diff.norm.pred$se.fit * 1.96

diffnorm.dat.pred <- cbind(newdat.diffnorm,diff.norm.pred)

ggplot(diffnorm.dat.pred,aes(x=diff.norm.ppt.min.z.highsev,y=fit,color=ppt.norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=low,ymax=high,fill=ppt.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip"))







#### Plot-level exploration ####


d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV > 3),]

sp <- "CONIF.ALLSP"
d.sp.curr <- d.sp[d.sp$species==sp,]
d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB")
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














#### Plot-level analysis ####


#d.sp
#d.plot

d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV > 3),]

sp <- "PINUS.ALLSP"
d.sp.curr <- d.sp[d.sp$species==sp,]
d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB")
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

#d.c <- d.c[complete.cases(d.c$HARDWOOD),]

library(lme4)

library(glmmADMB)
library(betareg)

d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)


#for pres/ab
m.full <- glmer(regen.presab.all ~ ppt.normal.sq_c +ppt.normal_c * diff.norm.ppt.z_c + (1|Fire),family="binomial",na.action="na.fail",data=d.c)
m.full <- glmmadmb(regen.presab.all.01 ~ ppt.normal.sq_c + ppt.normal_c * diff.norm.ppt.z_c + (1|Fire),family="binomial",data=d.c)

#for cover
#m.full <- betareg(SHRUB.pt ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal.sq_c + rad.march_c + seed_tree_distance_general_c,data=d.c)
#m.full <- glmmadmb(SHRUB.pt ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal.sq_c + rad.march_c + seed_tree_distance_general_c + (1|Fire),family="beta",data=d.c)

library(MuMIn)
library(snow)

cl <- makeCluster(4, type = "SOCK")
clusterExport(cl,c("d.c","glmer","glmmadmb","admbControl"))
a <- pdredge(m.full,cluster=cl, subset=(!(`ppt.normal.sq_c` & !`ppt.normal_c`)))
stopCluster(cl)

imp <- data.frame(var=names(importance(a)),imp=importance(a))
best.model.num <- as.character(row.names(a[1,]))
best.model <- get.models(a,best.model.num)[[1]]
summary(best.model)
best.model <- get.models(a,"15")[[1]]




### make predictions
library(arm)
library(merTools)

## for glmer, need to simulate coefs

m.sim <- sim(best.model)



diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100)


newdat <- data.frame(
  ppt.normal_c = c(rep(-1,100),rep(1,100)),
  ppt.normal.sq_c = c(rep(-1,100),rep(1,100)),
  ppt.norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.ppt.z_c = rep(diff.norm.seq,2),
  diff.norm.ppt.min.z_c = rep(diff.norm.seq,2),
  rad.march_c = 0,
  seed_tree_distance_general_c = -1)


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

mod.vars <- colnames(fixef(m.sim))

newdat.ordered <- newdat[,mod.vars]

# for each row of newdat.ordered, multiply data by all coef values, sum across, and compute median and conf int
preds <- matrix(nrow=nrow(newdat.ordered),ncol=3)
for(i in 1:nrow(newdat.ordered)) {
  
  newdat.row <- as.matrix(newdat.ordered[i,])
  
  prod <- sweep(fixef(m.sim),MARGIN=2,newdat.row,`*`)
  sums <- rowSums(prod)
  
  fit <- median(sums)
  upr <- quantile(sums,probs=0.975)
  lwr <- quantile(sums,probs=0.025)
  
  preds[i,] <- c(fit,lwr,upr)
  
}

preds <- inv.logit(preds)
colnames(preds) <- c("fit","lwr","upr")
preds <- as.data.frame(preds)


dat.preds <- cbind(newdat,preds)


ggplot(dat.preds,aes(x=diff.norm.ppt.z_c,y=fit,color=ppt.norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=ppt.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip"))



### Make predictions for glmmADMB


diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100)

newdat <- data.frame(
  ppt.normal_c = c(rep(-1,100),rep(1,100)),
  ppt.normal.sq_c = c(rep(-1,100),rep(1,100)),
  ppt.norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.ppt.z_c = rep(diff.norm.seq,2),
  diff.norm.ppt.min.z_c = rep(diff.norm.seq,2),
  rad.march_c = 0,
  seed_tree_distance_general_c = -1)








preds <- predict(best.model,newdat,interval="confidence",type="response")

dat.preds <- cbind(newdat,preds)

ggplot(dat.preds,aes(x=diff.norm.ppt.z_c,y=fit,color=ppt.norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=ppt.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip"))








#### Plot-level analysis without MuMIn (single basic model for all functional groups) ####


#d.sp
#d.plot

d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV > 3),]

sp <- "PINUS.ALLSP"
d.sp.curr <- d.sp[d.sp$species==sp,]
d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB")
d <- d[complete.cases(d[,vars.focal]),]
d.c <- center.df(d,vars.leave)
#d.c <- d

d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2

d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)


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

#d.c <- d.c[complete.cases(d.c$HARDWOOD),]

library(lme4)

library(glmmADMB)
library(betareg)

#for pres/ab
m.full <- glmer(regen.presab.all.01 ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal_c.sq + (1|Fire),family="binomial",na.action="na.fail",data=d.c)
m.full <- glmmadmb(regen.presab.all.01 ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal.sq_c + (1|Fire),family="binomial",data=d.c)

#for cover
#m.full <- betareg(SHRUB.pt ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal.sq_c + rad.march_c + seed_tree_distance_general_c,data=d.c)
#m.full <- glmmadmb(SHRUB.pt ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal.sq_c + rad.march_c + seed_tree_distance_general_c + (1|Fire),family="beta",data=d.c)

library(MuMIn)
library(snow)



### make predictions
library(arm)
library(merTools)

## for glmer, need to simulate coefs

m.sim <- sim(m.full)



diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100)


newdat <- data.frame(
  ppt.normal_c = c(rep(-1,100),rep(1,100)),
  ppt.normal.sq_c = c(rep(-1,100),rep(1,100)),
  ppt.norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.ppt.z_c = rep(diff.norm.seq,2),
  diff.norm.ppt.min.z_c = rep(diff.norm.seq,2),
  rad.march_c = 0,
  seed_tree_distance_general_c = -1)


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

mod.vars <- colnames(fixef(m.sim))

newdat.ordered <- newdat[,mod.vars]

# for each row of newdat.ordered, multiply data by all coef values, sum across, and compute median and conf int
preds <- matrix(nrow=nrow(newdat.ordered),ncol=3)
for(i in 1:nrow(newdat.ordered)) {
  
  newdat.row <- as.matrix(newdat.ordered[i,])
  
  prod <- sweep(fixef(m.sim),MARGIN=2,newdat.row,`*`)
  sums <- rowSums(prod)
  
  fit <- median(sums)
  upr <- quantile(sums,probs=0.975)
  lwr <- quantile(sums,probs=0.025)
  
  preds[i,] <- c(fit,lwr,upr)
  
}

preds <- inv.logit(preds)
colnames(preds) <- c("fit","lwr","upr")
preds <- as.data.frame(preds)


dat.preds <- cbind(newdat,preds)


ggplot(dat.preds,aes(x=diff.norm.ppt.z_c,y=fit,color=ppt.norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=ppt.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip"))



### Make predictions for glmmADMB


diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100)

newdat <- data.frame(
  ppt.normal_c = c(rep(-1,100),rep(1,100)),
  ppt.normal.sq_c = c(rep(-1,100),rep(1,100)),
  ppt.norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.ppt.z_c = rep(diff.norm.seq,2),
  diff.norm.ppt.min.z_c = rep(diff.norm.seq,2),
  rad.march_c = 0,
  seed_tree_distance_general_c = -1)

preds <- predict(m.full,newdat,interval="confidence",type="response")

dat.preds <- cbind(newdat,preds)

ggplot(dat.preds,aes(x=diff.norm.ppt.z_c,y=fit,color=ppt.norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=ppt.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip"))












#### trying to use brms ####


library(brms)

d.plot <- d.plot[(d.plot$survey.years.post %in% c(4,5)) & (d.plot$FIRE_SEV > 3),]


sp.opts <- c("PINUS.ALLSP","SHADE.ALLSP","HDWD.ALLSP","PIPO","ABCO","ABMA","CONIF.ALLSP","PSME","PILA","CADE27","PIJE","PIPJ")
cover.opts <- c("COV.SHRUB","COV.GRASS","COV.FORB","COV.HARDWOOD","COV.CONIFER")
cover.opts <- NULL
sp.opts <- c(cover.opts,sp.opts)

m.p <- list()
m <- list()

for(sp in sp.opts) {
  
      cat("Running model for: ",sp,"\n")
  
      if(sp %in% cover.opts) {
        d.sp.curr <- d.sp[d.sp$species=="PIPO",] #pick any species; for cover it doesn't matter; just need to thin to one row per plots
      } else {
        d.sp.curr <- d.sp[d.sp$species==sp,]
      }
      
      d <- merge(d.plot,d.sp.curr,all.x=TRUE,by="Regen_Plot")
      vars.leave <- c("Year.of.Fire","FORB","SHRUB","GRASS","CONIFER","HARDWOOD","FIRE_SEV","Year","firesev","fire.year","survey.years.post","regen.count.young","regen.count.old","regen.count.all","regen.presab.young","regen.presab.old","regen.presab.all")
      vars.focal <- c("ppt.normal","diff.norm.ppt.z","ppt.normal.sq","rad.march","seed_tree_distance_general","SHRUB","tmean.post","tmean.normal","diff.norm.tmean.z","diff.norm.tmean.max.z")
      d <- d[complete.cases(d[,vars.focal]),]
      d.c <- center.df(d,vars.leave)

      d.c$ppt.normal_c.sq <- d.c$ppt.normal_c^2
      d.c$tmean.normal_c.sq <- d.c$tmean.normal_c^2
      
      d.c$regen.presab.all.01 <- ifelse(d.c$regen.presab.all == TRUE,1,0)
      d.c$regen.presab.old.01 <- ifelse(d.c$regen.presab.old == TRUE,1,0)
      
      
      vars.focal.c <- paste0(vars.focal[-6],"_c")
      pairs(d.c[,vars.focal.c])
      
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
      
      
      if(sp %in% cover.opts) {
        
        sp.cov <- substr(sp,5,100)
        sp.cov <- paste0(sp.cov,".pt")
        
        d.c$cov.response <- d.c[,sp.cov]
        
        m[[sp]] <- brm(cov.response ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal_c.sq + tmean.normal_c * diff.norm.tmean.z_c + ppt.normal_c.sq,family="Beta",data=d.c,iter=2000,control = list(adapt_delta = 0.90),cores=3,chains=3)
        
        
      } else {
        #m[[sp]] <- brm(regen.presab.all.01 ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal_c.sq + tmean.normal_c * diff.norm.tmean.z_c + ppt.normal_c.sq,family="bernoulli",data=d.c,warmup=3000,iter=6000,control = list(adapt_delta = 0.90),cores=3,chains=3)
        m[[sp]] <- brm(regen.count.all.int ~ ppt.normal_c * diff.norm.ppt.z_c + ppt.normal_c.sq + tmean.normal_c * diff.norm.tmean.z_c + ppt.normal_c.sq + (1|Fire),family="zero_inflated_negbinomial",data=d.c,warmup=3000,iter=6000,control = list(adapt_delta = 0.90),cores=3,chains=3)
      }
}

for(sp in names(m)) {
  print("\n\n\n")
  print(sp)
  print("\n")
  print(summary(m[[sp]]))
  print("\n")
  print(loo(m[[sp]]))
  print("\n")
}


#### get loos
m.loos <- data.frame()
for(sp in names(m)) {
  a <- loo(m[[sp]])
  m.loo <- data.frame(sp=sp,LOOIC=a$looic,SE=a$se_looic)
  m.loos <- rbind(m.loos,m.loo)
}

m.loos


diff.norm.seq <- seq(from=-1.5,to=1.5,length.out=100)

newdat <- data.frame(
  ppt.normal_c = c(rep(-1,100),rep(1,100)),
  ppt.normal_c.sq = c(rep(1,100),rep(1,100)),
  ppt.norm.level = c(rep("low",100),rep("high",100)),
  diff.norm.ppt.z_c = rep(diff.norm.seq,2),
  diff.norm.ppt.min.z_c = rep(diff.norm.seq,2),
  tmean.normal_c = 0,
  diff.norm.tmean.z_c = 0,
  tmean.normal_c.sq = 0
  #diff.norm.ppt.min.z_c = rep(diff.norm.seq,2),
  #rad.march_c = 0,
  #seed_tree_distance_general_c = -1
)




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



### put all predictions in a list, with a column for species


sp.plots <- list()
dat.preds <- data.frame()

for(sp in names(m)) {
  
  
  m.sp <- m[[sp]]
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
  
  preds <- inv.logit(preds)
  #preds <- exp(preds)
  colnames(preds) <- c("fit","lwr","upr")
  preds <- as.data.frame(preds)
  
  sp.grp <- sp
  
  dat.preds.sp <- cbind(newdat,preds,sp.grp)
  
  dat.preds <- rbind(dat.preds,dat.preds.sp)


}

ggplot(dat.preds,aes(x=diff.norm.ppt.min.z_c,y=fit,color=ppt.norm.level)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr,fill=ppt.norm.level),alpha=0.3,color=NA) +
  guides(fill=guide_legend(title="Normal precip"),color=guide_legend(title="Normal precip")) +
  facet_wrap(~sp.grp,scales="free",ncol=5)




####testing
library(pROC)
predicted <- predict(m.bin.9)
observed <- d.2015$mort.bin
auc(observed,predicted)

