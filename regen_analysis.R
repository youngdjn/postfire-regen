setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

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

# if no data on seed tree distance (or it was recorded as being at/beyond the limit of the laser) used remote-sensing-based value
d.plot$seed.tree.any.comb <- ifelse(is.na(d.plot$seed.tree.any) | (d.plot$seed.tree.any >= 150),d.plot$dist.to.low,d.plot$seed.tree.any)

# quadratic climate terms
d.plot$ppt.normal.sq <- d.plot$ppt.normal^2
d.plot$tmean.normal.sq <- d.plot$tmean.normal^2

# variable for regen presence/absence
d.sp$regen.presab.young <- ifelse(d.sp$regen.count.young > 0,TRUE,FALSE)
d.sp$regen.presab.old <- ifelse(d.sp$regen.count.old > 0,TRUE,FALSE)


#### 1. Assign each plot a topoclimatic category ####

fires <- unique(d.plot$Fire)
d.plot$precip.category
d.plot$rad.category #radiation

for(fire in fires) {
  
  print(fire)
  
  ### Precipitation categories
  # determine what the precipitation breakpoints should be (here just using median)
  breaks <- quantile(d.plot[d.plot$Fire == fire,]$ppt.normal,probs=c(0.5),na.rm=TRUE)
  # categorize plots based on where they fall between the breakpoints  
  categories <- categorize(d.plot[d.plot$Fire==fire,]$ppt.normal,breaks,name="P")
  # store it into the plot data.frame
  d.plot[d.plot$Fire==fire,"precip.category"] <- categories
  
  ### Radiation categories #! note straylor cub do not have radiation yet
  # determine what the breakpoints should be (here just using median)
  breaks <- quantile(d.plot[d.plot$Fire == fire,]$rad.march,probs=c(0.5),na.rm=TRUE)
  # categorize plots based on where they fall between the breakpoints  
  categories <- categorize(d.plot[d.plot$Fire==fire,]$rad.march,breaks,name="R")
  # store it into the plot data.frame
  d.plot[d.plot$Fire==fire,"rad.category"] <- categories

}
  

### One variable reflecting the factorial combination of topoclimatic categories
d.plot$topoclim.cat <- paste(d.plot$precip.category,d.plot$rad.category,sep="_")
  
  
  
