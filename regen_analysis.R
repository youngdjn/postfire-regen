setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")


#### 0. Read in and clean data, thin to focal plots ####

# open intermediate data files
d.plot <- read.csv("data_intermediate/plot_level.csv",header=T,stringsAsFactors=FALSE)
d.sp <- read.csv("data_intermediate/speciesXplot_level.csv",header=T,stringsAsFactors=FALSE)

# thin to 5-year post fire plots #! need to allow 4-year fires once we add them and fix climate and regen summarization
d.plot <- d.plot[d.plot$survey.years.post == 5,]

# only Sierra Nevada fires #! need to add new fires to this list when they're in the dataset
sierra.fires <- c("STRAYLOR","CUB","RICH","DEEP","MOONLIGHT","ANTELOPE","BTU LIGHTENING","HARDING","BASSETTS","PENDOLA","AMERICAN RIVER","RALSTON","FREDS","SHOWERS","POWER")
d.plot <- d.plot[d.plot$Fire %in% sierra.fires,]

# fix incorrectly-named variable
d.plot$FORB <- d.plot$FORBE

# only medium-high and high severity (using field-assessed severity)
d.plot <- d.plot[d.plot$FIRE_SEV %in% c(4,5),]

# if no data on seed tree distance (or it was recorded as being at/beyond the limit of the laser) used remote-sensing-based value
d.plot$seed.tree.any.comb <- ifelse(is.na(d.plot$seed.tree.any) | (d.plot$seed.tree.any >= 150),d.plot$dist.to.low,d.plot$seed.tree.any)

# quadratic climate terms
d.plot$ppt.normal.sq <- d.plot$ppt.normal^2
d.plot$tmean.normal.sq <- d.plot$tmean.normal^2

# variable for regen presence/absence
d.sp$regen.presab.young <- ifelse(d.sp$regen.count.young > 0,TRUE,FALSE)
d.sp$regen.presab.old <- ifelse(d.sp$regen.count.old > 0,TRUE,FALSE)


