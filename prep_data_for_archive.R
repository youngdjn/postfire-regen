setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")


d.regen <- read.csv("../data_intermediate_processing_local/tree_summarized_sp.csv",header=TRUE,stringsAsFactors=FALSE)
d.plot <- read.csv("data_intermediate/plot_level.csv",header=TRUE,stringsAsFactors=FALSE)


resprout <- read.csv("../data_survey/Compiled/Resprouts.csv")
seedl <- read.csv("../data_survey/Compiled/tree_regen.csv")


## append seedling and resprout heights
seedl.ht <- seedl[,c("Regen_Plot","Species","tallest_ht_cm")]
respr.ht <- resprout[,c("Regen_Plot","Species","tallest_ht_cm")]
regen.ht <- rbind(seedl.ht,respr.ht)

# in case there are multiple seedling entries per species-plot combination, aggregate
regen.ht.agg <- aggregate(regen.ht$tallest_ht_cm,by=list(regen.ht$Regen_Plot,regen.ht$Species),FUN=max)
names(regen.ht.agg) <- c("Regen_Plot","species","tallest_ht_cm")

# merge in the seedling heights
d.regen <- merge(d.regen,regen.ht.agg,by=c("Regen_Plot","species"),all.x=TRUE)


#for both plot-leven and species-level data, remove the previously-surveyed plots
#! need to modify this to not filter out the Richter plots
d.plot <- d.plot[nchar(d.plot$Regen_Plot) < 10, ]
d.regen <- d.regen[nchar(d.regen$Regen_Plot) < 10, ]


#how many new burned plots? (not revisit)
sum((d.plot$new_or_revisit == "New") & (d.plot$FIRE_SEV >1),na.rm=TRUE)

#how many new control plots? (not revisit)
sum((d.plot$new_or_revisit == "New") & (d.plot$FIRE_SEV <2),na.rm=TRUE)

#how many revisits?
sum(d.plot$new_or_revisit == "Revisit",na.rm=TRUE)



#columns to keep in plot data
keep.cols <- c("Regen_Plot","Year","Fire","Year.of.Fire","Easting","Northing","ROCK","SHRUB","GRASS","FIRE_SEV","seed_tree_distance_general","seed_wall_distance","dominant_tree_1","dominant_tree_2","dominant_tree_3","dominant_shrub_1","dominant_shrub_2","dominant_shrub_ht_cm","new_or_revisit","burn_or_control","Year")
d.plot <- d.plot[,keep.cols]
new.col.names <- c("plot_ID","survey_year","fire_name","fire_year","easting","northing","cover_rock","cover_shrub","cover_grass","fire_severity","seed_tree_distance","seed_wall_distance","dominant_tree_1","dominant_tree_2","dominant_tree_3","dominant_shrub_1","dominant_shrub_2","dominant_shrub_ht","new_or_revisit","burn_or_control")
names(d.plot) <- new.col.names


#columns to keep in regen data
count.cols <- paste0("count.",0:12,"yr")
keep.cols <- c("Regen_Plot","species",count.cols,"unk_yr","tallest_ht_cm","surviving.trees.count","surviving.trees.ba")
d.regen <- d.regen[,keep.cols]
count.cols <- paste0("count_",0:12,"yr")
new.col.names <- c("plot_ID","species",count.cols,"count_unk_yr","tallest_ht","surviving_tree_count","surviving_trees_ba")
names(d.regen) <- new.col.names

#turn NAs in the count columns to 0



## write CSV files
write.csv(d.plot,"../data_archive/JFSP/plot_data.csv")
write.csv(d.regen,"../data_archive/JFSP/tree_data.csv")








