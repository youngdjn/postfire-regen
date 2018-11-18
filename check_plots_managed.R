setwd("~/UC Davis/Research Projects/Post-fire regen/Dev/postfire-regen")

library(sf)
library(tidyverse)


## For each plot, get all FACTS that overlap ##
# If fire is POWER, use Accmp date, else use Compl date
# Thin to FACTS within year of fire and year of survey

## Keep only FACTS that are planting, salvage (except helicopter)
## Concatenate into year and action

salvage = c("Salvage Cut (intermediate treatment, not regeneration)","Stand Clearcut (w/ leave trees) (EA/RH/FH)")
planting = c("Plant Trees")
release = c("Tree Release and Weed")

facts = st_read("../data_non_synced/CA_Activity_merged.shp")

facts = st_transform(facts,st_crs(plots))

facts_focal = facts %>%
  filter(ACTIV %in% c(salvage,planting,release))

facts_focal = facts_focal %>%
  mutate(yr_accomp = str_sub(DATE_A,1,4),
         yr_compl = str_sub(DATE_C,1,4)) %>%
  filter(!(is.null(yr_accomp)&is.null(yr_compl)))

plots = st_read("../shapefiles/plots_535.gpkg")
plots$rel_yrs = NA
plots$plant_yrs = NA
plots$salv_yrs = NA

intersection = st_intersects(plots,facts_focal)


for(i in 1:nrow(plots)) {
  
  plot = plots[i,]
  yr_fire = plot$Year.of.Fire
  yr_survey = yr_fire + plot$survey.years.post
  fire = plot$Fire
  
  facts_intersect_rows = intersection[[i]]
  facts_intersect = facts_focal[facts_intersect_rows,]
  
  facts_intersect = facts_intersect %>%
    mutate(yr_mgd = ifelse(fire == "Power" & is.na(yr_compl),yr_accomp,yr_compl)) %>%
    filter(!is.na(yr_mgd)) %>%
    filter(METHOD != "Helicopter") %>%
    filter((yr_mgd >= yr_fire) & (yr_mgd <= yr_survey))
  
  facts_int_rel = facts_intersect %>%
    filter(ACTIV %in% release)
  
  facts_int_plant = facts_intersect %>%
    filter(ACTIV %in% planting)
  
  facts_int_salv = facts_intersect %>%
    filter(ACTIV %in% salvage)
  
  rel_yrs = paste(facts_int_rel$yr_mgd,collapse = ", ")
  plant_yrs = paste(facts_int_plant$yr_mgd,collapse = ", ")
  salv_yrs = paste(facts_int_salv$yr_mgd,collapse = ", ")
  
  plots[i,"rel_yrs"] = rel_yrs
  plots[i,"plant_yrs"] = plant_yrs
  plots[i,"salv_yrs"] = salv_yrs
  
}


plots = plots %>%
  mutate(released = (rel_yrs != "") * 1,
         salvaged = (salv_yrs != "") * 1,
         planted = (plant_yrs != "") * 1) %>%
  mutate(managed = (released | salvaged | planted) * 1)

st_write(plots,"../shapefiles/plots_w_mgt.gpkg")


# export list of apparent facts-managed plots
plots_mgd = plots %>%
  filter(managed == 1)

write.csv(plots_mgd,"data_intermediate/plots_exclude_FACTS.csv")
