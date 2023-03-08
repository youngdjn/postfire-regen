library(tidyverse)
library(sf)
library(here)

excl = read_csv("data_intermediate/plots_exclude_FACTS.csv")
excl = unique(excl$Regen_Plot)

plot = read_csv("data_intermediate/plot_level_Science.csv") %>%
  select(Regen_Plot,Year,everything()) %>%
  filter(!(Regen_Plot %in% excl)) %>%
  #filter(Year > 2015) %>%
  filter(FIRE_SEV %in% c(4,5))
  

sp = read_csv("data_intermediate/speciesXplot_level_Science.csv") %>%
  filter(Regen_Plot %in% unique(plot$Regen_Plot))


### Prep species ###

# focal species
sp_foc = sp %>%
  filter(species %in% c("PILA","ABCO","ABMA","CADE","PICO","PSME","ABIES","CADE27", "CONIFER", "JUCA7","JUNIPERUS", "JUOC", "PICO", "PILA", "PIMO3", "PINUS", "PIPJ", "PISA2", "PMSE","PSMA","PSME","TAXUS","TOCA")) %>%
  select(-seed.tree.sp) %>%
  select(-adult.count,adult.ba) # for high-severity delivery only, don't include these

## scale back up the subsampled plots, then sum across species
sp_foc = sp_foc |>
  mutate(regen.dens = ifelse(subsampled, regen.count/15, regen.count/60)) |>
  group_by(Regen_Plot) |>
  summarize(conif_seedl_sqm = sum(regen.dens)) |>
  rename(plot_id = Regen_Plot)


### Prep plot ###

## convert albers to geo
p_sf = st_as_sf(plot,coords=c("Easting","Northing"),crs=26910) %>% st_transform(4326)
coords = st_coordinates(p_sf) %>% as.data.frame
p_sf$longitude = coords$X
p_sf$latitude = coords$Y
st_geometry(p_sf) = NULL

#st_write(p_sf %>% select(Regen_Plot),"test.gpkg")



p = p_sf %>%
  mutate(time_since_fire = Year - Year.of.Fire,
         fire_severity_category = "high") %>%
  select(plot_id = Regen_Plot,
         contributor,
         longitude,
         latitude,
         fire_year = Year.of.Fire,
         sample_year = Year,
         time_since_fire,
         distance_seed_source = seed_tree_distance_general,
         fire_severity_category,
         shrub_cover = SHRUB,
         grass_cover = GRASS,
         forb_cover = FORB,
         canopy_cover = LIVE_OCC
         )

# add species

d = left_join(p,sp_foc) %>%
  mutate(distance_seed_source = ifelse(distance_seed_source == 999,NA,distance_seed_source))

#write
write_csv(d,"data_intermediate/compiled_welch-young.csv")
