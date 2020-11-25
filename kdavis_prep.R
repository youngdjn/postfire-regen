library(tidyverse)
library(sf)

excl = read_csv("data_intermediate/plots_exclude_FACTS.csv")
excl = unique(excl$Regen_Plot)

plot = read_csv("data_intermediate/plot_level_Davis.csv") %>%
  select(Regen_Plot,Year,everything()) %>%
  filter(!(Regen_Plot %in% excl)) %>%
  filter(Year > 2015) %>%
  filter(FIRE_SEV %in% c(4,5))
  

sp = read_csv("data_intermediate/speciesXplot_level_Davis.csv") %>%
  filter(Regen_Plot %in% unique(plot$Regen_Plot))


### Prep species ###

# focal species
sp_foc = sp %>%
  filter(species %in% c("PIPJ","PILA","ABCO","ABMA","CADE","PICO","PSME")) %>%
  select(-seed.tree.sp) %>%
  select(-adult.count,adult.ba) # for high-severity delivery only, don't include these

# make wide
sp_wide = sp_foc %>%
  pivot_wider(id_cols = Regen_Plot,
              names_from = species,
              values_from = c(regen.count,subsampled)) %>% # ,adult.count,adult.ba
  mutate(across(starts_with("subsampled_"),  ~ifelse(.==TRUE,15,60) ))

names(sp_wide) = str_replace_all(names(sp_wide),"subsampled_","plot_size_")
names(sp_wide) = str_replace_all(names(sp_wide),"regen.","")

sp_wide = sp_wide %>%
  select(plot_id = Regen_Plot,everything())

### Prep plot ###

## convert albers to geo
p_sf = st_as_sf(plot,coords=c("Easting","Northing"),crs=26910) %>% st_transform(4326)
coords = st_coordinates(p_sf) %>% as.data.frame
p_sf$longitude = coords$X
p_sf$latitude = coords$Y
st_geometry(p_sf) = NULL

st_write(p_sf %>% select(Regen_Plot),"test.gpkg")



p = p_sf %>%
  mutate(contributor = "DJNY",
         time_since_fire = Year - Year.of.Fire,
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

d = left_join(p,sp_wide) %>%
  mutate(distance_seed_source = ifelse(distance_seed_source == 999,NA,distance_seed_source))

#write
write_csv(d,"data_intermediate/compiled_kdavis.csv")
