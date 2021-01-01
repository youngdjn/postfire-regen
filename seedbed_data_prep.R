library(tidyverse)
library(sf)

excl = read_csv("data_intermediate/plots_exclude_FACTS.csv")
excl = unique(excl$Regen_Plot)

excl2 = read_csv("data_intermediate/plots_exclude.csv") %>%
  filter(!is.na(Exclude)) %>%
  pull(Regen_Plot)

excl2 = unique(excl2)

excl3 = read_csv("data_intermediate/PlotsDeleted_Shive_etal_2018.csv") %>%
  pull(Regen_Plot)

excl3 = unique(excl3)

### Find plots with planted trees and exclude them
seedl <- read.csv("data_survey/Compiled/tree_regen_full.csv")
seedl.p <- seedl[seedl$seed_veg_plant == "P",]
plots.planted <- unique(seedl.p$Regen_Plot)

excl = c(excl, excl2, excl3, plots.planted) %>% unique

### Read plot-level data and filter to focal
plot = read_csv("data_intermediate/plot_level_full.csv") %>%
  select(Regen_Plot,Year,everything()) %>%
  filter(!(Regen_Plot %in% excl)) %>%
  #filter(Year > 2015) %>%
  filter(FIRE_SEV %in% c(4,5))


### Read and summarize seedling data
seedl <- read.csv("data_survey/Compiled/tree_regen_full.csv")

## filter to only kept plots
seedl = seedl %>%
  filter((Regen_Plot %in% plot$Regen_Plot))


seedl_genus = seedl %>%
  mutate(sp_coarse = recode(Species, "PIPO" = "Pine", "PSME" = "Dougfir", "ABCO" = "FirCedar", "CADE27" = "FirCedar", "PILA" = "Pine", "ABMA" = "FirCedar", "PIJE" = "Pine", "PINUS" = "Pine", "ABIES" = "FirCedar", "PICO" = "Pine", "PSMA" = "Dougfir", "PIMO3" = "Pine")) %>%
  filter(sp_coarse %in% c("Pine", "Dougfir", "FirCedar")) %>%
  mutate(sp_coarse = sp_coarse %>% droplevels()) %>%
  select(Regen_Plot, sp_coarse,
         age_00 = X0yr,
         age_01 = Ct_1yr,
         age_02 = Ct_2yr,
         age_03 = Ct_3yr,
         age_04 = Ct_4yr,
         age_05 = Ct_5yr,
         age_06 = X6yr,
         age_07 = X7yr,
         age_08 = X8yr,
         age_09 = X9yr) %>%
  group_by(Regen_Plot, sp_coarse) %>%
  summarize(across(age_00:age_09, sum)) %>%
  ungroup() %>%
  mutate(age_all = rowSums(select(.,age_00:age_09), na.rm=TRUE)) %>%
  
  ### thin to plot with n seedl > 5 by species
  filter(age_all > 4) %>%
  
  ## sum across fires
  mutate(fire = str_sub(Regen_Plot,1,3)) %>%
  mutate(fire = ifelse(fire == "CUB" & str_length(Regen_Plot) < 8, "CB2",fire)) %>%
  select(-Regen_Plot) %>%
  group_by(fire, sp_coarse) %>%
  summarize(across(age_00:age_09, sum, na.rm=TRUE),
            n_plots = n()) %>%
  ungroup() %>%
  
  ## only keep fires/species combos where data came from > 2 plots
  filter(n_plots > 2)


### Make long for plotting
seedl_long = seedl_genus %>%
  pivot_longer(c(-fire,-sp_coarse,-n_plots), names_to = "age", values_to = "count") %>%
  mutate(age = str_sub(age,5,7) %>% as.numeric) %>%
  mutate(fire = as.character(fire))


### Plot

ggplot(seedl_long %>% filter(sp_coarse == "Pine"), aes(x=age,y=count)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
  facet_wrap(vars(fire)) +
  theme_bw() +
  labs(title="Pines")

ggplot(seedl_long %>% filter(sp_coarse == "FirCedar"), aes(x=age,y=count)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
  facet_wrap(vars(fire)) +
  theme_bw() +
  labs(title="FirCedar")

ggplot(seedl_long %>% filter(sp_coarse == "Dougfir"), aes(x=age,y=count)) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
  facet_wrap(vars(fire)) +
  theme_bw() +
  labs(title="Dougfir")


## All species together
ggplot(seedl_long, aes(x=age,y=count, fill=sp_coarse)) +
  geom_bar(stat="identity", position = position_dodge()) +
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) +
  scale_fill_viridis_d(end=0.8, begin = 0.2) +
  facet_wrap(vars(fire), scales = "free_y") +
  theme_bw() +
  labs(x = "seedling age",
       y = "seedling count")



# ral, rus, elk have PIAT



