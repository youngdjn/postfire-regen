library(tidyverse)
library(sf)

modal <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

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
  filter(FIRE_SEV %in% c(4,5)) %>%
  # Pull in fire name and make unique
  mutate(fire = paste0(str_sub(Regen_Plot,1,3),Year))


### Read and summarize seedling data
seedl <- read.csv("data_intermediate_processing_local/regen_summarized_sp_full.csv")

## filter to only kept plots
seedl = seedl %>%
  filter((Regen_Plot %in% plot$Regen_Plot))


####!!!! Need a long table with seedling age counts by plot, so we can calc the year it germinated
####!!!! Also need to assign new 3-letter fire names for fires not sampled an adjacent years

seedl_genus = seedl %>%
  mutate(sp_coarse = recode(species, "PIPO" = "Pine", "PSME" = "Dougfir", "ABCO" = "FirCedar", "CADE27" = "FirCedar", "PILA" = "Pine", "ABMA" = "FirCedar", "PIJE" = "Pine", "PINUS" = "Pine", "ABIES" = "FirCedar", "PICO" = "Pine", "PSMA" = "Dougfir", "PIMO3" = "Pine")) %>%
  filter(sp_coarse %in% c("Pine", "Dougfir", "FirCedar")) %>%
  mutate(sp_coarse = sp_coarse %>% droplevels()) %>%
  select(Regen_Plot, sp_coarse,
         age_00 = count.0yr,
         age_01 = count.1yr,
         age_02 = count.2yr,
         age_03 = count.3yr,
         age_04 = count.4yr,
         age_05 = count.5yr,
         age_06 = count.6yr,
         age_07 = count.7yr,
         age_08 = count.8yr,
         age_09 = count.9yr,
         age_10 = count.10yr,
         age_11 = count.11yr,
         age_12 = count.12yr,
         age_13 = count.13yr,
         age_14 = count.14yr,
         age_15 = count.15yr) %>%

  ## bring in some plot data
  left_join(plot %>% select(Regen_Plot, SHRUB, Year, Year.of.Fire)) %>%
  
  # Pull in fire name and make unique
  mutate(fire = paste0(str_sub(Regen_Plot,1,3),Year)) %>%
  
  group_by(Regen_Plot, sp_coarse) %>%
  summarize(across(age_00:age_15, sum),
            across(c(fire, SHRUB, Year, Year.of.Fire), modal)) %>%
  ungroup() %>%
  mutate(age_all = rowSums(select(.,age_00:age_15), na.rm=TRUE))
  

### make it long so can convert ages to years
seedl_genus_long = seedl_genus %>%
  pivot_longer(age_00:age_15, names_to = "age", values_to = "seedl_count") %>%
  mutate(age = str_sub(age,5,7) %>% as.numeric) %>%
  # the 2016 and 2017 surveys aged 1 year high (e.g., a first-year was aged 1 instead of 0), so subtract one
  mutate(age = ifelse((age > 0) & (Year > 2015), age - 1, age)) %>%
  mutate(year_recruited = Year - age,
         seedl_count = ifelse(is.na(seedl_count),0,seedl_count)) %>%
  
  # remove plots without shrub (richter)
  filter(!is.na(SHRUB)) %>%
  
  # make a shrub low and high
  mutate(shrub_binary = ifelse(SHRUB > 50, "high", "low"))


### Aggregate across fire, species, shrub level
seedl_genus_firesumm = seedl_genus_long %>%
  group_by(sp_coarse, fire, shrub_binary, year_recruited) %>%
  summarize(seedl_count = sum(seedl_count),
            Year.of.Fire = modal(Year.of.Fire),
            n_plots = n()) %>%
  ungroup() %>%
  #remove years from > 2 prior to fire
  filter(year_recruited > (Year.of.Fire - 3))

### get the max number of plots for a given combo of fire, species, shrub_binary. Only keep the combos that had > 2
maxplots = seedl_genus_firesumm %>%
  group_by(fire, sp_coarse, shrub_binary) %>%
  summarize(max_nplots = max(n_plots))

seedl_genus_firesumm = seedl_genus_firesumm %>%
  left_join(maxplots) %>%
  filter(max_nplots > 2)





### Plot

d_plot = seedl_genus_firesumm %>%
  filter(shrub_binary == "low")

## All species together
p = ggplot(d_plot, aes(x=year_recruited,y=seedl_count, fill=sp_coarse)) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_vline(aes(xintercept = Year.of.Fire+.5), color="red") +
  scale_x_continuous(breaks=1986:2017) +
  scale_fill_viridis_d(end=0.8, begin = 0.2) +
  facet_wrap(vars(fire), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "year recruited",
       y = "seedling count",
       title = "Shrub cover < 50%")

png("seedbed/figures/recruitment_lowshrub.png", width=2000, height = 1500, res = 200)
p
dev.off()

d_plot = seedl_genus_firesumm %>%
  filter(shrub_binary == "high")

## All species together
p = ggplot(d_plot, aes(x=year_recruited,y=seedl_count, fill=sp_coarse)) +
  geom_bar(stat="identity", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_vline(aes(xintercept = Year.of.Fire+.5), color="red") +
  scale_x_continuous(breaks=1986:2017) +
  scale_fill_viridis_d(end=0.8, begin = 0.2) +
  facet_wrap(vars(fire), scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "year recruited",
       y = "seedling count",
       title = "Shrub cover > 50%")

png("seedbed/figures/recruitment_highshrub.png", width=1500, height = 1000, res = 200)
p
dev.off()




# ral, rus, elk have PIAT



