# postfire-regen  

## Scripts  
- regen_compile_data.R: take survey data and geospatial data (e.g. climate, elevation) to extract relevant data summarized by plot and by plotXspecies (the intermediate products--see below)  

- regen_compile_data_functions.R: functions used by the script above  

- regen_analysis.R: take intermediate data products, compute appropriate summaries and perform analyses  
  
- regen_analysis_functions.R: functions used by the script above  


## Initial inputs (used by regen_compile_data.R)  
- Access database tables (exported as csv)  
- Climatic variables   
- DEM  
- Solar radiation (March)  


## Intermediate data products (produced by regen_compile_data.R)  
1. Plot-level dataframe  
	- includes columns for: climate variables, shrub etc % cover, fire-level data repeated by plot, topographic variables, distance to seed tree,   
2. Species data by plot dataframe  
	- includes columns for: plot, species, number and basal area of adults (defined in a previous input script as greater than some diameter), total number of seedlings, number of “young” seedlings (0-2 years old, defined in previous script), number of “old” seedlings (germinated in first two years of fire, defined in previous script; accounts for fires of different ages), and number of “all” seedlings. All seedling count values are seedlings per year (over the years that apply--e.g., "Young seedlings" values are average number of seedlings that germinated per year in the last 3 years prior to plot survey)  
  
  
## Directories  
- root: principal scripts for data compilation and analysis  
- data_intermediate: intermediate data files (see header above)  
- code_old_reference: old scripts from previous iterations of analysis that could be useful for reference  
