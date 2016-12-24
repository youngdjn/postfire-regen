# postfire-regen  

## Scripts  
- regen_compile_data.R: take survey data and geospatial data (e.g. climate, elevation) to extract relevant data summarized by plot and by plotXspecies (the intermediate products--see below)  
  #### To-do:  
  - Update regen summarization to (a) specify old-young cutoff age on independent line, (b) have flexibility to deal with fires that are 4 years old.  
  - Seedlings counts should be counts per year  
  - Set change FORBE to FORB here rather than in regen_compile_data  
  
- regen_compile_data_functions.R: functions used by the script above  

- regen_analysis.R: take intermediate data products, compute appropriate summaries (via aggregation) and perform analyses  
  #### To-do:  
  - Improve topoclimate categorization to (a) allow more topoclimatic categories for fires that have more plots, to achieve a target of ~10(?) plots per category per fire, (b) make sure there are enough plots in all factorial combinations of categorization variables (e.g., what if there are no high-radiation plots in high-precipitation areas?), and (c) make sure there are enough control plots as well as high-severity plots  
  - Minimum of 5 highsev and 5 control plots in a topoclimatic category in order to include it in analysis

  
- regen_analysis_functions.R: functions used by the script above  


## Initial inputs (used by regen_compile_data.R)  
- Access database tables (exported as csv)  
- Climatic variables   
- DEM  
- Solar radiation (March)  


## Intermediate products (produced by regen_compile_data.R)  
1. Plot-level dataframe  
	- includes columns for: climate variables, shrub etc % cover, fire-level data repeated by plot, topographic variables, distance to seed tree,   
2. Species data by plot dataframe  
	- includes columns for: plot, species, number and basal area of adults (defined in a previous input script as greater than some diameter), total number of seedlings, number of “young” seedlings (1-3 years old, defined in previous script), number of “old” seedlings (4-5 years old, defined in previous script),   
  
  
## Functions  
FUNCTION for aggregation to topo-category of species data by plot-level topo data  
Input:   
Output: topo-aggregated dataframe with columns for fire, topo category, species, regen count, regen % of plots present, control (adult) count, control (adult) basal area, control (adult) % of plots present  
  - This is implemented in regen_analysis.R (steps #1 and 2; not a function--not sure that it has to be, but open to it)  


FUNCTION for melting aggregated data into a long form for plotting  
Returns long topo version for ggplot graphics with regen.(whatever variable) and control.(whatever variable) as separate rows in the variable column with each having their own values, using the melt() function. Keeps environmental predictors as columns, drops other regen and control data columns. Set up to do just one variable at a time for now.   
Input: topo-aggregated dataframe, variables to include  
Output: long topo-aggregated dataframe with just one metric of regen/control presence  
  - Not yet implemented  

## Directories  
- root: principal scripts for data compilation and analysis  
- data_intermediate: intermediate data files (see header above)
- code_old_reference: old scripts from previous iterations of analysis that could be useful for reference


## Final products  
1. Bayesian model by species/species group  
Inputs: topo-aggregated dataframe  
Outputs: parameter estimates, counterfactual plots  

2. Figure comparing regen to control  
Inputs: long format topo-aggregated dataframe  

3. Postfire climate on the relative abundance of old vs young seedlings  
Inputs: topo-aggregated dataframe  

