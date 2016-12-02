# postfire-regen

## To-do
Update regen summarization to (a) specify old-young cutoff age on independent line, (b) have flexibility to deal with fires that are 4 years old.


## Scripts
1. regen_compile_data.R: take survey data and geospatial data (e.g. climate, elevation) to extract relevant data summarized by plot and by plotXspecies (the intermediate products--see below)  
2. regen_compile_data_functions.R: functions used by the script above  


## Initial inputs  
Access database tables (exported)  
Climatic variables   
DEM  
Solar radiation (March)  


## Intermediate products  
1. Plot-level dataframe  
	- includes columns for: climate variables, shrub etc % cover, fire-level data repeated by plot, topographic variables, distance to seed tree,   
2. Species data by plot dataframe  
	- includes columns for: plot, species, number and basal area of adults (defined in a previous input script as greater than some diameter), total number of seedlings, number of “young” seedlings (1-3 years old, defined in previous script), number of “old” seedlings (4-5 years old, defined in previous script),   
  
  
## Functions  
FUNCTION for aggregation to topo-category of species data by plot-level topo data  
Input:   
Output: topo-aggregated dataframe with columns for fire, topo category, species, regen count, regen % of plots present, control (adult) count, control (adult) basal area, control (adult) % of plots present  


FUNCTION for melting aggregated data into a long form for plotting  
Returns long topo version for ggplot graphics with regen.(whatever variable) and control.(whatever variable) as separate rows in the variable column with each having their own values, using the melt() function. Keeps environmental predictors as columns, drops other regen and control data columns. Set up to do just one variable at a time for now.   
Input: topo-aggregated dataframe, variables to include  
Output: long topo-aggregated dataframe with just one metric of regen/control presence  


## Final products  
1. Bayesian model by species/species group  
Inputs: topo-aggregated dataframe  
Outputs: parameter estimates, counterfactual plots  

2. Figure comparing regen to control  
Inputs: long format topo-aggregated dataframe  

3. Postfire climate on the relative abundance of old vs young seedlings  
Inputs: topo-aggregated dataframe  

