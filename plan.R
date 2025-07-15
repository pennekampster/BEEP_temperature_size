#require(remotes)
#install_version("drake", version = "7.13.10", repos = "http://cran.us.r-project.org", dependencies=T, upgrade=T)

library(drake)
library(tidyverse)
library(lubridate)
library(data.table)

# model fitting
library(lme4)
library(lmerTest)
library(lmodel2)
library(broom.mixed)
library(optimx)
library(brms)
library(rstan)
library(ggfortify)

library(flextable)
library(modelsummary)
# fit reduced major axis regression
library(smatr)

# fit logistic growth
library(nls.multstart)  
library(nlstools)

# calculate AUC
library(zoo)

# output model tables
library(gtable)

# plotting
library(patchwork)
library(egg)
library("dotwhisker")
library(viridis)

# file operations
library(fs)

# libraries used in functions
library(AICcmodavg)
library(gt)
library(cowplot)



# source all functions
rm(list=ls())
source(here::here("R/funs.R"))

# build drake plan
plan <- drake_plan(

  # load & process the raw data
  raw_data = load_data(),
  processed_data = process_data(raw_data),
  
  # set up folder structure for outputs
  folder_ok = create_output_folder(),
  
  # aggregate size dataset for analysis
  mean_size_data = create_dataset_analysis_mean_size(processed_data),
   
  # fit logistic growth models as demand proxy and for time series alignment
  logistic_growth_parameters_for_aligning = 
  fit_logistic_growth(processed_data, response_var="density"),
   
  logistic_growth_parameters_as_demand_proxy = 
  fit_logistic_growth(processed_data, response_var="density"),
   
  aligned_timeseries_data = 
   align_mono_and_polyculture_time_series(
       processed_data, 
       logistic_growth_parameters_for_aligning, 
     response_var="biomass"),
   
  # calculate supply proxy and SD ratio
   supply_proxy_mono = calculate_supply_proxy_mono(aligned_timeseries_data, processed_data),
   mono_SD_data = create_mono_SD_data(logistic_growth_parameters_as_demand_proxy, supply_proxy_mono),
   
   # aggregate SD dataset 
   poly_SD_data = create_poly_SD_data(logistic_growth_parameters_as_demand_proxy, processed_data),
   
   
   
  # linear models on monocultures (associated with figure 1)
   
  #linear model analysis on size
  linear_models_size = linear_model_analysis_size(mean_size_data), 
  linear_models_supply = linear_model_analysis_supply(poly_SD_data, "AUC"), 
  linear_models_demand = linear_model_analysis_demand(logistic_growth_parameters_as_demand_proxy), 
  linear_model_SD = linear_model_analysis_SD(mono_SD_data),
  
  # Mixed model analyses across richness and temperature (associated with figure 3)
   
  # two-way mixed model analysis on size
  two_way_mixed_models_size = twoway_mixed_model_analysis_size(mean_size_data),
  twoway_parameters_size = visualize_twoway_mixed_model_parameters_size(two_way_mixed_models_size[[1]]),

  # investigate SD proxy across polycultures using mixed model
  mixed_model_SD = model_plot_SD_poly(poly_SD_data, poly_output_figure),

  two_way_mixed_models_SD_ratio = twoway_mixed_model_analysis_SD_ratio(poly_SD_data),
  twoway_parameters_SD_ratio = visualize_twoway_mixed_model_parameters_SD_ratio(two_way_mixed_models_SD_ratio[[1]]),


  # figures
  linear_composite_figure_SD = build_composite_figure_SD_model_mono(mono_SD_data),
  poly_output_figure = create_poly_output_figure(two_way_mixed_models_size[[2]]),

)

# visualise dependencies
#vis_drake_graph(plan)

# Re-run analysis from scratch
#clean()

# run drake plan
source(here::here("R/funs.R"))

# set lock environment to FALSE to run brms (not recommended by drake developer)
make(plan)

# look at intermediate result
# loadd(folder_ok)
# loadd(processed_data)


