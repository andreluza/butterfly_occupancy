### Evaluating multi-season occupancy models with autocorrelation fitted to heterogeneous datasets 
#### Andre Luza$^1$, Didier Alard$^1$, Frederic Barraquand$^2$
$^1$ UMR Biodiversité Gènes et Communautés, University of Bordeaux,
INRAE, Pessac, France\
$^2$ Institute of Mathematics of Bordeaux, University of Bordeaux,
CNRS, Bordeaux INP, Talence, France

# -----------------------------------------------

Brief description of the files and folders in this repository:

#### Root \
|- *SuppInfo-A-D.pdf*: Supporting Information A, B, C and D associated to the article\
|- *SuppInfo-E.pdf*: Supporting Information E associated to the article (analyses using the Random Walk Occupancy Model by Outwaithe et al. 2018)\
|- *SuppInfo-F.pdf*: Supporting Information F associated to the article (analyses using data of the four remaining species: the meadow brown, the false ringlet, the marsh fritillary, and the small copper)\
|- *SuppInfo-G.pdf*: Supporting Information G associated to the article (names/IDs of data contributors as appearing the FAUNA data sets)\

#### |- *Data*: Spatial data\
|-------- SpatialData: systematic grid (grid cells 1x1 km, French geographic system) from FAUNA\

#### |- *Processed_data*: Data obtained after processing raw data sets \
|-------- AltitudeHabitatStats.RData: Altitude data (From EUDEM) and habitat data (from CORINE land use data)\
|-------- SpeciesData: Occupancy data used in SDMs (from FAUNA data base)\
|-------- Water.RData: CORINE land cover data regarding water vs non-water sites\

#### |- *R*: R scripts (numbered from 1 to 20 to track different steps of data exploration and analyses)\
|--- *Simulations*\
|-------- R scripts are numbered sequentially, from 1 to 9. Thus they should be run sequentially:\
|-------- 1-study1-sc0-... .R: Scripts to reproduce original simulations (as per Doser & Stoudt 2024 = D&S) using Bernoulli sampling design\
|-------- 2-study1-sc1-... .R: Scripts to reproduce D&S simulations with higher spatial autocorrelation (lower values of \phi) and Bernoulli sampling design;\
|-------- 3-study2-sc0-... .R: Scripts to simulations with our Poisson sampling design ($\phi$ as in study 1-sc0)\
|-------- 4-study2-sc1-... .R: Scripts to simulations with our Poisson sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as occupancy covariate\
|-------- 5-study2-sc2-... .R: Scripts to simulations with our Poisson sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as occupancy and detection covariate (total overlap of covariates)\
|-------- 6-study2-sc3-... .R: Scripts to simulations with our Poisson sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as occupancy and detection covariate, and $v_{itj}$ (random covariate) in the detection model (partial overlap of covariates)\
|-------- 7-study3-sc1-... .R: Scripts to simulations with our Poisson sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as occupancy and detection covariate, and $v_{itj}$ (random covariate) in the detection model (partial overlap of covariates). Design including the effect of phenology + observer preferences on occupancy data\
|-------- 8-study3-sc2-... .R: Scripts to simulations with our Poisson sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as occupancy and detection covariate, and $v_{itj}$ (random covariate) in the detection model (partial overlap of covariates). Design including the effect of phenology + observer preferences + observation spot on occupancy data\
|-------- 9-study2-sc1-... .R: Scripts to simulations with our Poisson sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as occupancy covariate. Here the data was analyzed using the random walk model by Outhwaite et al. (2018) - https://doi.org/10.1016/j.ecolind.2018.05.010 \
|-------- 10-study-simulation-comparison-Identifiability: Script used to interpret simulations' output. Produce figures in the main text and SI\
|--- *Empirical analyses*\
|-------- 11-study-... .R: scripts used to organize butterfly data, extract covariates, and create the occupancy dataset available in the folder "Processed_data"\
|-------- 11-study-empirical-plot-occupancy-data.R: script used to plot butterfly data (Fig. 1)\
|-------- 12-study-empirical-buffer-DataAnalysis-... .R: script used to analyze data at Bordeaux + 10 km buffer scale. One script per species\
|-------- 13-study-empirical-buffer-Predictions-Maps-... .R: script used to make predictions and map. One script per species\
|-------- 14-study-empirical-buffer-Tables-PPO-... .R: script used to calculate Prior-Posterior overlap (Tables in Supporting Information). One script per species\
|-------- 15-study-empirical-buffer-DataAnalysis-... .R: script used to analyze data at Nouvelle Aquitaine scale. One script per species\
|-------- 16-study-empirical-buffer-Predictions-Maps-... .R: script used to make predictions and map at Nouvelle Aquitaine scale. One script per species\
|-------- 17-study-empirical-buffer-Tables-PPO--... .R: script used to calculate Prior-Posterior overlap (Tables in Supporting Information) at Nouvelle Aquitaine scale. One script per species\
|-------- 18-study-empirical-buffer-DataAnalysis-remaining-spp-stPGocc.R: script used to analyze data of the four remaining species - Bordeaux + 10 km buffer scale\
|-------- 19-study-empirical-buffer-Predictions-Maps-remaining-species.R: script used to make predictions and map distribution of the remaining species - Bordeaux + 10 km buffer scale\
|-------- 20-study-empirical-buffer-Tables-PPO-remaining-species.R: script used to calculate Prior-Posterior overlap (Tables in Supporting Information F) - Bordeaux + 10 km buffer scale\

#### |- *model_output*: results of the analyses of simulated and empirical data sets - these folders will be created while running scripts:
|--------------- *empirical*: results of empirical data analysis using spOccupancy; \
|--------------- *output_simulations*: rerun D&S simulations based on butterfly data: \
|---------------------- *"Xp_itj.rda"*: observation covariate $v_{itj}$ (generated once) \
|---------------------- *"sim-settings.rda"*: simulation constants and MCMC settings \
|---------------------- *"sampling_design_Poisson"*: sampling design used in study 2, scenarios 0 to 3 \
|---------------------- *"sampling_design_Poisson_phenology"*: sampling design used in study 3, scenario 1 (only phenology) \
|---------------------- *"sampling_design_Poisson_phenology_spot"*: sampling design used in study 3, scenario 2 (phenology and observation spot) \
|---------------------- *sims_D&S*: correct data and simulation results for D&S original sims - **study 1-scenario 0** \
|---------------------- *smooth_sims_D&S*: correct data and simulation results for D&S-like simulations with lower values of the spatial decay parameter $\phi$ - **study 1-scenario 1** \
|---------------------- *scenario_zero*: correct data and simulation results for scenario 0 - **study 2-scenario 0** \
|---------------------- *scenario_one*: correct data and simulation results for scenario 1 - **study 2-scenario 1** \
|---------------------- *scenario_two*: correct data and simulation results for scenario 2 - **study 2-scenario 2** \
|---------------------- *scenario_two_sparta*: correct data and simulation results for scenario 2 using the SPARTA approach (uncorrelated site random effects) (Outwhaite et al. 2018) - **study 2-scenario 2** \
|---------------------- *scenario_three*: correct data and simulation results for scenario 3 - **study 2-scenario 3** \
|---------------------- *scenario_phenology*: correct data and simulation results for the scenario incorporating phenology + observer preferences - **study 3 - scenario 1** \
|---------------------- *scenario_phenology_spot*: correct data and simulation results for the scenario incorporating phenology + observer preferences + observation spot - **study 3 - scenario 2** \
                      
#### End of the description                      




