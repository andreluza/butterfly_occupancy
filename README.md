Evaluating multi-season occupancy models with autocorrelation fitted to
heterogeneous datasets
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

André L Luza $^1$, Didier Alard $^{1,2}$, Frederic Barraquand $^3$

$^1$ UMR Biodiversité Gènes et Communautés, University of Bordeaux,
INRAE, Pessac, France  
$^2$ Observatoire FAUNA, University of Bordeaux, Pessac, France  
$^3$ Institute of Mathematics of Bordeaux, University of Bordeaux, CNRS,
Bordeaux INP, Talence, France

# ———————————————–

## Brief description of the files and folders in this repository:

#### Root 

\|- *SuppInfo-A-D.pdf*: Supporting Information A, B, C and D associated
to the article (also published as Online Supportiong Information)  
\|- *SuppInfo-E.pdf*: Supporting Information E: analyses using the
Random Walk Occupancy Model by Outwaithe et al. (2018)  
\|- *occ-model-rwalk-RE-simulations.txt*: Model of Outwaithe et
al. (2018) in BUGS language  
\|- *SuppInfo-F.pdf*: Supporting Information F: analyses using data of
the four remaining species: the meadow brown, the false ringlet, the
marsh fritillary, and the small copper  
\|- *SuppInfo-G.pdf*: Supporting Information G: names/IDs of data
contributors, as appearing in FAUNA data sets  
  

#### \|- *Data*: Spatial data 

\|——– SpatialData: systematic grid (grid cells 1x1 km, French geographic
system) from Observatoire FAUNA. The $1 \times 1$km spatial grid was
downloaded from: <a
href="https://observatoire-fauna.fr/ressources/publications?typePublication%5B%5D=fauna&amp;themesID%5B%5D=6\"
class="uri">https://observatoire-fauna.fr/ressources/publications?typePublication%5B%5D=fauna&amp;themesID%5B%5D=6\</a>

  

#### \|- *Processed_data*: Data obtained after processing raw data sets 

\|——– AltitudeHabitatStats.RData: Altitude data (from EUDEM) and habitat
data (from CORINE land use data)  
\|——– SpeciesData.RData: Occupancy data used in SDMs (from Observatoire
FAUNA data base)  
\|——– Water.RData: CORINE land cover data regarding water vs non-water
sites  
  

#### \|- *R*: R scripts 

<span style="color:blue"> \|— *Simulations*: R scripts used in
simulations. They are numbered sequentially, from 1 to 9. Thus, they
should be run respecting this sequence:  
</span> \|——– 1-study1-sc0-… .R: Scripts to reproduce original work of
Doser & Stoudt 2024 = D&S, using Bernoulli sampling design. “c1”:
generate simulated data; “c2”: apply the sampling design and fit the
multi-season occupancy model; “c3”: interpretation of results and
figures. Two separate folders that will host simulation results (path:
model_output/output_simulations) and figures are created here. Results
and figures of each study-scenario will be stored in one specific
folder.  
\|——– 2-study1-sc1-… .R: Scripts to simulate data with higher spatial
autocorrelation (lower values of $\phi$), under the original Bernoulli
sampling design of D&S. “c1”: generate simulated data; “c2”: apply the
sampling design and fit the multi-season occupancy model; “c3”:
interpretation of results and figures.  
\|——– 3-a-study2-Sim-Settings.R: Our simulation settings with initial
number of secondary sampling occasions $J=10$ (equivalent to sampling
months in our butterfly data). These settings are used in Study 2. The
code also create folders to store simulation results, figures, and the
theme for ggplot.  
\|——– 3-study2-sc0-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0). “c1”: creates the Poisson
sampling design - array $G_{itj}$; “c2”: generate simulated data, apply
the sampling design and fit the multi-season occupancy model; “c3”:
interpretation of results and figures  
\|——– 4-study2-sc1-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy covariate. “c1”: generate simulated data, apply the Poisson
sampling design already generated in the previous scenario, and fit the
multi-season occupancy model; “c2”: interpretation of results and
figures  
\|——– 5-study2-sc2-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate (total overlap of covariates). “c1”:
generate simulated data, apply the Poisson sampling design already
generated in the previous scenario, and fit the multi-season occupancy
model; “c2”: interpretation of results and figures  
\|——– 6-study2-sc3-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). “c1”: generate
simulated data, apply the Poisson sampling design already generated in
the previous scenario, and fit the multi-season occupancy model; “c2”:
interpretation of results and figures  
\|——– 7-study3-sc1-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). Design including
the effect of phenology + observer preferences on occupancy data. “c1”:
creates the Poisson sampling design - array $G_{itj}$ - with an effect
of phenology + observer preferences; “c2”: generate simulated data,
apply the sampling design and fit the multi-season occupancy model;
“c3”: interpretation of results and figures  
\|——– 8-study3-sc2-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). Design including
the effect of phenology + observer preferences + observation spot on
occupancy data. “c1”: creates the Poisson sampling design - array
$G_{itj}$ - with an effect of phenology + observer preferences +
observation spot in mid latitudes; “c2”: generate simulated data, apply
the sampling design and fit the multi-season occupancy model; “c3”:
interpretation of results and figures  
\|——– 9-study2-sc1-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and combination of covariates
as study2-sc1 - latitude $L_i$ as occupancy covariate, $v_{itj}$ as
detection covariate. Here the data was analyzed using the random walk
model by Outhwaite et al. (2018) -
<https://doi.org/10.1016/j.ecolind.2018.05.010>. “c1”: build the model
in BUGS language, generate simulated data, apply the Poisson sampling
design and fit the multi-season occupancy model; “c3”: interpretation of
results and figures  
\|——– 10-study-simulation-comparison-Identifiability: Script used to
interpret the output of simulations. This script will produce figures
shown in the main text and SI  
  

<span style="color:blue"> \|— *Empirical analyses*: R scripts used in
the analyses of empirical data. They are numbered sequentially, from 11
to 20. Thus, they should be run respecting this sequence:  
</span> \|——– 11-study-empirical-plot-occupancy-data.R: script used to
plot butterfly data (Fig. 1) with occupancy already available in the
folder “Processed_data” (raw data contain sensitive data. Access to raw
data must be requested from Observatoire FAUNA,
<https://observatoire-fauna.fr/>). The folder that will store figures
will be created here.  
\|——– 12-study-empirical-buffer-DataAnalysis-… .R: script used to
analyze data at Bordeaux + 10 km buffer scale. One script per species -
*P. icarus* and *L. dispar*.  
\|——– 13-study-empirical-buffer-Predictions-Maps-… .R: script used to
make predictions and map estimated distribution at buffer scale. One
script per species  
\|——– 14-study-empirical-buffer-Tables-PPO-… .R: script used to
calculate summary statistics for regression coefficients, including the
Prior-Posterior overlap (Tables in Supporting Information). One script
per species  
\|——– 15-study-empirical-buffer-DataAnalysis-… .R: script used to
analyze data at Nouvelle Aquitaine scale. One script per species  
\|——– 16-study-empirical-buffer-Predictions-Maps-… .R: script used to
make predictions and map estimated distribution at Nouvelle Aquitaine
scale. One script per species  
\|——– 17-study-empirical-buffer-Tables-PPO–… .R: script used to
calculate summary statistics for regression coefficients, including the
Prior-Posterior overlap (Tables in Supporting Information) at Nouvelle
Aquitaine scale. One script per species  
\|——– 18-study-empirical-buffer-DataAnalysis-remaining-spp-stPGocc.R:
script used to analyze data of the four remaining species - Bordeaux +
10 km buffer scale  
\|——– 19-study-empirical-buffer-Predictions-Maps-remaining-species.R:
script used to make predictions and map distribution of the remaining
species - Bordeaux + 10 km buffer scale  
\|——– 20-study-empirical-buffer-Tables-PPO-remaining-species.R: script
used to calculate summary statistics for regression coefficients,
including the Prior-Posterior overlap (Tables in Supporting Information
F). Results shown only at Bordeaux + 10 km buffer scale  
  

#### \|- *model_output*: results of the analyses of simulated and empirical data sets - these folders and files (RData) will be created while running scripts:

\|————— *empirical*: results of empirical data analysis (model output
and predictions) using the $spOccupancy$ package;  
\|————— *output_simulations*: simulation settings and results:  
\|———————- *“Xp_itj.rda”*: observation covariate $v_{itj}$ (generated
once)  
\|———————- *“sim-settings.rda”*: simulation constants and MCMC
settings  
\|———————- *“sampling_design_Poisson”*: sampling design used in study 2
from scenario 0 to 3  
\|———————- *“sampling_design_Poisson_phenology”*: sampling design used
in study 3 - scenario 1 (only phenology + observer preferences)  
\|———————- *“sampling_design_Poisson_phenology_spot”*: sampling design
used in study 3 - scenario 2 (phenology + observer preferences +
observation spot)  
\|———————- *sims_D&S*: correct data and simulation results as per Doser
& Stoudt (2024) - **study 1-scenario 0**  
\|———————- *smooth_sims_D&S*: correct data and simulation results for
analysis of the scenario of high spatial autocorrelation - low spatial
decay parameter $\phi$ - **study 1-scenario 1**  
\|———————- *scenario_zero*: simulation results for study 2 - scenario
0 - **study 2-scenario 0**  
\|———————- *scenario_one*: simulation results for study 2 - scenario 1 -
**study 2-scenario 1**  
\|———————- *scenario_two*: simulation results for study 2 - scenario 2 -
**study 2-scenario 2**  
\|———————- *scenario_two_sparta*: correct data, true $\psi_{it}$ from
**study 2-scenario 1**, and simulation results for study 2 - scenario 1
using the SPARTA approach (uncorrelated site random effects) (Outwhaite
et al. 2018) - **study 2-scenario 1**  
\|———————- *scenario_three*: simulation results for study 2 - scenario
3 - **study 2-scenario 3**  
\|———————- *scenario_phenology*: simulation results for the study 3 -
scenario 1 incorporating phenology + observer preferences - **study 3 -
scenario 1**  
\|———————- *scenario_phenology_spot*: correct data and simulation
results for study 3 - scenario 2 incorporating phenology + observer
preferences + observation spot - **study 3 - scenario 2**  
  

##### This article was produced using the following software and associated packages:

  

    ## R version 4.5.1 (2025-06-13)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 22.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Paris
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] MCMCvis_0.16.3          MASS_7.3-65             spBayes_0.4-8          
    ##  [4] spOccupancy_0.8.0       jagsUI_1.6.2            coda_0.19-4            
    ##  [7] tidyr_1.3.1             raster_3.6-32           sp_2.2-0               
    ## [10] terra_1.8-54            abind_1.4-5             spdep_1.3-11           
    ## [13] spData_2.3.4            kableExtra_1.4.0        knitr_1.50             
    ## [16] ggh4x_0.3.1             units_0.8-7             gridExtra_2.3          
    ## [19] rnaturalearthdata_1.0.0 rnaturalearth_1.0.1     reshape_0.8.9          
    ## [22] lubridate_1.9.4         ggbreak_0.1.4           dplyr_1.1.4            
    ## [25] ggplot2_3.5.2           sf_1.0-21               here_1.0.1             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rdpack_2.6.4       DBI_1.2.3          deldir_1.0-6       s2_1.0.7          
    ##  [5] rlang_1.1.6        magrittr_2.0.3     e1071_1.7-16       compiler_4.5.1    
    ##  [9] systemfonts_1.2.3  vctrs_0.6.5        stringr_1.5.1      pkgconfig_2.0.3   
    ## [13] wk_0.9.4           fastmap_1.2.0      magic_1.6-0        rmarkdown_2.29    
    ## [17] nloptr_2.0.0       purrr_1.0.4        xfun_0.52          aplot_0.2.5       
    ## [21] jsonlite_2.0.0     parallel_4.5.1     R6_2.6.1           stringi_1.8.7     
    ## [25] RColorBrewer_1.1-3 boot_1.3-32        Rcpp_1.0.14        iterators_1.0.14  
    ## [29] Matrix_1.7-4       splines_4.5.1      timechange_0.3.0   tidyselect_1.2.1  
    ## [33] rstudioapi_0.17.1  dichromat_2.0-0.1  yaml_2.3.10        doParallel_1.0.17 
    ## [37] codetools_0.2-19   lattice_0.22-5     tibble_3.3.0       plyr_1.8.9        
    ## [41] withr_3.0.2        evaluate_1.0.4     gridGraphics_0.5-1 proxy_0.4-27      
    ## [45] xml2_1.3.8         pillar_1.10.2      spAbundance_0.2.1  KernSmooth_2.23-26
    ## [49] foreach_1.5.2      reformulas_0.4.1   ggfun_0.1.8        generics_0.1.4    
    ## [53] rprojroot_2.0.2    scales_1.4.0       minqa_1.2.4        class_7.3-23      
    ## [57] glue_1.8.0         tools_4.5.1        lme4_1.1-37        RANN_2.6.1        
    ## [61] fs_1.6.6           rbibutils_2.3      nlme_3.1-168       patchwork_1.3.0   
    ## [65] Formula_1.2-4      cli_3.6.5          textshaping_0.3.6  viridisLite_0.4.2 
    ## [69] svglite_2.2.1      gtable_0.3.6       yulab.utils_0.2.0  digest_0.6.37     
    ## [73] classInt_0.4-11    ggplotify_0.1.2    farver_2.1.2       htmltools_0.5.8.1 
    ## [77] lifecycle_1.0.4    httr_1.4.7

  

#### End of the description
