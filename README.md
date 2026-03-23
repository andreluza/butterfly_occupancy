Evaluating multi-season occupancy models with autocorrelation fitted to
heterogeneous datasets
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

André L Luza $^1$, Didier Alard $^{1,2}$, Frédéric Barraquand $^3$

$^1$ UMR Biodiversité Gènes et Communautés, University of Bordeaux,
INRAE, Pessac, France  
$^2$ US Fauna, University of Bordeaux, Pessac, France  
$^3$ Institute of Mathematics of Bordeaux, University of Bordeaux, CNRS,
Bordeaux INP, Talence, France  

Correspondence to:  
André L Luza: <andre-luis.luza@u-bordeaux.fr>  
Frédéric Barraquand: <frederic.barraquand@u-bordeaux.fr>  

# ———————————————–

## Brief description of the files and folders in this repository:

#### Root 

\|- *SuppInfo-A-D.pdf*: Supporting Information A, B, C and D associated
to the article (also published as Online Supporting Information)  
\|- *SuppInfo-E.pdf*: Supporting Information E: Assessment of model
misspeficication, analysis of $\psi_{it} \times p_it$, and evaluation of
the study scenario 3-2 with NNGP=15 spatial neighbors.  
\|- *SuppInfo-F.pdf*: Supporting Information E: analyses using the
Random Walk Occupancy Model by Outwaithe et al. (2018)  
\|- *occ-model-rwalk-RE-simulations.txt*: Model of Outwaithe et
al. (2018) in BUGS language  
\|- *SuppInfo-G.pdf*: Supporting Information F: analyses using data of
the four remaining species: the meadow brown, the false ringlet, the
marsh fritillary, and the small copper  
\|- *SuppInfo-H.pdf*: Supporting Information G: names/IDs of data
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

<span style="color:blue"> \|— *SIMULATIONS*: R scripts used in
simulations. They are numbered sequentially, from 1 to 9. Thus, they
should be run respecting this sequence:  
</span> \|——– SIM-1-study1-sc0-… .R: Scripts to reproduce original work
of Doser & Stoudt 2024 (“D&S”) using Bernoulli sampling design. “c1”:
generate simulated data; “c2”: apply the sampling design and fit the
multi-season occupancy model; “c3”: interpretation of results and
figures. Two separate folders that will host simulation results (path:
model_output/output_simulations) and figures are created here. Results
and figures of this study-scenario are stored in one specific folder in
the path: “model_output –\> output_simulations –\> sims_D&S”.  
\|——– SIM-2-study1-sc1-… .R: Scripts to simulate data with higher
spatial autocorrelation (lower values of $\phi$), under the original
Bernoulli sampling design of D&S. “c1”: generate simulated data; “c2”:
apply the sampling design and fit the multi-season occupancy model;
“c3”: interpretation of results and figures. Results and figures of this
study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\> smooths_sims_D&S”.  
\|——– SIM-3-a-study2-Sim-Settings.R: Our simulation settings with
initial number of secondary sampling occasions $J=10$ (equivalent to
sampling months in our butterfly data). These settings are used in Study
2. The code also create folders to store simulation results, figures,
and the theme for ggplot (as per D&S).  
\|——– SIM-3-study2-sc0-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0). “c1”: creates the Poisson
sampling design - array $G_{itj}$; “c2”: generate simulated data, apply
the sampling design and fit the multi-season occupancy model; “c3”:
interpretation of results and figures. Results and figures of this
study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\> scenario_zero”.  
\|——– SIM-4-study2-sc1-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy covariate. “c1”: generate simulated data, apply the Poisson
sampling design already generated in the previous scenario, and fit the
multi-season occupancy model; “c2”: interpretation of results and
figures. Results and figures of this study-scenario are stored in one
specific folder in the path: “model_output –\> output_simulations –\>
scenario_one”.  
\|——– SIM-5-study2-sc2-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate (total overlap of covariates). “c1”:
generate simulated data, apply the Poisson sampling design already
generated in the previous scenario, and fit the multi-season occupancy
model; “c2”: interpretation of results and figures. Results and figures
of this study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\> scenario_two”.  
\|——– SIM-6-study2-sc3-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). “c1”: generate
simulated data, apply the Poisson sampling design already generated in
the previous scenario, and fit the multi-season occupancy model; “c2”:
interpretation of results and figures. Results and figures of this
study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\> scenario_three”.  
\|——– SIM-7-study3-sc1-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). Design including
the effect of phenology + observer preferences on occupancy data. “c1”:
creates the Poisson sampling design - array $G_{itj}$ - with an effect
of phenology + observer preferences; “c2”: generate simulated data,
apply the sampling design and fit the multi-season occupancy model;
“c3”: interpretation of results and figures. Results and figures of this
study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\> scenario_phenology”.  
\|——– SIM-8-study3-sc2-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). Design including
the effect of phenology + observer preferences + observation spot on
occupancy data. “c1”: creates the Poisson sampling design - array
$G_{itj}$ - with an effect of phenology + observer preferences +
observation spot in mid latitudes; “c2”: generate simulated data, apply
the sampling design and fit the multi-season occupancy model; “c3”:
interpretation of results and figures. Results and figures of this
study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\> scenario_phenology_spot”.  
\|——– SIM-9-study3-sc2-…NNGP-15.R: Scripts used to analyze simulated
data, as per study scenario 3-2, using NNGP=15 spatial neighbors in the
Gaussian Process approximation. Results and figures of this
study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\>
scenario_phenology_spot-NNGP-15-review”.  
\|——– SIM-10-study3-sc3-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). Design including
the effect of phenology + observer preferences + observation spot of
I\*50% sites on occupancy data. “c1”: creates the Poisson sampling
design - array $G_{itj}$ - with an effect of phenology + observer
preferences + observation spot in mid latitudes; “c2”: generate
simulated data, apply the sampling design and fit the multi-season
occupancy model; “c3”: interpretation of results and figures. Results
and figures of this study-scenario are stored in one specific folder in
the path: “model_output –\> output_simulations –\>
scenario_phenology_spot_review”.  
\|——– SIM-11-study2-sc1-… .R: Scripts to analyze data from study
scenario 2-1 with the random walk model by Outhwaite et al. (2018) -
<https://doi.org/10.1016/j.ecolind.2018.05.010>. “c1”: build the model
in BUGS language, generate simulated data, apply the Poisson sampling
design and fit the multi-season occupancy model; “c2”: interpretation of
results and figures. Results and figures of this study-scenario are
stored in one specific folder in the path: “model_output –\>
output_simulations –\> scenario_two_sparta”.  
\|——– SIM-12-study-simulation-comparison-Identifiability: Script used to
interpret the output of simulations. This script will produce figures
shown in the main text and SI  
\|——– SIM-12-simulation-psi-p-product: Evaluation of the product of
$\psi_{it} \times p_{it}$ for each scenario.  

  

<span style="color:blue"> \|— *EMPIRICAL DATA ANALYSES*: R scripts used
in the analyses of empirical data. They have the prefix “EMP” and are
numbered sequentially, from 4 to 11. Code 1 to 3 involve the handling of
sensitive information about threatened species and are not shown in the
repository; code 4 from to 11 handle processed data used in our
analyzes. Code should be run respecting this sequence:  
</span> \|——– EMP-4-plot-occupancy-data.R: script used to plot butterfly
data (Fig. 1) with occupancy already available in the folder
“Processed_data” (raw data contain sensitive data. Access to raw data
must be requested from Observatoire FAUNA,
<https://observatoire-fauna.fr/>). The folder that will store figures
will be created here.  
\|——– EMP-5-model-fit-buffer-L-dispar-stPGOcc.R: script used to fit
models to large copper *Lycaena dispar* data at Bordeaux + 10 km buffer
scale.  
\|——– EMP-5-model-fit-buffer-L-dispar-stPGOcc.R: script used to fit
models to common blue *Polyommatus icarus* data at Bordeaux + 10 km
buffer scale.  
\|——– EMP-6-buffer-Predictions-Maps-L-dispar.R: script used to make
predictions and map the estimated distribution of the large copper at
buffer scale.  
\|——– EMP-6-buffer-Predictions-Maps-P-icarus.R: script used to make
predictions and map the estimated distribution of the common blue at
buffer scale.  
\|——– EMP-7-buffer-Tables-PPO-L-dispar.R: script used to calculate
summary statistics of regression coefficients, including the
Prior-Posterior overlap (Tables in Supporting Information), from the
model fitted to large copper data at buffer scale.  
\|——– EMP-7-buffer-Tables-PPO-P-icarus.R: script used to calculate
summary statistics of regression coefficients, including the
Prior-Posterior overlap (Tables in Supporting Information), from the
model fitted to common blue data at buffer scale.  
\|——– EMP-8-model-fit-NouvAquit-L-dispar-stPGOcc.R: script used to fit
models to large copper *Lycaena dispar* data at Nouvelle Aquitaine
scale.  
\|——– EMP-8-model-fit-NouvAquit-P-icarus-stPGOcc.R: script used to fit
models to common blue *Polyommatus icarus* data at Nouvelle Aquitaine
scale.  
\|——– EMP-9-NouvAquit-Predictions-Maps-L-dispar.R: script used to make
predictions and map the estimated distribution of the large copper at
the Nouvelle Aquitaine scale.  
\|——– EMP-9-NouvAquit-Predictions-Maps-P-icarus.R: script used to make
predictions and map the estimated distribution of the common blue at the
Nouvelle Aquitaine scale.  
\|——– EMP-10-NouvAquit-Tables-PPO-L-dispar.R: script used to calculate
summary statistics for regression coefficients, including the
Prior-Posterior overlap (Tables in Supporting Information), from the
model fitted to large copper data at the Nouvelle Aquitaine scale.  
\|——– EMP-10-NouvAquit-Tables-PPO-P-icarus.R: script used to calculate
summary statistics for regression coefficients, including the
Prior-Posterior overlap (Tables in Supporting Information), from the
model fitted to common blue data at the Nouvelle Aquitaine scale.  
\|——– EMP-11-buffer-model-fit-remaining-species.R: script used to fit
models to data of *Coenonympha oedippus*, *Euphydryas aurinia*, *Lycaena
phlaeas*, and *Maniola jurtina* at Bordeaux + 10 km buffer scale.  
\|——– EMP-11-buffer-Predictions-Maps-remaining-species.R: script used to
make predictions and map distribution of *Coenonympha oedippus*,
*Euphydryas aurinia*, *Lycaena phlaeas*, and *Maniola jurtina* at
Bordeaux + 10 km buffer scale  
\|——– EMP-11-buffer-Tables-PPO-remaining-species.R: script used to
calculate summary statistics for regression coefficients, including the
Prior-Posterior overlap (Tables in Supporting Information F), from
models fitted to *Coenonympha oedippus*, *Euphydryas aurinia*, *Lycaena
phlaeas*, and *Maniola jurtina* at Bordeaux + 10 km buffer scale.  
  

#### \|- *model_output*: General simulation settings, sampling designs, and results of the analyses of simulated and empirical data sets. These folders and files (RData) will be created while running scripts:

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
in study 3 - scenario 1 (only phenology + observer preferences).  
\|———————- *“sampling_design_Poisson_phenology_spot”*: sampling design
used in study 3 - scenario 2 (phenology + observer preferences +
observation spot I*25% of the sites).  
\|———————- *”sampling_design_Poisson_phenology_spot_st3_sc3”*: sampling
design used in study 3 - scenario 2 (phenology + observer preferences +
observation spot I*50% of the sites).  
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
\|———————- *scenario_phenology_spot-NNGP-15-review*: analysis of
simulated data for study 3 - scenario 2 with NNGP=15 spatial
neighbors.  
\|———————- *scenario_phenology_spot_review*: correct data and simulation
results for study 3 - scenario 3 incorporating phenology + observer
preferences + observation spot of 50% of the sites - **study 3 -
scenario 3**  

  

#### \|- *R-mis-specification*: Code used in the assessment of model misspecification (three scenarios analyzed):

\|——– SIM-1-study1-sc0-… .R: Scripts to reproduce original work of Doser
& Stoudt 2024 (“D&S”) using Bernoulli sampling design. “c1”: generate
simulated data with mis.spec.type = “probit”; “c2”: apply the sampling
design and fit the multi-season occupancy model (with logit link
function); “c3”: interpretation of results and figures. Two separate
folders that will host simulation results (path:
model_output/output_simulations) and figures are created here. Results
and figures of this study-scenario are stored in one specific folder in
the path: “model_output –\> output_simulations –\> sims_D&S”.  
\|——– SIM-3-a-study2-Sim-Settings.R: Our simulation settings with
initial number of secondary sampling occasions $J=10$ (equivalent to
sampling months in our butterfly data). These settings are used in Study
2. The code also create folders to store simulation results, figures,
and the theme for ggplot (as per D&S).  
\|——– SIM-5-study2-sc2-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate (total overlap of covariates). “c1”:
generate simulated data with pnorm link for $\psi_{it}$ and $p_{itj}$,
apply the Poisson sampling design already generated in the previous
scenario, and fit the multi-season occupancy model (with logit link
function); “c2”: interpretation of results and figures. Results and
figures of this study-scenario are stored in one specific folder in the
path: “model_output –\> output_simulations –\> scenario_two”.  
\|——– SIM-8-study3-sc2-… .R: Scripts to simulations with our Poisson
sampling design ($\phi$ as in study 1-sc0) and latitude $L_i$ as
occupancy and detection covariate, and $v_{itj}$ (random covariate) in
the detection model (partial overlap of covariates). Design including
the effect of phenology + observer preferences + observation spot on
occupancy data. “c1”: creates the Poisson sampling design - array
$G_{itj}$ - with an effect of phenology + observer preferences +
observation spot in mid latitudes; “c2”: generate simulated data with
pnorm link for $\psi_{it}$ and $p_{itj}$, apply the sampling design and
fit the multi-season occupancy model (with logit link function); “c3”:
interpretation of results and figures. Results and figures of this
study-scenario are stored in one specific folder in the path:
“model_output –\> output_simulations –\> scenario_phenology_spot”.  

#### \|- *model_output_mis-specification*: Results of the assessment of model misspecification (three scenarios analyzed):

\|————— *figures*: Scatter plots showing the relationship between true
and estimated $\psi_{it}$ under a misspecified model;  
\|————— *output_simulations*: simulation settings and misspecification
results for the following studies - scenarios:  
\|———————- *“sim-settings.rda”*: simulation constants and MCMC
settings  
\|———————- *sims_D&S*: correct data and simulation results as per Doser
& Stoudt (2024) - **study 1-scenario 0**  
\|———————- *scenario_two*: simulation results for study 2 - scenario 2 -
**study 2-scenario 2**  
\|———————- *scenario_phenology_spot*: correct data and simulation
results for study 3 - scenario 2 incorporating phenology + observer
preferences + observation spot of I\*25% sites.  

##### This article was produced using the following software and associated packages:

  

    ## R version 4.3.3 (2024-02-29)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 24.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Paris
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] MCMCvis_0.16.5          MASS_7.3-60.0.1         spBayes_0.4-8          
    ##  [4] spOccupancy_0.8.0       coda_0.19-4.1           tidyr_1.3.2            
    ##  [7] raster_3.6-32           sp_2.2-1                terra_1.9-1            
    ## [10] abind_1.4-8             spdep_1.4-2             spData_2.3.4           
    ## [13] kableExtra_1.4.0        knitr_1.51              ggh4x_0.3.1            
    ## [16] units_1.0-1             gridExtra_2.3           rnaturalearthdata_1.0.0
    ## [19] rnaturalearth_1.2.0     reshape_0.8.10          lubridate_1.9.5        
    ## [22] ggbreak_0.1.6           dplyr_1.2.0             ggplot2_4.0.2          
    ## [25] sf_1.1-0                here_1.0.2             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rdpack_2.6.6       DBI_1.3.0          deldir_2.0-4       s2_1.1.9          
    ##  [5] rlang_1.1.7        magrittr_2.0.4     otel_0.2.0         e1071_1.7-17      
    ##  [9] compiler_4.3.3     systemfonts_1.0.5  vctrs_0.7.1        stringr_1.6.0     
    ## [13] pkgconfig_2.0.3    wk_0.9.5           fastmap_1.2.0      magic_1.6-1       
    ## [17] rmarkdown_2.30     nloptr_2.2.1       purrr_1.2.1        xfun_0.56         
    ## [21] aplot_0.2.9        parallel_4.3.3     R6_2.6.1           stringi_1.8.7     
    ## [25] RColorBrewer_1.1-3 boot_1.3-30        Rcpp_1.1.1         iterators_1.0.14  
    ## [29] Matrix_1.6-5       splines_4.3.3      timechange_0.4.0   tidyselect_1.2.1  
    ## [33] rstudioapi_0.18.0  dichromat_2.0-0.1  yaml_2.3.12        doParallel_1.0.17 
    ## [37] codetools_0.2-19   lattice_0.22-5     tibble_3.3.1       plyr_1.8.9        
    ## [41] withr_3.0.2        S7_0.2.1           evaluate_1.0.5     gridGraphics_0.5-1
    ## [45] proxy_0.4-29       xml2_1.5.2         pillar_1.11.1      spAbundance_0.2.1 
    ## [49] KernSmooth_2.23-22 foreach_1.5.2      reformulas_0.4.4   ggfun_0.2.0       
    ## [53] generics_0.1.4     rprojroot_2.1.1    scales_1.4.0       minqa_1.2.8       
    ## [57] class_7.3-22       glue_1.8.0         tools_4.3.3        lme4_2.0-1        
    ## [61] RANN_2.6.2         fs_1.6.7           rbibutils_2.4.1    nlme_3.1-164      
    ## [65] patchwork_1.3.2    Formula_1.2-5      cli_3.6.5          rappdirs_0.3.4    
    ## [69] viridisLite_0.4.3  svglite_2.1.3      gtable_0.3.6       yulab.utils_0.2.4 
    ## [73] digest_0.6.39      classInt_0.4-11    ggplotify_0.1.3    farver_2.1.2      
    ## [77] htmltools_0.5.9    lifecycle_1.0.5

  

#### End of the description
