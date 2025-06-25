Auxiliary information
================
2025-05-22













### Occupancy model with spatially uncorrelated/unstructured random effects

Below we describe the multi-season site occupancy model of (Outhwaite et
al. 2018) which we used in the sensitivity analyses. This model
represents an improvement from the occupancy model used in (Isaac et al.
2014) in the sense that it uses a random walk (normally distributed)
prior for the temporal random effect parameter $b_t$, thus enabling the
sharing of occupancy information between consecutive years to estimate
occupancy trends. The model can be described in the following manner
(respecting the hierarchical structure we used in our main text, and
using parameter notation as (Outhwaite et al. 2018)).

The first part of the model is an occupancy state process, where the
latent occupancy state $z_{it}$ of a single species at $i=1, ..., I$
sites and $t=1, ..., T$ primary occasions is modeled. The state $z_{it}$
is drawn from a Bernoulli distribution depending on the probability of
occupancy of the site $\psi_{it}$. The $\psi_{it}$ is modeled in
function of covariates varying over sites and / or primary occasions
$X_{it}$, and $\boldsymbol{\beta}$ is a vector of regression
coefficients including the intercept and slope (eq. 1):

$$
\tag{1}
\begin{equation}
\begin{split}
z_{it}|\psi_{it} \sim \mathcal{B}(\psi_{it}), \\
\text{logit}(\psi_{it}) = X_{it}^T\boldsymbol{\beta} + u_{i} + b_{t}
\end{split}
\end{equation}
$$

The original model formulation has no covariates. Some formulations
include one intercept coming from a Beta ($\mathcal{B}e$) distribution.
We also changed the formulation of the spatially uncorrelated random
effects $u_{i}$ , which here were defined through two Normal
distributions with $\mu=0$ and variance $\sigma^2_u$ for each site $i$
to ensure the identifiability of the model intercept. Thus we used eq.
2:

$$
\tag{2}
\begin{equation}
\begin{split}
u_i \sim 
\begin{cases} 
\mathcal{N}(0,10^4) & \text{for } i=1 \\
\mathcal{N}(u_{i},\tau_u) & \text{for } i>1
\end{cases}
\end{split}
\end{equation}
$$

where $\tau_u$ is precision parameter calculated as $1/\sigma^2_u$. This
parameter controls the variance of occupancy probability across sites.
The prior for $u_1$ is Normal with $\mu_b=0$ and variance $10^4$, and
priors for successive sites are determined via the second term in eq. 2.
The prior for $\sigma^2_u$ was $\sigma_u\sim U(0,5)$.

The temporal random effect $b_{t}$ follows a random walk model defined
in eq. 3:

$$
\tag{3}
\begin{equation}
\begin{split}
b_t \sim 
\begin{cases} 
\mathcal{N}(\mu_b,10^3) & \text{for } t=1 \\
\mathcal{N}(b_{t-1},\tau_b) & \text{for } t>1
\end{cases}
\end{split}
\end{equation}
$$

where $\tau_b$ is precision parameter calculated as $1/\sigma^2_b$, this
parameter controls the variance of species occupancy probability in time
$t$ relative to the previous time $t-1$. The prior for $b_1$ came from a
Normal distribution with $\mu_b=0$ and variance $10^3$, and priors for
successive time points are determined via the second term in eq. 2. The
prior for $\sigma^2_b$ was $\sigma_b\sim U(0,5)$.

The detection model of Outhwaite et al. (2018) includes the year effect
in the model intercept. Since we want something similar to Doser and
Stoudt (2024) model, the detection model we implemented eq. 4:

$$
\tag{4}
\begin{equation}
\begin{split}
y_{itj} | z_{it} \sim \mathcal{B}(z_{it} \times p_{itj}) ,\\ 
\text{logit}(p_{itj}) = v_{itj}^T\boldsymbol{\alpha}
\end{split}
\end{equation}
$$

where $\boldsymbol{\alpha}$ is a vector of regression coefficients,
including the intercept and slopes that represent the effect of the
covariates $v_{itj}$ on $p_{itj}$ (eq. 4). The species encounter history
$y_{itj}$ is then conditionally related to the latent occupancy
$z_{it}$, meaning that for a truly occupied site and primary occasion
$z_{it}=1$, the species will be detected in one individual secondary
occasion with probability $p_{j}$ (eq. 4). If unoccupied $z_{it}=0$ then
the species can not be detected.

Prior choice has a large influence on the performance of Outhwaite’s
model (see supplementary material in (Outhwaite et al. 2018)). The
original model formulation had
$\beta_0 \sim \mathcal{B}e(\alpha,\beta)$, with $\alpha$ and $\beta$
controlling the shape of the distribution (not to be confused with model
coefficients presented above). In our analyzes, we followed the authors
and set both parameters to $\alpha, \beta = 2$. The same prior was set
to the detection intercept $\alpha_0$. The inverse logit function was
then applied to this quantity. Priors for regression slopes (from both
occupancy and detection models) were taken from a Normal distribution
with average and precision as $\mathcal{N} (0,0.001)$.

This model was fitted to 42 simulated data sets under the design of
study 2 - scenario 1 (Poisson sampling design, effect of latitude $L_i$
on occupancy, effect of observation-level covariate $v_{itj}$ on
detection). The model was written in BUGS language and analyzed in JAGS
through the R package jagsUI (Kellner 2024). Parameter estimates using
the real data set was done with three parallel MCMC chains of 10,000
iterations each, burn-in of 5,000 iterations, and thinning each 10
iterations (yielding 1,500 posterior distribution draws in total). It
took 18 days to run the model for 42 data sets x 16 autocorrelation
scenarios per data set using a Dell Inc. Precision Tower 3620 with
processor Intel(R)Xeon(R)E3-1240 v5x8 with 32GB RAM and 8 cores.

<img src="figures/TuningD&S_sims/Scenario1_sparta/Figure-3-st2-sc1-sparta.png" width="100%" height="100%" />
Figure 1: Relationship between the true site occupancy probability
$\psi_{it}$ (X-axis) and the estimated site occupancy probability
$\hat{\psi_{it}}$ (Y-axis). Spatial and temporal autocorrelation sub
scenarios are represented along columns and rows, respectively, with the
following levels: high $\phi=15$, low $\phi=3.75$, high $\rho=0.9$, low
$\rho=0.5$, high $\sigma^2_T=1.5$, low $\sigma^2_T=0.3$, high
$\sigma^2=1.5$, low $\sigma^2=0.3$. These simulations used data from
study 1-scenario 1. Each gray line represents the locally estimated
scatter plot smoothing (LOESS) relationship between $\psi_{it}$ and
$\hat{\psi_{it}}$ per simulated data set. The black line depicts the
averaged relationship across the 42 simulated data sets, and the red
dashed line depicts an 1:1 relationship.

<img src="figures/TuningD&S_sims/Scenario1_sparta/maps_comparison_sparta.png" width="100%" height="100%" />
Figure 2: Spatial map of true site occupancy probability $\psi_{it}$
(top), estimated site occupancy probability $\hat{\psi_{it}}$ (middle),
and the difference between them (bottom) across scenarios of high (left)
and low spatial and temporal autocorrelation (right). Estimates of
$\hat{\psi_{it}}$ were produced by the model with spatially uncorrelated
site random effects. Positive differences (in blue) indicate
overestimation relative to the true occupancy, whereas negative
differences (in red) indicate underestimate relative to the true
occupancy. Results for the first simulated data set and the first
primary period of estimation ($t=1$). Data shown in the left plot were
simulated with $\phi=3.75$, $\sigma^2=0.3$, $\rho=0.5$, $\sigma^2_T=0.3$
(strong spatial autocorrelation and weak temporal autocorrelation),
whereas in the right the values of these parameters were $\phi=15$,
$\sigma^2=1.5$, $\rho=0.9$, $\sigma^2_T=1.5$ (weak spatial
autocorrelation and strong temporal autocorrelation).

<img src="figures/TuningD&S_sims/Scenario1_sparta/coeff_densities.png" width="100%" height="100%" />

Figure 3: Density plots showing the model coefficients estimated by the
model with spatially uncorrelated random effects. The black vertical
line depicts the truth of each parameter. Spatial and temporal
autocorrelation sub scenarios are represented along columns and rows,
respectively. The density was created using the output of the analysis
of 42 data sets simulated with $\phi=15$, $\sigma^2=1.5$, $\rho=0.9$,
$\sigma^2_T=1.5$.

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-doser2024fractional" class="csl-entry">

Doser, Jeffrey W., and Sara Stoudt. 2024. “‘Fractional Replication’ in
Single‐visit Multi‐season Occupancy Models: Impacts of Spatiotemporal
Autocorrelation on Identifiability.” *Methods in Ecology and Evolution*
15 (2): 358–72. <https://doi.org/10.1111/2041-210X.14275>.

</div>

<div id="ref-isaac2014statistics" class="csl-entry">

Isaac, Nick J. B., Arco J. Van Strien, Tom A. August, Marnix P. De
Zeeuw, and David B. Roy. 2014. “Statistics for Citizen Science:
Extracting Signals of Change from Noisy Ecological Data.” Edited by
Barbara Anderson. *Methods in Ecology and Evolution* 5 (10): 1052–60.
<https://doi.org/10.1111/2041-210X.12254>.

</div>

<div id="ref-kellner2024jagsUI" class="csl-entry">

Kellner, Ken. 2024. *jagsUI: A Wrapper Around ’Rjags’ to Streamline
’JAGS’ Analyses*. <https://CRAN.R-project.org/package=jagsUI>.

</div>

<div id="ref-outhwaite2018prior" class="csl-entry">

Outhwaite, Charlotte L, Richard E Chandler, Gary D Powney, Ben Collen,
Richard D Gregory, and Nick JB Isaac. 2018. “Prior Specification in
Bayesian Occupancy Modelling Improves Analysis of Species Occurrence
Data.” *Ecological Indicators* 93: 333–43.
<https://www.sciencedirect.com/science/article/pii/S1470160X18303443?via%3Dihub>.

</div>

</div>
