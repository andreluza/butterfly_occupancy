# ---------------------------------------------------------------

# Settings for simulations - used in studies 2 to 3
# on study 2 - scenario 3 we added one coefficient in alpha object (coefficients of the detection model)

# ---------------------------------------------------------------

# load packages & functions
rm(list=ls())
source ("R/packages.R")
source ("R/functions.R")


# create directory to receive results of the scenarios ------------------
dir.create(here ("model_output", "output_simulations","scenario_zero"))
dir.create(here ("model_output", "output_simulations","scenario_one"))
dir.create(here ("model_output", "output_simulations","scenario_two"))
dir.create(here ("model_output", "output_simulations","scenario_three"))
dir.create(here ("model_output", "output_simulations","scenario_phenology"))
dir.create(here ("model_output", "output_simulations","scenario_phenology_spot"))


# ggplot theme -------------------------------------
my_theme <- theme(legend.position = 'bottom', 
                  strip.text = element_text(size=12),
                  strip.text.y = element_text(color = 'black'),
                  strip.text.x = element_text(color = 'black'), 
                  text = element_text(family="LM Roman 10"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
                  axis.text.y = element_text(size = 10),
                  axis.title = element_text(size=15))


# MCMC settings -----------------------------------------------
n.samples <- 25000
batch.length <- 25
n.burn <- 15000
n.thin <- 10
n.chains <- 3 # changed (previously was 1 chain)
accept.rate <- 0.43
n.batch <- (n.samples / batch.length)*n.chains

# short MCMC settings -----------------------------------------------
#n.samples <- 100
#batch.length <- 25
#n.burn <- 50
#n.thin <- 1
#n.chains <- 3 # changed (previously was 1 chain)
#accept.rate <- 0.43
#n.batch <- (n.samples / batch.length)*n.chains

# Parameters for simulation -----------------------------------------------
# Number of data sets for each scenario
n.sims <- 100
set.seed(10110)
my.seeds <- sample(1:100000, n.sims, replace = FALSE)

# Spatial locations
I.x <- 30
I.y <- 40
I <- I.x * I.y

# Matrix of spatial locations
s.x <- seq(0, 1, length.out = I.x)
s.y <- seq(0, 1, length.out = I.y)
coords <- as.data.frame(expand.grid(s.x, s.y))
coords$Site <- seq(1,nrow(coords))
# Number of years (keep constant for now, but will likely want to explore 
#                  this as well). 
n.time <- 10

# Number of secondary occasions
J<-10

# Occurrence coefficients ---------------------
# Generate a single covariate
beta <- c(0, 0.5)
p.occ <- length(beta)
# Detection coefficient ---------------
# A single covariate on detection
alpha <- c(0, -0.5)
# Spatial parameters ------------------
sp <- TRUE
# Assume an exponential correlation model for now. 
cov.model <- 'exponential'
# Spatial variances
sigma.sq.vals <- c(0.3, 1.5)
# Spatial decay
phi.vals <- c(3 / .2, 3 / .8)
# Temporal parameters -----------------
rho.vals <- c(0.5, 0.9)
sigma.sq.t.vals <- c(0.3, 1.5)
# Total number of simulation scenarios
n.scenarios <- length(sigma.sq.vals) * length(phi.vals) * 
  length(rho.vals) * length(sigma.sq.t.vals)
# Different combinations of all the four parameters that vary
scenario.vals <- expand.grid(sigma.sq = sigma.sq.vals, phi = phi.vals, 
                             rho = rho.vals, sigma.sq.t = sigma.sq.t.vals)

# save
save.image(here("model_output", "output_simulations", "sim-settings.RData"))

# end
rm(list=ls())
