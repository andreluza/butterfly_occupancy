# -----------------------------------------------------------------------------

# Simulation study 1 - scenario 0: replication of D&S Simulations

# -------------------------------------------------------------------------------
# Original description
# main-sim-mixed.R: script to run a set of simulations to assess identifiability of a "mixed design" multi-season occupancy model under varying levels of spatial and temporal autocorrelation. 
# Authors: Jeffrey W. Doser and Sara Stoudt
# Approximate run time: 10 days
# -------------------------------------------------------------------------------

rm(list = ls())
library(spOccupancy)
library(here)

# MCMC settings
n.samples <- 25000
batch.length <- 25
n.burn <- 15000
n.thin <- 10
n.chains <- 3 # changed (previously was 1 chain)
accept.rate <- 0.43
n.batch <- (n.samples / batch.length)*n.chains

# Parameters for simulation -----------------------------------------------
# Number of data sets for each scenario
n.sims <- 100
# Spatial locations (1200 total)
J.x <- 30
J.y <- 40
I <- J.x * J.y
# Number of years (keep constant for now, but will likely want to explore 
#                  this as well). 
n.time <- 10
# Occurrence coefficient --------------
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
# NOTE: simTOcc generates data across a unit square. When using an exponential 
#       correlation function, the effective spatial range (distance at which 
#       correlation between sites drops to 0.05) is 3 / phi. Thus, the following
#       values correspond to effective spatial ranges of 20% and 80% of the 
#       study region. 
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
# psi
psi.true <- array(NA, dim = c(I, n.time, n.sims, n.scenarios))
psi.mean.samples <- array(NA, dim = c(I, n.time, n.sims, n.scenarios))
psi.low.samples <- array(NA, dim = c(I, n.time, n.sims, n.scenarios))
psi.high.samples <- array(NA, dim = c(I, n.time, n.sims, n.scenarios))
# reg coeff occupancy
beta.mean.samples <- array(NA, dim = c(p.occ, n.sims, n.scenarios))
beta.low.samples <- array(NA, dim = c(p.occ, n.sims, n.scenarios))
beta.high.samples <- array(NA, dim = c(p.occ, n.sims, n.scenarios))
# reg coeff detection
alpha.mean.samples <- array(NA, dim = c(length(alpha), n.sims, n.scenarios))
alpha.low.samples <- array(NA, dim = c(length(alpha), n.sims, n.scenarios))
alpha.high.samples <- array(NA, dim = c(length(alpha), n.sims, n.scenarios))
# spatial random effect
w.mean.samples <- array(NA, dim = c(n.sims, I, n.scenarios))
w.low.samples <- array(NA, dim = c(n.sims, I, n.scenarios))
w.high.samples <- array(NA, dim = c(n.sims, I, n.scenarios))
# temporal autocorrelation
eta.mean.samples <- array(NA, dim = c(n.sims, n.time, n.scenarios))
eta.low.samples <- array(NA, dim = c(n.sims, n.time, n.scenarios))
eta.high.samples <- array(NA, dim = c(n.sims, n.time, n.scenarios))
# spatial autocorrelation
theta.mean.samples <- array(NA, dim = c(n.sims, ncol(scenario.vals), n.scenarios))
theta.low.samples <- array(NA, dim = c(n.sims, ncol(scenario.vals), n.scenarios))
theta.high.samples <- array(NA, dim = c(n.sims, ncol(scenario.vals), n.scenarios))
# all samples
theta.samples <- array(NA, dim = c(n.batch,ncol(scenario.vals),n.sims,  n.scenarios))
eta.samples <- array(NA, dim = c(n.batch, n.time,n.sims, n.scenarios))

# Load the data sets
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# Run simulations ---------------------------------------------------------
for (j in 1:n.sims) {
  print(paste("Currently on simulation set ", j, " out of ", n.sims, sep = ''))
  for (i in 1:n.scenarios) {
    print(paste("Currently on scenario ", i, " out of ", n.scenarios, sep = ''))
    
    # index
    curr.indx <- (j - 1) * n.scenarios + i
    dat <- dat.full[[curr.indx]]
    rep.indx <- matrix(ifelse(rbinom(I * n.time, 1, 0.1) == 1, 2, 1), I, n.time)
    psi.true[, , j, i] <- dat$psi
    phi.tune <- 0.5
    # Prep the data for spOccupancy -------------------------------------------
    # Site x Replicate
    y <- dat$y[, , 1:2, drop = FALSE]
    # Occurrence Covariates
    X <- dat$X
    # Detection Covariates
    X.p <- dat$X.p[, , 1:2, , drop = FALSE]
    # Mixed design
    for (t in 1:n.time) {
      y[which(rep.indx[, t] == 1), t, 2] <- NA
      X.p[which(rep.indx[, t] == 1), t, 2, ] <- NA
    }
    
    # Coordinates
    coords <- dat$coords
    # Package all data into a list
    occ.covs <- list(int = X[, , 1],
                     occ.cov.1 = X[, , 2])
    det.covs <- list(det.cov.1 = X.p[, , , 2])
    data.list <- list(y = y, occ.covs = occ.covs, det.covs = det.covs, coords = coords)
    # Priors
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                       alpha.normal = list(mean = 0, var = 2.72),
                       sigma.sq.ig = c(a = 2, b = 1.5),
                       sigma.sq.t.ig = c(2, 1),
                       phi.unif = c(a = 3 / 1, b = 3 / 0.05))
    # Starting values
    z.init <- apply(y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    inits.list <- list(beta = 0, alpha = 0, sigma.sq.t = 0.5, phi = 3 / .5, 
                       sigma.sq = 1, rho = 0, z = z.init)
    # Tuning
    tuning.list <- list(phi = phi.tune, rho = 0.5)
    
    # Fit the model with stPGOcc
    out <- stPGOcc(occ.formula = ~ occ.cov.1,
                   det.formula = ~ det.cov.1,
                   data = data.list,
                   n.batch = n.batch/n.chains,
                   batch.length = batch.length,
                   inits = inits.list,
                   priors = prior.list,
                   accept.rate = 0.43,
                   cov.model = cov.model,
                   tuning = tuning.list,
                   n.omp.threads = 1,
                   verbose = TRUE,
                   ar1 = TRUE,
                   NNGP = TRUE,
                   n.neighbors = 5,
                   n.report = 25,
                   n.burn = n.burn,
                   n.thin = n.thin,
                   n.chains = n.chains)
    
    # occupancy probability
    psi.mean.samples[, , j, i] <- apply(out$psi.samples, c(2, 3), mean)
    psi.low.samples[, , j, i] <- apply(out$psi.samples, c(2, 3), quantile, 0.025)
    psi.high.samples[, , j, i] <- apply(out$psi.samples, c(2, 3), quantile, 0.975)
    # regression coeff occupancy
    beta.mean.samples[, j, i] <- apply(out$beta.samples, 2, mean)
    beta.low.samples[, j, i] <- apply(out$beta.samples, 2, quantile, 0.025)
    beta.high.samples[, j, i] <- apply(out$beta.samples, 2, quantile, 0.975)
    # regression coeff detection
    alpha.mean.samples[, j, i] <- apply(out$alpha.samples, 2, mean)
    alpha.low.samples[, j, i] <- apply(out$alpha.samples, 2, quantile, 0.025)
    alpha.high.samples[, j, i] <- apply(out$alpha.samples, 2, quantile, 0.975)
    # spatial random effect
    w.mean.samples[j, , i] <- apply(out$w.samples, 2, mean)
    w.low.samples[j, , i] <- apply(out$w.samples, 2, quantile, 0.025)
    w.high.samples[j, , i] <- apply(out$w.samples, 2, quantile, 0.975)
    # temporal random effect
    eta.mean.samples[j, , i] <- apply(out$eta.samples, 2, mean)
    eta.low.samples[j, , i] <- apply(out$eta.samples, 2, quantile, 0.025)
    eta.high.samples[j, , i] <- apply(out$eta.samples, 2, quantile, 0.975)
    # theta (phi and sigma^2)
    theta.mean.samples[j, , i] <- apply(out$theta.samples, 2, mean)
    theta.low.samples[j, , i] <- apply(out$theta.samples, 2, quantile, 0.025)
    theta.high.samples[j, , i] <- apply(out$theta.samples, 2, quantile, 0.975)
    # all draws of theta, eta and w
    theta.samples[,,j,i] <- out$theta.samples
    eta.samples [,,j,i]<- out$eta.samples
    
    # save each run
    save(theta.samples, 
     eta.samples,
     psi.mean.samples, 
     psi.low.samples, 
     psi.high.samples,
     beta.mean.samples, 
     beta.low.samples, 
     beta.high.samples,
     alpha.mean.samples, 
     alpha.low.samples, 
     alpha.high.samples,
     w.mean.samples,
     w.low.samples,
     w.high.samples,
     eta.mean.samples,
     eta.low.samples, 
     eta.high.samples,
     theta.mean.samples, 
     theta.low.samples, 
     theta.high.samples,
     psi.true, 
     beta, 
     alpha,
     scenario.vals, file = here ("model_output", 
                                 "output_simulations", 
                                 "sims_D&S",
                                 paste0 ("sim-mixed-stPGOcc-results-SimSce",curr.indx,".rda"))
    )
    
    # the loop makes the array each time larger, so
    # remove files to avoid reach the upper memory bound
    # remove three output before the present output
    unlink (here ("model_output", 
                  "output_simulations", 
                  "sims_D&S",
                  paste0("sim-mixed-stPGOcc-results-SimSce",curr.indx-3,".rda")),recursive=F)

  } # i (n.scenarios)
} # j (n.sims)
