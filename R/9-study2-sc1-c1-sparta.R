
# -------------------------------------------------------

# Simulation study 2 - 1

# SCENARIO 0 - changing from Bernoulli to a Poisson sampling
# SCENARIO 1 - SCENARIO 0 + X_it = psi in function of L_i

# Using the approach of sparta (spatially uncorrelated random effects)
#  and a Random Walk model for temporal autocorrelation

# ---------------------------------------------------------------

rm(list = ls())
library(spOccupancy)
library(here)
library(abind)
library(jagsUI)

# write the model in BUGS language
# JAGS code for SPARTA model plus random walk prior
# on the year effect of the state model + intercept + halfcauchy hyperpriors
# Source: https://github.com/BiologicalRecordsCentre/sparta/blob/master/inst/models/SPARTA_ranwalk_halfcauchy.txt
# or inverse gamma hyperpriors: https://github.com/BiologicalRecordsCentre/sparta/blob/master/inst/models/SPARTA_ranwalk_intercept_inversegamma.txt

sink ("occ-model-rwalk-RE-simulations.txt")
cat ("

  model{

  # State model -----------------------------------------------
  for (i in 1:nsite){ 
    for (t in 1:nyear){   
      z[i,t] ~ dbern(muZ[i,t]) 
      logit(muZ[i,t])<- beta0 + 
                        beta1 * lat[i,t]+
                        a[t] +  # temporal random walk
                        eta[i]  # site (uncorrelated) random effects
    }}   
  
  # State model priors (temporal random walk)
  a[1] ~ dnorm(0, 0.001)    
  for(t in 2:nyear){
    a[t] ~ dnorm(a[t-1], tau.a)
  }
  
  tau.a <- 1/(sd.a * sd.a) #dgamma(0.001,0.001)    ## decide on this
  # sd.a ~ dt(0, 1, 1)T(0,) #pow(tau.a, -0.5)       ## decide on this
  # or still
  sd.a ~ dunif (0,5) # half uniform prior
  
  # set eta[1] to zero to enable beta0 to be estimated
  eta[1]~dnorm(0,1e-04) # or eta[1]<-0
  # sample from the normal distribution the remaining effects
  for (i in 2:nsite) {
    eta[i] ~ dnorm(0, tau2)       
  }
  
  # precision
  tau2 <- 1/(sigma2 * sigma2) # dgamma(0.001,0.001)   ## decide on this
  #sigma2 ~ dt(0, 1, 1)T(0,) # pow(tau2, -0.5)         ## decide on this 
  sigma2 ~ dunif (0,5) # half uniform prior
  
  # regression coefficients
  beta0 <- log(beta0.u/(1-beta0.u)) ## decide on this
  beta0.u ~ dbeta(2,2)     ## decide on this
  
  # slopes
  beta1 ~ dnorm(0, 0.001)
  
  # Observation model priors  -----------------------------------
  # make it equivalent to other implementations we did based on Doser et al. (2022) spOccupancy R package
  # priors for detection regression coefs
  # regression coefficients
  alpha0 <- log(alpha0.u/(1-alpha0.u)) ## decide on this
  alpha0.u ~ dbeta(2,2)     ## decide on this
  
  # slopes
  alpha1 ~ dnorm(0, 0.001)
  
  # Likelihood
  for (i in 1:nsite){ 
    for (t in 1:nyear){   
      for (j in 1:nrep) {
        
        # detection model
        logit (p[i,t,j]) <- alpha0+
                            alpha1*v_itj[i,t,j]
        
        # observation conditional on occurrence
        y[i,t,j] ~ dbern (z[i,t]*p[i,t,j])
  
      }
    }
  }
      
  # Derived parameters -----------------------------------
  for (t in 1:nyear) {  
    psi.fs[t] <- sum(z[1:nsite, t])/nsite
  }

     }", fill=T)
sink()

# create dir to receive the results -------------------------
dir.create (here ("model_output", "output_simulations", "scenario_two_sparta"))
      
# load simulation settings ------------------------------------------------
load(here("model_output", "output_simulations", "sim-settings.RData"))

# load sampling design
load (file=here ("model_output", 
                 "output_simulations", 
                 "sampling_design_Poisson.rda")
)

# intermediary objects
# psi
# n.sims <- 10
psi.true <- array(NA, dim = c(I, n.time, n.sims, n.scenarios))

# Load the data sets
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# load our new X.p with more secondary occasions ---------------------
load (here ("model_output", "output_simulations", "Xp_itj.rda"))

# Run simulations ---------------------------------------------------------
for (s in 1:n.sims) {
  print(paste("Currently on simulation set ", s, " out of ", n.sims, sep = ''))
  for (sc in 1:n.scenarios) {
    print(paste("Currently on scenario ", sc, " out of ", n.scenarios, sep = ''))
    set.seed(my.seeds[s])
    
    # index
    curr.indx <- (s - 1) * n.scenarios + sc
    dat <- dat.full[[curr.indx]]
    
    # sim data for spOccupancy -------------------------------------------
    # L_i will replace X_it ----------------------
    L_i <- scale(dat.full[[curr.indx]]$coords[, 2])[,1]
    
    # create design matrix for covariates
    X <- array(1, dim = c(I, n.time, length(beta)))
    X[, , 2] <- L_i # set latitude here
    
    # Linear predictor for occupancy (psi)
    logit_psi <- array(0, dim = c(I, n.time))
    
    # simulate occupancy         
    for (t in 1:n.time) {
      for (i in 1:I) {
        logit_psi[i, t] <- X[i, t, ] %*% beta + 
          dat.full[[curr.indx]]$w[i,] + 
          dat.full[[curr.indx]]$eta[t,]
      }
    }
    
    # Logistic function to get occupancy probabilities
    psi <- 1 / (1 + exp(-logit_psi))
    psi.true[, , s, sc] <- psi
    
    # Generate binary occupancy data (sample from a Binomial)
    Z <- matrix(rbinom(I * n.time, 
                       size = 1, 
                       prob = psi), 
                nrow = I, 
                ncol = n.time)
    
    # Linear predictor for detection (p)
    logit_p <- array(0, dim = c(I, n.time, J))
    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p[i, t, j] <-  X.p[i, t, j, ] %*% alpha
          
        }
      }
    }
    
    # Logistic function to get detection probabilities
    p <- 1 / (1 + exp(-logit_p))
    
    # Sample detection non-detection data
    y <- array(dim = c(I, n.time, J))
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          y[i, t, j] <- rbinom(1, size = Z[i, t], prob = p[i, t, j])
          
        }
      }
    }
    
    # Imposing our mixed design --------------
    rep.indx <- ifelse(G_itj >= 1, 2, 1)
    
    # Prep the data for spOccupancy -------------------------------------------
    # Site x Replicate
    y_mixed <-y # [, , 1:2, drop = FALSE]
    
    # Our Mixed Design (only in the second survey)
    for (t in 1:n.time) {
      for (j in 1:dim(y_mixed)[3]) {
        y_mixed[which(rep.indx[,t,j] == 1), t, j] <- NA # set NAs
      }
    }
    
    # Coordinates
    coords <- dat$coords
    
    # bundle data into a list
    jags.data <- list( y = y_mixed,
                       nsite = nrow(y_mixed),
                       nyear = ncol(y_mixed),
                       nrep = dim(y_mixed)[3],
                       
                       # site covs
                       lat = X[ , , 2],
                       
                       # det covs
                       v_itj=X.p[, , , 2] 
                       )
    str(jags.data)
    
    # Starting values
    z.init <- apply(jags.data$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    inits <- function () {list(z = z.init, 
                               beta0 = rnorm (0,sqrt(2.72)),
                               alpha0 = rnorm (0,sqrt(2.72)),
                               sigma2 = runif(1,0,5),
                               sd.a = runif(1,0,5))}
    
    # parameters to track
    params <- c("beta0","beta1",
                # precision 
                "tau2","sigma2","tau.a","sd.a",
                # detection
                "alpha0","alpha1",
                # derived
                "psi.fs",
                # RE
                "eta","a",
                # psi
                "muZ",
                # p
                "p"
    )
    
    # Fit the model with jags
    out <- jags (jags.data, 
                 inits = inits, 
                 parameters.to.save = params,
                 model.file = "occ-model-rwalk-RE-simulations.txt",
                 n.chains = 3, 
                 n.thin = 10, 
                 n.iter = 10000, 
                 n.burnin = 5000,
                 parallel=T)
    
    # save
    save(out$mean,
         out$summary,
         psi.true, beta, scenario.vals, file = here ("model_output", 
                                                     "output_simulations", 
                                                     "scenario_two_sparta",
                                                     paste0("sim-mixed-stPGOcc-results_",curr.indx,".rda"))
    )
    
    
    
  } # i (n.scenarios)
} # j (n.sims)

# end
