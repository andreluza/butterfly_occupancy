
# -------------------------------------------------------

# Simulation study 3 - scenario 2

# PS: For review only -- number of neighbors increased from 5 to 15

# -------------------------------------------------------

# Scenario phenology & observer preferences for mid season + Observation spot

## Study 2
# SCENARIO 0 - changing from Bernoulli to a Poisson sampling
# SCENARIO 1 - SCENARIO 0 + X_it = psi in function of L_i
# SCENARIO 2 - SCENARIO 0 + SCENARIO 1 + p in function of L_i and V_itj
# SCENARIO 3 - SCENARIO 0 + SCENARIO 1 + SCENARIO 2 + p in function of L_i

## Study 3
# SCENARIO 1 - SCENARIO 0 + SCENARIO 1 + SCENARIO 2 + SCENARIO 3 + observer preferences for mid season & J=10
# SCENARIO 2 - SCENARIO 0 + SCENARIO 1 + SCENARIO 2 + SCENARIO 3 + phenology & observer preferences for mid season  & J=10 + spot

# ---------------------------------------------------------------

rm(list = ls())
library(spOccupancy)
library(here)
library(abind)
library(reshape)
library(gridExtra)
library(ggplot2)
library(dplyr)

# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# load sampling design
load (file=here ("model_output", 
                 "output_simulations", 
                 "sampling_design_Poisson_phenology_spot.rda")
)

# Detection coefficient ---------------
# A single covariate on detection
alpha <- c(0, -0.5, -0.5) # added one coefficient for latitude

# intermediary objects
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
theta.samples <- array(NA, dim = c(3000,ncol(scenario.vals),n.sims,  n.scenarios))
eta.samples <- array(NA, dim = c(3000, n.time,n.sims, n.scenarios))

# Load the data sets
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# load our new X.p with more secondary occasions ---------------------
load (here ("model_output", "output_simulations", "Xp_itj.rda"))

# create folder to host results
dir.create (here ("model_output", 
      "output_simulations", 
      "scenario_phenology_spot-NNGP-15-review"))

# Run simulations ---------------------------------------------------------
for (s in 1:20) { # use 20
  print(paste("Currently on simulation set ", s, " out of ", n.sims, sep = ''))
  for (sc in 1:n.scenarios) {
    print(paste("Currently on scenario ", sc, " out of ", n.scenarios, sep = ''))
    set.seed(my.seeds[s])
    # index
    curr.indx <- (s - 1) * n.scenarios + sc
    dat <- dat.full[[curr.indx]]
    #psi.true[, , s, sc] <- dat$psi
    phi.tune <- 0.5
    
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
    
    # bind latitude
    lat_bind <- replicate (J,X[, , 2])
    V_itj_L_i <- abind(X.p,lat_bind)
    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p[i, t, j] <- V_itj_L_i [i, t, j, ] %*% alpha
          
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
    y_mixed <-y# [, , 1:2, drop = FALSE]
    
    # Our Mixed Design (only in the second survey)
    for (t in 1:n.time) {
      for (j in 1:dim(y_mixed)[3]) {
        y_mixed[which(rep.indx[,t,j] == 1), t, j] <- NA # set NAs
      }
    }
    
    # identify missing data
    missing_sites <- apply (y_mixed,1,sum,na.rm=T)
    
    # Package all data into a list
    occ.covs <- list(int = X[which(missing_sites>0), , 1],
                     occ.cov.1 = X[which(missing_sites>0), , 2])
    det.covs <- list(int = V_itj_L_i[which(missing_sites>0), , , 1],
                     det.cov.1=V_itj_L_i[which(missing_sites>0), , , 2],
                     det.cov.2=V_itj_L_i[which(missing_sites>0), , , 3]
    )
    str(data.list <- list(y = y_mixed[which(missing_sites>0),,], 
                          occ.covs = occ.covs, det.covs = det.covs, 
                          coords = dat$coords[which(missing_sites>0),]))
    
    # Priors
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                       alpha.normal = list(mean = 0, var = 2.72),
                       sigma.sq.ig = c(a = 2, b = 1.5),
                       sigma.sq.t.ig = c(2, 1),
                       phi.unif = c(a = 3 / 1, b = 3 / 0.05))
    
    # Starting values
    z.init <- apply(y_mixed[which(missing_sites>0),,], c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
    inits.list <- list(beta = 0, alpha = 0, sigma.sq.t = 0.5, phi = 3 / .5, 
                       sigma.sq = 1, rho = 0, z = z.init)
    
    # Tuning
    tuning.list <- list(phi = phi.tune, rho = 0.5)
    
    # Fit the model with stPGOcc
    out <- stPGOcc(occ.formula = ~ occ.cov.1,
                   det.formula = ~ det.cov.1+det.cov.2,
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
                   n.neighbors = 15, # number of neighbors increase from 5 to 15
                   n.report = 25,
                   n.burn = n.burn,
                   n.thin = n.thin,
                   n.chains = n.chains)
    
    
    # make predictions
    out.pred.occ <- predict(out, X, dat$coords, verbose = FALSE,t.cols=seq(1,n.time),
                        type = 'occupancy')
    
    # occupancy probability
    psi.mean.samples[, , s, sc] <- apply(out.pred.occ$psi.0.samples, c(2, 3), mean,na.rm=T)
    psi.low.samples[, , s, sc] <- apply(out.pred.occ$psi.0.samples, c(2, 3), quantile, 0.025,na.rm=T)
    psi.high.samples[, , s, sc] <- apply(out.pred.occ$psi.0.samples, c(2, 3), quantile, 0.975,na.rm=T)
    # regression coeff occupancy
    beta.mean.samples[, s, sc] <- apply(out$beta.samples, 2, mean,na.rm=T)
    beta.low.samples[, s, sc] <- apply(out$beta.samples, 2, quantile, 0.025,na.rm=T)
    beta.high.samples[, s, sc] <- apply(out$beta.samples, 2, quantile, 0.975,na.rm=T)
    # regression coeff detection
    alpha.mean.samples[, s, sc] <- apply(out$alpha.samples, 2, mean,na.rm=T)
    alpha.low.samples[, s, sc] <- apply(out$alpha.samples, 2, quantile, 0.025,na.rm=T)
    alpha.high.samples[, s, sc] <- apply(out$alpha.samples, 2, quantile, 0.975,na.rm=T)
    # spatial random effect
    w.mean.samples[s, , sc] <- apply(out.pred.occ$w.0.samples, 2, mean,na.rm=T)
    w.low.samples[s, , sc] <- apply(out.pred.occ$w.0.samples, 2, quantile, 0.025,na.rm=T)
    w.high.samples[s, , sc] <- apply(out.pred.occ$w.0.samples, 2, quantile, 0.975,na.rm=T)
    # temporal random effect
    eta.mean.samples[s, , sc] <- apply(out$eta.samples, 2, mean,na.rm=T)
    eta.low.samples[s, , sc] <- apply(out$eta.samples, 2, quantile, 0.025,na.rm=T)
    eta.high.samples[s, , sc] <- apply(out$eta.samples, 2, quantile, 0.975,na.rm=T)
    # theta (phi and sigma^2)
    theta.mean.samples[s, , sc] <- apply(out$theta.samples, 2, mean,na.rm=T)
    theta.low.samples[s, , sc] <- apply(out$theta.samples, 2, quantile, 0.025,na.rm=T)
    theta.high.samples[s, , sc] <- apply(out$theta.samples, 2, quantile, 0.975,na.rm=T)
    # all draws of theta, eta and w
    theta.samples[,, s, sc] <- out$theta.samples
    eta.samples [,, s, sc]<- out$eta.samples
    
    # save
    save(psi.mean.samples, psi.low.samples, psi.high.samples,
         beta.mean.samples, beta.low.samples, beta.high.samples,
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
         theta.samples,
         eta.samples,
         # true
         psi.true, beta, scenario.vals, file = here ("model_output", 
                                                     "output_simulations", 
                                                     "scenario_phenology_spot-NNGP-15-review",
                                                     paste0("sim-mixed-stPGOcc-results_",curr.indx,".rda"))
    )
    
    # the loop makes the array each time larger, so
    # remove files to avoid reach the upper memory bound
    # remove three output before the present output
    unlink (here ("model_output", 
                  "output_simulations", 
                  "scenario_phenology_spot-NNGP-15-review",
                  paste0("sim-mixed-stPGOcc-results_",curr.indx-3,".rda")),recursive=F)
    
    
  } # i (n.scenarios)
} # j (n.sims)


# ------------------------------------------------------
# illustrate the scenario

s=1
sc=1

# index
curr.indx <- (s - 1) * n.scenarios + sc
dat <- dat.full[[curr.indx]]
psi.true[, , s, sc] <- dat$psi
phi.tune <- 0.5

# sim data for spOccupancy -------------------------------------------
# L_i will replace X_it ----------------------
L_i <- scale(dat$coords[, 2])[,1]

# create design matrix for covariates
X <- array(1, dim = c(I, n.time, length(beta)))
X[, , 2] <- L_i # set latitude here

# Linear predictor for occupancy (psi)
logit_psi <- array(0, dim = c(I, n.time))

# simulate occupancy         
for (t in 1:n.time) {
  for (i in 1:I) {
    logit_psi[i, t] <- X[i, t, ] %*% beta + 
      dat$w[i,] + 
      dat$eta[t,]
  }
}

# Logistic function to get occupancy probabilities
psi <- 1 / (1 + exp(-logit_psi))

# Generate binary occupancy data (sample from a Binomial)
Z <- matrix(rbinom(I * n.time, 
                   size = 1, 
                   prob = psi), 
            nrow = I, 
            ncol = n.time)

# Reshape for plotting
occupancy_data <- data.frame(
  X = rep(dat$coords[,1], n.time),
  Y = rep(dat$coords[,2], n.time),
  Time = rep(1:n.time, each = I),
  Occupied = as.vector(Z)
)

# Linear predictor for detection (p)
logit_p <- array(0, dim = c(I, n.time, J))

# bind latitude
lat_bind <- replicate (J,X[, , 2])
V_itj_L_i <- abind(X.p,lat_bind)

# simulate detection    
for (t in 1:n.time) {
  for (i in 1:I) {
    for (j in 1:J) {
      
      # 0 in the logit scale will produce p = 0.5
      logit_p[i, t, j] <- V_itj_L_i [i, t, j, ] %*% alpha
      
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
y_mixed <-y# [, , 1:2, drop = FALSE]

# Our Mixed Design (only in the second survey)
for (t in 1:n.time) {
  for (j in 1:dim(y_mixed)[3]) {
    y_mixed[which(rep.indx[,t,j] == 1), t, j] <- NA # set NAs
  }
}

# Reshape for plotting
observation_data_full <- data.frame(
  
  X = rep(dat$coords[,1],20),
  Y = rep(dat$coords[,2],20),
  Site=melt(y_mixed)[,1],
  Time = melt(y_mixed)[,2],#rep(1:nT, each = M * J),
  Survey = melt(y_mixed)[,3], #rep(1:J, each = M, times = nT),
  Observation = melt(y_mixed)[,4], #as.vector(y)
  ObsCov = melt(X.p[, , , 2])
  
)

# --------------------------------
# plot variables (one dataset)
# Plot the true occupancy data for the first year
proj_Z <- grid.arrange (
  
  ggplot(occupancy_data %>%
           filter (Time==1), 
         aes(x = X, y = Y, fill = factor(Occupied))) +
    geom_tile() +
    scale_fill_manual(values = c("black", "violet"),na.value = "gray90", 
                      labels = c("Not occupied", "Occupied")) +
    coord_fixed() +
    
    labs(title = "True occupancy (for t=1)",
         x = "Longitude", y = "Latitude", fill = "Detection")+
   # my_theme+
    theme(legend.position = "right",
          legend.text = element_text(angle=0))
  
  ,
  
  
  # fluctuation in occupancy over time
  data.frame (year=seq(1,n.time),
              Zt=colSums(Z)/nrow(Z)) %>%
    ggplot ()+
    geom_line(aes(year, Zt))+
    geom_point(aes(year, Zt))+
    #my_theme+
    labs (x="Year Index",y="Proportion of\n truly occuppied sites")
  ,
  
  ncol=2
)

# plot
proj_X <- grid.arrange (
  
  # Plot the latitude
  ggplot(X[, 1, 2] %>%
           data.frame (val=.,
                       coords), aes(x = Var1, y = Var2, fill = val)) +
    geom_tile() +
    coord_fixed() +
    labs(title = expression(Latitude[i]),
         x = "", y = "Latitude", fill = expression(X[i]))+
    scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
    scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
    #my_theme +
    theme(legend.position = "right",
          legend.text = element_text(angle=0),
          legend.key.height = unit(0.5,"cm"))
  
  ,
  
  # plot spatial random effect
  ggplot(cbind (dat$coords %>%
                  data.frame(.),
                w= dat$w[,1]),
         aes(x = X1, y = X2, fill = w)) +
    geom_tile() +
    coord_fixed() +
    labs(title = "Spatial Random effect",
         x = "", y = "", fill = expression(omega[i]))+
    scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
    scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
    #my_theme +
    theme(legend.position = "right",
          legend.text = element_text(angle=0),
          legend.key.height = unit(0.5,"cm"))
  
  , 
  
  table(apply (is.na(y_mixed[,1,])!=T,1,sum)) %>%
    data.frame () %>%
    ggplot(aes(x=Var1,y=Freq)) +
    geom_col()+
    #my_theme+
    xlab("Secondary occasions")+
    ylab("Cells")
  
  
  ,
  ncol=3)

# observation covariate and observation data
proj_Xp_Y <- grid.arrange (
  
  
  ggplot(observation_data_full %>%
           filter (Time ==1),
         aes(x = X, y = Y, fill = Y)) +
    geom_tile() +
    coord_fixed() +
    facet_wrap(~ Survey , nrow = 2) +
    labs(title = "Latitude (t=1)",
         x = "Longitude", y = "Latitude", fill = expression(Lat[i])) +
    scale_colour_viridis_c(option = "magma",
                           direction=1,na.value = "gray90") +
    scale_fill_viridis_c(option = "magma",direction=1,
                         na.value = "gray90")+
    #my_theme +
    
    theme(legend.position = "right",
          legend.text = element_text(angle=0),
          legend.key.height = unit(0.5,"cm"),
          axis.title = element_text(size=8),
          axis.text.x = element_text(size=4),
          axis.text.y = element_text(size=4)
    )
  
  ,
  
  
  ggplot(observation_data_full %>%
           filter (Time ==1),
         aes(x = X, y = Y, fill = ObsCov.value)) +
    geom_tile() +
    coord_fixed() +
    facet_wrap(~ Survey , nrow = 2) +
    labs(title = "Observation-level Covariate (t=1)",
         x = "Longitude", y = "Latitude", fill = expression(v[ij])) +
    scale_colour_viridis_c(option = "magma",
                           direction=1,na.value = "gray90") +
    scale_fill_viridis_c(option = "magma",direction=1,
                         na.value = "gray90")+
    #my_theme +
    
    theme(legend.position = "right",
          legend.text = element_text(angle=0),
          legend.key.height = unit(0.5,"cm"),
          axis.title = element_text(size=8),
          axis.text.x = element_text(size=4),
          axis.text.y = element_text(size=4)
    )
  
  ,
  
  ggplot(observation_data_full %>%
           filter (Time==1), 
         aes(x = X, y = Y, fill = factor(Observation))) +
    geom_tile() +
    scale_fill_manual(values = c("black", "violet"),na.value = "gray90", 
                      labels = c("Not detected", "Detected")) +
    facet_wrap(~ Survey , nrow = 2) +
    coord_fixed() +
    
    labs(title = "Detection / Non Detection data",
         x = "Longitude", y = "", fill = "Detection")+
    #my_theme+
    theme(legend.position = "right",
          legend.text = element_text(angle=0),
          axis.title = element_text(size=8),
          axis.text.x = element_text(size=4),
          axis.text.y = element_text(size=4)
    )
  
  ,
  
  nrow=1,ncol=3)


# arrange plot
png(here("figures", "TuningD&S_sims", "Scenario_phenology_spot.png"),
    width = 1000, height = 1000,units = "px")
  
  grid.arrange (proj_Z,
                proj_X,
                proj_Xp_Y,ncol=1)

dev.off()


