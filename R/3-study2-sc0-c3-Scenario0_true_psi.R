
# -------------------------------------------------------

# Simulation study 2

# ---------------------------------------------------------------

# SCENARIO 0 - changing from Bernoulli to a Poisson sampling
# using the sampling scheme created in the code '3-study2-sc0-c1-Poisson-Sampling-Design.R'
# source code as reference: https://rdrr.io/cran/spOccupancy/src/R/simTOcc.R


# ---------------------------------------------------------------

rm(list = ls())
library(spOccupancy)
library(here)

# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# load sampling design
load (file=here ("model_output", 
                 "output_simulations", 
                 "sampling_design_Poisson.rda")
)

# intermediary objects
# psi
psi.true <- array(NA, dim = c(I, n.time, n.sims, n.scenarios))

# Load the data sets
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# simulate a new X.p with more secondary occasions ---------------------
# from the simTocc source code
p.det <- length(alpha)
X.p <- array(NA, dim = c(I, n.time, J, p.det))
X.p[, , , 1] <- 1
if (p.det > 1) {
  for (i in 1:I) {
    for (t in 1:n.time) {
      for (j in 1:J) {
        X.p[i, t, j, 2:p.det] <- rnorm(p.det - 1)
      } # j
    } # t
  } # i  
}

# save 
save (X.p, file=here ("model_output", "output_simulations", "Xp_itj.rda"))

# Run simulations ---------------------------------------------------------
for (s in 1:n.sims) {
  print(paste("Currently on simulation set ", s, " out of ", n.sims, sep = ''))
  for (sc in 1:n.scenarios) {
    print(paste("Currently on scenario ", sc, " out of ", n.scenarios, sep = ''))out.rwalk.RE
    set.seed(my.seeds[s])
    
    # index
    curr.indx <- (s - 1) * n.scenarios + sc
    dat <- dat.full[[curr.indx]]
    phi.tune <- 0.5
    
    # sim data for spOccupancy -------------------------------------------
    # Linear predictor for occupancy (psi)
    logit_psi <- array(0, dim = c(I, n.time))
    
    # simulate occupancy         
    for (t in 1:n.time) {
      for (i in 1:I) {
        logit_psi[i, t] <- beta %*% dat.full[[curr.indx]]$X[i, t, ]  + 
          dat.full[[curr.indx]]$w[i,] + 
          dat.full[[curr.indx]]$eta[t,]
      }
    }
    
    # Logistic function to get occupancy probabilities
    psi <- 1 / (1 + exp(-logit_psi))
    psi.true[, , s, sc] <- psi
    
    # save
    save(psi.true, file = here ("model_output", 
                                "output_simulations", 
                                "scenario_zero",
                                paste0("sim-mixed-stPGOcc-psi-true_",curr.indx,".rda"))
    )
    
    # the loop makes the array each time larger, so
    # remove files to avoid reach the upper memory bound
    # remove three output before the present output
    unlink (here ("model_output", 
                  "output_simulations", 
                  "scenario_zero",
                  paste0("sim-mixed-stPGOcc-psi-true_",curr.indx-3,".rda")),recursive=F)
    
    
  } # i (n.scenarios)
} # j (n.sims)

# end
