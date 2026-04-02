
# -------------------------------------------------------

#  COMPARISON OF STUDIES, SCENARIOS AND SUB SCENARIOS
#  
#  Code to produce most figures used in the main text

# TODO --- check if it is better to use the product or the realization of this product Bernoulli (psi * p)


# -------------------------------------------------------

rm(list = ls())
gc()
require(spOccupancy) # package for data analyses
require(here) # transit between folders
require(ggplot2) # plots
require(dplyr)# organize data
require(tidyr)# organize data
require(reshape) # organize data
require(gridExtra) # organize plots
require(ggbreak) # organize plots
require(abind) # construct arrays

# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# load model output
sceDS <- new.env()
sceDSsmooth <- new.env()
sce0 <- new.env()
sce1 <- new.env()
sce2 <- new.env()
sce3 <- new.env()
scePhen <- new.env()
scePhenSpot <- new.env()
scePhenSpot2 <- new.env()

# Scenario 1 - 0 (baseline scenario) ------------------------------------------

# load in the environments
load(file = here ("model_output", "output_simulations", "sims_D&S",
                  "sim-mixed-stPGOcc-results-merged-B.rda"),sceDS)
sceDS$study <- 1
sceDS$sc <- 0

# Relating the true psi_it x p_it with estimated psi_it x p_it
# Load the data sets
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# Organize labels for plots ------------------------------
# spatial params
scenario.vals <- scenario.vals %>%
  mutate (scenario = seq(1,nrow(scenario.vals)))

# spatial labels
scenario.vals$spatial <- rep(c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                               expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                               expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                               expression(paste(sigma, " "^2, " High, ", phi, "  Low"))),
                             4)
# temporal labels
scenario.vals$time <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                        expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                        expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                        expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  High")),
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  High")),
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  High")),
                        expression(paste(sigma, " "[T]^2, " High, ", rho, "  High")))

# choose the simulation run and scenario to evaluate
s = 1
sc = 13

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
        
    # index -- this will be used to extract data from the simulated data
    curr.indx <- (s - 1) * n.scenarios + sc
    
    # true
    prod_occ_det <- (dat.full[[curr.indx]]$psi) * apply(dat.full[[curr.indx]]$p[ , , 1:2], c(1,2), mean) # first two survey occasions
    
    # estimated
    
    # Linear predictor for detection (p)
    logit_p <- array(0, dim = c(I, n.time, 2)) # 2 survey occasions
    
    # recover detection using the detection covariates and estimated coefficients    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:2) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p[i, t, j] <- dat.full[[curr.indx]]$X.p[i, t, j, ] %*% sceDS$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_true <- 1 / (1 + exp(-logit_p))
    
    # estimated product
    est_occ_det <-  sceDS$psi.mean.samples [ , , s, sc] * apply (p_true, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (true = as.vector (prod_occ_det),
           estimated = as.vector (est_occ_det),
           sim = s,
           scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
    
  }
  ) # close scenario
  }
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)

# save
png(here("figures","sims_D&S","Figure-psi-p-st1-sc0.png"),
    width = 600, height = 600,units = "px")
  
    fig.psi.p.prod

dev.off()
  
# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(sceDS)
rm(fig.psi.p.prod)
gc()

# Scenario 1 - 1 ------------------------------------------
load(file = here ("model_output", "output_simulations", "smooth_sims_D&S",
                  "sim-mixed-stPGOcc-results-merged-B.rda"),sceDSsmooth)
sceDSsmooth$study <- 1
sceDSsmooth$sc <- 1

# Load the simulated data sets - higher autocorrelation
load(here ("model_output", "output_simulations", "smooth_sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # index -- this will be used to extract data from the simulated data
    curr.indx <- (s - 1) * n.scenarios + sc
    
    # true
    prod_occ_det <- sceDSsmooth$psi.true [ , , s, sc] * apply(dat.full[[curr.indx]]$p[ , , 1:2], c(1,2), mean) # first two survey occasions
    
    # estimated
    # Linear predictor for detection (p)
    logit_p <- array(0, dim = c(I, n.time, 2)) # 2 survey occasions
    
    # recover detection using the detection covariates and estimated coefficients    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:2) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p[i, t, j] <- dat.full[[curr.indx]]$X.p[i, t, j, ] %*% sceDSsmooth$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_true <- 1 / (1 + exp(-logit_p))
    
    # estimated product
    est_occ_det <-  sceDSsmooth$psi.mean.samples [ , , s, sc] * apply (p_true, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (#true.psi = 
                                true = as.vector (prod_occ_det),
                                estimated = as.vector (est_occ_det),
                                sim = s,
                                scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
    
  }
  ) # close scenario
}
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","sims_D&S_smooth","Figure-psi-p-st1-sc1.png"),
    width = 600, height = 600,units = "px")

  fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(sceDSsmooth)
rm(fig.psi.p.prod)
gc()


# Scenario 2 - 0 ------------------------------------------
# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# load covariate
load (file=here ("model_output", "output_simulations", "Xp_itj.rda"))

load(file = here ("model_output", "output_simulations", "scenario_zero",
                  "sim-mixed-stPGOcc-results-merged-B.rda"),sce0)
sce0$study <- 2
sce0$sc <- 0

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # Recover true detection probability
    logit_p_true <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_true[i, t, j] <- X.p[i, t, j, ] %*% alpha
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_true <- 1 / (1 + exp(-logit_p_true))
    
    # true product
    prod_occ_det <- sce0$psi.true [ , , s, sc] * apply(p_true, c(1,2), mean) # first two survey occasions
    
    # estimated detection
    logit_p_hat <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_hat[i, t, j] <- X.p[i, t, j, ] %*% sce0$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_hat <- 1 / (1 + exp(-logit_p_hat))
    
    # estimated product
    est_occ_det <-  sce0$psi.mean.samples [ , , s, sc] * apply (p_hat, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (
      true = as.vector (prod_occ_det),
      estimated = as.vector (est_occ_det),
      sim = s,
      scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
    
  }
  ) # close scenario
}
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","TuningD&S_sims","Scenario0","Figure-psi-p-st2-sc0.png"),
    width = 600, height = 600,units = "px")

  fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(avg.df)
rm(sce0)
rm(fig.psi.p.prod)
gc()

# Scenario 2 - 1 ------------------------------------------
# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

load(file = here ("model_output", "output_simulations", "scenario_one",
                  "sim-mixed-stPGOcc-results-merged-B.rda"),sce1)

sce1$study <- 2
sce1$sc <- 1

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # Recover true detection probability
    logit_p_true <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_true[i, t, j] <- X.p[i, t, j, ] %*% alpha
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_true <- 1 / (1 + exp(-logit_p_true))
    
    # true product
    prod_occ_det <- sce1$psi.true [ , , s, sc] * apply(p_true, c(1,2), mean) # first two survey occasions
    
    # estimated detection
    logit_p_hat <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_hat[i, t, j] <- X.p[i, t, j, ] %*% sce1$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_hat <- 1 / (1 + exp(-logit_p_hat))
    
    # estimated product
    est_occ_det <-  sce1$psi.mean.samples [ , , s, sc] * apply (p_hat, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (
      true = as.vector (prod_occ_det),
      estimated = as.vector (est_occ_det),
      sim = s,
      scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
    
    }
    ) # close scenario
  }
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","TuningD&S_sims","Scenario1","Figure-psi-p-st2-sc1.png"),
    width = 600, height = 600,units = "px")

  fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(sce1)
rm(fig.psi.p.prod)
gc()

# Scenario 2 - 2 ------------------------------------------
load(file = here ("model_output", "output_simulations", "scenario_two",
                  "sim-mixed-stPGOcc-results-merged-B.rda"),sce2)
# load true psi 
#load(file = here ("model_output", "output_simulations", "scenario_two",
#                  "sim-mixed-stPGOcc-psi-true_1600.rda"),sce2)
sce2$study <- 2
sce2$sc <- 2

# Load the data sets to have latitude
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # index -- this will be used to extract data from the simulated data
    curr.indx <- (s - 1) * n.scenarios + sc
    
    # L_i will replace X_it ----------------------
    L_i <- scale(dat.full[[curr.indx]]$coords[, 2])[,1]
    
    # create design matrix for covariates
    X <- array(1, dim = c(I, n.time, length(beta)))
    X[, , 2] <- L_i # set latitude here
    
    # bind latitude
    lat_bind <- replicate (J,X[, , 2])
    L_i_det <- abind(X.p,lat_bind)
    
    # remove V_itj \ X.p
    L_i_det <- L_i_det [,,,-2]
    
    # Linear predictor for detection (p)
    logit_p <- array(0, dim = c(I, n.time, J))
    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p[i, t, j] <- L_i_det [i, t, j, ] %*% alpha
          
        }
      }
    }
    
    # Logistic function to get detection probabilities
    p_true <- 1 / (1 + exp(-logit_p))
    
    # true product
    prod_occ_det <- sce2$psi.true [ , , s, sc] * apply(p_true, c(1,2), mean) # first two survey occasions
    
    # estimated detection
    logit_p_hat <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and the estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_hat[i, t, j] <- L_i_det [i, t, j, ] %*% sce2$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_hat <- 1 / (1 + exp(-logit_p_hat))
    
    # estimated product
    est_occ_det <-  sce2$psi.mean.samples [ , , s, sc] * apply (p_hat, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (
      true = as.vector (prod_occ_det),
      estimated = as.vector (est_occ_det),
      sim = s,
      scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
  }
  ) # close scenario
}
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data [,-9]%>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","TuningD&S_sims","Scenario2","Figure-psi-p-st2-sc2.png"),
    width = 600, height = 600,units = "px")

  fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(sce2)
rm(fig.psi.p.prod)
gc()

# Scenario 2 - 3 ------------------------------------------
load(file = here ("model_output", "output_simulations", "scenario_three",
                  "sim-mixed-stPGOcc-results_1600.rda"),sce3)
sce3$study <- 2
sce3$sc <- 3

# Load the data sets to have latitude
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# Detection coefficient ---------------
# Now alpha has three coefficients
alpha <- c(0, -0.5, -0.5) # added one coefficient for latitude

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # index -- this will be used to extract data from the simulated data
    curr.indx <- (s - 1) * n.scenarios + sc
    
    # L_i will replace X_it ----------------------
    L_i <- scale(dat.full[[curr.indx]]$coords[, 2])[,1]
    
    # create design matrix for covariates
    X <- array(1, dim = c(I, n.time, length(beta)))
    X[, , 2] <- L_i # set latitude here
    
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
    p_true <- 1 / (1 + exp(-logit_p))
    
    # true product
    prod_occ_det <- sce3$psi.true [ , , s, sc] * apply(p_true, c(1,2), mean) # first two survey occasions
    
    # estimated detection
    logit_p_hat <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and the estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_hat[i, t, j] <- V_itj_L_i [i, t, j, ] %*% sce3$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_hat <- 1 / (1 + exp(-logit_p_hat))
    
    # estimated product
    est_occ_det <-  sce3$psi.mean.samples [ , , s, sc] * apply (p_hat, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (
      true = as.vector (prod_occ_det),
      estimated = as.vector (est_occ_det),
      sim = s,
      scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
  }
  ) # close scenario
}
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","TuningD&S_sims","Scenario3","Figure-psi-p-st2-sc3.png"),
    width = 600, height = 600,units = "px")

fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(sce3)
rm(fig.psi.p.prod)
gc()


# Scenario 3 - 1 ------------------------------------------
load(file = here ("model_output", "output_simulations", "scenario_phenology",
                  "sim-mixed-stPGOcc-results_1600.rda"),scePhen)
# load true psi 
#load(file = here ("model_output", "output_simulations", "scenario_phenology",
#                  "sim-mixed-stPGOcc-psi-true_1600.rda"),scePhen)

scePhen$study <- 3
scePhen$sc <- 1

# Load the data sets to have latitude
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# Detection coefficient ---------------
# Now alpha has three coefficients
alpha <- c(0, -0.5, -0.5) # added one coefficient for latitude

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # index -- this will be used to extract data from the simulated data
    curr.indx <- (s - 1) * n.scenarios + sc
    
    # L_i will replace X_it ----------------------
    L_i <- scale(dat.full[[curr.indx]]$coords[, 2])[,1]
    
    # create design matrix for covariates
    X <- array(1, dim = c(I, n.time, length(beta)))
    X[, , 2] <- L_i # set latitude here
    
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
    p_true <- 1 / (1 + exp(-logit_p))
    
    # true product
    prod_occ_det <- scePhen$psi.true [ , , s, sc] * apply(p_true, c(1,2), mean) # first two survey occasions
    
    # estimated detection
    logit_p_hat <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and the estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_hat[i, t, j] <- V_itj_L_i [i, t, j, ] %*% scePhen$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_hat <- 1 / (1 + exp(-logit_p_hat))
    
    # estimated product
    est_occ_det <-  scePhen$psi.mean.samples [ , , s, sc] * apply (p_hat, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (
      true = as.vector (prod_occ_det),
      estimated = as.vector (est_occ_det),
      #true = unlist(lapply (as.vector (prod_occ_det), rbinom, n = 1, size = 1)),
      #estimated = unlist(lapply (as.vector (est_occ_det), rbinom, n = 1, size = 1)),
      
      sim = s,
      scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
   
   # c(sum(plot.df$true), sum(plot.df$estimated))
    plot((plot.df$true), (plot.df$estimated),col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
    }
   ) # close scenario
  }
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","TuningD&S_sims","ScenarioPhenology","Figure-psi-p-st3-sc1.png"),
    width = 700, height = 600,units = "px")

  fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(scePhen)
rm(fig.psi.p.prod)
gc()


# Scenario 3 - 2 ------------------------------------------

load(file = here ("model_output", "output_simulations", "scenario_phenology_spot",
                  "sim-mixed-stPGOcc-results_1600.rda"),scePhenSpot)
# load true psi 
#load(file = here ("model_output", "output_simulations", "scenario_phenology_spot",
#                  "sim-mixed-stPGOcc-psi-true_1600.rda"),scePhenSpot)

scePhenSpot$study <- 3
scePhenSpot$sc <- 2

# Load the data sets to have latitude
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# Detection coefficient ---------------
# Now alpha has three coefficients
alpha <- c(0, -0.5, -0.5) # added one coefficient for latitude

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # index -- this will be used to extract data from the simulated data
    curr.indx <- (s - 1) * n.scenarios + sc
    
    # L_i will replace X_it ----------------------
    L_i <- scale(dat.full[[curr.indx]]$coords[, 2])[,1]
    
    # create design matrix for covariates
    X <- array(1, dim = c(I, n.time, length(beta)))
    X[, , 2] <- L_i # set latitude here
    
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
    p_true <- 1 / (1 + exp(-logit_p))
    
    # true product
    prod_occ_det <- scePhenSpot$psi.true [ , , s, sc] * apply(p_true, c(1,2), mean) # first two survey occasions
    
    # estimated detection
    logit_p_hat <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and the estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_hat[i, t, j] <- V_itj_L_i [i, t, j, ] %*% scePhenSpot$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_hat <- 1 / (1 + exp(-logit_p_hat))
    
    # estimated product
    est_occ_det <-  scePhenSpot$psi.mean.samples [ , , s, sc] * apply (p_hat, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (
      true = as.vector (prod_occ_det),
      estimated = as.vector (est_occ_det),
      sim = s,
      scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
  }
  ) # close scenario
}
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","TuningD&S_sims","ScenarioSpot","Figure-psi-p-st3-sc2.png"),
    width = 700, height = 600,units = "px")

  fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(scePhenSpot)
rm(fig.psi.p.prod)
gc()

# Scenario 3 - 2 ------------------------------------------

load(file = here ("model_output", "output_simulations", "scenario_phenology_spot_2",
                  "sim-mixed-stPGOcc-results_1600.rda"),scePhenSpot2)
scePhenSpot2$study <- 3
scePhenSpot2$sc <- 3

# Load the data sets to have latitude
load(here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
dat.full <- dat

# Detection coefficient ---------------
# Now alpha has three coefficients
alpha <- c(0, -0.5, -0.5) # added one coefficient for latitude

# calculate the product of psi_it and p_it (aggregated over J) for all sim runs and scenarios
prod_psi_det <- lapply (seq(1,n.sims), function (s) {
  lapply (seq (1, n.scenarios), function (sc) {
    
    # index -- this will be used to extract data from the simulated data
    curr.indx <- (s - 1) * n.scenarios + sc
    
    # L_i will replace X_it ----------------------
    L_i <- scale(dat.full[[curr.indx]]$coords[, 2])[,1]
    
    # create design matrix for covariates
    X <- array(1, dim = c(I, n.time, length(beta)))
    X[, , 2] <- L_i # set latitude here
    
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
    p_true <- 1 / (1 + exp(-logit_p))
    
    # true product
    prod_occ_det <- scePhenSpot2$psi.true [ , , s, sc] * apply(p_true, c(1,2), mean) # first two survey occasions
    
    # estimated detection
    logit_p_hat <- array(0, dim = c(I, n.time, J)) # J survey occasions
    
    # recover detection using the detection covariates and the estimated coefficients    
    # simulate detection    
    for (t in 1:n.time) {
      for (i in 1:I) {
        for (j in 1:J) {
          
          # 0 in the logit scale will produce p = 0.5
          logit_p_hat[i, t, j] <- V_itj_L_i [i, t, j, ] %*% scePhenSpot2$alpha.mean.samples[, s, sc]
          
        }
      }
    }
    # Logistic function to get detection probabilities
    p_hat <- 1 / (1 + exp(-logit_p_hat))
    
    # estimated product
    est_occ_det <-  scePhenSpot2$psi.mean.samples [ , , s, sc] * apply (p_hat, c(1,2), mean) 
    
    # organize data frame
    prod_psi_det <- data.frame (
      true = as.vector (prod_occ_det),
      estimated = as.vector (est_occ_det),
      sim = s,
      scenario = sc) 
    
    # right join the scenarios
    plot.df <- cbind (prod_psi_det,
                      scenario.vals [match(prod_psi_det$scenario, scenario.vals$scenario),])
    
    # edit labels to facets
    time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
                   expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
    spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                      expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                      expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                      expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
    plot.df$scenario <- as.numeric(plot.df$scenario)
    plot.df$spatial <- ifelse(plot.df$scenario %in% c(1, 5, 9, 13), 'A', 
                              ifelse(plot.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                     ifelse(plot.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
    plot.df$time <- ifelse(plot.df$scenario %in% c(1, 2, 3, 4), 'A', 
                           ifelse(plot.df$scenario %in% c(5, 6, 7, 8), 'B', 
                                  ifelse(plot.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
    plot.df$spatial <- factor(plot.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                              labels = spatial.labs)
    plot.df$time <- factor(plot.df$time, levels = c('A', 'B', 'C', 'D'), 
                           labels = time.labs)
    
    plot(plot.df$true, plot.df$estimated,col=rgb(0,0,1,alpha=0.2))
    abline(0,1)
    abline (lm (plot.df$estimated ~ plot.df$true))
    plot.df
    
  }
  ) # close scenario
}
) # close simulation

# dissolve lists into df
plot.data <- do.call( rbind, # dissolve scenarios inside simulations
                      lapply (prod_psi_det, function (s)  # dissolve simulations
                        do.call(rbind, s))  
) 
# obtain the average relationship
# data frame with average
avg.df <- plot.data[,-9] %>%
  group_by(scenario, sim) %>%
  arrange(true, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(true),
            est.avg = mean(estimated)) %>%
  ungroup()

# plot data
fig.psi.p.prod <- plot.data[,-9] %>% # remove duplicated column
  #filter (scenario == 1) %>%
  ggplot() + 
  #geom_point(aes(x=true, y=estimated), alpha = 0.2) +
  geom_abline() + 
  geom_smooth(aes(x=true, y=estimated,group=(sim)),
              col = "gray", alpha=.2,linewidth=0.2,se=F) +
  labs ( x = bquote ("True "*psi[it]*" x "*p[it]*""),
         y = bquote ("Estimated "*psi[it]*" x "*p[it]*"")) + 
  theme_light(base_size = 16) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), 
              se = FALSE, col = 'blue4', lineend = 'round', 
              lwd = 0.5) +
  facet_grid(time~spatial, 
             labeller = label_parsed)
# save
png(here("figures","TuningD&S_sims","ScenarioSpot2","Figure-psi-p-st3-sc3.png"),
    width = 700, height = 600,units = "px")

  fig.psi.p.prod

dev.off()

# free space
rm(plot.data)
rm(prod_psi_det)
rm(dat.full)
rm(dat)
rm(avg.df)
rm(scePhenSpot2)
rm(fig.psi.p.prod)
gc()

# end







