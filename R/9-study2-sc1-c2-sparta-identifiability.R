
# -------------------------------------------------------

# Interpretation study 2 scenario 1
# Using the approach of sparta (spatially uncorrelated random effects)


# SCENARIO 0 - changing from Bernoulli to a Poisson sampling
# SCENARIO 1 - SCENARIO 0 + X_it = psi in function of L_i

# ---------------------------------------------------------------

rm(list = ls())
gc()
library(spOccupancy)
require(here)
require(ggplot2)
require(dplyr)
require(reshape)
require(gridExtra)
require(tidyr)

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


# create a folder to host the figures
dir.create(here ("figures", "sims_present_paper","Scenario1_sparta"))

# load model output
list_output <- list.files (here ("model_output", "output_simulations", "scenario_two_sparta"),
                           pattern = "sim-mixed-stPGOcc-results")

# Load the data sets
load(file = here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))

# load true psi 
load(file = here ("model_output", "output_simulations", "scenario_one",
                  "sim-mixed-stPGOcc-psi-true_1600.rda"))

# ordering the output
order_output <- as.numeric(
  gsub (".rda", "",
        gsub ("sim-mixed-stPGOcc-results_","", list_output)
  ))
list_output <- list_output[order(order_output)]  # reorder

# load and extract estimated psi_it
hat_psi_it <- lapply (list_output, function (i) {
  # load 
  load(file = here ("model_output", 
                    "output_simulations", 
                    "scenario_two_sparta",
                  i))
  # extract
  hat_psi_it <- output_save$mean$muZ
  
  # rm output_save
  rm(output_save)
  ;
  # return
  hat_psi_it
  
  }
)

# plot one single simulation run + the truth
plot.df <- lapply (seq (1:42), function (sim) 
  do.call(rbind, lapply (seq (1:16), function (sce) {
    
    curr.indx <- (sim - 1) * n.scenarios + sce
    plot.df <- cbind (Estimated = melt(hat_psi_it[[curr.indx]])[,3],
                      melt(psi.true[,,sim,sce]),
                      sim=sim,
                      scenario=sce) 
    
    plot.df
})))
plot.df <- do.call(rbind,plot.df)

# list to array
#arr_hat_psi_it<-array(unlist(hat_psi_it),dim=c(1200,10,42,16)) # nsites x nyears x nsims x nscenarios
#arr_hat_psi_it <- lapply (hat_psi_it,melt)
#arr_hat_psi_it <- do.call(rbind , arr_hat_psi_it)

# Generate the data -------------------------------------------------------
# Number of data sets for each scenario
n.sims <- 42

# Spatial locations
J.x <- 30
J.y <- 40
J <- J.x * J.y

# Number of years
n.time <- 10

# Five replicates
n.rep <- matrix(5, J, n.time)

# Occurrence coefficient -----------------
# Generate a single covariate
beta <- c(0, 0.5)
p.occ <- length(beta)

# Detection coefficient ---------------
# A single covariate on detection
alpha <- c(0, -0.5)
# Spatial parameters ------------------
sp <- TRUE
# Assume an exponential correlation model
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

# ----------------------------------------------------

# SUMMARIZE AS DOSER & STOUDT
colnames(plot.df) <- c('Estimated','site', 'year', 'True', 'sim', 'scenario')
plot.df$diff <- plot.df$True - plot.df$Estimated # calculate difference

# remove after all runs
#plot.df <- plot.df %>%
#  filter(is.na(val) !=T)

# data frame with average
avg.df <- plot.df %>%
  group_by(scenario, sim) %>%
  arrange(True, .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(val.avg = mean(True),
            est.avg = mean(Estimated)) %>%
  ungroup()

# Summary plot for single-visit logit
avg.df$scenario <- as.numeric(avg.df$scenario)
avg.df$spatial <- ifelse(avg.df$scenario %in% c(1, 5, 9, 13), 'A', 
                         ifelse(avg.df$scenario %in% c(2, 6, 10, 14), 'B', 
                                ifelse(avg.df$scenario %in% c(3, 7, 11, 15), 'C', 'D')))
avg.df$time <- ifelse(avg.df$scenario %in% c(1, 2, 3, 4), 'A', 
                      ifelse(avg.df$scenario %in% c(5, 6, 7, 8), 'B', 
                             ifelse(avg.df$scenario %in% c(9, 10, 11, 12), 'C', 'D')))
spatial.labs <- c(expression(paste(sigma, " "^2, " Low, ", phi, "  High")), 
                  expression(paste(sigma, " "^2, " High, ", phi, "  High")),
                  expression(paste(sigma, " "^2, " Low, ", phi, "  Low")),
                  expression(paste(sigma, " "^2, " High, ", phi, "  Low")))
avg.df$spatial <- factor(avg.df$spatial, levels = c('A', 'B', 'C', 'D'), 
                         labels = spatial.labs)
time.labs <- c(expression(paste(sigma, " "[T]^2, " Low, ", rho, "  Low")), 
               expression(paste(sigma, " "[T]^2, " Low, ", rho, "  High")), 
               expression(paste(sigma, " "[T]^2, " High, ", rho, "  Low")), 
               expression(paste(sigma, " "[T]^2, " High, ", rho, "  High"))) 
avg.df$time <- factor(avg.df$time, levels = c('A', 'B', 'C', 'D'), 
                      labels = time.labs)
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

# summary plot
# Figure 3 of associated manuscript. 
fig.3.plot <- plot.df %>%
  ggplot() +
  theme_light(base_size = 16) +
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  #geom_point(aes(x = val, y = est, group = sim),alpha=0.01,col="gray90")+
  geom_smooth(aes(x = True, y = Estimated, group = sim), col = 'gray', alpha = 0.1, se = FALSE, 
              lineend = 'round', lwd = 0.25) +
  geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), se = FALSE, col = 'black', lineend = 'round', 
              lwd = 0.5) + 
  labs(x = 'True Occupancy Probability', y = 'Estimated Occupancy Probability') + 
  my_theme

# save
png(here("figures","sims_present_paper", "Scenario1_sparta", "Figure-3-st2-sc1-sparta.png"),
    width = 600, height = 600,units = "px")

  fig.3.plot

dev.off()


# -----------------------------------------------------
# map of true and estimated psi_it

# Reshape for plotting  ---- do for scenario 14, high correlations
occupancy_data<- data.frame(
  X = dat[[1]]$coords[,1],
  Y = dat[[1]]$coords[,2],
  True = psi.true[,1,2,14],
  Estimated = (hat_psi_it[[30]])[,1],
  Study="2",
  sc="1"
) %>% 
  mutate (Difference =  Estimated - True)

# Plot the true occupancy data for the first year
map_true <- ggplot(occupancy_data %>%
                                  reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                                  filter (variable == "True"), 
                                aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma",limits=c(0,1)) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Truth",
    x = "", y = "Latitude", fill = "Occupancy")+
  my_theme+
  theme(legend.position = "right",
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.text = element_text(angle=0,size=6),
        legend.title = element_text(size=7),
        strip.text.y = element_text(size=6),
        strip.text.x = element_text(size=6),
        legend.key.height = unit(0.5,"cm"),
        strip.background=element_rect(colour="black",
                                      fill="gray90"))

# estimated 
map_estimated <- ggplot(occupancy_data %>%
                     reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                     filter (variable == "Estimated"), 
                   aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma",limits=c(0,1)) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Estimated",
       x = "", y = "Latitude", fill = "Occupancy")+
  my_theme+
  theme(legend.position = "right",
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.text = element_text(angle=0,size=6),
        legend.title = element_text(size=7),
        strip.text.y = element_text(size=6),
        strip.text.x = element_text(size=6),
        legend.key.height = unit(0.5,"cm"),
        strip.background=element_rect(colour="black",
                                      fill="gray90"))

# Plot the true occupancy data for the first year
map_difference <- ggplot(occupancy_data%>%
                           reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                           filter (variable == "Difference"), 
                         aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(limits=c(-1,1))+
  coord_fixed() +
  theme_minimal() +
  labs(title = "Difference",
    x = "Longitude", y = "", fill = "Difference")+
  my_theme+
  theme(legend.position = "right",
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.text = element_text(angle=0,size=6),
        legend.title = element_text(size=7),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=6),
        legend.key.height = unit(0.5,"cm"),
        strip.background=element_rect(colour="black",
                                      fill="gray90"))

# do also for scenario 3 - low correlations overall ------------------
# Reshape for plotting  
occupancy_data<- data.frame(
  X = dat[[1]]$coords[,1],
  Y = dat[[1]]$coords[,2],
  True = psi.true[,1,2,3],
  Estimated = (hat_psi_it[[19]])[,1],
  Study="2",
  sc="1"
) %>% 
  mutate (Difference = Estimated - True )

# Plot the true occupancy data for the first year
map_true_l <- ggplot(occupancy_data %>%
                     reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                     filter (variable == "True"), 
                   aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma",limits=c(0,1)) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Truth",
       x = "", y = "Latitude", fill = "Occupancy")+
  my_theme+
  theme(legend.position = "right",
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.text = element_text(angle=0,size=6),
        legend.title = element_text(size=7),
        strip.text.y = element_text(size=6),
        strip.text.x = element_text(size=6),
        legend.key.height = unit(0.5,"cm"),
        strip.background=element_rect(colour="black",
                                      fill="gray90"))

# estimated 
map_estimated_l <- ggplot(occupancy_data %>%
                          reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                          filter (variable == "Estimated"), 
                        aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma",limits=c(0,1)) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Estimated",
       x = "", y = "Latitude", fill = "Occupancy")+
  my_theme+
  theme(legend.position = "right",
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.text = element_text(angle=0,size=6),
        legend.title = element_text(size=7),
        strip.text.y = element_text(size=6),
        strip.text.x = element_text(size=6),
        legend.key.height = unit(0.5,"cm"),
        strip.background=element_rect(colour="black",
                                      fill="gray90"))


# Plot the true occupancy data for the first year
map_difference_l <- ggplot(occupancy_data%>%
                           reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                           filter (variable == "Difference"), 
                         aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(limits=c(-1,1))+
  coord_fixed() +
  theme_minimal() +
  labs(title = "Difference",
       x = "Longitude", y = "", fill = "Difference")+
  my_theme+
  theme(legend.position = "right",
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        legend.text = element_text(angle=0,size=6),
        legend.title = element_text(size=7),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=6),
        legend.key.height = unit(0.5,"cm"),
        strip.background=element_rect(colour="black",
                                      fill="gray90"))
# arrange plots
p1<-grid.arrange(map_true_l, 
                 map_estimated_l,
                 map_difference_l,
                 nrow=3)
p2<- grid.arrange(map_true, 
             map_estimated,
             map_difference,
             nrow=3)

# save all them
png (here("figures","sims_present_paper" ,"Scenario1_sparta","maps_comparison_sparta.png"),
     width = 13, height = 17, units = "cm",res=150)

  grid.arrange (
    p1
    ,
    p2,
    ncol=2

  )

dev.off()


# evaluating regression coefficients -------------------------------------------
# load and extract estimated betas
hat_betas <- lapply (list_output, function (i) {
  # load 
  load(file = here ("model_output", 
                    "output_simulations", 
                    "scenario_two_sparta",
                    i))
  # extract
  hat_betas <- output_save$mean[grep('beta',names(output_save$mean))]
  # rm output_save
  rm(output_save)
  ;
  # return
  hat_betas
  
}
)


# load  and extract estimated alphas -------------------------------------------
hat_alphas <- lapply (list_output, function (i) {
  # load 
  load(file = here ("model_output", 
                    "output_simulations", 
                    "scenario_two_sparta",
                    i))
  # extract
  hat_alphas <- output_save$mean[grep('alpha',names(output_save$mean))]
  # rm output_save
  rm(output_save)
  ;
  # return
  hat_alphas
  
}
)

# save
png(here("figures","sims_present_paper", "Scenario1_sparta", "coeff_densities.png"),
    width = 500, height = 500,units = "px")

# plot (densities)
grid.arrange(
  
  # occupancy
  
  do.call(rbind, lapply (hat_betas,unlist)) %>%
    ggplot(aes(x=beta0))+
    geom_density(fill="white")+
    geom_vline(xintercept = beta[1])+
    labs(x = bquote(""*hat(beta)[0]*""))+
    my_theme
  
  ,
  do.call(rbind, lapply (hat_betas,unlist)) %>%
    ggplot(aes(x=beta1))+
    geom_density(fill="white")+
    geom_vline(xintercept = beta[2])+
    labs(x = bquote(""*hat(beta)[1]*""))+
    my_theme+
    theme(axis.title.y = element_blank())

  , # detection
  
  do.call(rbind, lapply (hat_alphas,unlist)) %>%
    ggplot(aes(x=alpha0))+
    geom_density(fill="white")+
    geom_vline(xintercept = alpha[1])+
    labs(x = bquote(""*hat(alpha)[0]*""))+
    my_theme
  ,
  do.call(rbind, lapply (hat_alphas,unlist)) %>%
    ggplot(aes(x=alpha1))+
    geom_density(fill="white")+
    geom_vline(xintercept = alpha[2])+
    labs(x = bquote(""*hat(alpha)[1]*""))+
    my_theme+
    theme(axis.title.y = element_blank())
  ,
  ncol=2,nrow=2
  
  )

dev.off()


