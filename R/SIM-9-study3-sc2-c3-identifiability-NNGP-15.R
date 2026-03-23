
# -------------------------------------------------------

# Simulation study 3 - scenario 2 with phenology & observer preferences for midseason & observation spot effect

# manuscript revision NGPP = 15


# -------------------------------------------------------

rm(list = ls())
library(spOccupancy)
require(here)
require(ggplot2)
require(dplyr)
require(reshape)
require(gridExtra)
require(tidyr)

# load D&S model output -------------------------------------
#load(file = here ("model_output", "output_simulations", "sims_D&S",
#                  "sim-mixed-stPGOcc-results-merged.rda"))

# keep only psi.true
#rm(list=setdiff(ls(), c("psi.true")))

# create a folder to host the figures
dir.create(here ("figures", "sims_present_paper","ScenarioSpot-NGPP-15"))

# Load the data sets
load(file = here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))

# load model output
load(file = here ("model_output", "output_simulations", "scenario_phenology_spot-NNGP-15-review",
                  "sim-mixed-stPGOcc-results_320.rda"))


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

# Generate the data -------------------------------------------------------
# Number of data sets for each scenario
n.sims <- 100

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
plot.df <- as.data.frame.table(psi.true)
colnames(plot.df) <- c('site', 'year', 'sim', 'scenario', 'val')
psi.mean.df <- as.data.frame.table(psi.mean.samples)
colnames(psi.mean.df) <- c('site', 'year', 'sim', 'scenario', 'val')
plot.df$est <- psi.mean.df$val
plot.df$diff <- plot.df$val - plot.df$est

# remove after all runs
plot.df <- plot.df %>%
  filter(is.na(val) !=T)

# data frame with average
avg.df <- plot.df %>%
 group_by(scenario, sim) %>%
 arrange(val, .by_group = TRUE) %>%
 mutate(fake.id = 1:n()) %>%
 ungroup() %>%
 group_by(scenario, fake.id) %>%
 summarize(val.avg = mean(val),
	    est.avg = mean(est)) %>%
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
  geom_smooth(aes(x = val, y = est, group = sim), col = 'gray', alpha = 0.1, se = FALSE, 
 	      lineend = 'round', lwd = 0.25) +
   geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2) +
   geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), se = FALSE, col = 'black', lineend = 'round', 
 	      lwd = 0.5) + 
   labs(x = 'True Occupancy Probability', y = 'Estimated Occupancy Probability') + 
   my_theme

# save
png(here("figures","sims_present_paper", "ScenarioSpot-NGPP-15", "Figure-3-st3-sc2.png"),
    width = 20, height = 20, units = "cm", res=200)

  fig.3.plot

dev.off()


# plot of points (explore why the most overestimate rather than underestimate psi_it)
# Figure 3 of associated manuscript. 
plot.df %>%
  filter (sim == "A") %>%
  ggplot() +
  theme_light(base_size = 16) +
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  geom_point(aes(x = val, y = est, group = sim, col = as.numeric(site)),alpha=0.01)+
  scale_colour_gradient2(low="blue", mid= "gray", high="red",midpoint=600)+
  #geom_smooth(aes(x = val, y = est, group = sim), col = 'gray', alpha = 0.1, se = FALSE, 
  #            lineend = 'round', lwd = 0.25) +
  geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2) +
  geom_smooth(data = avg.df, aes(x = val.avg, y = est.avg), se = FALSE, col = 'black', lineend = 'round', 
              lwd = 0.5) + 
  labs(x = 'True Occupancy Probability', y = 'Estimated Occupancy Probability') + 
  my_theme


# data to plot spatial params ----------------------------------------
dat.plot <- melt(theta.mean.samples)%>% 
  dplyr::rename(#"nbatch"="X1",
                "param"="X2" ,
                "sims"= "X1",
                "scenario" = "X3",
                "val.theta"="value") %>%
  filter(is.na(val.theta) !=T) %>%
  right_join(.,scenario.vals,
             by="scenario") %>%
  data.frame()

# calculate the averages to plot
av.dat.plot <-  dat.plot%>%
  group_by (param , scenario) %>%
  summarise(mean_par = mean(val.theta))

# match
dat.plot <- cbind (dat.plot,
                   av.dat.plot[match(
                     paste (dat.plot$param, av.dat.plot$scenario),
                     paste (av.dat.plot$param, av.dat.plot$scenario)
                   ),"mean_par"])

# Ensure time and spatial columns are character vectors
dat.plot$spatial <- as.character(dat.plot$spatial)
dat.plot$time <- as.character(dat.plot$time)

# Convert time and spatial columns to factors for proper faceting
dat.plot$spatial <- factor(dat.plot$spatial, levels = unique(dat.plot$spatial))
dat.plot$time <- factor(dat.plot$time, levels = unique(dat.plot$time))

# parameter
dat.plot$paramLab <- as.factor(dat.plot$param)
levels ( dat.plot$paramLab ) <- c(expression(paste(sigma, " "^2)),
                                  expression(paste(phi)),
                                  expression(paste(sigma, " "[T]^2)),
                                  expression(paste(rho))
                                   
                               )



# Select the paramLab for the current param (assuming i is defined and valid)
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 1])

png(here("figures","sims_present_paper", "ScenarioSpot-NGPP-15", "Density_sigma2-st3-sc2.png"),
    width = 20, height = 20, units = "cm", res=200)

# plot
# spatial sigma
dat.plot %>%
  filter (param == 1) %>%
  ggplot()+
  geom_density(aes(x=val.theta,
                     group=scenario),fill="white")+
  geom_vline(aes(xintercept=sigma.sq))+
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  labs(x=bquote(""*hat(sigma^2)*""))

dev.off()


# spatial phi
# Select the paramLab for the current param (assuming i is defined and valid)
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 2])

png(here("figures","sims_present_paper", "ScenarioSpot-NGPP-15", "Density_phi-st3-sc2.png"),
    width = 20, height = 20, units = "cm", res=200)

dat.plot %>%
  filter (param == 2) %>%
  ggplot()+
  geom_density(aes(x=val.theta,
                   group=scenario),fill="white")+
  geom_vline(aes(xintercept=phi))+
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  labs(x=bquote(""*hat(phi)*""))


dev.off()

# temporal rho
# Select the paramLab for the current param (assuming i is defined and valid)
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 3])

png(here("figures","sims_present_paper", "ScenarioSpot-NGPP-15", "Density_sigma2_t-st3-sc2.png"),
    width = 20, height = 20, units = "cm", res=200)

dat.plot %>%
  filter (param == 3) %>%
  ggplot()+
  geom_density(aes(x=val.theta,
                   group=scenario),fill="white")+
  geom_vline(aes(xintercept=sigma.sq.t))+
  xlim(c(0,7))+
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  labs(x=bquote(""*hat(sigma[T]^2)*""))

dev.off()

# temporal sigma
# Select the paramLab for the current param (assuming i is defined and valid)
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 4])

png(here("figures","sims_present_paper", "ScenarioSpot-NGPP-15", "Density_rho-st3-sc2.png"),
    width = 20, height = 20, units = "cm", res=200)

dat.plot %>%
  filter (param == 4) %>%
  
  ggplot()+
  geom_density(aes(x=val.theta,
                   group=scenario),fill="white")+
  geom_vline(aes(xintercept=rho))+
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  labs(x=bquote(""*hat(rho)*""))

dev.off()

# end