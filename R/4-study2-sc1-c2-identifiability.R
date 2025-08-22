
# -------------------------------------------------------

# Simulation study 2 - scenario 1

# Produce output/figures used in the main text and SI

# -------------------------------------------------------

rm(list = ls())
library(spOccupancy)
require(here)
require(ggplot2)
require(dplyr)
require(reshape)
require(gridExtra)
require(tidyr)

# create a folder to host the figures
dir.create(here ("figures", "sims_present_paper","Scenario1"))

# Load the data sets
load(file = here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))

# load model output
load(file = here ("model_output", "output_simulations", "scenario_one",
                  "sim-mixed-stPGOcc-results_1600.rda"))

# load true psi 
#load(file = here ("model_output", "output_simulations", "scenario_one",
#                  "sim-mixed-stPGOcc-psi-true_1600.rda"))


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
png(here("figures","sims_present_paper", "Scenario1", "Figure-3-st1-sc1.png"),
    width = 600, height = 600,units = "px")

fig.3.plot

dev.off()


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

png(here("figures","sims_present_paper", "Scenario1", "Density_sigma2-st1-sc1.png"),
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

png(here("figures","sims_present_paper", "Scenario1", "Density_phi-st1-sc1.png"),
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

png(here("figures","sims_present_paper", "Scenario1", "Density_sigma2_t-st1-sc1.png"),
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

png(here("figures","sims_present_paper", "Scenario1", "Density_rho-st1-sc1.png"),
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

# -------------------------------------------

# bias-variance trade-off
# sigma 2
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 1])
bias.sigma.sq <- dat.plot %>%
  filter (param == 1) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = ((mean(val.theta)-sigma.sq)^2),
             Variance = (val.theta-mean(val.theta))^2, 
             #VarianceFUN=(var(val.theta)),
             MSE=(val.theta-sigma.sq)^2) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = mean(Bias),
             Variance = mean(Variance), 
             #VarianceFUN=mean(VarianceFUN),
             MSE=mean(MSE)) %>%
  pivot_longer(-c(scenario,spatial,time,paramLab,param), names_to = "Metric", values_to = "Value") %>%
  # plot
  ggplot()+
  geom_histogram(aes(y=Value,
                     x=Metric,
                     group=Metric),
                 stat = "identity")+
  
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  theme(#axis.text.x = element_blank(),
        axis.ticks =  element_blank(),
        axis.title.x = element_blank())+
  ylab("Error")

# check
test <- dat.plot %>%
  filter (param == 1) %>%
  filter (scenario ==8)
mean((mean(test$val.theta) - test$sigma.sq)^2)


# phi
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 2])
bias.phi<-dat.plot %>%
  filter (param == 2) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = ((mean(val.theta)-phi)^2),
             Variance = (val.theta-mean(val.theta))^2, 
             #VarianceFUN=(var(val.theta)),
             MSE=(val.theta-phi)^2) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = mean(Bias),
             Variance = mean(Variance), 
             #VarianceFUN=mean(VarianceFUN),
             MSE=mean(MSE)) %>%
  pivot_longer(-c(scenario,spatial,time,paramLab,param), names_to = "Metric", values_to = "Value") %>%
  # plot
  ggplot()+
  geom_histogram(aes(y=Value,
                     x=Metric,
                     group=Metric),
                 stat = "identity")+
  
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  theme(axis.ticks =  element_blank(),
        axis.title.x = element_blank())+
  ylab("Error")

# sigma 2 T
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 3])
bias.sigma.sqT <- dat.plot %>%
  filter (param == 3) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = ((mean(val.theta)-sigma.sq.t)^2),
             Variance = (val.theta-mean(val.theta))^2, 
             #VarianceFUN=(var(val.theta)),
             MSE=(val.theta-sigma.sq.t)^2) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = mean(Bias),
             Variance = mean(Variance), 
             #VarianceFUN=mean(VarianceFUN),
             MSE=mean(MSE)) %>%
  pivot_longer(-c(scenario,spatial,time,paramLab,param), names_to = "Metric", values_to = "Value") %>%
  # plot
  ggplot()+
  geom_histogram(aes(y=Value,
                     x=Metric,
                     group=Metric),
                 stat = "identity")+
  
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  theme(axis.ticks =  element_blank(),
        axis.title.x = element_blank())+
  ylab("Error")

# rho
current_param_lab <- unique(dat.plot$paramLab[dat.plot$param == 4])
bias.rho<-dat.plot %>%
  filter (param == 4) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = ((mean(val.theta)-rho)^2),
             Variance = (val.theta-mean(val.theta))^2, 
             #VarianceFUN=(var(val.theta)),
             MSE=(val.theta-rho)^2) %>%
  group_by(spatial,time,paramLab,param,scenario) %>%
  summarize (Bias = mean(Bias),
             Variance = mean(Variance), 
             #VarianceFUN=mean(VarianceFUN),
             MSE=mean(MSE)) %>%
  pivot_longer(-c(scenario,spatial,time,paramLab,param), names_to = "Metric", values_to = "Value") %>%
  # plot
  ggplot()+
  geom_histogram(aes(y=Value,
                     x=Metric,
                     group=Metric),
                 stat = "identity")+
  
  facet_grid(time~spatial,scales="free",labeller = label_parsed) +
  ggtitle(parse(text = as.character(current_param_lab)) )+
  my_theme+
  theme(axis.ticks =  element_blank(),
        axis.title.x = element_blank())+
  ylab("Error")

# arrange plots to save
# sigma sq
png(here("figures","sims_present_paper", "Scenario1", "biasVatT_sigma-sq-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

  bias.sigma.sq

dev.off()

# phi
png(here("figures","sims_present_paper", "Scenario1", "biasVatT_phi-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

  bias.phi

dev.off()

# rho
png(here("figures","sims_present_paper", "Scenario1", "biasVatT_rho-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

  bias.rho
  
dev.off()

# sigma sq T
png(here("figures","sims_present_paper", "Scenario1", "biasVatT_sigma-sqT-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

  bias.sigma.sqT

dev.off()


# -------------------------------------------

# project psi and w in space

# add coordinates
table(dat[[1]]$coords[, 1] == dat[[2]]$coords[, 1] )
table(dat[[1]]$coords[, 1] == dat[[100]]$coords[, 1] )
# all coordinates are the same
plot.df$x <- dat[[1]]$coords[, 1] 
plot.df$y <- dat[[1]]$coords[, 2] 

# but not the spatial effect
# test
# scenario <- 1 # up to 16
# sim <- 3 # up to n.sim
year <- 1 # year don't matter here
# scenario.sim <- 33 # do scenario 1, 17, 33 ... resemble to each other?

# do this to repeat the 16 scenarios n.sims times
scenario.vals.rep <- lapply (seq (1,n.sims), function (i) 
    cbind (scenario.vals,
           sim = i)) # bind a column to depict the combination of simulations
scenario.vals.rep <- do.call(rbind, scenario.vals.rep)# melt
# now this scenario.vals.rep object has the same dim as the simulated data 'dat'

# loop over scenarios and sims and evaluate correlation table and plot
check_spatial_patterns <- lapply (seq(1,nrow(scenario.vals.rep)), function (i) {
  
       # build the dataframe with spatial coords, true and estimated psi and w 
      sp_data <- data.frame (x=dat[[i]]$coords[, 1] ,
                             y=dat[[i]]$coords[, 2] ,
                             scenario = scenario.vals.rep[i,"scenario"],
                             sim = scenario.vals.rep[i,"sim"],
                             psi.true = psi.true[,year,
                                                 scenario.vals.rep[i,"sim"],
                                                 scenario.vals.rep[i,"scenario"]],
                             psi.est = psi.mean.samples[,year,
                                                        scenario.vals.rep[i,"sim"],
                                                        scenario.vals.rep[i,"scenario"]],
                             w.true=dat[[i]]$w,
                             w.est=w.mean.samples[scenario.vals.rep[i,"sim"],,scenario.vals.rep[i,"scenario"]]
      ) 
      
      #  (melt to plot)
      sp_data_melt <- sp_data %>% 
        melt (id.vars = c("x", "y")) 
      
      # plot
      sp_params <- sp_data_melt %>%
        filter(variable %in% c("scenario","sim") ==F) %>%
        ggplot(aes(x = x, y = y, fill = value)) + 
        facet_wrap(~variable,nrow=2) +
        scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
                             na.value = NA, limits = c(min(sp_data_melt$value), max(sp_data_melt$value))) + 
        scale_x_continuous(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0)) + 
        geom_raster() + 
        theme_light(base_size = 22) + 
        labs(x = 'Easting', y = 'Northing', fill = "value") +
        my_theme+
        theme(legend.position = "none")
        
      
      # plot psi
      rel_psi <- sp_data %>%
        dplyr::select (c("psi.true","psi.est")) %>%
        ggplot(aes(x=psi.true, y=psi.est))+
        geom_point(alpha=0.3,col="gray40")+# loop over scenarios and sims and evaluate correlation table and plot
        geom_abline(col="red", lty = 2)+
        geom_smooth(method= "lm",se=F,col="black")+
        labs(y=bquote("Estimated occupancy "*(hat(psi[i]))*""),
             x=bquote("True occupancy "*(psi[i])*""))+
        my_theme
      
      # plot w
      rel_w<-sp_data %>%
        dplyr::select (c("w.true","w.est")) %>%
        ggplot(aes(x=w.true, y=w.est))+
        geom_point(alpha=0.3,col="gray40")+
        geom_abline(col="red", lty = 2)+
        geom_smooth(method= "lm",se=F,col="black")+
        labs(y=bquote("Estimated spatial random effect "*(hat(w[i]))*""),
             x=bquote("True spatial random effect "*(w[i])*""))+
        my_theme  
      
      # correlation
      corr_table <- sp_data %>%
        dplyr::select(c("psi.true","psi.est","w.true","w.est"))%>%
        cor() 
      
      # grid arrange
      arranged_plot <- grid.arrange (sp_params,
                                     rel_psi,
                                      rel_w,
                            layout_matrix=rbind(c(1,1,2),
                                        c(1,1,3)))
      # organize output
      output<- list (dataset = sp_data,
                     corr_table=corr_table,
                     arranged_plot=arranged_plot)
      ;
      # return
      output
})

# see the correlation matrix
sapply (check_spatial_patterns, "[","corr_table")

# transform the output in an array and calculate stats
output_array_sp<-array(unlist(
          sapply (check_spatial_patterns, "[","corr_table")),
      dim=c(4,4,length(check_spatial_patterns)))

# set names
dimnames(output_array_sp) <- list(c("psi.true","psi.est","w.true","w.est"),
                   c("psi.true","psi.est","w.true","w.est"),
                   seq(1,length(check_spatial_patterns)))
apply (output_array_sp,c(1,2),mean)
apply (output_array_sp,c(1,2),sd)

# plot correlations
corr.psi<-data.frame (corr=output_array_sp[1,2,]) %>%
  ggplot()+
  geom_histogram(aes(x=corr))+
  geom_density(aes(x=corr),fill="white",alpha=0.5)+
  geom_vline(xintercept=apply (output_array_sp,c(1,2),mean,na.rm=T)[1,2],lty = 2)+
  labs (x=  bquote("Pearson's correlation coefficient"*(rho(psi[i],hat(psi[i])))*""))+
  my_theme

#spatial random effect
corr.w<-data.frame (corr=output_array_sp[3,4,]) %>%
  ggplot()+
  geom_histogram(aes(x=corr))+
  geom_density(aes(x=corr),fill="white",alpha=0.5)+
  geom_vline(xintercept=apply (output_array_sp,c(1,2),mean,na.rm=T)[3,4],lty = 2)+
  labs (x=  bquote("Pearson's correlation coefficient"*(rho(w[i],hat(w[i])))*""))+
  my_theme

# arrange plot
png(here ("figures","sims_present_paper", "Scenario1", "panel_sp-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

grid.arrange(
  check_spatial_patterns[[1]]$arranged_plot,
  check_spatial_patterns[[3]]$arranged_plot,
  corr.psi,
  corr.w,
  nrow=2,ncol=2
)

dev.off()

# check groups in psi
corr_data <- cbind (melt(output_array_sp[1,2,]),
       scenario=seq(1,16)) %>%
  group_by(scenario) %>%
  summarise(corr=mean(value,na.rm=T)) %>%
  right_join(scenario.vals, by="scenario") 

# Ensure time and spatial columns are character vectors
corr_data$spatial <- as.character(corr_data$spatial)
corr_data$time <- as.character(corr_data$time)

# Convert time and spatial columns to factors for proper faceting
corr_data$spatial <- factor(corr_data$spatial, levels = unique(corr_data$spatial))
corr_data$time <- factor(corr_data$time, levels = unique(corr_data$time))

# the two first spatial scenarios had lower correlation between true and estimated psi
corr_data %>%
  ggplot()+
  geom_histogram(aes(x=corr))+
  facet_grid(time~spatial,scales="fixed",labeller = label_parsed) +
  labs (x=  bquote("Average Pearson's correlation coefficient"*(rho(psi[i],hat(psi[i])))*""))+
  my_theme +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())
  
# plot with all estimates
df.spatial <- do.call(rbind, sapply (check_spatial_patterns, "[","dataset")) 
df.spatial <- cbind(df.spatial,
                    scenario.vals[match(df.spatial$scenario,scenario.vals$scenario),c("spatial","time")] )

# Ensure time and spatial columns are character vectors
df.spatial$spatial <- as.character(df.spatial$spatial)
df.spatial$time <- as.character(df.spatial$time)

# Convert time and spatial columns to factors for proper faceting
df.spatial$spatial <- factor(df.spatial$spatial, levels = unique(df.spatial$spatial))
df.spatial$time <- factor(df.spatial$time, levels = unique(df.spatial$time))

# averaged df
av.df.spatial <-df.spatial %>%
  group_by(scenario) %>%
  arrange(psi.est, 
          psi.true,
          w.est,
          w.true,
          .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(psi.est = mean(psi.est),
            psi.true = mean(psi.true),
            w.est = mean(w.est),
            w.true = mean(w.true)
            ) %>%
  ungroup()

# match scneario data
av.df.spatial <- cbind(av.df.spatial,
                    scenario.vals[match(av.df.spatial$scenario,scenario.vals$scenario),c("spatial","time")] )

# Ensure time and spatial columns are character vectors
av.df.spatial$spatial <- as.character(av.df.spatial$spatial)
av.df.spatial$time <- as.character(av.df.spatial$time)

# Convert time and spatial columns to factors for proper faceting
av.df.spatial$spatial <- factor(av.df.spatial$spatial, levels = unique(av.df.spatial$spatial))
av.df.spatial$time <- factor(av.df.spatial$time, levels = unique(av.df.spatial$time))

# plot
fig.3.plot.w <- df.spatial%>%
  
  ggplot() +
  theme_light(base_size = 16) +
  geom_point(aes(x = w.true, y = w.est, group = sim),alpha=0.3,col="gray60")+
  geom_smooth(aes(x = w.true, y = w.est, group = sim), col = 'white', alpha = 0.1, se = FALSE, 
              lineend = 'round', lwd = 0.25) +
  geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2) +
  geom_smooth(data = av.df.spatial, aes(x = w.true, y = w.est), se = FALSE, 
              col = 'black', 
              lwd = 0.5,
              lineend = 'round') + 
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  labs(y=bquote("Estimated spatial random effect "*(hat(w[i]))*""),
               x=bquote("True spatial random effect "*(w[i])*""))+
  my_theme

# plot 
png(here("figures","sims_present_paper", "Scenario1", "Figure-3-spatial-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

fig.3.plot.w

dev.off()

# -------------------------------------

# temporal autocorrelation
# loop over scenarios and sims and evaluate correlation table and plot
check_time_patterns <- lapply (seq(1,nrow(scenario.vals.rep)), function (i) {
  
      # build the dataframe with spatial coords, true and estimated psi and w 
      t_data <- data.frame (year=seq(1,10),
                             eta.true = dat[[i]]$eta,
                             eta.est = eta.mean.samples[scenario.vals.rep[i,"sim"],
                                                        ,
                                                        scenario.vals.rep[i,"scenario"]       
                                                        ],
                            scenario=scenario.vals.rep[i,"scenario"],
                            sim=scenario.vals.rep[i,"sim"]) 
      # plot psi
      rel_eta <- t_data %>%
        dplyr::select (c("eta.true","eta.est")) %>%
        ggplot(aes(x=eta.true, y=eta.est))+
        geom_point()+
        geom_abline(col="red", )+
        geom_smooth(method= "lm",se=F,col="black")+
        labs(y=bquote("Estimated time random effect "*(hat(eta[t]))*""),
             x=bquote("True time random effect "*(eta[t])*""))+
        my_theme
      
      
      # organize output
      output<- list (dataset = t_data,
                     corr_table=cor(t_data%>%
                                      dplyr::select (c("eta.true","eta.est"))), # correlation
                     arranged_plot=rel_eta)
      ;
      # return
      output
      
})

# see the correlation matrix
sapply (check_time_patterns, "[","corr_table")

# transform the output in an array and calculate stats
output_array_t<-array(unlist(
  sapply (check_time_patterns, "[","corr_table")),
  dim=c(2,2,length(check_time_patterns)))

# set names
dimnames(output_array_t) <- list(c("eta.true","eta.est"),
                                  c("eta.true","eta.est"),
                                  seq(1,length(check_time_patterns)))
apply (output_array_t,c(1,2),mean,na.rm=T)
apply (output_array_t,c(1,2),sd,na.rm=T)

# check
plot(check_time_patterns[[1]]$arranged_plot)


# plot with all estimates
df.temporal <- do.call(rbind, sapply (check_time_patterns, "[","dataset")) 
df.temporal <- cbind(df.temporal,
                    scenario.vals[match(df.temporal$scenario,scenario.vals$scenario),c("spatial","time")] )

# Ensure time and spatial columns are character vectors
df.temporal$spatial <- as.character(df.temporal$spatial)
df.temporal$time <- as.character(df.temporal$time)

# Convert time and spatial columns to factors for proper faceting
df.temporal$spatial <- factor(df.temporal$spatial, levels = unique(df.temporal$spatial))
df.temporal$time <- factor(df.temporal$time, levels = unique(df.temporal$time))

# averaged df
av.df.temporal <-df.temporal %>%
  group_by(scenario) %>%
  arrange(eta.true, 
          eta.est,
          .by_group = TRUE) %>%
  mutate(fake.id = 1:n()) %>%
  ungroup() %>%
  group_by(scenario, fake.id) %>%
  summarize(eta.est = mean(eta.est),
            eta.true = mean(eta.true)
  ) %>%
  ungroup()

# match scneario data
av.df.temporal <- cbind(av.df.temporal,
                       scenario.vals[match(av.df.temporal$scenario,scenario.vals$scenario),c("spatial","time")] )

# Ensure time and spatial columns are character vectors
av.df.temporal$spatial <- as.character(av.df.temporal$spatial)
av.df.temporal$time <- as.character(av.df.temporal$time)

# Convert time and spatial columns to factors for proper faceting
av.df.temporal$spatial <- factor(av.df.temporal$spatial, levels = unique(av.df.temporal$spatial))
av.df.temporal$time <- factor(av.df.temporal$time, levels = unique(av.df.temporal$time))

# plot
fig.3.plot.eta <- df.temporal%>%
  
  ggplot() +
  theme_light(base_size = 16) +
  geom_point(aes(x = eta.true, y = eta.est, group = sim),alpha=0.5,col="gray20")+
  geom_smooth(aes(x = eta.true, y = eta.est, group = sim), col = 'white', alpha = 1, se = FALSE, 
              lineend = 'round', lwd = 0.25) +
  geom_abline(slope = 1, intercept = 0, col = 'red', lty = 2) +
  geom_smooth(data = av.df.temporal, aes(x = eta.true, y = eta.est), se = FALSE, 
              col = 'black', 
              lwd = 0.5,
              lineend = 'round') + 
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  labs(y=bquote("Estimated temporal random effect "*(hat(eta[t]))*""),
       x=bquote("True temporal random effect "*(eta[t])*""))+
  my_theme

png(here("figures","sims_present_paper", "Scenario1", "Figure-3-temporal-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

fig.3.plot.eta

dev.off()

# -----------------------------------------------------

psi_p_data <- data.frame (sim=melt(beta.mean.samples[1,,])[,1],
            scenario=melt(beta.mean.samples[1,,])[,2],
            beta0=melt(beta.mean.samples[1,,])[,3],
            beta0.lwr=melt(beta.low.samples[1,,])[,3],
            beta0.upr=melt(beta.high.samples[1,,])[,3],
            beta1=melt(beta.mean.samples[2,,])[,3],
            beta1.lwr=melt(beta.low.samples[2,,])[,3],
            beta1.upr=melt(beta.high.samples[2,,])[,3],
            alpha0= melt(alpha.mean.samples[1,,])[,3],
            alpha0.lwr= melt(alpha.low.samples[1,,])[,3],
            alpha0.upr= melt(alpha.high.samples[1,,])[,3],
            alpha1= melt(alpha.mean.samples[2,,])[,3],
            alpha1.lwr= melt(alpha.low.samples[2,,])[,3],
            alpha1.upr= melt(alpha.high.samples[2,,])[,3]
            ) %>%
  
  right_join(scenario.vals,by="scenario") 

# Ensure time and spatial columns are character vectors
psi_p_data$spatial <- as.character(psi_p_data$spatial)
psi_p_data$time <- as.character(psi_p_data$time)

# Convert time and spatial columns to factors for proper faceting
psi_p_data$spatial <- factor(psi_p_data$spatial, levels = unique(psi_p_data$spatial))
psi_p_data$time <- factor(psi_p_data$time, levels = unique(psi_p_data$time))

# density
density.psi.p <- psi_p_data %>%
  mutate (spatial=as.factor(spatial),
          time = as.factor(time)) %>%
  
      ggplot( aes(x = plogis(beta0), y = plogis(alpha0))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+   
               geom_point(alpha=0.1) + 
               scale_fill_viridis_d(direction = -1,option="magma")+
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote("logit-1"*(hat(alpha)[0])*""),
       x = bquote("logit-1"*(hat(beta)[0])*""))+
    geom_point (aes(x=plogis(beta[1]),y=plogis(alpha[1])),col="white",shape=4)+
  my_theme
                         
# save
png(here("figures","sims_present_paper", "Scenario1", "density.psi.p-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

  density.psi.p

dev.off()


# try this
# https://stackoverflow.com/questions/13613157/plot-3d-density
#require(plotly)
#library(MASS)
#
#test_sc1 <- psi_p_data %>%
#  filter (is.na(alpha0)!=T) %>%
#  filter (scenario == 1) 
#den3d <- kde2d(test_sc1$beta0, 
#               test_sc1$alpha0)
#
## plot 
#library(plotly)
#plot_surface1<-plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()
#
#kaleido(plot_surface1, here("figures","sims_present_paper", "Scenario1", #"surface_sce1.png"))
#
#
#test_sc16 <- psi_p_data %>%
#  filter (is.na(alpha0)!=T) %>%
#  filter (scenario == 16) 
#den3d_sc16 <- kde2d(test_sc16$beta0, 
#               test_sc16$alpha0)
#
## the new part:
#require(patchwork)
#
#plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()
#plot_ly(x=den3d_sc16$x, y=den3d_sc16$y, z=den3d_sc16$z) %>% add_surface()
#

# add bars of lower and upper CI

# density
density.psi.p.bars <- psi_p_data %>%
  mutate (spatial=as.factor(spatial),
          time = as.factor(time)) %>%
  
  ggplot( aes(x = plogis(beta0), y = plogis(alpha0))) +
  
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+   
  geom_errorbar(
    aes(ymin = plogis(alpha0.lwr), 
        ymax = plogis(alpha0.upr)), 
    color = "gray",alpha=0.5
  )+
  geom_errorbar(
    aes(xmin= plogis(beta0.lwr),
        xmax= plogis(beta0.upr)), 
    color = "gray",alpha=0.5
  )+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma")+
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote("logit-1"*(hat(alpha)[0])*""),
       x = bquote("logit-1"*(hat(beta)[0])*""))+
  geom_point (aes(x=plogis(beta[1]),y=plogis(alpha[1])),col="white",shape=4)+
  my_theme

# save
png(here("figures","sims_present_paper", "Scenario1", "density.psi.p.bars-st1-sc1.png"),
    width = 800, height = 800, units = "px")

  density.psi.p.bars

dev.off()

# relationship between intercept and regression coef
png(here("figures","sims_present_paper", "Scenario1", "rel.beta0.beta1-st1-sc1.png"),
    width = 20, height = 20, units = "cm", res=200)

grid.arrange (
  
  # occupancy intercept and reg coeff
  psi_p_data %>%
    
  ggplot(aes (x=beta1,y=plogis(beta0)))+
    geom_errorbar(
      aes(ymin = plogis(beta0.lwr), 
          ymax = plogis(beta0.upr)), 
      color = "gray",alpha=0.5
    )+
    geom_errorbar(
      aes(xmin= (beta1.lwr),
          xmax= (beta1.upr)), 
      color = "gray",alpha=0.5
    ) +
    
  geom_point()+
  geom_smooth(col="orange2")+
  my_theme+
  # show extreme values
  geom_point(data=psi_p_data[(which(plogis(psi_p_data$beta0) < 0.25 | 
                                      plogis(psi_p_data$beta0) > 0.75)),], 
             (aes (x=beta1,y=plogis(beta0))),
             col="red",size=0.8)+
  labs (y=bquote("logit-1"*(hat(beta[0]))*""),
        x=bquote(""*hat(beta[1])*""),
        title="A")+
    xlim(c(0.3,0.8))+
    ylim(c(0,1))+
    
    
  geom_point(aes(x=0.5,y=0.5),col="orange",size=4,shape=4)
  
  ,
  
  # detection intercept and reg coeff
  psi_p_data %>%
      ggplot(aes (x=alpha1,y=plogis(alpha0)))+
    geom_errorbar(
      aes(ymin = plogis(alpha0.lwr), 
          ymax = plogis(alpha0.upr)), 
      color = "gray",alpha=0.5
    )+
    geom_errorbar(
      aes(xmin= (alpha1.lwr),
          xmax= (alpha1.upr)), 
      color = "gray",alpha=0.5
    ) +
      geom_point()+
      geom_smooth(col="orange2")+
      my_theme+
      # show extreme values
      geom_point(data=psi_p_data[(which(plogis(psi_p_data$beta0) < 0.25 | 
                                          plogis(psi_p_data$beta0) > 0.75)),], 
                 (aes (x=alpha1,y=plogis(alpha0))),
                 col="red",size=0.8)+
      labs (y=bquote("logit-1"*(hat(alpha[0]))*""),
            x=bquote(""*hat(alpha[1])*""),
            title="B")+
    ylim(c(0.33,0.72))+
    xlim(c(-.8,-.2))+
      geom_point(aes(x=-0.5,y=0.5),col="orange",size=4,shape=4)
    
  
  
  , ncol=2 
)

dev.off()

# end