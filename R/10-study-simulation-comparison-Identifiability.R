
# -------------------------------------------------------

#  COMPARISON OF STUDIES, SCENARIOS AND SUB SCENARIOS
#  
#  Code to produce most figures used in the main text

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

# ggplot theme
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
dir.create(here("figures", "Scenario-Comparisons"))

# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# load model output
sceDSsmooth <- new.env()
sceDS <- new.env()
sce0 <- new.env()
sce1 <- new.env()
sce2 <- new.env()
sce3 <- new.env()
scePhen <- new.env()
scePhenSpot <- new.env()

# load in the environments
# # 
load(file = here ("model_output", "output_simulations", "sims_D&S",
                  "sim-mixed-stPGOcc-results_1600.rda"),sceDS)
sceDS$study <- 1
sceDS$sc <- 0

load(file = here ("model_output", "output_simulations", "smooth_sims_D&S",
                  "sim-mixed-stPGOcc-results_1600.rda"),sceDSsmooth)
sceDSsmooth$study <- 1
sceDSsmooth$sc <- 1


load(file = here ("model_output", "output_simulations", "scenario_zero",
                  "sim-mixed-stPGOcc-results_1600.rda"),sce0)
sce0$study <- 2
sce0$sc <- 0

load(file = here ("model_output", "output_simulations", "scenario_one",
                  "sim-mixed-stPGOcc-results_1600.rda"),sce1)

sce1$study <- 2
sce1$sc <- 1

load(file = here ("model_output", "output_simulations", "scenario_two",
                  "sim-mixed-stPGOcc-results_1600.rda"),sce2)
sce2$study <- 2
sce2$sc <- 2

load(file = here ("model_output", "output_simulations", "scenario_three",
                  "sim-mixed-stPGOcc-results_1600.rda"),sce3)
sce3$study <- 2
sce3$sc <- 3

load(file = here ("model_output", "output_simulations", "scenario_phenology",
                  "sim-mixed-stPGOcc-results_1600.rda"),scePhen)

scePhen$study <- 3
scePhen$sc <- 1

load(file = here ("model_output", "output_simulations", "scenario_phenology_spot",
                  "sim-mixed-stPGOcc-results_1600.rda"),scePhenSpot)

scePhenSpot$study <- 3
scePhenSpot$sc <- 2

# -------------------------------------------

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

# error plot ------------------------------------------
# psi_it ----------------

# estimated values
BVDec_psi_it <- lapply (list (sceDS,
                              sceDSsmooth,
                              sce0,
                              sce1,
                              sce2,
                              sce3,
                              scePhen,
                              scePhenSpot), function (i) {
  
  # estimated values
  BVDec_psi_it <- melt(i$psi.mean.samples) %>% 
        # bind the true
        cbind ( melt(i$psi.true) %>% 
                  dplyr::select("value") %>%
                  mutate (true = value) %>%
                  dplyr::select(-"value")) %>%
        # filter any NA
        filter(is.na(value) !=T) %>%
        # names variables
        dplyr::rename("sites"="X1",
                      "years"="X2" ,
                      "sims"= "X3",
                      "scenario" = "X4",
                      "value"="value") %>%
        #  bind sub scenarios IDs
        right_join(.,scenario.vals,
                   by="scenario") %>%
        # bind scenario IDs
        mutate(Study=i$study,
               sc=i$sc)%>%
        data.frame() %>%
        
        # calculate the difference of each value relative to the true
        mutate (diff = ((value) - true),
                # calculate variance
                var = (value - mean(value))^2) %>%
          
        # agreggate and make the bias variance decomposition
        group_by(spatial,time,scenario,Study,sc) %>%
        summarize (Bias = (mean(diff)),
                   Variance = mean(var),
                   MSE=(Bias)^2
                   ) %>%
        
        # long format
        pivot_longer(-c(scenario,spatial,time,Study,sc), names_to = "Metric", values_to = "Value")
        
      ;
  
  BVDec_psi_it

  }
  
)

# melt to plot
BVDec_psi_it <- do.call(rbind, BVDec_psi_it)
  
# Ensure time and spatial columns are character vectors
BVDec_psi_it$spatial <- as.character(BVDec_psi_it$spatial)
BVDec_psi_it$time <- as.character(BVDec_psi_it$time)

# Convert time and spatial columns to factors for proper faceting
BVDec_psi_it$spatial <- factor(BVDec_psi_it$spatial, levels = unique(BVDec_psi_it$spatial))
BVDec_psi_it$time <- factor(BVDec_psi_it$time, levels = unique(BVDec_psi_it$time))

# order metric order
BVDec_psi_it$Metric <- factor(BVDec_psi_it$Metric,
                               levels = c("Bias","Variance", "MSE"))

# plot
fig_error_psi <- ggplot(BVDec_psi_it %>%
                          filter (Metric == "MSE"))+
  
  # bars
  geom_col(aes(y=Value,
               x=paste(Study,sc,sep="-"),
               fill=Metric,
               group=sc),
               position="stack") +

  # grid
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +

  
  # color for study
  scale_fill_manual(values=c("MSE"="gray40","Bias"="gray60", "Variance"="gray20"))+
  scale_colour_manual(values=c("MSE"="gray80","Bias"="gray60", "Variance"="black"))+
  ggtitle(parse(text = as.character(expression(psi[it]))) )+
  my_theme+
  theme(#axis.text.x = element_blank(),
    legend.position = "none",
    axis.ticks =  element_blank()#,
    #axis.title.x = element_blank()
    )+
  labs(y="Mean Squared Error", x="Study-scenario")

# save
png (here("figures", "Scenario-Comparisons","Error_psi_it.png"),
     width = 600, height = 600,units = "px")

  fig_error_psi

dev.off()


# phi -------------------------------------------------

# estimated values
BVDec_phi <- lapply (list (sceDS,
                              sceDSsmooth,
                              sce0,
                              sce1,
                              sce2,
                              sce3,
                              scePhen,
                              scePhenSpot), function (i) {
                                
                                # estimated values
                                BVDec_phi <- melt(i$theta.mean.samples[,2,]) %>% 
                                  filter(is.na(value) !=T) %>%
                                  dplyr::rename("sims"= "X1",
                                                "scenario" = "X2",
                                                "value"="value") %>%
                                  # bind the true
                                  right_join(.,scenario.vals,
                                             by="scenario") %>%
                                  # bind scenario IDs
                                  mutate(Study=i$study,
                                         sc=i$sc)%>%
                                  data.frame() %>%
                                  
                                  # calculate the difference of each value relative to the true
                                  mutate (diff = ((value) - phi),
                                          # calculate variance
                                          var = (value - mean(value))^2) %>%
                                  
                                  # agreggate and make the bias variance decomposition
                                  group_by(spatial,time,scenario,Study,sc) %>%
                                  summarize (Mean = mean (value),
                                             LCI = quantile(value,0.025),
                                             UCI=quantile(value,0.975),
                                             Bias = (mean(diff)),
                                             Variance = mean(var),
                                             MSE=(Bias)^2
                                  ) %>%
                                  
                                  # long format
                                  pivot_longer(-c(scenario,spatial,time,Study,sc), names_to = "Metric", values_to = "Value")
                                
                                ;
                                
                                BVDec_phi
                                
                              }
                        
)

# melt to plot
BVDec_phi <- do.call(rbind, BVDec_phi)

# Ensure time and spatial columns are character vectors
BVDec_phi$spatial <- as.character(BVDec_phi$spatial)
BVDec_phi$time <- as.character(BVDec_phi$time)

# Convert time and spatial columns to factors for proper faceting
BVDec_phi$spatial <- factor(BVDec_phi$spatial, levels = unique(BVDec_phi$spatial))
BVDec_phi$time <- factor(BVDec_phi$time, levels = unique(BVDec_phi$time))

# order metric order
BVDec_phi$Metric <- factor(BVDec_phi$Metric,
                              levels = c("Mean","LCI","UCI", "Bias","Variance", "MSE"))

# plot
fig_error_phi <- ggplot(BVDec_phi %>%
                          filter (Metric == "MSE"))+
  
  # bars
  geom_col(aes(y=Value,
               x=paste(Study,sc,sep="-"),
               fill=Metric,
               group=sc),
           position="stack") +
  
  # grid 
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  
  
  # color for study
  scale_fill_manual(values=c("MSE"="gray40","Bias"="gray60", "Variance"="gray20"))+
  scale_colour_manual(values=c("MSE"="gray80","Bias"="gray60", "Variance"="black"))+
  ggtitle(parse(text = as.character(expression(phi))) )+
  my_theme+
  theme(#axis.text.x = element_blank(),
    legend.position = "none",
    axis.ticks =  element_blank()
  )+
  labs(y="Mean Squared Error", x="Study-scenario") 

# save
png (here("figures", "Scenario-Comparisons","Error_phi.png"),
     width = 600, height = 600,units = "px")

  fig_error_phi

dev.off()


# phi estimates across scenarios --------------------
phi_decay <- ggplot(BVDec_phi %>%
         filter (Metric == "Mean"))+
  
  # bars
   geom_line(aes(y=Value,
                 x=paste(Study,sc,sep="-"),
                 colour=Metric,
                 group=Metric),
             stat = "identity",
             linewidth=2,
             inherit.aes = T) +
  
  # lower CI
  geom_line(data=BVDec_phi %>%
              filter (Metric == "LCI"),
            aes(y=Value,
                x=paste(Study,sc,sep="-"),
                colour=Metric,
                group=Metric),
            stat = "identity",
            linewidth=1,
            inherit.aes = F) +
  # upper CI
  geom_line(data=BVDec_phi %>%
              filter (Metric == "UCI"),
            aes(y=Value,
                x=paste(Study,sc,sep="-"),
                colour=Metric,
                group=Metric),
            stat = "identity",
            linewidth=1,
            inherit.aes = F) +
  
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  
  
  # color for study
  scale_fill_manual(values=c("MSE"="gray40","Bias"="gray60", "Variance"="gray20"))+
  scale_colour_manual(values=c("MSE"="gray80","Bias"="gray60", "Variance"="black"))+
  ggtitle(parse(text = as.character(expression(phi))) )+
  my_theme+
  theme(#axis.text.x = element_blank(),
    legend.position = "none",
    axis.ticks =  element_blank()#,
    #axis.title.x = element_blank()
  )+
  labs(y=bquote ("Spatial decay "*hat(phi)*""), x="Study-scenario")

# save
png (here("figures", "Scenario-Comparisons","Phi_decay_comparison.png"),
     width = 600, height = 600,units = "px")

  phi_decay

dev.off()


# rho -------------------------------------------------

# estimated values
BVDec_rho <- lapply (list (sceDS,
                              sceDSsmooth,
                              sce0,
                              sce1,
                              sce2,
                              sce3,
                              scePhen,
                              scePhenSpot), function (i) {
                                
                                # estimated values
                                BVDec_rho <- melt(i$theta.mean.samples[,4,]) %>% 
                                  filter(is.na(value) !=T) %>%
                                  dplyr::rename("sims"= "X1",
                                                "scenario" = "X2",
                                                "value"="value") %>%
                                  # bind the true
                                  right_join(.,scenario.vals,
                                             by="scenario") %>%
                                  # bind scenario IDs
                                  mutate(Study=i$study,
                                         sc=i$sc)%>%
                                  data.frame() %>%
                                  
                                  # calculate the difference of each value relative to the true
                                  mutate (diff = ((value) - rho),
                                          # calculate variance
                                          var = (value - mean(value))^2) %>%
                                  
                                  # agreggate and make the bias variance decomposition
                                  group_by(spatial,time,scenario,Study,sc) %>%
                                  summarize (Bias = (mean(diff)),
                                             Variance = mean(var),
                                             MSE=(Bias)^2
                                  ) %>%
                                  
                                  # long format
                                  pivot_longer(-c(scenario,spatial,time,Study,sc), names_to = "Metric", values_to = "Value")
                                
                                ;
                                
                                BVDec_rho
                                
                              }
                        
)

# melt to plot
BVDec_rho <- do.call(rbind, BVDec_rho)

# Ensure time and spatial columns are character vectors
BVDec_rho$spatial <- as.character(BVDec_rho$spatial)
BVDec_rho$time <- as.character(BVDec_rho$time)

# Convert time and spatial columns to factors for proper faceting
BVDec_rho$spatial <- factor(BVDec_rho$spatial, levels = unique(BVDec_rho$spatial))
BVDec_rho$time <- factor(BVDec_rho$time, levels = unique(BVDec_rho$time))

# order metric order
BVDec_rho$Metric <- factor(BVDec_rho$Metric,
                              levels = c("Bias","Variance", "MSE"))

# plot
fig_error_rho <- ggplot(BVDec_rho%>%
                          filter (Metric == "MSE"))+
  
  # bars
  geom_col(aes(y=Value,
               x=paste(Study,sc,sep="-"),
               fill=Metric,
               group=sc),
           position="stack") +
  
  # grid
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  
  
  # color for study
  scale_fill_manual(values=c("MSE"="gray40","Bias"="gray60", "Variance"="gray20"))+
  scale_colour_manual(values=c("MSE"="gray80","Bias"="gray60", "Variance"="black"))+
  ggtitle(parse(text = as.character(expression(rho))) )+
  my_theme+
  theme(#axis.text.x = element_blank(),
    legend.position = "none",
    axis.ticks =  element_blank()
  )+
  labs(y="Mean Squared Error", x="Study-scenario")

# save
png (here("figures", "Scenario-Comparisons","Error_rho.png"),
     width = 600, height = 600,units = "px")

  fig_error_rho

dev.off()

# sigma^2_T -------------------------------------------------

# estimated values
BVDec_sigma2_T <- lapply (list (sceDS,
                           sceDSsmooth,
                           sce0,
                           sce1,
                           sce2,
                           sce3,
                           scePhen,
                           scePhenSpot), function (i) {
                             
                             # estimated values
                             BVDec_sigma2_T <- melt(i$theta.mean.samples[,3,]) %>% 
                               filter(is.na(value) !=T) %>%
                               dplyr::rename("sims"= "X1",
                                             "scenario" = "X2",
                                             "value"="value") %>%
                               # bind the true
                               right_join(.,scenario.vals,
                                          by="scenario") %>%
                               # bind scenario IDs
                               mutate(Study=i$study,
                                      sc=i$sc)%>%
                               data.frame() %>%
                               
                               # calculate the difference of each value relative to the true
                               mutate (diff = ((value) - sigma.sq.t),
                                       # calculate variance
                                       var = (value - mean(value))^2) %>%
                               
                               # agreggate and make the bias variance decomposition
                               group_by(spatial,time,scenario,Study,sc) %>%
                               summarize (Bias = (mean(diff)),
                                          Variance = mean(var),
                                          MSE=(Bias)^2
                               ) %>%
                               
                               # long format
                               pivot_longer(-c(scenario,spatial,time,Study,sc), names_to = "Metric", values_to = "Value")
                             
                             ;
                             
                             BVDec_sigma2_T
                             
                           }
                     
)

# melt to plot
BVDec_sigma2_T <- do.call(rbind, BVDec_sigma2_T)

# Ensure time and spatial columns are character vectors
BVDec_sigma2_T$spatial <- as.character(BVDec_sigma2_T$spatial)
BVDec_sigma2_T$time <- as.character(BVDec_sigma2_T$time)

# Convert time and spatial columns to factors for proper faceting
BVDec_sigma2_T$spatial <- factor(BVDec_sigma2_T$spatial, levels = unique(BVDec_sigma2_T$spatial))
BVDec_sigma2_T$time <- factor(BVDec_sigma2_T$time, levels = unique(BVDec_sigma2_T$time))

# order metric order
BVDec_sigma2_T$Metric <- factor(BVDec_sigma2_T$Metric,
                           levels = c("Bias","Variance", "MSE"))

# plot
fig_error_sigma2_T <- ggplot(BVDec_sigma2_T%>%
                          filter (Metric == "MSE"))+
  
  # bars
  geom_col(aes(y=Value,
               x=paste(Study,sc,sep="-"),
               fill=Metric,
               group=sc),
           position="stack") +
  
  # grid
  facet_grid(time ~ spatial, 
             labeller = label_parsed) +
  
  
  # color for study
  scale_fill_manual(values=c("MSE"="gray40","Bias"="gray60", "Variance"="gray20"))+
  scale_colour_manual(values=c("MSE"="gray80","Bias"="gray60", "Variance"="black"))+
  ggtitle(parse(text = as.character(expression(sigma^2[T]))))+
  my_theme+
  theme(#axis.text.x = element_blank(),
    legend.position = "none",
    axis.ticks =  element_blank()
  )+
  labs(y="Mean Squared Error", x="Study-scenario")

# save
png (here("figures", "Scenario-Comparisons","Error_sigma2_T.png"),
     width = 600, height = 600,units = "px")

  fig_error_sigma2_T

dev.off()


# -----------------------------------------------------
# map of true and estimated psi_it
# load data 
load(file = here ("model_output", "output_simulations", "sims_D&S", "sim-data-correct.rda"))
gc()

# choose a sim run to show

chosen_run <- 2

# Reshape for plotting  ---- do for scenario 14, high correlations
occupancy_data <- lapply (list (sceDS,
                             sceDSsmooth,
                             sce0,
                             sce1,
                             sce2,
                             sce3,
                             scePhen,
                             scePhenSpot), function (i) { 
                               
              data.frame(
                X = dat[[1]]$coords[,1],
                Y = dat[[1]]$coords[,2],
                True = i$psi.true[,1,chosen_run,14],
                Estimated = i$psi.mean.samples[,1,chosen_run,14],
                Study=i$study,
                sc=i$sc
                )

                               })
# melt data
occupancy_data <- do.call(rbind,occupancy_data)

# Plot the true occupancy data for the first year
map_true_vs_estimated <- ggplot(occupancy_data %>%
         mutate (Difference =  Estimated - True) %>%
           reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
         filter (variable != "Difference"), 
         aes(x = X, y = Y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option="magma") +
    coord_fixed() +
    theme_minimal() +
    facet_grid(variable~Study+sc)+
    labs(#title = "Site Occupancy Probability (for t=5)",
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
map_difference <- ggplot(occupancy_data %>%
                                  mutate (Difference = Estimated - True ) %>%
                                  reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                                  filter (variable == "Difference"), 
                                aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2()+
  coord_fixed() +
  theme_minimal() +
  facet_grid(variable~Study+sc)+
  labs(#title = "Site Occupancy Probability (for t=5)",
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


png (here("figures", "Scenario-Comparisons","maps_comparison_sc14.png"),
     width = 20, height = 17,units = "cm",res=200)

  grid.arrange(map_true_vs_estimated+theme(plot.margin=unit(c(0,1,1,0.5), "cm")),
               map_difference+theme(plot.margin=unit(c(-9,1,1,0.5), "cm")),
               nrow=2)

dev.off()

# do also for scenario 3 - low correlations overall ------------------
occupancy_data_low <- lapply (list (sceDS,
                                sceDSsmooth,
                                sce0,
                                sce1,
                                sce2,
                                sce3,
                                scePhen,
                                scePhenSpot), function (i) { 
                                  
                                  data.frame(
                                    X = dat[[1]]$coords[,1],
                                    Y = dat[[1]]$coords[,2],
                                    True = i$psi.true[,1,chosen_run,3],
                                    Estimated = i$psi.mean.samples[,1,chosen_run,3],
                                    Study=i$study,
                                    sc=i$sc
                                  )
                                  
                                })
# melt data
occupancy_data_low <- do.call(rbind,occupancy_data_low)

# Plot the true occupancy data for the first year
map_true_vs_estimated_low <- ggplot(occupancy_data_low %>%
                                  mutate (Difference = Estimated - True) %>%
                                  reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                                  filter (variable != "Difference"), 
                                aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma") +
  coord_fixed() +
  theme_minimal() +
  facet_grid(variable~Study+sc)+
  labs(#title = "Site Occupancy Probability (for t=5)",
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
map_difference_low <- ggplot(occupancy_data_low %>%
                           mutate (Difference = Estimated-True) %>%
                           reshape::melt (id.vars = c("X","Y","Study","sc")) %>%
                           filter (variable == "Difference"), 
                         aes(x = X, y = Y, fill = value)) +
  geom_tile() +
  scale_fill_gradient2()+
  coord_fixed() +
  theme_minimal() +
  facet_grid(variable~Study+sc)+
  labs(#title = "Site Occupancy Probability (for t=5)",
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


png (here("figures", "Scenario-Comparisons","maps_comparison_sc3.png"),
     width = 20, height = 17,units = "cm",res=200)
  
  grid.arrange(map_true_vs_estimated_low+theme(plot.margin=unit(c(0,1,1,0.5), "cm")),
               map_difference_low+theme(plot.margin=unit(c(-9,1,1,0.5), "cm")),
               nrow=2)

dev.off()

# -----------------------------------------------------
# relationship between spatial and temporal random effects (theta)

theta_data <- lapply (list (sceDS,
                            sceDSsmooth,
                            sce0,
                            sce1,
                            sce2,
                            sce3,
                            scePhen,
                            scePhenSpot), function (i) {

                            # bind data
                            theta_data <-  rbind (
                                  
                              data.frame (# occupancy intercept and coefficients
                                    sim=melt(i$theta.mean.samples)[,1],
                                    scenario=melt(i$theta.mean.samples)[,3],
                                    "sigma_sq"=melt(i$theta.mean.samples[,1,])[,3],
                                    "phi"=melt(i$theta.mean.samples[,2,])[,3],
                                    "sigma_sqT"=melt(i$theta.mean.samples[,3,])[,3],
                                    "rho"=melt(i$theta.mean.samples[,4,])[,3],
                                    
                                    study=i$study,
                                    sc=i$sc
                                  )) %>%
                                  
                                  right_join(scenario.vals,by="scenario")
                                
                              ;
                              
                              theta_data
                              
                            })

# melt to a df
theta_data <- do.call(rbind,theta_data)

# Ensure time and spatial columns are character vectors
theta_data$spatial <- as.character(theta_data$spatial)
theta_data$time <- as.character(theta_data$time)

# Convert time and spatial columns to factors for proper faceting
theta_data$spatial <- factor(theta_data$spatial, levels = unique(theta_data$spatial))
theta_data$time <- factor(theta_data$time, levels = unique(theta_data$time))

# adjust the true for the scenario 1-1 (lower values of phi)
to_replace <- theta_data [which(theta_data$study == 1 & theta_data$sc == 1),"phi.y"]
theta_data [which(theta_data$study == 1 & theta_data$sc == 1),"phi.y"] <- ifelse (to_replace == 15,1,0.5)

# density
density_phi_sigma <- theta_data %>%
  filter(scenario == c(1:4) &
           is.na(phi.x)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = phi.x, y = sigma_sq)) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  #xlim(c(0,1))+
  ylim(c(0,4))+
  theme_bw()+
  labs(y = bquote(""*hat(sigma^2)*""),
       x = bquote(""*hat(phi)*""))+
  geom_point (aes(x=phi.y,y=sigma.sq),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_phi_sigma1-4.png"),
     width = 17, height = 17,units = "cm",res=150)

  density_phi_sigma

dev.off()

# density
density_phi_sigma2 <- theta_data %>%
  filter(scenario == c(13:16) &
           is.na(phi.x)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = phi.x, y = sigma_sq)) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  #xlim(c(0,1))+
  ylim(c(0,4))+
  theme_bw()+
  labs(y = bquote(""*hat(sigma^2)*""),
       x = bquote(""*hat(phi)*""))+
  geom_point (aes(x=(phi.y),y=(sigma.sq)),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_phi_sigma13-16.png"),
     width = 17, height = 17,units = "cm",res=150)

  density_phi_sigma2

dev.off()

# density rho and sigma_sqT ---------------------
density_rho_sigmaT <- theta_data %>%
  filter(scenario == c(1:4) &
           is.na(phi.x)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = rho.x, y = sigma_sqT)) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(-1,1))+
  #ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote(""*hat(sigma^2)[T]*""),
       x = bquote(""*hat(rho)*""))+
  geom_point (aes(x=(rho.y),y=(sigma.sq.t)),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_rho_sigmaT1-4.png"),
     width = 17, height = 17,units = "cm",res=150)

  density_rho_sigmaT

dev.off()

# density rho and sigma_sqT ---------------------
density_rho_sigmaT2 <- theta_data %>%
  filter(scenario == c(13:16) &
           is.na(phi.x)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = rho.x, y = sigma_sqT)) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(-1,1))+
  #ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote(""*hat(sigma^2)[T]*""),
       x = bquote(""*hat(rho)*""))+
  geom_point (aes(x=(rho.y),y=(sigma.sq.t)),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_rho_sigmaT13-16.png"),
     width = 17, height = 17,units = "cm",res=150)

  density_rho_sigmaT2

dev.off()

# relationship between intercepts and regression coefficients-----------------------------------------------------

psi_p_data <- lapply (list (sceDS,
                            sceDSsmooth,
                            sce0,
                            sce1,
                            sce2,
                            sce3,
                            scePhen,
                            scePhenSpot), function (i) {
                              
    
    # if three regression coefficients in detection
    if (nrow(i$alpha.mean.samples) ==3) {
    
    # bind data
    psi_p_data <-  rbind (
              
            data.frame (# occupancy intercept and coefficients
                        sim=melt(i$beta.mean.samples[1,,])[,1],
                        scenario=melt(i$beta.mean.samples[1,,])[,2],
                        beta0=melt(i$beta.mean.samples[1,,])[,3],
                        beta0.lwr=melt(i$beta.low.samples[1,,])[,3],
                        beta0.upr=melt(i$beta.high.samples[1,,])[,3],
                        beta1=melt(i$beta.mean.samples[2,,])[,3],
                        beta1.lwr= melt(i$beta.low.samples[2,,])[,3],
                        beta1.upr=melt(i$beta.high.samples[2,,])[,3],
                        
                        # detection intercept and coefficients
                        alpha0= melt(i$alpha.mean.samples[1,,])[,3],
                        alpha0.lwr= melt(i$alpha.low.samples[1,,])[,3],
                        alpha0.upr= melt(i$alpha.high.samples[1,,])[,3],
                        alpha1= melt(i$alpha.mean.samples[2,,])[,3],
                        alpha1.lwr= melt(i$alpha.low.samples[2,,])[,3],
                        alpha1.upr= melt(i$alpha.high.samples[2,,])[,3],
                        alpha2= melt(i$alpha.mean.samples[3,,])[,3],
                        alpha2.lwr= melt(i$alpha.low.samples[3,,])[,3],
                        alpha2.upr= melt(i$alpha.high.samples[3,,])[,3],
                        
                        study=i$study,
                        sc=i$sc
            )) %>%
  
  right_join(scenario.vals,by="scenario")
    
    } else {
      
      # bind data
      psi_p_data <-  rbind (
        
        data.frame (# occupancy intercept and coefficients
          sim=melt(i$beta.mean.samples[1,,])[,1],
          scenario=melt(i$beta.mean.samples[1,,])[,2],
          beta0=melt(i$beta.mean.samples[1,,])[,3],
          beta0.lwr=melt(i$beta.low.samples[1,,])[,3],
          beta0.upr=melt(i$beta.high.samples[1,,])[,3],
          beta1=melt(i$beta.mean.samples[2,,])[,3],
          beta1.lwr= melt(i$beta.low.samples[2,,])[,3],
          beta1.upr=melt(i$beta.high.samples[2,,])[,3],
          
          # detection intercept and coefficients
          alpha0= melt(i$alpha.mean.samples[1,,])[,3],
          alpha0.lwr= melt(i$alpha.low.samples[1,,])[,3],
          alpha0.upr= melt(i$alpha.high.samples[1,,])[,3],
          alpha1= melt(i$alpha.mean.samples[2,,])[,3],
          alpha1.lwr= melt(i$alpha.low.samples[2,,])[,3],
          alpha1.upr= melt(i$alpha.high.samples[2,,])[,3],
          alpha2= NA,
          alpha2.lwr= NA,
          alpha2.upr= NA,
          
          study=i$study,
          sc=i$sc
        )) %>%
        
        right_join(scenario.vals,by="scenario")
      
    }
    
    ;
    
    psi_p_data

})

# melt to a df
psi_p_data <- do.call(rbind,psi_p_data)
  
# Ensure time and spatial columns are character vectors
psi_p_data$spatial <- as.character(psi_p_data$spatial)
psi_p_data$time <- as.character(psi_p_data$time)

# Convert time and spatial columns to factors for proper faceting
psi_p_data$spatial <- factor(psi_p_data$spatial, levels = unique(psi_p_data$spatial))
psi_p_data$time <- factor(psi_p_data$time, levels = unique(psi_p_data$time))

# density
density_intercepts1 <- psi_p_data %>%
  filter(scenario == c(1:4) &
           is.na(alpha0)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
      ggplot( aes(x = plogis(beta0), y = plogis(alpha0))) +
      geom_density_2d_filled(contour_var = "ndensity",bins=6)+
               geom_point(alpha=0.1) + 
               scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote("logistic"*(hat(alpha)[0])*""),
       x = bquote("logistic"*(hat(beta)[0])*""))+
    geom_point (aes(x=plogis(beta[1]),y=plogis(alpha[1])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_intercepts_sc1-4.png"),
     width =18, height = 18,units = "cm",res=150)

  density_intercepts1

dev.off()

# density
density_intercepts2 <- psi_p_data %>%
  filter(scenario == c(13:16) &
           is.na(alpha0)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = plogis(beta0), y = plogis(alpha0))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote("logistic"*(hat(alpha)[0])*""),
       x = bquote("logistic"*(hat(beta)[0])*""))+
  geom_point (aes(x=plogis(beta[1]),y=plogis(alpha[1])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_intercepts_sc13-16.png"),
     width =18, height = 18,units = "cm",res=150)

  density_intercepts2

dev.off()

# regression coeff --------------------------------------
# density
density_betas1 <- psi_p_data %>%
  filter(scenario == c(1:4) &
           is.na(beta0)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = plogis(beta0), y = (beta1))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  #ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote(""*hat(beta)[1]*""),
       x = bquote("logistic"*(hat(beta)[0])*""))+
  geom_point (aes(x=plogis(beta[1]),y=(beta[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_betas_sc1-4.png"),
     width =18, height = 18,units = "cm",res=150)


  density_betas1

dev.off()

# other scenarios
density_betas2 <- psi_p_data %>%
  filter(scenario == c(13:16) &
           is.na(alpha0)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = plogis(beta0), y = (beta1))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  #ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote(""*hat(beta)[1]*""),
       x = bquote("logistic"*(hat(beta)[0])*""))+
  geom_point (aes(x=plogis(beta[1]),y=(beta[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_betas_sc13-16.png"),
     width = 18, height = 18, units = "cm",res=150)

  density_betas2

dev.off()

# alphas ------------------------------
# density
density_alphas1 <- psi_p_data %>%
  filter(scenario == c(1:4) &
           is.na(beta0)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = plogis(alpha0), y = (alpha1))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  #ylim(c(0,1))+
  theme_bw()+
  labs(y = bquote(""*hat(alpha)[1]*""),
       x = bquote("logistic"*(hat(alpha)[0])*""))+
  geom_point (aes(x=plogis(alpha[1]),y=(alpha[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_alphas_sc1-4.png"),
     width = 18, height = 18, units = "cm",res=150)

  density_alphas1

dev.off()

# other sncearios
density_alphas2 <- psi_p_data %>%
  filter(scenario == c(13:16) &
           is.na(alpha0)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = plogis(alpha0), y = (alpha1))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  #ylim(c(-1,0))+
  theme_bw()+
  labs(y = bquote(""*hat(alpha)[1]*""),
       x = bquote("logistic"*(hat(alpha)[0])*""))+
  geom_point (aes(x=plogis(alpha[1]),y=(alpha[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_alphas_sc13-16.png"),
     width = 18, height = 18, units = "cm",res=150)

  density_alphas2

dev.off()

# second alpha coefficient ---------------------
# density
density_alphas2 <- psi_p_data %>%
  filter(scenario == c(1:4) &
           is.na(alpha2)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = plogis(alpha0), y = (alpha2))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  ylim(c(-1,0))+
  theme_bw()+
  labs(y = bquote(""*hat(alpha)[2]*""),
       x = bquote("logistic"*(hat(alpha)[0])*""))+
  geom_point (aes(x=plogis(alpha[1]),y=(alpha[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_alphas2_sc1-4.png"),
     width = 18, height = 12, units = "cm",res=150)

  density_alphas2

dev.off()

# other sncearios
density_alphas2_b <- psi_p_data %>%
  filter(scenario == c(13:16) &
           is.na(alpha2)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = plogis(alpha0), y = (alpha2))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  xlim(c(0,1))+
  ylim(c(-1,0))+
  theme_bw()+
  labs(y = bquote(""*hat(alpha)[2]*""),
       x = bquote("logistic"*(hat(alpha)[0])*""))+
  geom_point (aes(x=plogis(alpha[1]),y=(alpha[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_alphas2_sc13-16.png"),
     width = 18, height = 12, units = "cm",res=150)

  density_alphas2_b

dev.off()

# alpha 1 and alpha2 ---- low temporal autocorrelation
# other sncearios
density_alphas <- psi_p_data %>%
  filter(scenario == c(1:4) &
           is.na(alpha2)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = (alpha1), y = (alpha2))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  #xlim(c(0,1))+
  ylim(c(-1,0))+
  theme_bw()+
  labs(y = bquote(""*hat(alpha)[2]*""),
       x = bquote(""*hat(alpha)[1]*""))+
  geom_point (aes(x=(alpha[2]),y=(alpha[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_alphas_pairs_sc1-4.png"),
     width = 18, height = 12, units = "cm",res=150)

  density_alphas

dev.off()


# alpha 1 and alpha2 ---- high temporal autocorrelation
# other sncearios
density_alphas <- psi_p_data %>%
  filter(scenario == c(13:16) &
           is.na(alpha2)!=T ) %>%
  mutate (spatial=as.factor(spatial)) %>%
  
  ggplot( aes(x = (alpha1), y = (alpha2))) +
  geom_density_2d_filled(contour_var = "ndensity",bins=6)+
  geom_point(alpha=0.1) + 
  scale_fill_viridis_d(direction = -1,option="magma") +
  facet_grid(study+sc ~ spatial, 
             labeller = label_parsed) +
  #xlim(c(0,1))+
  ylim(c(-1,0))+
  theme_bw()+
  labs(y = bquote(""*hat(alpha)[2]*""),
       x = bquote(""*hat(alpha)[1]*""))+
  geom_point (aes(x=(alpha[2]),y=(alpha[2])),col="gray80",shape=4)+
  theme(axis.text.x = element_text(angle=45))

# save
png (here("figures", "Scenario-Comparisons","density_alphas_pairs_sc13-16.png"),
     width = 18, height = 12, units = "cm",res=150)

  density_alphas

dev.off()

# build the likelihood surfaces  --------------------

# https://stackoverflow.com/questions/13613157/plot-3d-density
require(plotly)
library(MASS)

# all random effects high
test_sc2 <- theta_data %>%
  #filter (is.na(alpha0)!=T) %>%
  filter (study == 1 & sc == 0 & scenario == 2) 
den3d <- kde2d((test_sc2$sigma_sq), 
               (test_sc2$phi.x))

# plot 
plot_surface2<-plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% 
  add_surface()
plot_surface2
htmlwidgets::saveWidget(as_widget(plot_surface2), here("figures", "Scenario-Comparisons","st1-sc0_all_high.html"))
