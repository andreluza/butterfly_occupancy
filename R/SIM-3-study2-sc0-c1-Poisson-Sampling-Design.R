# ---------------------------------------------------------------

# Sampling design for simulations (with J=10)

# Using Poisson distribution to find the number of secondary occasions per site and year; then the function 'sample' with Uniform probabilities will be used to distribute J secondary occasions to sites and years. The issues with sample function used without replacement and with different probability inclusion (see https://larryzhangnz.netlify.app/post/note-on-sample/ ; https://freerangestats.info/blog/2024/08/31/ppswor) don't apply here. Thus, with uniform probabilities we might be better off using sample() without replacement, if it is well-defined for uniform probabilities. This seems to be the case:

### Testing uniform sample() without replacement if you want

# sample_mat = matrix (NA,nrow=100000,ncol=5)
# for (i in 1:100000){sample_mat[i,]=sample(1:10,5)}
# table(sample_mat)

# ---------------------------------------------------------------

# load packages & functions
rm(list=ls())
source ("R/packages.R")
source ("R/functions.R")


# directory to store figures
dir.create(here ("figures", "sims_present_paper"))

# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# ----------------------------------------------------------------------
# implement our Poisson-based design
# Poisson distribution to set surveys to sites
# representing the Poisson sampling design
require(ggplot2)
require(dplyr)

# set lambda
J<-10
lambda<-1+0.1 # average sampling intensity

# our design
p1 <- data.frame (p_SO = dpois(seq(0,J,1), lambda=lambda),
                  SO = seq(0,J,1)) %>%
  ggplot(aes(x=SO,y=p_SO))+
  geom_col()+
  geom_smooth(aes(x=SO,y=p_SO),method = "glm", formula=y ~ splines::ns(x, 2),
              se = F, 
              method.args = list(family = "poisson"),
              linetype = "dashed",col="orange")+
  theme_minimal()+
  labs (x="Secondary periods",y="Probability")+
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=11))+
  annotate(geom="text",x=4,y=0.3,label = bquote(""*lambda*"=1.1"),size=8)

# save
png(here ("figures", "sims_present_paper",  "poisson_data_design.png"))
p1
dev.off()

# built our empirically motivated data design ---------------------------

set.seed(1234)
D_it <- matrix (NA, nrow=I,ncol=n.time)
for (t in 1:n.time)  {
  
  D_it[,t] <- rpois(n=I, lambda=1.1) # sites in D&S could receive either 1 or 2 visits (p_J = 1 + 0,1)
  
}

# sample : distribute these D_it surveys across J secondary occasions
# create an array of sites, years, secondary occasions
G_itj <- array (0, 
                dim = c(I,n.time,J)
)

# vector of uniform probabilities
p_J <-  rep(1,J)/J

# for each site and year
for (i in 1:I) {
  for (t in 1:n.time) {
    
    # Uses the Brewer's method to select a sample of units (unequal probabilities, without replacement, fixed sample size). 
    # set 0 if the site did not receive any survey month in Poisson distribution
    if (D_it[i,t] == 0) {
      
      G_itj[i,t,] <- 0
      
    } else {
      
      # run sampling when D_it>0        
      SecOcc <- sample (x=seq(1,J),
                               size=D_it[i,t],
                               replace = F,
                               prob = p_J)
      
      #SecOcc<-SecOcc[is.na(SecOcc)!=T]
      # set 1 if the site was sampled in each J
      G_itj[i,t,SecOcc] <- 1 
      
    } # close ifelse
  }
}

# distribution of monthly surveys
require(dplyr)
require(ggplot2)
data.frame (table(D_it[,1])) %>%
  ggplot(aes(x=Var1,y=Freq))+
  geom_col()+
  theme_minimal()+
  labs (x="Secondary occasions",y="Probability")
# any spatial gap?
range(rowSums(D_it>0))

# plot
data.frame (table(G_itj[,1,1])) %>%
  ggplot(aes(x=Var1,y=Freq))+
  geom_col()+
  theme_minimal()+
  labs (x="Secondary occasions",y="Probability")

# check
table(rowSums(D_it) == apply(G_itj>0,1,sum,na.rm=T))

# plot
dist_surveys_sites <-  reshape::melt(D_it) %>%
  dplyr::rename("Site" = "X1",
                "Year" = "X2",
                "Nsurveys" = "value") %>%
  mutate("Year" = paste0("t=",Year))%>%
  right_join (coords, by="Site")  %>%
  # and plot
  filter(Year == "t=1") %>% 
  ggplot(aes(Var1, Var2)) +
  geom_tile(aes(fill = Nsurveys)) +
  labs(x = "X", y = "Y")+
  #facet_wrap(~Year, ncol = 1, scales = "fixed") +
  my_theme+
  theme(legend.position = "none",
        strip.background = element_rect(colour = "orange4"),
        panel.background =  element_rect(colour = "orange4")) +
  scale_fill_gradient(low = "white", high = "orange4") +
  geom_text(aes(label = ifelse(Nsurveys==0,NA,Nsurveys)), size = 2)

dist_surveys_sites

# reshape to plot
dist_surveys_sites_occ <- reshape::melt(G_itj)%>%
  dplyr::rename("Site" = "X1",
                "Year" = "X2",
                "SecOcc" = "X3",
                "Nsurveys" = "value") %>%
  mutate("Year" = paste0("t=",Year))%>%
  mutate("SecOcc" = paste0("j=",SecOcc))%>%
  right_join (coords, by="Site")  %>% 
  # and plot
  filter(Year == "t=1") %>% 
  ggplot(aes(Var1, Var2)) +
  geom_tile(aes(fill = Nsurveys)) +
  labs(x = "X", y = "Y")+
  geom_text(aes(label = ifelse(Nsurveys==0,NA,Nsurveys)), size = 2)+
  
  facet_wrap(~SecOcc, nrow = 2, scales = "fixed") +
  my_theme+
  theme(legend.position = "none",
        strip.background = element_rect(colour = "orange4"),
        panel.background =  element_rect(colour = "orange4")) +
  scale_fill_gradient(low = "white", high = "orange4") 

dist_surveys_sites_occ

# save the design to be used in the other scenarios
save (D_it,G_itj,
      file=here ("model_output", 
                 "output_simulations", 
      "sampling_design_Poisson.rda")
      )
      
# end
rm(list=ls())





