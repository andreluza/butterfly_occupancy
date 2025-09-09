# ---------------------------------------------------------------

# Sampling design with phenology + observer sampling effects for simulations (with J=10)

# The issue about the 'sample' function using non-uniform probabilities applies here (see https://larryzhangnz.netlify.app/post/note-on-sample/ ; https://freerangestats.info/blog/2024/08/31/ppswor). As no established technique is available up to this date to solve the problem, we used a noisy Gaussian curve and imposed a threshold on it to mimic phenology + observer sampling effects (mid season preferences) on the observation data.

# Nota bene: the key to the phenology + observer sampling effects that we propose is that you first render the phenology of the year/site a bit noisy (some times you will select values outside the center ones, though not often), through a vector of "seasonal intensity of sampling", and then you just take the indices corresponding to the k largest values of that sample, where k is the number of values that you want to sample.
# Here k is the value that you randomly draw with D_it, which follows the Poisson distribution. 

# ---------------------------------------------------------------

# load packages & functions
rm(list=ls())
source ("R/packages.R")
source ("R/functions.R")

# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# load the previous design to have the same D_it
load (file=here ("model_output", 
                 "output_simulations", 
                 "sampling_design_Poisson.rda")
)
rm(G_itj) # rm as we will create a new one

# ----------------------------------------------------------------------
# Redesign our Poisson-based design
# Poisson distribution to set surveys to sites

set.seed(1234)

# sample : distribute the D_it surveys across J secondary occasions
# create an array of sites, years, secondary occasions
G_itj <- array (0, 
                dim = c(I,n.time,J)
)

# phenology  ------------------------------------
# function to create a subtly noisy Gaussian-like curve (so that sometimes the core months won't be sampled) 
seasonal_effect <- function(j, m, s) { 
  
  p<-exp(-((j-m)/s)^2) # gaussian curve
  q<-p*exp(rnorm(J,0,0.33)) # noisy gaussian-like curve
  res<-list(p=p,q=q)
  return(res)
  
}

# illustrate the effect (4 runs)
par(mfrow=c(2,2))
lapply (seq(1,4), function (i){
  
  plot(seasonal_effect(j=1:J,m=J/2,s=J/4)[[2]],type="b",xlab="Secondary occasion",ylab="Sampling probability",ylim=c(0,2))
  points(seasonal_effect(j=1:J,m=J/2,s=J/4)[[1]],pch=19,col="gray")
  lines(seasonal_effect(j=1:J,m=J/2,s=J/4)[[1]],col="gray")
  

})

# for each site and year
for (i in 1:I) {
  for (t in 1:n.time) {
    
    # set 0 if the site did not receive any survey month in Poisson distribution
    if (D_it[i,t] == 0) {
      
      G_itj[i,t,] <- 0
      
    } else {
    # run sampling\thresholding when D_it>0        
      
      # create the vector of probs for that year
      p_J <- seasonal_effect(j=1:J,m=J/2,s=J/4)
      
      # thresholding q to select the sampled months
      k <- D_it[i,t]
      selected_months<-order(p_J$q,decreasing=TRUE)[1:k]
      
      # set 1 if the site was sampled in each J
      G_itj[i,t,selected_months] <- 1 
    
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

# plot
data.frame (table(G_itj[,1,7])) %>%
  ggplot(aes(x=Var1,y=Freq))+
  geom_col()+
  theme_minimal()+
  labs (x="Secondary occasions",y="Probability")

# check
table(rowSums(D_it) == apply(G_itj>0,1,sum,na.rm=T)) 

# spatial gaps?
range(rowSums(D_it>0)) # no

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
  mutate("SecOccB" = as.character(paste0("j=",SecOcc)))%>%
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
      file=here ("model_output", "output_simulations", 
                 "sampling_design_Poisson_phenology.rda")
)

# end
rm(list=ls())





