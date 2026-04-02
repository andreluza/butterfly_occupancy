# ---------------------------------------------------------------

# 3 -3 restablishing the amount of data

# Since the amount of data declined from 3.1 to 3.2, we need to restablish it

# Thus, we will spread surveys using Poisson distribution with \lambda=1.1 * (I/(I*0.25)), so that we will have 4 times more surveys in the spot to have similar amount of data between 3 - 1 and 3 - 3

# ---------------------------------------------------------------

# load packages & functions
rm(list=ls())
source ("R/packages.R")
source ("R/functions.R")

# load simulation settings
load(here("model_output", "output_simulations", "sim-settings.RData"))

# phenology function  ------------------------------------
# function to create a subtly noisy Gaussian-like curve (so that sometimes the core months won't be sampled) 
seasonal_effect <- function(j, m, s) { 
  
  p<-exp(-((j-m)/s)^2) # gaussian curve
  q<-p*exp(rnorm(J,0,0.33)) # noisy gaussian-like curve
  res<-list(p=p,q=q)
  return(res)
  
}

# ----------------------------------------------------------------------
# implement our Poisson-based design
# Poisson distribution to set surveys to sites

set.seed(1234)
D_it <- matrix (NA, nrow=I,ncol=n.time)
for (t in 1:n.time)  {

  # create the vector of site probs for that year
  # site probs
  p_T <- seasonal_effect(j=1:I,m=I/2,s=I/2)
    
  # thresholding q to select 50% of the sampled sites
  selected_sites<-order(p_T$q,decreasing=TRUE)[1:(I*0.25)]
  unsampled_sites <- seq(1,I)[which(seq(1,I) %in% selected_sites == F)]
    
  # set 0 if the site is not among the sampled sites in each year
  lambda_vector <- rep(1.1 * (I/(I*0.25)), I) # value of lambda for sampled sites # sites in D&S could receive 1 to J visits (p_J = 1 + 0,1)
  lambda_vector [unsampled_sites] <- 0 # value of lambda for unsampled sites
  
  # distributing visits across the sampled sites
  D_it[,t] <- rpois(n=I, lambda=lambda_vector) 
  
}
sum(D_it)

# Now we can define which sites were sampled using the phenology function
# illustrate the effect (4 runs)
par(mfrow=c(2,2))
lapply (seq(1,4), function (i){
  
  plot(seasonal_effect(j=1:I, m=I/2, s=I/2)[[2]],type="b",xlab="Site",ylab="Sampling probability",ylim=c(0,2))
  points(seasonal_effect(j=1:I, m=I/2, s=I/2)[[1]],pch=19,col="gray")
  lines(seasonal_effect(j=1:I, m=I/2, s=I/2)[[1]],col="gray")
  
  
})

# then we include the seasonal effect -----------------------------------

# sample : distribute these D_it surveys across J secondary occasions
# create an array of sites, years, secondary occasions
G_itj <- array (0, 
                dim = c(I,n.time,J)
)

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
sum(G_itj>0)
prod(dim(G_itj))

# spatial gaps?
range(rowSums(D_it>0)) # yes

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
  #my_theme+
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
  #my_theme+
  theme(legend.position = "none",
        strip.background = element_rect(colour = "orange4"),
        panel.background =  element_rect(colour = "orange4")) +
  scale_fill_gradient(low = "white", high = "orange4") 

dist_surveys_sites_occ

# save the design to be used in the other scenarios
save (D_it,G_itj,
      file=here ("model_output", "output_simulations", 
                 "sampling_design_Poisson_phenology_spot_st3_sc3.rda")
)

# end
rm(list=ls())
