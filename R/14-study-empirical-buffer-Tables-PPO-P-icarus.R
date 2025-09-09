

# -----------------------------------------------


#     Interpretation Empirical data analyses - P icarus
#     Analysis within the buffer of 10 km^2 around Bordeaux

#     Prior-posterior overlap (Results section)
#     Table with coefficients (Supporting Information)


# -----------------------------------------------

# load packages
rm(list=ls())
gc()
source ("R/packages.R")

# load output (batch) ----------------------------------
list_output <- list.files (here ("model_output", "empirical"),pattern = "Polyommatus*") 
list_output <- list_output[grep("urban", list_output)]

# load water  -------------------------------
load (file = here ("Processed_data",
                   "Water.RData"))

# Create Table of coefficients / prior posterior overlap assessment
#load package
require(MCMCvis) # prior posterior overlap
require(posterior)
table_summ <- lapply(list_output, function (ng){
      
      # load output
      load (file=here ("model_output", "empirical",ng))
  
      res_output <- rbind (
          # regression coefs
          t(rbind (apply(out$beta.samples,2,mean), 
                   apply(out$beta.samples,2,quantile,c(0.025,0.975)),
                   out$rhat$beta 
              )
            ),
        
          # random effects 
          t (rbind (apply(out$theta.samples,2,mean),
                    apply(out$theta.samples,2,quantile,c(0.025,0.975)),
                    out$rhat$theta
                    )
             ),
          # detection coefs
          t (rbind (apply(out$alpha.samples,2,mean),
                    apply(out$alpha.samples,2,quantile,c(0.025,0.975)),
                    out$rhat$alpha
                    )
                    )
             )
        
      # identify output
      res_output <- data.frame (res_output)
      res_output <-cbind (res_output,
                          model=ng)
      colnames(res_output) <- c("Average", "LCI","UCI","RHat","Model")
      
      # prior posterior overlap ------------------------------------
      
      #run the function for just beta parameters
      # occupancy coeffs
      beta_overlap <- MCMCtrace(out$beta.samples,
                priors = rnorm(n = 1200,0,sqrt(2.72)), 
                PPO_out = TRUE,
                plot = T,
                pdf=F,
                type = "density")
      
      # detection coeffs
      alpha_overlap <- MCMCtrace(out$alpha.samples,
                                priors = rnorm(n = 1200,0,sqrt(2.72)), 
                                PPO_out = TRUE,
                                plot = T,
                                pdf=F,
                                type = "density")
      
      # theta coeffs
      # sigma_sq
      sigma_sq_overlap <- MCMCtrace(out$theta.samples,
                                 priors = 1/rgamma (n = 1200,shape = 2, rate = 1.5), # inverse gamma prior
                                 PPO_out = TRUE,
                                 plot = T,
                                 pdf=F,
                                 type = "density")
      sigma_sq_overlap <- sigma_sq_overlap[which(sigma_sq_overlap$param == "sigma.sq"),] # select the parameter (this prior don't apply to all theta)
      # sigma_sq_T
      sigma_sq_T_overlap <- MCMCtrace(out$theta.samples,
                                    priors = 1/rgamma (n = 1200,shape = 2, rate = 1), # inverse gamma prior 
                                    PPO_out = TRUE,
                                    plot = T,
                                    pdf=F,
                                    type = "density")
      sigma_sq_T_overlap <- sigma_sq_T_overlap[which(sigma_sq_T_overlap$param == "sigma.sq.t"),] # select the parameter (this prior don't apply to all theta)
      
      
      # phi  -  be flexible to consider different priors (weakly informative and informative)
      if (length(grep("InfPrior", ng))==0) {
        
        phi_overlap <- MCMCtrace(out$theta.samples,
                                        priors = runif(1200, 3 / 1, 3 / 0.005),
                                        PPO_out = TRUE,
                                        plot = T,
                                        pdf=F,
                                        type = "density")
        phi_overlap <- phi_overlap[which(phi_overlap$param == "phi"),] # select the parameter (this prior don't apply to all theta)
        
      } else {
        
        phi_overlap <- MCMCtrace(out$theta.samples,
                                 priors = runif(1200, 3 / 6, 3 / 1),
                                 PPO_out = TRUE,
                                 plot = T,
                                 pdf=F,
                                 type = "density")
        phi_overlap <- phi_overlap[which(phi_overlap$param == "phi"),] # select the parameter (this prior don't apply to all theta)
        
        
        
      }
      
      # rho
      rho_overlap <- MCMCtrace(out$theta.samples,
                               priors = runif(1200, -1, 1),
                               PPO_out = TRUE,
                               plot = T,
                               pdf=F,
                               type = "density")
      rho_overlap <- rho_overlap[which(rho_overlap$param == "rho"),] # select the parameter (this prior don't apply to all theta)
      
      # bind results in a single data frame
      res_output <- cbind (res_output,
             
             rbind (beta_overlap,
                   sigma_sq_overlap,
                   phi_overlap,
                   sigma_sq_T_overlap,
                   rho_overlap,
                   alpha_overlap)
      )
      
      # set expression names
      res_output$param <- c("$\\beta_0$",
                             "$\\beta_{LAT}$",
                             "$\\beta_{LAT^2}$",
                             "$\\beta_{LON}$",
                             "$\\beta_{ALT}$",
                              "$\\beta_{NON-WATER}$",
                            "$\\beta_{NON-WATER^2}$",
                            "$\\beta_{URBAN}$",
                            "$\\sigma^2$", 
                             "$\\phi$",
                             "$\\sigma^2_T$",
                             "$\\rho$", 
                             "$\\alpha_0$",
                             "$\\alpha_{LAT}$",
                             "$\\alpha_{OBS}$",
                            "$\\alpha_{MON}$",
                            "$\\alpha_{MON^2}$",
                            "$\\alpha_{NON-WATER}$"
      )
      
      
      rm(out) # save memory
      ; # return
      
      res_output
      
}
)

# melt
table_summ <- do.call(rbind, table_summ)
rownames(table_summ) <- NULL

# format table  to show in the main text
tab_kable<-table_summ %>%
  dplyr::select (param, Average, LCI,UCI,RHat,percent_PPO, Model) %>%
  # results with 15 neighbors
  filter (Model == "buffer_outputNNG15_urban_Polyommatus .Rdata") %>%
  # bind results with informative priors
  cbind (table_summ %>%
           dplyr::select (Average, LCI,UCI,RHat,percent_PPO, Model) %>%
           filter (Model == "buffer_outputNNG15InfPrior_urban_Polyommatus .Rdata"))

# do the table
knitr::kable(tab_kable [,-which(colnames(tab_kable) == "Model")],
             align="c",
             format="latex",
             digits = 3, 
             col.names = c("Parameter","Mean","LCI","UCI","RHat","Overlap","Mean","LCI","UCI","RHat","Overlap"),
             caption = "Estimated parameters of the multi-season site-occupancy model fitted to P. icarus data.", booktabs = T) %>%
  add_header_above(c("-"=1,"Weak prior, NG=15" = 5, "Informative prior, NG=15" = 5))


# end


