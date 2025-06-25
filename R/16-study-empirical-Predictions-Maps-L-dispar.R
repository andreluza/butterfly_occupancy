

# --------------------------------

# Obtaining predictions from the model - Lycaena dispar
# Predictions for the whole Nouvelle-Aquitaine region, from 2009 to 2023 (average was used in the map)

# using spOccupancy package
# models with nngp=15
# weak prior for phi
# informative prior for phi

# helpful for make predictions in large-data settings
# https://groups.google.com/g/spocc-spabund-users/c/uoRUji5vkUQ
# https://groups.google.com/g/spocc-spabund-users/c/Q_A4V15NYSk/m/AzfwhRZ1AwAJ?utm_medium=email&utm_source=footer
# https://01503778298842392382.googlegroups.com/attach/1de27213a332e/prediction-split-example.R?part=0.1&view=1&vt=ANaJVrGTulV6MlIJKkIE1r2GFaI5SaE1HhV0DaDa5e9Ioaj07sEiwzeZ0sOpegw2B6860tciX_uGKPoViSTuwQYw-vUcEcnYHW8TEU57r6gQsSU-vWKxpfY

# --------------------------------------

# load packages
rm(list=ls())
gc()
source ("R/packages.R")

# ggplot theme
my_theme <- theme(legend.position = 'bottom', 
                  strip.text = element_text(size=6),
                  strip.text.y = element_text(color = 'black'),
                  strip.text.x = element_text(color = 'black'), 
                  legend.text = element_text(size=6),
                  text = element_text(family="LM Roman 10"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle = 45, hjust = 1, size = 5), 
                  axis.text.y = element_text(size = 5),
                  axis.title = element_text(size=10))

# load processed data ------------------------
load (file = here("Processed_data", 
                  "Occupancy_data_spOccupancy.RData"))

# load altitude and habitat  -------------------------------
load (file = here ("Processed_data",
                   "AltitudeHabitatStats.RData"))

# load water  -------------------------------
load (file = here ("Processed_data",
                   "Water.RData"))

# load spatial data -------------------------------
cells_NAquitane <- st_read(dsn=here ("Data", "SpatialData","Maillage_1x1km"),
                           layer="1x1km_n-a")

# select 15 years of data
sp_table_years <- sp_table_years[,(ncol(sp_table_years)-14):ncol(sp_table_years),,]
base_table_years <- base_table_years[,(ncol(base_table_years)-14):ncol(base_table_years),]

# where NAs in table of spp, there is NAs in total records
table(is.na(sp_table_years[,,,6]) == (base_table_years==0))

# Organize Data for spOcc analysis ---------------
# Package all data into a list
# occ
X <- array(1, dim = c(nrow(sp_table_years), ncol(sp_table_years), 9)) # intercept + n covariates (lat, lat2, & long, etc )
X[, , 2] <-  scale(cell_centroid_df[, "Y"])[,1]
X[, , 3] <-  (scale(cell_centroid_df[, "Y"]^2)[,1])
X[, , 4] <-  scale(cell_centroid_df[, "X"])[,1]
X[, , 5] <-  scale(altitude_stats$altitude_mean)[,1]
X[, , 6] <-  scale(1-extract.water.summary[,"1"])[,1] # non water = 1 minus water(extract.water.summary)
X[, , 7] <-  (scale(1-extract.water.summary[,"1"]^2)[,1])
X[, , 8] <-  scale(rowSums(extract.hab.NA.summary [,c("35","36", "37","40","41","42","43")]))[,1] # codes below
X[, , 9] <-  (scale(rowSums(extract.hab.NA.summary [,c("35","36", "37","40","41","42","43")])^2)[,1]) # codes below
#35 166 166 255 255 411 - Inland marshes  <-----------------------
#36 77 77 255 255 412 - Peat bogs <-----------------------
#37 204 204 255 255 421 - Salt marshes  <-----------------------
#40 0 204 242 255 511 - Water courses
#41 128 242 230 255 512 - Water bodies
#42 0 255 166 255 521 - Coastal lagoons
#43 166 255 230 255 522 - Estuaries

# imput altitude because there are 10 NAs
# set the minimum of the scaled values because NAs are in the lowlands/coastline
X[, , 5] [is.na(X[, , 5])] <- min(X[, , 5],na.rm=T)

# detection covariates
X.p <- array (1, dim=c(nrow(sp_table_years),
                       ncol(sp_table_years),
                       dim(sp_table_years)[3],
                       6)) # intercept + 4 covariates
X.p[,,,2] <- X[, , 2] # latitude
X.p[,,,3] <- scale((obs_table_years[,(ncol(obs_table_years)-14):ncol(obs_table_years),]))# observers 

# bind months (phenology)
months_bind  <- replicate (ncol(sp_table_years),
                           do.call(cbind,lapply (seq(1,10,1), function (i) rep(i,nrow(sp_table_years))))
)

# change the order to fit the right format
months_bind <- aperm(months_bind, c(1,3,2))
X.p[,,,4] <- scale(months_bind)
X.p[,,,5] <- scale(months_bind^2)

# non-water
X.p[,,,6] <- X[, , 6] # non-water habitats

# select the species
sp<- grep("Lycaena dispar", sp_list) # L. dispar

# identify missing data - zeroed encounter histories
missing_sites <- rowSums(is.na(sp_table_years[,,,sp]))==max(rowSums(is.na(sp_table_years[,,,sp]))) # 10 years x 10 months with zeroes 

# compare the missing sites between 
table(missing_sites) # sp detection table
table(rowSums(base_table_years[])==0) # and the effort table

gc()

# naive occupancy
mean(apply (sp_table_years[missing_sites==F,,,sp],2,sum,na.rm=T)/nrow(sp_table_years[missing_sites==F,,,]))*100

# load the models to  make predictions ----------------------------------------
NG15weak <- new.env()
NG15inf <- new.env()
load (file=here ("model_output", 
                 "empirical",
                 paste0 ("outputNNG15_", substr(sp_list[sp],1,12),".Rdata")),NG15weak)
NG15weak$lab <- "Ng15weak"
gc()
load (file=here ("model_output", 
                 "empirical",
                 "outputNNG15InfPrior_Lycaena disp.Rdata"),NG15inf)
NG15inf$lab <- "Ng15inf"

# make predictions 
gc()
list_output <- list (NG15weak,
                     NG15inf)

# generate predictions for all sites (predictions in groups of sites)
# save chunks
generate_predictions <- lapply (seq(1,length(list_output)), function (out) {
                             
                             # Do the predictions in chunks of sites #-------------------------------------------------------
                             vals <- split(1:nrow(cell_centroid_df), ceiling(seq_along(1:nrow(cell_centroid_df)) / 1500))
                             psi.quants <- array(NA, dim = c(nrow(cell_centroid_df), ncol(sp_table_years)))  
                             w.quants <- array(NA, dim = c(nrow(cell_centroid_df), 1)) #
                            for (j in 1:length(vals)) {
                               print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
                               curr.indx <- vals[[j]]
                               out.pred <- predict(list_output[[out]]$out, X[curr.indx, , -c(6,7), drop = FALSE],
                                                   cell_centroid_df[curr.indx, ],
                                                   t.cols = seq(1,ncol(sp_table_years)), 
                                                   n.omp.threads = 10, 
                                                   verbose = T)
                               psi.quants[curr.indx, ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
                               w.quants[curr.indx, ] <- apply(out.pred$w.0.samples, 2, mean)
                               
                               rm(out.pred)
                               gc()
                             }
                             
                             # save these predictions as they took some time to run
                             save (psi.quants,w.quants, 
                                   file = here ("model_output", "empirical", paste ("predictions-", list_output[[out]]$lab,"-", substr(sp_list[sp],1,12),".Rdata")))
                             
                             
                             }
                           )
                             

# load predictions
generate_predictions <- lapply (list("Ng15weak", "Ng15inf"), function (i) {

  load (file = here ("model_output", "empirical", paste ("predictions-", i,"-", substr(sp_list[sp],1,12),".Rdata")))

  #  return what interests
  res <- list (psi_it = psi.quants,
               w_i = w.quants)
  res
  }
)


# Produce plots for all models (weak and informative priors)
res_plots <- lapply (seq(1,length(generate_predictions)), function (out) {
  
    # estimates for non-sampled sites
    data_non_sampled <- data.frame (cells = rownames(cell_centroid_df [missing_sites==T,]),
                                    psi_i = apply (generate_predictions[[out]]$psi_it[missing_sites==T,],1, # expected values (average of predictions) of psi_i
                                                   mean),
                                    w_i = generate_predictions[[out]]$w_i [missing_sites==T,] # expected values (average of predictions) of w_i
                                                 )
    # estimates for sampled sites
    data_sampled <- data.frame (cells = rownames(cell_centroid_df [missing_sites==F,]),
                                psi_i = apply (list_output[[out]]$out$psi.samples,2,mean), # point estimates of psi_i
                                w_i = apply (list_output[[out]]$out$w.samples,2,mean)) # point estimates of omega_i
    # table(data_sampled$cells %in% data_non_sampled$cells) # check if they are all different
    gc()
    # sampled within Bordeaux's buffer
    #data_sampled <- data_sampled[which(data_sampled$cells %in% gironde_cells$CODE_10KM),]
    
    # bind the results
    res_data <- rbind (data_sampled,
                       data_non_sampled)
    
    table(data_sampled$cells %in% data_non_sampled$cells)
    
    # match cell ID
    res_data <- res_data[match (rownames(cell_centroid_df),res_data$cells),]
    # table(res_data$cells == gironde_cells$CODE_10KM) # check
    gc() # clean
    
    # map of non-sampled areas ----- out-of-sample predictions needed -------------
    # occupancy
    p1_missing<-ggplot(data = cells_NAquitane[missing_sites==T,]) +
      geom_sf(fill="white")+
      geom_sf(data=cbind (cells_NAquitane[missing_sites==T,],
                          bar_psi_i = res_data$psi_i[missing_sites==T]),
              aes (fill=bar_psi_i,
                   col=bar_psi_i),
              alpha=0.75)+
      scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(0,1))+
      scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(0,1))+
      ggtitle("")+
      theme(plot.title = element_text(size=10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 3), 
            axis.text.y = element_text(size = 3),
            strip.text.x = element_text(face = "italic",size=7),
            legend.key.size = unit(0.4, 'cm'), #change legend key size
            legend.background = element_blank(),
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=8), #change legend title font size
            legend.text = element_text(size=8,angle=45)) +
      labs (col = expression(paste('E(', hat(psi[i]),')')),
            fill = expression(paste('E(', hat(psi[i]),')')))+
      my_theme
    
    # spatial random effect
    p2_missing<-ggplot(data = cells_NAquitane[missing_sites==T,]) +
      geom_sf(fill="white")+
      geom_sf(data=cbind (cells_NAquitane[missing_sites==T,],
                          bar_w_i = res_data$w_i[missing_sites==T]),
              aes (fill=bar_w_i,
                   col=bar_w_i),
              alpha=0.75)+
      scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(min(res_data$w_i),max(res_data$w_i)))+
      scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(min(res_data$w_i),max(res_data$w_i)))+
      ggtitle("")+
      theme(plot.title = element_text(size=10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 3), 
            axis.text.y = element_text(size = 3),
            strip.text.x = element_text(face = "italic",size=7),
            legend.key.size = unit(0.4, 'cm'), #change legend key size
            legend.background = element_blank(),
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=8), #change legend title font size
            legend.text = element_text(size=8,angle=45)) +
      labs (col = expression(paste('E(', hat(omega[i]),')')),
            fill = expression(paste('E(', hat(omega[i]),')')))+
      my_theme
    
    # show predictions for all sites --------------------------
    # E(psi_i)
    p1 <- ggplot(data = cells_NAquitane) +
      geom_sf(fill="white")+
      geom_sf(data=cbind (cells_NAquitane,
                          bar_psi_i = res_data$psi_i), # what have data
              aes (fill=bar_psi_i,
                   col=bar_psi_i),
              alpha=0.75)+
      scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(0,1))+
      scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(0,1))+
      theme(plot.title = element_text(size=10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 3), 
            axis.text.y = element_text(size = 3),
            strip.text.x = element_text(face = "italic",size=7),
            legend.key.size = unit(0.4, 'cm'), #change legend key size
            legend.background = element_blank(),
            legend.key.height = unit(0.35, 'cm'), #change legend key height
            legend.key.width = unit(0.4, 'cm'), #change legend key width
            legend.title = element_text(size=8), #change legend title font size
            legend.text = element_text(size=8,angle=45))+
      labs (col = expression(paste('', y[i],'')),
            fill = expression(paste('E(', hat(psi[i]),')')))+
      my_theme+
      guides(col="none")
    
    # the spatial random effect
    p2 <- ggplot(data = cells_NAquitane) +
      geom_sf(fill="white")+
      geom_sf(data=cbind (cells_NAquitane,
                          bar_w_i = res_data$w_i), # what have data),
              aes (fill=bar_w_i,
                   col=bar_w_i),
              alpha=0.75)+
      scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(min(res_data$w_i),max(res_data$w_i)))+
      scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90",limits=c(min(res_data$w_i),max(res_data$w_i)))+
      ggtitle("")+
      theme(plot.title = element_text(size=10),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 3), 
            axis.text.y = element_text(size = 3),
            strip.text.x = element_text(face = "italic",size=7),
            legend.key.size = unit(0.4, 'cm'), #change legend key size
            legend.background = element_blank(),
            legend.key.height = unit(0.35, 'cm'), #change legend key height
            legend.key.width = unit(0.4, 'cm'), #change legend key width
            legend.title = element_text(size=8), #change legend title font size
            legend.text = element_text(size=8,angle=45))+
      labs (col = expression(paste('', y[i],'')),
            fill = expression(paste('E(', hat(omega[i]),')')))+
      my_theme+
      guides(col="none")
    
    # summarize in-sample estimates (point estimates)
    gc()
    summ_tab_psi <-apply(list_output[[out]]$out$psi.samples,c(1,3),mean)
    summ_tab_psi <- data.frame (psi = apply(summ_tab_psi,2,mean),
                     uci = apply(summ_tab_psi,2,quantile, 0.975),
                     lci = apply(summ_tab_psi,2,quantile, 0.025),
                     year = seq(2009,2023))
    gc()
    
    # plot
    p3<- summ_tab_psi %>%
      ggplot(aes(x=year,psi)) +
      geom_ribbon(aes(ymin=lci, ymax=uci),fill="white")+
      geom_line(linewidth=1,col="black")+
      geom_line(data = data.frame (psi = apply (sp_table_years[,,,sp],2,sum,na.rm=T)/nrow(sp_table_years[missing_sites==F,,,]),year = seq(2009,2023)), 
                                             aes (x=year, y=psi))+
        ggtitle("")+
        labs(x="Year", y = expression(paste('E(', hat(psi[t]),')')))+
        my_theme+
        ylim(c(0,1))
      
      # phenology
      # summarize
      par(mfrow=c(1,1),mar=c(4,4,4,4))
      point_est <- apply (list_output[[out]]$out$alpha.samples,2,mean)
      lci<-apply (list_output[[out]]$out$alpha.samples,2,quantile,0.025)
      uci<-apply (list_output[[out]]$out$alpha.samples,2,quantile,0.975)
      
      
      # organize data to plot
      phen_plot <- data.frame (p = plogis(point_est[1] + 
                                            point_est[5]*seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10)+
                                            point_est[6]*(seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10)^2)),
                               lci = plogis(lci[1] + 
                                              lci[5]*seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10)+
                                              lci[6]*(seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10)^2)),
                               uci = plogis(uci[1] + 
                                              uci[5]*seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10)+
                                              uci[6]*(seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10)^2)),
                               x = seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10),
                               month=c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")
                               
      ) 
      phen_plot$month <-factor(phen_plot$month,
                               levels = c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
      
      p4 <- ggplot(phen_plot, aes(x=x, y = p)) +
        geom_ribbon(aes(ymin=lci, ymax=uci),fill="white")+
        #geom_point()+
        geom_line(linewidth=1)+
        my_theme+ #theme(axis.ticks.x = element_blank())+
        labs(y = expression(paste('E(', hat(p[j]),')')),
             x="Month")+
        ggtitle("")+
        scale_x_continuous(
          breaks=seq(min(X.p[,,,5]),max(X.p[,,,5]),length.out=10),
          labels=c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")
        )
      # create a list of resulting plots
      res_plots <- list (p1=p1,
                         p2=p2,
                         p3=p3,
                         p4=p4,
                         p1_missing = p1_missing,
                         p2_missing = p2_missing
      )
      ;
      res_plots
      
})

# plot and save results
png (here ("figures", "empirical","Fig_summ_results_Ldispar_Aquitaine.png"),width = 25,height = 20,units = "cm",res=400)
  
  grid.arrange(# weak prior
    res_plots[[1]]$p1+ggtitle("A"),
    res_plots[[1]]$p2,
    res_plots[[1]]$p3,
    res_plots[[1]]$p4,
    
    # 15 neightbors informative
    res_plots[[2]]$p1+ggtitle("B"),
    res_plots[[2]]$p2,
    res_plots[[2]]$p3,
    res_plots[[2]]$p4,
    
  nrow=2,ncol=4
)

dev.off()

# plot occupancy and random effects for out of sample cells/sites
png (here ("figures", "empirical","estimates_missing_cells_Ldispar_Aquitaine.png"),width = 12,height = 20,units = "cm",res=400)
  
  grid.arrange(# weak prior
    res_plots[[1]]$p1_missing+ggtitle("A"),
    res_plots[[1]]$p2_missing,
    
    # 15 neighbors informative
    res_plots[[2]]$p1_missing+ggtitle("B"),
    res_plots[[2]]$p2_missing,
    
    nrow=2,ncol=2
  )

dev.off()

# end




