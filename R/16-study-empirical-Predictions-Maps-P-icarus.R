
# --------------------------------

# Obtaining predictions from the model - Polyommatus icarus
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

# create a dir to receive the results
dir.create(here("model_output", "empirical"))

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
# occ: lat+lat2+lon+elev+non_water+non_water2+ grasslands+grasslands2
X <- array(1, dim = c(nrow(sp_table_years), ncol(sp_table_years), 9)) # intercept + n covariates (lat, lat2, long, etc)
X[, , 2] <-  scale(cell_centroid_df[, "Y"])[,1]
X[, , 3] <-  (scale(cell_centroid_df[, "Y"])[,1])^2
X[, , 4] <-  scale(cell_centroid_df[, "X"])[,1]
X[, , 5] <-  scale(altitude_stats$altitude_mean)[,1]
X[, , 6] <-  scale(1-extract.water.summary[,"1"])[,1] # non water = 1 minus water(extract.water.summary)
X[, , 7] <-  (scale((1-(extract.water.summary[,"1"])^2))[,1])
X[, , 8] <-  scale((extract.hab.NA.summary [,c("26")]))[,1] # codes below
X[, , 9] <-  (scale((extract.hab.NA.summary [,c("26")]))[,1])^2 # codes below
#26 204 242 77 255 321 - Natural grasslands     <-----------------------
#27 166 255 128 255 322 - Moors and heathland   <-----------------------
#28 166 230 77 255 323 - Sclerophyllous vegetation 
#30 230 230 230 255 331 - Beaches - dunes - sands
#32 204 255 204 255 333 - Sparsely vegetated areas

# plot 
ggplot() +
  geom_sf(fill="white")+
  geom_sf(data= cells_NAquitane)+
  geom_sf(data = cbind (cells_NAquitane,
                        alt = X[, 1, 5]),
          aes(fill=alt,col=alt))+
  scale_fill_viridis_c(na.value = "red")+
  scale_colour_viridis_c(na.value = "red")

# plot vars
ggplot() +
  geom_sf(fill="white")+
  geom_sf(data= cells_NAquitane)+
  geom_sf(data = cbind (cells_NAquitane,
                        land = X[, 1, 6]),
          aes(fill=land,col=land))+
  scale_fill_viridis_c(na.value = "red")+
  scale_colour_viridis_c(na.value = "red")

# input altitude because there are 10 NAs
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

# change the order to fit the correct format
months_bind <- aperm(months_bind, c(1,3,2))
X.p[,,,4] <- scale(months_bind)[,1]
X.p[,,,5] <- (scale(months_bind)[,1])^2

# non-water
X.p[,,,6] <- X[, , 6] # non-water habitats

# select species ------------------------------
sp<- grep("Polyommatus icarus", sp_list)  # P. icarus

# identify missing data - zeroed encounter histories
missing_sites <- rowSums(is.na(sp_table_years[,,,sp]))==max(rowSums(is.na(sp_table_years[,,,sp]))) # 10 years x 10 months with zeroes 
table(missing_sites == (rowSums(base_table_years)==0)) # 10 years x 10 months with zeroes 

# naive occupancy
# calculate and save naive occupancy
sites_det_year <- apply (sp_table_years[,,,sp],c(1,2), sum,na.rm=T)>0 # if site was sampled in year t
sites_det_year <- colSums(sites_det_year) 
# year effort
eff_year <- colSums(apply (base_table_years,c(1,2), sum,na.rm=T)>0) # total number of sites sampled in each year
naive_occ <- data.frame (sites_det_year = sites_det_year,
                         total_sites_year = eff_year,
                         naive_occ = sites_det_year/eff_year)

# average naive occupancy
mean(naive_occ$naive_occ)*100
mean(naive_occ$sites_det_year)

gc()

# load the models to  make predictions ----------------------------------------
NG15weak <- new.env()
NG15inf <- new.env()
load (file=here ("model_output", 
                 "empirical",
                 paste0 ("outputNNG15_", substr(sp_list[sp],1,12),".Rdata")),NG15weak)
NG15weak$lab <- "Ng15weak"
load (file=here ("model_output", 
                 "empirical",
                 "outputNNG15InfPrior_Polyommatus .Rdata"),NG15inf)
NG15inf$lab <- "Ng15inf"

# make predictions 
gc()
list_output <- list (NG15weak,
                     NG15inf)

# generate predictions for all sites (predictions in groups of sites)
# save chunks
generate_predictions <- lapply (seq(1,length(list_output)), function (out) {
                             
                             # Do the predictions in chunks of sites -------------------------------------------------------
                             vals <- split(1:nrow(cell_centroid_df), ceiling(seq_along(1:nrow(cell_centroid_df)) / 1500))
                             psi.quants <- array(NA, dim = c(nrow(cell_centroid_df), ncol(sp_table_years))) # 
                             psi.quants.lci <- array(NA, dim = c(nrow(cell_centroid_df), ncol(sp_table_years)))  
                             psi.quants.uci <- array(NA, dim = c(nrow(cell_centroid_df), ncol(sp_table_years)))  
                             w.quants <- array(NA, dim = c(nrow(cell_centroid_df), 1)) #
                             for (j in 1:length(vals)) {
                               print(paste("Currently on set ", j, " out of ", length(vals), sep = ''))
                               curr.indx <- vals[[j]]
                               out.pred <- predict(list_output[[out]]$out, X[curr.indx, , , drop = FALSE],
                                                   cell_centroid_df[curr.indx, ],
                                                   t.cols = seq(1,ncol(sp_table_years)), 
                                                   n.omp.threads = 10, 
                                                   verbose = T)
                               psi.quants[curr.indx, ] <- apply(out.pred$psi.0.samples, c(2, 3), mean)
                               psi.quants.lci[curr.indx, ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.025)
                               psi.quants.uci[curr.indx, ] <- apply(out.pred$psi.0.samples, c(2, 3), quantile, 0.975)
                               w.quants[curr.indx, ] <- apply(out.pred$w.0.samples, 2, mean)
                               
                               rm(out.pred)
                               gc()
                             }
                             
                             # save these predictions as they took some time to run
                             save (psi.quants,psi.quants.uci,psi.quants.lci,w.quants, 
                                   file = here ("model_output", "empirical", paste ("predictions-", list_output[[out]]$lab,"-", substr(sp_list[sp],1,12),"#.Rdata")))
                             
                             }
                           )
                             
# load predictions
generate_predictions <- lapply (list("Ng15weak", "Ng15inf"), function (i) {

  load (file = here ("model_output", "empirical", paste ("predictions-", i,"-", substr(sp_list[sp],1,12),"#.Rdata")))

  #  return what interests
  res <- list (psi_it = psi.quants,
               psi_it_lci = psi.quants.lci,
               psi_it_uci = psi.quants.uci,
               
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
                              psi_i = apply (generate_predictions[[out]]$psi_it[missing_sites==F,],1,mean), # point estimates of psi_i
                              w_i = (generate_predictions[[out]]$w_i[missing_sites==F,])) # point estimates of omega_i
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
    scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
    scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
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
                        bar_psi_i = ifelse (X[, 1, 5] == min(X[, 1, 5] ), NA,res_data$psi_i)), # what have data
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
                        bar_w_i = ifelse (X[, 1, 5] == min(X[, 1, 5] ), NA, res_data$w_i)), # what have data),
            aes (fill=bar_w_i,
                 col=bar_w_i),
            alpha=0.75)+
    scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
    scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
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
  #summ_tab_psi <-apply(out.pred.occ$psi.0.samples,c(1,3),function (x) sum(x)/ncol(out$out$psi.samples))
  # summarize in-sample estimates (point estimates)
  #summ_tab_psib <-apply(list_output[[out]]$out$psi.samples,c(1,3),function (x) sum(x)/ncol(list_output[[out]]$out$psi.samples))
  #table(round(summ_tab_psib,3) == round(summ_tab_psi,3))
  
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
    geom_line(data = cbind (naive_occ, 
                            year = seq(2009,2023)), 
              aes (x=year, y=naive_occ))+
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
                                        point_est["phen"]*seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)+
                                        point_est["phen2"]*(seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)^2)),
                           lci = plogis(lci[1] + 
                                          lci["phen"]*seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)+
                                          lci["phen2"]*(seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)^2)),
                           uci = plogis(uci[1] + 
                                          uci["phen"]*seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)+
                                          uci["phen2"]*(seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)^2)),
                           x = seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100),
                           month=c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov")
                           
  ) 
  phen_plot$month <-factor(phen_plot$month,
                           levels = c("Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
  # plot phenology and observer effect
  p4 <- ggplot(phen_plot, aes(x=x, y = p)) +
    geom_ribbon(aes(ymin=lci, ymax=uci),fill="white")+
    #geom_point()+
    geom_line(linewidth=1)+
    my_theme+ #theme(axis.ticks.x = element_blank())+
    labs(y = expression(paste('E(', hat(p[j]),')')),
         x="Month")+
    ggtitle("")+
    scale_x_continuous(
      breaks=seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=10),
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
png (here ("figures", "empirical","Fig_summ_results_Picarus_Aquitaine.png"),width = 25,height = 20,units = "cm",res=400)

  grid.arrange(# weak prior
    res_plots[[1]]$p1+ggtitle("A"),
    res_plots[[1]]$p2,
    res_plots[[1]]$p3,
    res_plots[[1]]$p4,
    
    # 15 inf priors sigma
    res_plots[[2]]$p1+ggtitle("B"),
    res_plots[[2]]$p2,
    res_plots[[2]]$p3,
    res_plots[[2]]$p4,
    
    
  nrow=2,ncol=4
)

dev.off()

# plot occupancy and random effects for out of sample cells/sites
# plot and save results
png (here ("figures", "empirical","estimates_missing_cells_Picarus_Aquitaine.png"),width = 12,height = 20,units = "cm",res=400)

grid.arrange(# weak prior
  res_plots[[1]]$p1_missing+ggtitle("A"),
  res_plots[[1]]$p2_missing,
  
  # 30 neightbors
  res_plots[[2]]$p1_missing+ggtitle("B"),
  res_plots[[2]]$p2_missing,
  
  nrow=2,ncol=2
)

dev.off()

# end




