

# --------------------------------

#     Obtaining Predictions from the Model
#     Analysis within the buffer of 10 km^2 around Bordeaux
#     Modeling the distribution of P. icarus using data from 2000-2023
#     Code will produce plots shown in the main text and SI

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

# choosing a smaller scale for prediction --------------------------------------
# Gironde department
# source: https://www.actualitix.com/blog/shapefiles-des-departements-de-france.html
cells_Gironde <- st_read(dsn=here ("Data", "SpatialData","33-gironde"),
                         layer="33-")

a <- ggplot() +
  geom_sf(fill="white")+
  geom_sf(data= cells_NAquitane)

# bordeaux distance
bordeaux_distance <- st_distance (cells_Gironde,
                                  cells_Gironde %>%
                                    filter (NOM_COMM == "BORDEAUX"))

# intersection Gironde - NAq
cells_buffer_bordeaux <- (st_intersection(cells_Gironde %>%
                                            cbind (dist=as.numeric(bordeaux_distance)) %>%
                                            filter (dist < 10000)
                                          ,
                                          cells_NAquitane))

# cells in the region of Gironde
gironde_cells <- cells_NAquitane [which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),]
# altitude_stats <- altitude_stats[which(cells_NAquitane$CODE_10KM %in% inter_gironde_NAq$CODE_10KM),]

# select sites
# select sites (cells from Nouvelle Aquitaine within the Bordeaux + 10 km buffer)
sel_sites <- which(cells_NAquitane$CODE_10KM %in% unique(cells_buffer_bordeaux$CODE_10KM))
sp_table_years <- (sp_table_years[sel_sites,,,])
base_table_years <- (base_table_years[sel_sites,,])

# where NAs in table of spp, there is NAs in total records
table(is.na(sp_table_years[,,,6]) == (base_table_years==0))

# Organize Data for spOcc analysis ---------------
# Package all data into a list
# occ: lat+lat2+lon+elev+non_water+non_water2+urban
X <- array(1, dim = c(nrow(sp_table_years), ncol(sp_table_years), 8)) # intercept + n covariates (lat, lat2, long, etc)
X[, , 2] <-  scale(cell_centroid_df[sel_sites, "Y"])[,1]
X[, , 3] <-  (scale(cell_centroid_df[sel_sites, "Y"])[,1])^2
X[, , 4] <-  scale(cell_centroid_df[sel_sites, "X"])[,1]
X[, , 5] <-  scale(altitude_stats$altitude_mean[sel_sites])[,1]
X[, , 6] <-  scale(1-extract.water.summary[sel_sites,"1"])[,1] # non water = 1 minus water(extract.water.summary)
X[, , 7] <-  (scale((1-extract.water.summary[sel_sites,"1"]))[,1])^2
X[, , 8] <-  scale((extract.hab.NA.summary [sel_sites,c("1")]))[,1] # 
#1 230 0 77 255 111 - Continuous urban fabric
#2 255 0 0 255 112 - Discontinuous urban fabric

# plot urban areas
ggplot() +
  geom_sf(fill="white")+
  #geom_sf(data= cells_buffer_bordeaux)+
  geom_sf(data = cbind (gironde_cells,
                        urban = X[, 1, 8]),
          aes(fill=urban,col=urban))+
  scale_fill_viridis_c(na.value = "red")+
  scale_colour_viridis_c(na.value = "red")

# imput altitude because there are 10 NAs
# where are the NAs
ggplot() +
  geom_sf(fill="white")+
  #geom_sf(data= cells_NAquitane)+
  geom_sf(data = cbind (gironde_cells,
                        alt = X[, 1, 5]),
          aes(fill=alt,col=alt))+
  scale_fill_viridis_c(na.value = "red")+
  scale_colour_viridis_c(na.value = "red")

# set the minimum of the scaled values because NAs are in the lowlands/coastline
X[, , 5] [is.na(X[, , 5])] <- min(X[, , 5],na.rm=T)

# detection covariates
X.p <- array (1, dim=c(nrow(sp_table_years),
                       ncol(sp_table_years),
                       dim(sp_table_years)[3],
                       6)) # intercept + 4 covariates
X.p[,,,2] <- X[, , 2] # latitude
X.p[,,,3] <- scale((obs_table_years[sel_sites,,]))#[,(ncol(obs_table_years)-9):ncol(obs_table_years),]))# observers

# bind months (phenology)
months_bind  <- replicate (ncol(sp_table_years),
                           do.call(cbind,lapply (seq(1,10,1), function (i) rep(i,nrow(sp_table_years))))
)

# change the order to fit the right format
months_bind <- aperm(months_bind, c(1,3,2))
X.p[,,,4] <- scale(months_bind)[,1]# linear effect of month
X.p[,,,5] <- (scale(months_bind)[,1])^2# quadratic effect of month

# non-water
X.p[,,,6] <- X[, , 6] # non-water habitats

# select species -------------------------------
sp<- grep("Polyommatus icarus", sp_list)

# identify missing data - zeroed encounter histories
missing_sites <- rowSums(is.na(sp_table_years[,,,sp]))==max(rowSums(is.na(sp_table_years[,,,sp]))) # 10 years x 10 months with zeroes 
table(missing_sites == (rowSums(base_table_years)==0)) # 10 years x 10 months with zeroes 

# missing cells in the buffer
table (rownames (base_table[which(rowSums(base_table_years) == 0),]) %in% unique(cells_buffer_bordeaux$CODE_10KM))
table(rowSums(base_table_years) == 0)/sum(table(rowSums(base_table_years) == 0)) # ~ 20% of the sites were sampled in 24 years

# Coordinates
coords <- cell_centroid_df[,c("X","Y")][which(missing_sites==F),]
cell_centroid_df_buffer <- cell_centroid_df[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),]

# Pack all data into a list
# 
occ.covs <- list(int = X[which(missing_sites==F), , 1],
                 lat = X[which(missing_sites==F), , 2],
                 lat2 = X[which(missing_sites==F), , 3],
                 lon=X[which(missing_sites==F), , 4],
                 elev=X[which(missing_sites==F), , 5],
                 non_water=X[which(missing_sites==F), , 6],
                 non_water2=X[which(missing_sites==F), , 7],
                 urban=X[which(missing_sites==F), , 8],
                 siteID = (seq(1,nrow(X[which(missing_sites==F), , 6]))
                           
                 ))
det.covs <- list(int = X.p[which(missing_sites==F), , , 1],
                 lat=X.p[which(missing_sites==F), , , 2],
                 obs=X.p[which(missing_sites==F), , , 3],
                 phen=X.p[which(missing_sites==F), , , 4],
                 phen2=X.p[which(missing_sites==F), , , 5],
                 non_water=X.p[which(missing_sites==F), , ,6]
                 
)

# bundle data
str(data.list.full <- list(y = sp_table_years [which(missing_sites==F),,,sp], # select common blue 
                           occ.covs = occ.covs, 
                           det.covs = det.covs, 
                           coords = coords))

# number of detections
table(data.list.full$y) 
sum(apply (data.list.full$y,1,max,na.rm=T)) # sites with detection

# calculate and save naive occupancy
sites_det_year <- apply (data.list.full$y,c(1,2), sum,na.rm=T)>0 # if site was sampled in year t
sites_det_year <- colSums(sites_det_year) 
# year effort
eff_year <- colSums(apply (base_table_years,c(1,2), sum,na.rm=T)>0) # total number of sites sampled in each year
naive_occ <- data.frame (sites_det_year = sites_det_year,
                         total_sites_year = eff_year,
                         naive_occ = sites_det_year/eff_year)

# average naive occupancy
mean(naive_occ$naive_occ)*100
mean(naive_occ$sites_det_year)

# missing data in the buffer
table(gironde_cells$CODE_10KM %in% rownames(coords))

# missing data in the buffer
missing_buffer <- gironde_cells$CODE_10KM [which(missing_sites==T)]

# load the models to  make predictions ----------------------------------------
NG15weak <- new.env()
NG15inf <- new.env()
#NG15infSigma <- new.env()
load (file=here ("model_output", 
                 "empirical",
                 "buffer_outputNNG15_urban_Polyommatus .Rdata"),NG15weak)
NG15weak$lab <- "Ng15weak"
load (file=here ("model_output", 
                 "empirical",
                 "buffer_outputNNG15InfPrior_urban_Polyommatus .Rdata"),NG15inf)
NG15inf$lab <- "Ng15inf"

# make predictions and produce plots for all models
res_plots <- lapply (list (NG15weak,
                           NG15inf), function (out) {
      
                             # predict
                             out.pred.occ <- predict(out$out, 
                                                     X, 
                                                     cell_centroid_df_buffer,
                                                     verbose = T,
                                                     t.cols=seq(1,ncol(data.list.full$y)),
                                                     type = 'occupancy'
                             )
                             
                             # map of all cells -----------------------------------------------------------
                             # estimates for non-sampled sites
                             data_non_sampled <- data.frame (cells = gironde_cells$CODE_10KM[which(missing_sites ==T)],
                                                             psi_i = apply (out.pred.occ$psi.0.samples[,which(missing_sites ==T),],2, # expected values (average of predictions) of psi_i
                                                                            mean),
                                                             w_i = apply (out.pred.occ$w.0.samples[,which(missing_sites ==T)],2, # expected values (average of predictions) of w_i
                                                                          mean))
                             # estimates for sampled sites
                             data_sampled <- data.frame (cells = gironde_cells$CODE_10KM[which(missing_sites ==F)],
                                                         psi_i = apply (out.pred.occ$psi.0.samples[,which(missing_sites ==F),],2,mean), # point estimates of psi_i
                                                         w_i = apply (out.pred.occ$w.0.samples[,which(missing_sites ==F)],2,mean)) # point estimates of omega_i
                             # table(data_sampled$cells %in% data_non_sampled$cells) # check if they are all different
                             
                             # sampled within Bordeaux's buffer
                             #data_sampled <- data_sampled[which(data_sampled$cells %in% gironde_cells$CODE_10KM),]
                             
                             # bind the results
                             res_data <- rbind (data_sampled,
                                                data_non_sampled)
                             
                             data_sampled$cells %in% data_non_sampled$cells
                             table(res_data$cells %in% gironde_cells$CODE_10KM)
                             
                             # match cell ID
                             res_data <- res_data[match (gironde_cells$CODE_10KM,res_data$cells),]
                             # table(res_data$cells == gironde_cells$CODE_10KM) # check
                             
                             
                             # map of non-sampled areas ----- out-of-sample predictions needed -------------
                             # occupancy
                             p1_missing<-ggplot(data = gironde_cells[missing_sites==T,]) +
                               geom_sf(fill="white")+
                               geom_sf(data=cbind (gironde_cells[missing_sites==T,],
                                                   bar_psi_i = apply (out.pred.occ$psi.0.samples[,which(missing_sites ==T),],2,mean)),
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
                             p2_missing<-ggplot(data = gironde_cells[which(missing_sites==T),]) +
                               geom_sf(fill="white")+
                               geom_sf(data=cbind (gironde_cells[missing_sites==T,],
                                                   bar_psi_i = apply (out.pred.occ$w.0.samples[,which(missing_sites==T)],2,mean)),
                                       aes (fill=bar_psi_i,
                                            col=bar_psi_i),
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
                             
                             # E(psi_i)
                             p1 <- ggplot(data = gironde_cells) +
                               geom_sf(fill="white")+
                               geom_sf(data=cbind (gironde_cells,
                                                   bar_psi_i = apply (out.pred.occ$psi.0.samples,2,mean),
                                                   #ifelse (extract.water.summary[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),"1"] >0.8, NA, res_data$psi_i),
                                                   det_data_i = !missing_sites), # what have data
                                       aes (fill=bar_psi_i,
                                            col=as.factor(ifelse (det_data_i == T,1,0))),
                                       alpha=0.75)+
                               scale_colour_viridis_d(option = "magma",direction=1,na.value = "gray90",begin=0.75)+
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
                             p2 <- ggplot(data = gironde_cells) +
                               geom_sf(fill="white")+
                               geom_sf(data=cbind (gironde_cells,
                                                   bar_w_i = apply (out.pred.occ$w.0.samples,2,mean),#ifelse (extract.water.summary[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),"1"] >0.8, NA, res_data$w_i),
                                                   det_data_i = !missing_sites), # what have data),
                                       aes (fill=bar_w_i,
                                            col=as.factor(ifelse (det_data_i == T,1,0))),
                                       alpha=0.75)+
                               scale_colour_viridis_d(option = "magma",direction=1,na.value = "gray90",begin=0.75)+
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
                             summ_tab_psi <-apply(out.pred.occ$psi.0.samples,c(1,3),mean)
                             # do by hand to see if it correspond
                             #summ_tab_psib<-lapply (seq(1,nrow(out.pred.occ$psi.0.samples)), function (i)
                             
                             #  colSums(out.pred.occ$psi.0.samples[i,,])/ncol(out$out$psi.samples)
                             
                             #  )
                             #table(do.call(rbind,summ_tab_psib) == summ_tab_psi)
                             
                             # calculate mean yearly occupancy
                             print(data.frame (psi = apply(summ_tab_psi,2,mean),
                                               uci = apply(summ_tab_psi,2,quantile, 0.975),
                                               lci = apply(summ_tab_psi,2,quantile, 0.025)) %>%
                                     colMeans()*100)
                             
                             # plot
                             p3<- data.frame (psi = apply(summ_tab_psi,2,mean),
                                              uci = apply(summ_tab_psi,2,quantile, 0.975),
                                              lci = apply(summ_tab_psi,2,quantile, 0.025),
                                              year = seq(2000,2023)) %>%
                               ggplot(aes(x=year,psi)) +
                               geom_ribbon(aes(ymin=lci, ymax=uci),fill="white")+
                               geom_line(linewidth=1,col="black")+
                               geom_line(data = cbind (naive_occ, 
                                                       year = seq(2000,2023)), 
                                         aes (x=year, y=naive_occ))+
                               ggtitle("")+
                               labs(x="Year", y = expression(paste('E(', hat(psi[t]),')')))+
                               my_theme+
                               ylim(c(0,1))
                             
                             # phenology
                             # summarize
                             par(mfrow=c(1,1),mar=c(4,4,4,4))
                             point_est <- apply (out$out$alpha.samples,2,mean)
                             lci<-apply (out$out$alpha.samples,2,quantile,0.025)
                             uci<-apply (out$out$alpha.samples,2,quantile,0.975)
                             
                             # organize data to plot
                             phen_plot <- data.frame (p = plogis(point_est[1] + 
                                                                   point_est[4]*seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)+
                                                                   point_est[5]*(seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)^2)),
                                                      lci = plogis(lci[1] + 
                                                                     lci[4]*seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)+
                                                                     lci[5]*(seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)^2)),
                                                      uci = plogis(uci[1] + 
                                                                     uci[4]*seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)+
                                                                     uci[5]*(seq(min(X.p[,,,4]),max(X.p[,,,4]),length.out=100)^2)),
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
png (here ("figures", "empirical","Fig_summ_results_Picarus_buffer.png"),width = 27,height = 18,units = "cm",res=200)
  
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
png (here ("figures", "empirical","estimates_missing_cells_Picarus_buffer.png"),width = 12,height = 14,units = "cm",res=150)

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




