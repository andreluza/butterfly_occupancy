
# --------------------------------

# Empirical data analyses - species not in the main text and supp info
# Using the data collected within the buffer of 10km^2 around Bordeaux

# using spOccupancy package
# models with nngp=15
# weak prior for phi
# informative prior for phi

# --------------------------------------

# load packages
rm(list=ls())
source ("R/packages.R")

# create a dir to receive the results
dir.create(here("model_output", "empirical"))

# load processed data ------------------------
load (file = here("Processed_data", 
                  "Occupancy_data_spOccupancy.RData"))

# number of sites sampled - more sites were sampled in 2021 than 2018 (which had + bfly records)
(apply (sp_table_years, c(2), function (x) sum(is.na(x)==F))/nrow(sp_table_years))
mean(apply (sp_table_years, c(2), function (x) sum(is.na(x)==F))/nrow(sp_table_years)) # the averaged % sampled cells across all 24 years

table(apply (is.na(sp_table_years[,,,1]),1,sum)==apply (is.na(sp_table_years[,,,2]),1,sum))
table(apply (is.na(sp_table_years[,,,1]),1,sum)==apply (is.na(sp_table_years[,,,6]),1,sum))
table(apply (is.na(sp_table_years[,,,1]),1,sum)==apply(base_table_years==0,1,sum))

# load altitude and habitat  -------------------------------
load (file = here ("Processed_data",
                   "AltitudeHabitatStats.RData"))

# load water  -------------------------------
load (file = here ("Processed_data",
                   "Water.RData"))

# check covariates - water
cbind(cells_NAquitane,
      altitude_stats) %>%
  ggplot()+
  geom_sf(aes(fill = altitude_mean ,col=altitude_mean ))

# check covariates - land
cbind(cells_NAquitane,
      land = 1-extract.water.summary[,"1"]) %>% # the column "1" represents water extracted from CORINE habitat 
  ggplot()+
  geom_sf(aes(fill = land ,col=land ))

# load spatial data -------------------------------
cells_NAquitane <- st_read(dsn=here ("Data", "SpatialData","Maillage_1x1km"),
                           layer="1x1km_n-a")

# choosing a smaller scale for prediction --------------------------------------
# Gironde department
# source: https://www.actualitix.com/blog/shapefiles-des-departements-de-france.html
cells_Gironde <- st_read(dsn=here ("Data", "SpatialData","33-gironde"),
                           layer="33-")

# bordeaux distance
bordeaux_distance <- st_distance (cells_Gironde,
                                  cells_Gironde %>%
                                    filter (NOM_COMM == "BORDEAUX"))

# intersection Gironde - NAq
cells_buffer_bordeaux <- (st_intersection(cells_Gironde %>%
                                            cbind (dist=as.numeric(bordeaux_distance)) %>%
                                            filter (dist < 10000) # here the distance from Bordeaux is defined
                                          ,
                                          cells_NAquitane))

# cells in the region of Gironde
gironde_cells <- cells_NAquitane [which(cells_NAquitane$CODE_10KM %in% unique(cells_buffer_bordeaux$CODE_10KM)),]

# select sites
sp_table_years <- (sp_table_years[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),,,])
base_table_years <- (base_table_years[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),,])

# where NAs in table of spp, there is NAs in total records
table(is.na(sp_table_years[,,,6]) == (base_table_years==0))

# Organize Data for spOcc analysis ---------------
# Package all data into a list
# occ
X <- array(1, dim = c(nrow(sp_table_years), ncol(sp_table_years), 9)) # intercept + n covariates (lat, lat2, long, etc)
X[, , 2] <-  scale(cell_centroid_df[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM), "Y"])[,1]
X[, , 3] <-  (scale(cell_centroid_df[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM), "Y"]^2)[,1])
X[, , 4] <-  scale(cell_centroid_df[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM), "X"])[,1]
X[, , 5] <-  scale(altitude_stats$altitude_mean[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM)])[,1]
X[, , 6] <-  scale(1-extract.water.summary[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),"1"])[,1] # non water = 1 minus water(extract.water.summary)
X[, , 7] <-  (scale(1-(extract.water.summary[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),"1"]^2))[,1])
X[, , 8] <-  scale((extract.hab.NA.summary [which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),c("1")]))[,1] # 
X[, , 9] <-  (scale((extract.hab.NA.summary [which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),c("2")]))[,1]) # 
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
X.p[,,,3] <- scale((obs_table_years[which(cells_NAquitane$CODE_10KM %in% cells_buffer_bordeaux$CODE_10KM),,]))#[,(ncol(obs_table_years)-9):ncol(obs_table_years),]))# observers

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

#  fit the model ----------------------------------------------------------
# MCMC settings 
n.samples <- 100000
batch.length <- 100
n.burn <- 98000
n.thin <- 5
n.chains <- 3
accept.rate <- 0.43
n.batch <- (n.samples / batch.length)*n.chains

# select species -------------------------------
# run for all the four species at once ------------------------------
sp <- which( sp_list %in% 
               c ("Coenonympha oedippus (Fabricius, 1787)",  
                  "Euphydryas aurinia (Rottemburg, 1775)",
                  "Lycaena phlaeas (Linnaeus, 1761)",      
                  "Maniola jurtina (Linnaeus, 1758)")
)  # select species

# apply data organization and model run for all four species
lapply (sp, function (sel_sp) {
  
      # identify missing data - zeroed encounter histories
      missing_sites <- rowSums(is.na(sp_table_years[,,,sel_sp]))==max(rowSums(is.na(sp_table_years[,,,sel_sp]))) # 10 years x 10 months with zeroes 
      table(missing_sites == (rowSums(base_table_years)==0)) # 10 years x 10 months with zeroes 
      
      # missing cells in the buffer
      table (rownames (base_table[which(rowSums(base_table_years) == 0),]) %in% unique(cells_buffer_bordeaux$CODE_10KM))
      table(rowSums(base_table_years) == 0)/sum(table(rowSums(base_table_years) == 0)) # ~ 20% of the sites were sampled in 24 years
      
      # Coordinates
      coords <- cell_centroid_df[,c("X","Y")][which(missing_sites==F),]
      
      # Pack all data into a list
      occ.covs <- list(int = X[which(missing_sites==F), , 1],
                       lat = X[which(missing_sites==F), , 2],
                       lat2 = X[which(missing_sites==F), , 3],
                       lon=X[which(missing_sites==F), , 4],
                       elev=X[which(missing_sites==F), , 5],
                       non_water=X[which(missing_sites==F), , 6],
                       non_water2=X[which(missing_sites==F), , 7],
                       urban=X[which(missing_sites==F), , 8]#,
                       #siteID = (seq(1,nrow(X[which(missing_sites==F), , 6])))
                                 
                       )
      det.covs <- list(int = X.p[which(missing_sites==F), , , 1],
                       lat=X.p[which(missing_sites==F), , , 2],
                       obs=X.p[which(missing_sites==F), , , 3],
                       phen=X.p[which(missing_sites==F), , , 4],
                       phen2=X.p[which(missing_sites==F), , , 5],
                       non_water=X.p[which(missing_sites==F), , ,6]
                       
      )
      
      # bundle data
      str(data.list.full <- list(y = sp_table_years [which(missing_sites==F),,,sel_sp], # select common blue 
                                 occ.covs = occ.covs, 
                                 det.covs = det.covs, 
                                 coords = coords))
      table(data.list.full$y>=0)
      sum(data.list.full$y,na.rm=T) # number of detections
      
      # naive occupancy
      table(apply (data.list.full$y,1,sum,na.rm=T)>0)/nrow (data.list.full$y)
      
      # missing data in the buffer
      table(gironde_cells$CODE_10KM %in% rownames(coords))
      
      # Analysis under weakly informative priors ----------------------------------
      prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                         alpha.normal = list(mean = 0, var = 2.72),
                         sigma.sq.ig = c(a = 2, b = 1.5),
                         sigma.sq.t.ig = c(2, 1),
                         phi.unif = c(a = 3 / 1, b = 3 / 0.05))
      
      # Starting values
      z.init <- apply(data.list.full$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
      inits.list <- list(beta = 0, alpha = 0, sigma.sq.t = 0.5, phi = 3 / .5, 
                         sigma.sq = 1, rho = 0, z = z.init)
      
      # Tuning
      tuning.list <- list(phi = 0.5, rho = 0.5)
      
      # RUN THE ANALYSIS WITH NNG=15 ----------------------------------------
      # Fit the model with stPGOcc
      out <- stPGOcc(occ.formula = ~ lat+lat2+lon+elev+non_water+non_water2+urban,
                     det.formula = ~ non_water+lat+obs+phen+phen2,
                     data = data.list.full,
                     n.batch = n.batch/n.chains,
                     batch.length = batch.length,
                     inits = inits.list,
                     priors = prior.list,
                     accept.rate = 0.43,
                     cov.model = "exponential",
                     tuning = tuning.list,
                     n.omp.threads = 3, # TODO: change as necessary. 
                     verbose = TRUE,
                     ar1 = TRUE,
                     NNGP = TRUE,
                     n.neighbors = 15,
                     n.report = 25,
                     n.burn = n.burn,
                     n.thin = n.thin,
                     n.chains = n.chains)
      
      # save model output
      save (out,
            file=here ("model_output", 
                       "empirical",
                       paste0 ("buffer_outputNNG15_urban_", substr(sp_list[sel_sp],1,12),".Rdata")))
      # remove heavy objects
      rm(out)
      rm(data.list.full)
      rm(occ.covs)
      rm(det.covs)
      gc()
      
      
})

# end