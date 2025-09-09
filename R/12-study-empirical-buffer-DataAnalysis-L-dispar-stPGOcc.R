
# --------------------------------

# Empirical data analyses - Lycaena dispar
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
# Bordeaux distance
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
# select sites (cells from Nouvelle Aquitaine within the Bordeaux + 10 km buffer)
sel_sites <- which(cells_NAquitane$CODE_10KM %in% unique(cells_buffer_bordeaux$CODE_10KM))
# subset table
sp_table_years <- (sp_table_years[sel_sites,,,])
base_table_years <- (base_table_years[sel_sites,,])

# where NAs in table of spp, there is NAs in total records
table(is.na(sp_table_years[,,,6]) == (base_table_years==0))

# Organize Data for spOcc analysis ---------------
# Package all data into a list
# Occupancy covariates
X <- array(1, dim = c(nrow(sp_table_years), ncol(sp_table_years), 10)) # intercept + n covariates (lat, lat2, & long, etc )
X[, , 2] <-  scale(cell_centroid_df[sel_sites, "Y"])[,1]
X[, , 3] <-  (scale(cell_centroid_df[sel_sites, "Y"])[,1])^2
X[, , 4] <-  scale(cell_centroid_df[sel_sites, "X"])[,1]
X[, , 5] <-  scale(altitude_stats$altitude_mean[sel_sites])[,1]
X[, , 6] <-  scale(1-extract.water.summary[sel_sites,"1"])[,1] # non water = 1 minus water(extract.water.summary)
X[, , 7] <-  (scale((1-extract.water.summary[sel_sites,"1"]))[,1])^2
X[, , 8] <-  scale(rowSums(extract.hab.NA.summary [sel_sites,c("35","36", "37","40","41","42","43")]))[,1] # codes below
X[, , 9] <-  (scale(rowSums(extract.hab.NA.summary [sel_sites,c("35","36", "37","40","41","42","43")]))[,1])^2 # codes below
#35 166 166 255 255 411 - Inland marshes  <-----------------------
#36 77 77 255 255 412 - Peat bogs <-----------------------
#37 204 204 255 255 421 - Salt marshes  <-----------------------
#40 0 204 242 255 511 - Water courses
#41 128 242 230 255 512 - Water bodies
#42 0 255 166 255 521 - Coastal lagoons
#43 166 255 230 255 522 - Estuaries
X[, , 10] <-  scale((extract.hab.NA.summary [sel_sites,c("1")]))[,1] # 
#1 230 0 77 255 111 - Continuous urban fabric

# Input altitude because there are 10 NAs
# # where are the NAs
(a <- ggplot() +
    geom_sf(fill="white")+
    geom_sf(data= cells_NAquitane)+
    geom_sf(data = cbind (cells_NAquitane,
                          alt = altitude_stats$altitude_mean),
            aes(fill=alt,col=alt))+
    scale_fill_viridis_c(na.value = "red")+
    scale_colour_viridis_c(na.value = "red"))

# Set the minimum of the scaled values because NAs are in the lowlands/coastline
X[, , 5] [is.na(X[, , 5])] <- min(X[, , 5],na.rm=T)

# detection covariates
X.p <- array (1, dim=c(nrow(sp_table_years),
                       ncol(sp_table_years),
                       dim(sp_table_years)[3],
                       6)) # intercept + 5 covariates
X.p[,,,2] <- X[, , 2] # latitude
X.p[,,,3] <- scale((obs_table_years[sel_sites,,]))# n observers

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

# select species ------------------------------
sp<- grep("Lycaena dispar", sp_list)

# identify missing data - zeroed encounter histories
missing_sites <- rowSums(is.na(sp_table_years[,,,sp]))==max(rowSums(is.na(sp_table_years[,,,sp]))) # 10 years x 10 months with zeroes 

# compare the missing sites between 
table(missing_sites) # sp detection table
table(rowSums(base_table_years[])==0) # and the effort table

# missing cells in the buffer
table (rownames (base_table[which(rowSums(base_table_years) == 0),]) %in% unique(cells_buffer_bordeaux$CODE_10KM))
table(rowSums(base_table_years) == 0)/sum(table(rowSums(base_table_years) == 0))

# Coordinates
coords <- cell_centroid_df[,c("X","Y")][which(missing_sites==F),]

# Pack all data into lists
occ.covs <- list(int = X[which(missing_sites==F), , 1],
                 lat = X[which(missing_sites==F), , 2],
                 lat2 = X[which(missing_sites==F), , 3],
                 lon=X[which(missing_sites==F), , 4],
                 elev=X[which(missing_sites==F), , 5],
                 non_water=X[which(missing_sites==F), , 6],
                 non_water2=X[which(missing_sites==F), , 7],
                 marshes=X[which(missing_sites==F), , 8],
                 marshes2=X[which(missing_sites==F), , 9],
                 urban=X[which(missing_sites==F), , 10],
                 siteID = (seq(1,nrow(X[which(missing_sites==F), , 6]))
                 
))
det.covs <- list(int = X.p[which(missing_sites==F), , , 1],
                 lat=X.p[which(missing_sites==F), , , 2],
                 obs=X.p[which(missing_sites==F), , , 3],
                 phen=X.p[which(missing_sites==F), , , 4],
                 phen2=X.p[which(missing_sites==F), , , 5],
                 non_water=X.p[which(missing_sites==F), , , 6]
                 
)

# bundle data
str(data.list.full <- list(y = sp_table_years [which(missing_sites==F),,,sp], # select common blue 
                           occ.covs = occ.covs, 
                           det.covs = det.covs, 
                           coords = coords))
table(data.list.full$y>=0)
sum(data.list.full$y,na.rm=T)

# naive occupancy
table(apply (data.list.full$y,1,sum,na.rm=T)>0)/nrow(data.list.full$y)

# missing data in the buffer
table(gironde_cells$CODE_10KM %in% rownames(coords))

# missing data in the buffer
missing_buffer <- gironde_cells$CODE_10KM [which(gironde_cells$CODE_10KM %in% rownames(coords)==F)]

# plot observed data within the buffer
ggplot(data = gironde_cells) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (gironde_cells,
                      detect = as.factor(apply (sp_table_years[,,,sp],
                                      1,max,na.rm=T))),
          aes (fill=detect,
               col=detect),
          alpha=0.75)+
  scale_fill_viridis_d(option="magma")

#  fit the model ----------------------------------------------------------
# MCMC settings 
n.samples <- 100000
batch.length <- 100
n.burn <- 98000
n.thin <- 5
n.chains <- 3
accept.rate <- 0.43
n.batch <- (n.samples / batch.length)*n.chains

# Analysis under weakly informative priors ----------------------------------
# May want to explore sensitivity of results to these a bit more.  
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
out <- stPGOcc(occ.formula = ~ lat+lat2+lon+elev+marshes+marshes2+urban,
               det.formula = ~ lat+obs+phen+phen2+non_water,
               data = data.list.full,
               n.batch = n.batch/n.chains,
               batch.length = batch.length,
               inits = inits.list,
               priors = prior.list,
               accept.rate = 0.43,
               cov.model = "exponential",
               tuning = tuning.list,
               n.omp.threads = 3, 
               verbose = TRUE,
               ar1 = TRUE,
               NNGP = TRUE,
               n.neighbors = 15,
               n.report = 25,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains)

summary(out)

# save model output
save (out,
      file=here ("model_output", 
                 "empirical",
                 paste0 ("buffer_outputNNG15_urban_", substr(sp_list[sp],1,12),".Rdata")))

rm(out)
gc()

# Informative Priors ----------------------------------------------------------

# May want to explore sensitivity of results to these a bit more.  
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72),
                   sigma.sq.ig = c(a = 2, b = 1.5),
                   sigma.sq.t.ig = c(2, 1),
                   phi.unif = c(a = 3 / 6, b = 3 / 1)) # interpreted as somehow 6 km up to 1 km distance neighbors

# Starting values
z.init <- apply(data.list.full$y, c(1, 2), function(a) as.numeric(sum(a, na.rm = TRUE) > 0))
inits.list <- list(beta = 0, alpha = 0, sigma.sq.t = 0.5, phi = mean (c(a = 3 / 6, b = 3 / 1)), 
                   sigma.sq = 1, rho = 0, z = z.init)

# Tuning
tuning.list <- list(phi = mean (c(a = 3 / 6, b = 3 / 1)), rho = 0.5)

# RUN THE ANALYSIS WITH NNG=15 ----------------------------------------
# Fit the model with stPGOcc
out <- stPGOcc(occ.formula = ~ lat+lat2+lon+elev+marshes+marshes2+urban,
               det.formula = ~ lat+obs+phen+phen2+non_water,
               data = data.list.full,
               n.batch = n.batch/n.chains,
               batch.length = batch.length,
               inits = inits.list,
               priors = prior.list,
               accept.rate = 0.43,
               cov.model = "exponential",
               tuning = tuning.list,
               n.omp.threads = 3, 
               verbose = TRUE,
               ar1 = TRUE,
               NNGP = TRUE,
               n.neighbors = 15,
               n.report = 25,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains)

summary(out)

# save model output
save (out, 
      file=here ("model_output", 
                 "empirical",
                 paste0 ("buffer_outputNNG15InfPrior_urban_", substr(sp_list[sp],1,12),".Rdata")))

# end