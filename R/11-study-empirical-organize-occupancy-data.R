

# ----------------------------------------------------------------

#  Creating an occupancy data (array of sites x years x months x species) to be analyzed and released.

#  This array consider butterfly community data to obtain absences

# --------------------------------------

# load packages
rm(list=ls())
gc()
source ("R/packages.R")

# ggplot theme (D&S theme - nice plots)
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

# load data
load(here ("Processed_data",
           "ProcessedDetectionData.RData"))

# ----------------------------------------------------

# Effort table (total of records)
# Community data used as a basis to obtain non-detections
# all sampled cells (with the detection of at least one species)
base_table <- (tapply (dataPointsPolygonClean$DiscDet, # detections
                       list (dataPointsPolygonClean$Maille1, # sites
                             dataPointsPolygonClean$YearMonth # survey
                       ),
                       FUN= function(x) sum(x,na.rm=T)
))
range(base_table,na.rm=T)
base_table[is.na(base_table)] <- 0 # if sum resulted in NA, replace by zero (ie no species was sampled in this cell, year and month)

# add missing cells to the visit table
base_table <- rbind (base_table,
                              matrix(0, 
                                     nrow = length(unique(cells_NAquitane$CODE_10KM)) - nrow(base_table),
                                     ncol = ncol (base_table),
                                     dimnames = list(
                                        unique(cells_NAquitane$CODE_10KM) [which(!unique(cells_NAquitane$CODE_10KM) %in% rownames(base_table))],
                                                     colnames(base_table))
                              )
)
# match names
base_table <- base_table [match(unique(cells_NAquitane$CODE_10KM), rownames(base_table)),]
table(rownames(base_table) == unique(cells_NAquitane$CODE_10KM)) # check site names

# values reported in the MS regarding surveys per site
range(rowSums(base_table))

# create a vector as basis to match secondary occasion names
months <- unique(dataPointsPolygonClean$Month)[order(unique(dataPointsPolygonClean$Month))]
years <- unique(dataPointsPolygonClean$Year); years <- years[order(years,decreasing=F)]
(match_sec_occ <- lapply (months , function (i)
  
  paste (i, years,sep=" ") # bind months and yers
  ))

# melt
match_sec_occ <- do.call(rbind,match_sec_occ)
match_sec_occ<- melt(match_sec_occ)

# add missing columns () - secondary occasions
base_table <- cbind (base_table,
                     
                     matrix(0, nrow = nrow(base_table),
                            ncol = length(match_sec_occ$value [which(match_sec_occ$value %in% colnames(base_table)==F)]),
                            dimnames = list(
                              rownames(base_table),
                              match_sec_occ$value [which(match_sec_occ$value %in% colnames(base_table)==F)])
                     )
)

# match base table  (sampling events)
(base_table <- base_table [,match(match_sec_occ$value, colnames(base_table))]) # correct temporal ordering
table(base_table>0) 

# range of the number of observations across cells (reported in the MS)
range(apply(base_table,1,sum))
mean(apply(base_table,1,sum)) # average N records
mean(apply(base_table,1,sd)) # sd

# count the number of spp per site x year x month ------------------------------------------------------------
spp_table <- (tapply (dataPointsPolygonClean$TaxNomVal, # detections
                       list (dataPointsPolygonClean$Maille1, # sites
                             dataPointsPolygonClean$YearMonth # survey
                       ),
                       FUN= function(x) length(unique(x))
))
spp_table[is.na(spp_table)] <- 0 # if sum resulted in NA, replace by zero (ie none species was sampled in this cell and date)

# add missing cells to the community table
spp_table <- rbind (spp_table,
                     matrix(0, 
                            nrow = length(unique(cells_NAquitane$CODE_10KM)) - nrow(spp_table),
                            ncol = ncol (spp_table),
                            dimnames = list(
                              unique(cells_NAquitane$CODE_10KM) [which(!unique(cells_NAquitane$CODE_10KM) %in% rownames(spp_table))],
                              colnames(spp_table))
                     )
)

# match names
spp_table <- spp_table [match(unique(cells_NAquitane$CODE_10KM), rownames(spp_table)),]
table(rownames(spp_table) == unique(cells_NAquitane$CODE_10KM))

# add missing columns ()
spp_table <- cbind (spp_table,
                     
                    matrix(0, nrow = nrow(spp_table),
                           ncol = length(match_sec_occ$value [which(match_sec_occ$value %in% colnames(spp_table)==F)]),
                           dimnames = list(
                             rownames(spp_table),
                             match_sec_occ$value [which(match_sec_occ$value %in% colnames(spp_table)==F)])
                    )
)

# match spp_table  (sampling events)
spp_table <- spp_table [,match(match_sec_occ$value, colnames(spp_table))]
colnames(spp_table) == match_sec_occ$value

# summary
range(spp_table)

# number of observers per site x year x month ------------------------------------------------------------
obs_table <- (tapply (dataPointsPolygonClean$AnonObserverID, # observers per sampling event
                      list (dataPointsPolygonClean$Maille1, # sites
                            dataPointsPolygonClean$YearMonth # survey
                      ),
                      FUN= function(x) length(unique(x)))
)

obs_table[is.na(obs_table)] <- 0 # if sum resulted in NA, replace by zero (ie no observer was out in this cell and date)

# add missing cells to the # observer table
obs_table <- rbind (obs_table,
                    matrix(0, 
                           nrow = length(unique(cells_NAquitane$CODE_10KM)) - nrow(obs_table),
                           ncol = ncol (obs_table),
                           dimnames = list(
                             unique(cells_NAquitane$CODE_10KM) [which(!unique(cells_NAquitane$CODE_10KM) %in% rownames(obs_table))],
                             colnames(obs_table))
                    )
)
# match names
obs_table <- obs_table [match(unique(cells_NAquitane$CODE_10KM), rownames(obs_table)),]
table(rownames(obs_table) == unique(cells_NAquitane$CODE_10KM))

# add missing columns ()
obs_table <- cbind (obs_table,
                    
                    matrix(0, nrow = nrow(obs_table),
                           ncol = length(match_sec_occ$value [which(match_sec_occ$value %in% colnames(obs_table)==F)]),
                           dimnames = list(
                             rownames(obs_table),
                             match_sec_occ$value [which(match_sec_occ$value %in% colnames(obs_table)==F)])
                    )
)

# match obs_table colnames (sampling events)
obs_table <- obs_table [,match(match_sec_occ$value, colnames(obs_table))]
colnames(obs_table) == match_sec_occ$value

# summary
range(obs_table)

# number of dates per site x year x month -------------------------------------------------------
dates_table <- (tapply (dataPointsPolygonClean$DateDebut, # detections
                      list (dataPointsPolygonClean$Maille1, # sites
                            dataPointsPolygonClean$YearMonth # survey
                      ),
                      FUN= function(x) length(unique(x))
))
dates_table[is.na(dates_table)] <- 0 # if sum resulted in NA, replace by zero (ie none species was sampled in this cell and date)

# add missing cells to the # observer table
dates_table <- rbind (dates_table,
                    matrix(0, 
                           nrow = length(unique(cells_NAquitane$CODE_10KM)) - nrow(dates_table),
                           ncol = ncol (dates_table),
                           dimnames = list(
                             unique(cells_NAquitane$CODE_10KM) [which(!unique(cells_NAquitane$CODE_10KM) %in% rownames(dates_table))],
                             colnames(dates_table))
                    )
)

# match names  (sampling events)
dates_table <- dates_table [match(unique(cells_NAquitane$CODE_10KM), rownames(dates_table)),]
table(rownames(dates_table) == unique(cells_NAquitane$CODE_10KM))

# add missing columns ()
dates_table <- cbind (dates_table,
                    
                      matrix(0, nrow = nrow(dates_table),
                             ncol = length(match_sec_occ$value [which(match_sec_occ$value %in% colnames(dates_table)==F)]),
                             dimnames = list(
                               rownames(dates_table),
                               match_sec_occ$value [which(match_sec_occ$value %in% colnames(dates_table)==F)])
                      )
)

# match dates_table
dates_table <- dates_table [,match(match_sec_occ$value, colnames(dates_table))]
colnames(dates_table) == match_sec_occ$value

# detection of the target species ----------------------------------------------
# create one table per species
# select focal spp
sp_list <- c("Coenonympha oedippus (Fabricius, 1787)",
             "Euphydryas aurinia (Rottemburg, 1775)",
             "Lycaena dispar (Haworth, 1802)",
             "Lycaena phlaeas (Linnaeus, 1761)",
             "Maniola jurtina (Linnaeus, 1758)",
             "Polyommatus icarus (Rottemburg, 1775)" )

# subset the dataset to have data for focal spp  
dataPointsPolygonClean_subsetSpp <- dataPointsPolygonClean %>%
  
  filter (TaxNomVal %in% sp_list)

# detection table
species_table <- (tapply (dataPointsPolygonClean_subsetSpp$DiscDet, # detections
                          list (dataPointsPolygonClean_subsetSpp$Maille1, # sites
                                dataPointsPolygonClean_subsetSpp$YearMonth, # months
                                dataPointsPolygonClean_subsetSpp$TaxNomVal # species
                          ),
                          FUN= function(x) sum(x,na.rm=T)
))
range(species_table,na.rm=T) # rang N records
species_table[species_table>=1] <- 1 # if any observation, set 1

# build the matrix to input surveys (as these six focal species were not sighted in all combinations of surveys)
input_matrix_surveys <- array(NA,dim = c(
  nrow(species_table), # cells not in the detection table
  sum(!(colnames(base_table)) %in% colnames(species_table)), # sampling events (years x months) not in the detection table
  dim(species_table)[3]) , # N spp
  
  # set names
  dimnames= list (rownames(species_table),
                  colnames(base_table) [(!(colnames(base_table)) %in% colnames(species_table))],
                  dimnames(species_table)[[3]])
)

# bind missing surveys (second dimension of species_table)
species_table <- (abind(species_table,
                        input_matrix_surveys,
                        along=2
))

# match colnames (surveys)
species_table <- species_table[,match (colnames(base_table),colnames(species_table)),]

# build the matrix to input sites (as the six focal species were not sighted in all sites with samples)
input_matrix_sites <- array(NA,
                              dim = c(
                                  sum(!(rownames(base_table)) %in% rownames(species_table)), # cells not in the detection table
                                  ncol(base_table), # cells not in the detection table
                                  dim(species_table)[3]) , # n visits
                                  
                              dimnames= list (rownames(base_table) [!(rownames(base_table)) %in% rownames(species_table)],
                                              colnames(base_table),
                                              dimnames(species_table)[[3]])
)

# bind missing cells (1st dimension of species_table)
species_table <- (abind(species_table,
                        input_matrix_sites,
                        along=1
))

# match rownames (sites)
species_table <- species_table[match (rownames(base_table),rownames(species_table)),,]
table(rownames(species_table) == (cells_NAquitane$CODE_10KM)) # check if names correspond
dim(species_table)

# --------------------------------
# Considering community data (base_table) to obtain absences
# see the number of 1s and NAs - note that we still don't have zeros so we need to use the base_table to have them
lapply (seq(1,dim(species_table)[3]), function (i)
  
  c("1"=sum(species_table[,,i]>=1,na.rm=T), # the number of 1s
    "0"=sum(species_table[,,i]==0,na.rm=T), # zeros
    "NA"=sum(is.na(species_table[,,i])) # the number of NAs
    
  ))

# using community data (base_table) to replace NAs by 0 when any species was detected
for (i in 1:dim(species_table)[3]) {
  # keep the array structure using the for loop
  species_table[,,i] <- ifelse (is.na(species_table[,,i]) & base_table >= 1, # if NA in the detection of a specific species and there is detection of at least one other species  
                                0, # set zero (a survey occurred (base_table >= 1), but the focal spp was not found)
                                species_table[,,i]) # otherwise, keep NAs or 1s
  
}

# check after the replacement (the number of NAs should decrease, and the zeros should appear)
lapply (seq(1,dim(species_table)[3]), function (i)
  
  c("1"=sum(species_table[,,i]>=1,na.rm=T), # the number of 1s
    "0"=sum(species_table[,,i]==0,na.rm=T), # zeros
    "NA"=sum(is.na(species_table[,,i])) # the number of NAs
  ))
# do some checks
sum((base_table)==0) # the number of zeros in the base table is equal to the number of NAs in the detection tables
table(is.na(species_table[,,1])) # test for the first species 
table(is.na(species_table[,,6])) # test for the 6th species 
# fine!

## geographic coordinates (cell centroid)
cell_centroid <- cells_NAquitane %>%
  st_centroid() %>%
  filter (CODE_10KM %in% rownames(base_table)) 

# subset
cell_centroid <- cell_centroid[match(rownames(base_table), cell_centroid$CODE_10KM),]

#table(cell_centroid$id_unique == rownames(base_table))
cell_centroid_df <- cell_centroid %>%
  st_coordinates()
rownames(cell_centroid_df) <- cell_centroid$CODE_10KM

# -------------------------
# 
# Now transform these data into arrays

# sequence of j months and t years
df_years <- cbind (seq(1,ncol(base_table),10),
                   seq(10,ncol(base_table),10))

# species_table is organized with sampling events (month x year combinations) across columns. We need to 
# create a third dimension to represent the months
# separate the columns 
sp_table_years <-  lapply (seq(1,nrow(df_years)), function (i)
  
  species_table [,seq (df_years[i,1],df_years[i,2]),]
  
)
length(sp_table_years) # length equal to the number of years
lapply(sp_table_years,colnames) # check colnames
dim(sp_table_years[[1]]) # each table now has site x year x sp dims 

# shell_array
shell_array <- array (NA, dim = c(dim(sp_table_years[[1]])[1],
                                  dim(sp_table_years[[1]])[2],
                                  dim(sp_table_years[[1]])[3],
                                  nrow(df_years)
))

# fill with each individual table
for (i in 1:length(sp_table_years)) {
  
  shell_array [,,,i] <- sp_table_years[[i]]
  
}

# make some checks
table(shell_array[,,,1] == species_table[,1:10,]) # data for the first year match?
table(shell_array[,,,2] == species_table[,11:20,]) # data for the 2nd year match?
table(shell_array[,,,5] == species_table[,41:50,]) # data for the 5nd year match?
table(shell_array[,,,10] == species_table[,91:100,]) # data for the 10nd year match?
table(shell_array[,,,24] == species_table[,231:240,]) # data for the last 24th year match?

# rename
sp_table_years <- shell_array

# change the order to fit the right format: site x year x month x species
sp_table_years <- aperm(sp_table_years, c(1,4,2,3))

# still match?
table(sp_table_years[,1,,] == species_table[,1:10,]) # check 1st year
table(sp_table_years[,24,,] == species_table[,231:240,]) # check 24th year

# do the same edit for the number of observers ----------------------
obs_table_years <-  lapply (seq(1,nrow(df_years)), function (i)
  
  obs_table [,seq (df_years[i,1],df_years[i,2])]
  
)

# shell_array
shell_array <- array (NA, dim = c(dim(obs_table_years[[1]])[1],
                                  dim(obs_table_years[[1]])[2],
                                  nrow(df_years)
))

# fill
for (i in 1:length(obs_table_years)) {
  
  shell_array [,,i] <- obs_table_years[[i]]
  
}

# check if every data is equal
table(shell_array[,,1] == obs_table[,1:10])
table(shell_array[,,2] == obs_table[,11:20])
table(shell_array[,,1] == obs_table[,1:10])
table(shell_array[,,10] == obs_table[,91:100])
table(shell_array[,,24] == obs_table[,231:240])

# rename
obs_table_years <- shell_array

# change the order to fit the right format (site x year x month)
obs_table_years <- aperm(obs_table_years, c(1,3,2))

# existence of values correspond between the arrays
table((obs_table_years >= 0) == (sp_table_years[,,,1] >= 0))

# do the same edit for base_table (effort table) ------------------------------------
base_table_years <-  lapply (seq(1,nrow(df_years)), function (i)
  
  base_table [,seq (df_years[i,1],df_years[i,2])]
  
)

# shell_array
shell_array <- array (NA, dim = c(dim(base_table_years[[1]])[1],
                                  dim(base_table_years[[1]])[2],
                                  nrow(df_years)
))

# fill
for (i in 1:length(base_table_years)) {
  
  shell_array [,,i] <- base_table_years[[i]]
  
}

# check
table(shell_array[,,1] == base_table[,1:10])
table(shell_array[,,2] == base_table[,11:20])
table(shell_array[,,1] == base_table[,1:10])
table(shell_array[,,10] == base_table[,91:100])
table(shell_array[,,24] == base_table[,231:240])

# rename
base_table_years <- shell_array

# change the order to fit the right format
base_table_years <- aperm(base_table_years, c(1,3,2))

# existence of values correspond between arrays
table((obs_table_years >= 0) == (base_table_years >= 0))
table((sp_table_years[,,,1] >= 0) == (base_table_years >= 0)) # sp1
table((sp_table_years[,,,6] >= 0) == (base_table_years >= 0)) # sp6

# always used the same order of sites
table(rownames(species_table) == cells_NAquitane$CODE_10KM)
table(cell_centroid$CODE_10KM == cells_NAquitane$CODE_10KM)

# Save the objects to be used as input of the multi-season site occupancy models
save (dataPointsPolygonClean_subsetSpp,
      species_table,
      sp_table_years,
      base_table,
      base_table_years,
      cell_centroid_df,
      obs_table_years,
      sp_list,
      cells_NAquitane,
      file = here("Processed_data", 
                  "Occupancy_data_spOccupancy.RData"))

# clean work space
rm(list=ls())
# end

