

# ----------------------------------------------------------------

#  Summarizing the data set motivating the simulations
#  
#  Butterfly data from the whole Nouvelle-Aquitaine region ,at the scale of 1x1 km (systematic grid). 

#  Here we created intermediary files for creating occupancy (occurrence) data. The observations were done through points, linestrings and polygons (geometries did not follow the TypeGeoms so we merged/pooled all data for analyzes).
#  We did some editing in the name of observers to minimize redundancy/ortographic errors, and also filtered data considering the years (2000-2023) and month of analysis (February to November). 

# ------------------------------------------------------------------

# load packages
rm(list=ls())
source ("R/packages.R")

# create a dir to receive the figures
dir.create(here("figures", "empirical"))

# load spatial data (maille) -------------------------------
cells_NAquitane <- st_read(dsn=here ("Data", "SpatialData","Maillage_1x1km"),
                           layer="1x1km_n-a")

# Gironde department ----------------------------
# source: https://www.actualitix.com/blog/shapefiles-des-departements-de-france.html
cells_Gironde <- st_read(dsn=here ("Data", "SpatialData","33-gironde"),
                           layer="33-")

# communes of Bordeaux Metropole
communes <- c("BORDEAUX",
              "AMBARES-ET-LAGRAVE",
              "AMBES",
              "ARTIGUES-PRES-BORDEAUX"  ,
              "BASSENS",
              "BEGLES",
              "BLANQUEFORT",
              "BOULIAC",
              "LE BOUSCAT",
              "BRUGES",
              "CARBON-BLANC",
              "CENON",
              "EYSINES",
              "FLOIRAC",
              "GRADIGNAN",
              "LE HAILLAN",
              "LORMONT",
              "MARTIGNAS-SUR-JALLE" ,
              "MERIGNAC",
              "PAREMPUYRE",
              "PESSAC",
              "SAINT-AUBIN-DE-MEDOC"  ,
              "SAINT-LOUIS-DE-MONTFERRAND" ,
              "SAINT-MEDARD-EN-JALLES"  ,
              "SAINT-VINCENT-DE-PAUL" ,
              "LE TAILLAN-MEDOC"  ,
              "TALENCE",
              "VILLENAVE-D'ORNON" )
              
# plot              
a <- ggplot() +
  geom_sf(fill="white")+
  geom_sf(data= cells_NAquitane)

# bordeaux distance
#bordeaux_metropole <- cells_Gironde %>%
#             filter (NOM_COMM %in% communes)

# bordeaux distance
bordeaux_distance <- st_distance (cells_Gironde,
                                  cells_Gironde %>%
                                    filter (NOM_COMM == "BORDEAUX"))

#  plot
(a <- ggplot()  + geom_sf(data=cells_Gironde %>%
                           cbind (dist=as.numeric(bordeaux_distance)) %>%
                           filter (dist < 10000),
                         aes(fill= NOM_COMM))+
  theme(legend.position = "none"))

# plot
#(a <- a  + geom_sf(data=bordeaux_metropole,
#            aes(fill= NOM_COMM))+
#  theme(legend.position = "none"))

# intersection Gironde - NAquitaine
#cells_metropole <- st_intersection(cells_NAquitane,
#                                   bordeaux_metropole)

# intersection Buffer - Gironde -- buffer of 10 km around Bordeaux
cells_buffer_bordeaux <- (st_intersection(cells_Gironde %>%
                                        cbind (dist=as.numeric(bordeaux_distance)) %>%
                                        filter (dist < 10000)
                                      ,
                                      cells_NAquitane))

# matching between Bordeaux Metropole and an area of 10 km around the city
#a + geom_sf(data=cells_buffer_bordeaux,
#            aes(fill=NOM_REGION))+ 
#  geom_sf(data=cells_metropole,
#          aes(fill=NOM_COMM),col= "black")+
#  theme(legend.position = "none")

#   load butterfly (community) data ------------------------------
#  families
families <- list.dirs(here::here("Data","SpeciesData", "Community Data"),full.names = F)[-1]
families <- families[-grep("2025", families)]

# families in the atlas data
# families_atlas <- list.dirs(here::here("Data","SpeciesData", "Atlas Data"),full.names = F)[-1]

#  occurrence data summary (points) -------------------------------------
#k = species[1]
#i= method [2]
dataPoints<- lapply (families, function (k) 
  
  tryCatch({
    
    # load
    dat<-read.csv (here ("Data","SpeciesData",
                         "Community Data",
                         k,
                         "Point.csv"), 
                   sep= "\t")
    
    # spatial object
    datSpatPts <- (st_sf (cbind (dat,
                                 st_as_sfc(dat[,"GeomWkt"])),
                          crs = 2154))
    
    # remove the third dimension to avoid issues
    datSpatPts <- st_zm(datSpatPts)
    
    # changing counts
    # if presence then must (at least) have min and max = 1  individual
    datSpatPts[which (is.na(datSpatPts$DenbrMin) & datSpatPts$StatPresen=="présent"),"DenbrMin"] <- 1
    datSpatPts[which (is.na(datSpatPts$DenbrMax) & datSpatPts$StatPresen=="présent"),"DenbrMax"] <- 1
    datSpatPts$DenbrMean <- (datSpatPts$DenbrMax+datSpatPts$DenbrMin)/2
    
    # absence
    datSpatPts[datSpatPts$StatPresen=="absent","DenbrMin"] <- 0
    datSpatPts[datSpatPts$StatPresen=="absent","DenbrMax"] <- 0
    datSpatPts[datSpatPts$StatPresen=="absent","DenbrMean"] <- 0
    
    # discrete detection
    datSpatPts$DiscDet <- ifelse (datSpatPts$StatPresen == "présent", 1,0)
    
    # data over time
    # Convert to class date
    # https://epirhandbook.com/en/working-with-dates.html
    datSpatPts <- datSpatPts %>% 
      mutate(DateDebut = as.Date(DateDebut, format = "%Y-%m-%d"))%>% 
      mutate(Year = year(DateDebut),
             Month=month(DateDebut,label = T,abbr=T),
             
             YearMonth = zoo::as.yearmon(DateDebut),
             # add the information that this is occurrence data
             dataType = "Occurrence")
    
    
    datSpatPts
    
  },  error = function(e) return ("NULL"))
  
)

# bind all data
dataPoints <- do.call(rbind, dataPoints)

# data linestrings ---------------------------------------------
#i= method [2]
dataLinestrings<- lapply (families, function (k) 
  
  tryCatch({
    
    # load
    dat<-read.csv (here ("Data","SpeciesData",
                         "Community Data",
                         k,
                         "Linestring.csv"), 
                   sep= "\t")
    
    # spatial object
    datSpatPts <- (st_sf (cbind (dat,
                                 st_as_sfc(dat[,"GeomWkt"])),
                          crs = 2154))
    
    # remove the third dimension to avoid issues
    datSpatPts <- st_zm(datSpatPts)
    
    # changing counts
    # if presence then must (at least) have min and max = 1  individual
    datSpatPts[which (is.na(datSpatPts$DenbrMin) & datSpatPts$StatPresen=="présent"),"DenbrMin"] <- 1
    datSpatPts[which (is.na(datSpatPts$DenbrMax) & datSpatPts$StatPresen=="présent"),"DenbrMax"] <- 1
    datSpatPts$DenbrMean <- (datSpatPts$DenbrMax+datSpatPts$DenbrMin)/2
    
    # absence
    datSpatPts[datSpatPts$StatPresen=="absent","DenbrMin"] <- 0
    datSpatPts[datSpatPts$StatPresen=="absent","DenbrMax"] <- 0
    datSpatPts[datSpatPts$StatPresen=="absent","DenbrMean"] <- 0
    
    # discrete detection
    datSpatPts$DiscDet <- ifelse (datSpatPts$StatPresen == "présent", 1,0)
    
    # data over time
    # Convert to class date
    # https://epirhandbook.com/en/working-with-dates.html
    datSpatPts <- datSpatPts %>% 
      mutate(DateDebut = as.Date(DateDebut, format = "%Y-%m-%d"))%>% 
      mutate(Year = year(DateDebut),
             Month=month(DateDebut,label = T,abbr=T),
             
             YearMonth = zoo::as.yearmon(DateDebut),
             # add the information that this is occurrence data
             dataType = "Occurrence")
    
    
    datSpatPts
    
  },  error = function(e) return ("NULL"))
  
)

# bind all data
dataLinestrings <- do.call(rbind, dataLinestrings)

# load polygons ---------------------------------------------------
dataPolygons<- lapply (families, function (k) 
  
  tryCatch({
    
    # load
    dat<-read.csv (here ("Data","SpeciesData",
                         "Community Data",
                         k,
                         "Polygon.csv"), 
                   sep= "\t")
    
    # spatial object
    datSpatPol <- (st_sf (cbind (dat,
                                 st_as_sfc(dat[,"GeomWkt"])),
                          crs = 2154))
    
    # remove the third dimension to avoid issues
    datSpatPol <- st_zm(datSpatPol)
    
    # changing counts
    # if presence then must (at least) have min and max = 1  individual
    datSpatPol[which (is.na(datSpatPol$DenbrMin) & datSpatPol$StatPresen=="présent"),"DenbrMin"] <- 1
    datSpatPol[which (is.na(datSpatPol$DenbrMax) & datSpatPol$StatPresen=="présent"),"DenbrMax"] <- 1
    datSpatPol$DenbrMean <- (datSpatPol$DenbrMax+datSpatPol$DenbrMin)/2
    
    # absence
    datSpatPol[datSpatPol$StatPresen=="absent","DenbrMin"] <- 0
    datSpatPol[datSpatPol$StatPresen=="absent","DenbrMax"] <- 0
    datSpatPol[datSpatPol$StatPresen=="absent","DenbrMean"] <- 0
    
    # discrete detection
    datSpatPol$DiscDet <- ifelse (datSpatPol$StatPresen == "présent", 1,0)
    
    # data over time
    # Convert to class date
    # https://epirhandbook.com/en/working-with-dates.html
    datSpatPol <- datSpatPol %>% 
      mutate(DateDebut = as.Date(DateDebut, format = "%Y-%m-%d"))%>% 
      mutate(Year = year(DateDebut),
             Month=month(DateDebut,label = T,abbr=T),
             YearMonth = zoo::as.yearmon(DateDebut),
             # add the information that this is occurrence data
             dataType = "Occurrence")
    
    datSpatPol
    
  },  error = function(e) return ("NULL"))
  
)

# bind all data
dataPolygons <- do.call(rbind, dataPolygons)

# -------------------------------------------------

# cleaning of point data
dataPointsClean <- dataPoints %>%
  filter (StatPresen == "présent") %>% # only presence data 
  filter (Maille1 != "") %>% #choose only observations in the scale of 1x1km 
  mutate(#scientific_name = gsub("\\s*\\([^\\)]+\\)","",TaxNomVal), # remove words inside parenthesis
         Observer = tolower (Observer),# lower case
         ObserverMod = gsub("\\s*\\([^\\)]+\\)","",Observer),# remove words inside parenthesis) # lower case
         ObserverMod = gsub("\\s*\\(\\)\\s*", "", ObserverMod),# remove parenthesis and space
         ObserverMod = gsub("�", "é", ObserverMod, fixed = TRUE), # remove encoding issue
         ObserverMod = stringi::stri_trans_general(ObserverMod, "Latin-ASCII"),   # removing diacritical marks (accents)
         ObserverMod = gsub("-", " ", ObserverMod ), # remove dash
         ObserverMod = gsub("_", " ", ObserverMod ), # remove underline
         ObserverMod = gsub("[()]", "", ObserverMod), # remove unique parenthesis
         ObserverMod = gsub("/", ",", ObserverMod), # replace /
         ObserverMod = gsub("&", ",", ObserverMod), # replace &
         ObserverMod = gsub("mathilde  lasfargue", "mathilde lasfargue", ObserverMod),
         ObserverMod = gsub("\\?", "e", ObserverMod),# adjust one name
         ObserverMod = gsub("chazaud, p.", "chazaud p.", ObserverMod),# adjust one name
         ObserverMod = gsub("chretien, p.", "chretien p.", ObserverMod),# adjust one name
         ObserverMod = gsub("corradini, p.", "corradini p.", ObserverMod),# adjust one name
         ObserverMod = gsub("davy costaille et michael pasquet", "davy costaille, michael pasquet", ObserverMod),# adjust one name
         ObserverMod = gsub("dufay, c.", "dufay c.", ObserverMod),# adjust one name
         ObserverMod = gsub("durand, g.", "durand g.", ObserverMod),# adjust one name
         ObserverMod = gsub("gaubet, p.", "gaubet p.", ObserverMod),# adjust one name
         ObserverMod = gsub(";", ",", ObserverMod),# adjust one name
         ObserverMod = gsub("s.vrignaud", "s. vrignaud", ObserverMod),# adjust one name
         ObserverMod = gsub("yannig bernard (eliomys), yves suffran (bordeaux métropole)", "yves suffran (bordeaux métropole), yannig bernard (eliomys)", ObserverMod),# adjust one name
         ObserverMod = gsub("pierre grsivard", "pierre grisvard", ObserverMod)) # adjust one name

# cleaning of linestring data
dataLinestringsClean <- dataLinestrings %>%
  filter (StatPresen == "présent") %>% # only presence data 
  filter (Maille1 != "") %>% #choose only observations in the scale of 1x1km 
  mutate(#scientific_name = gsub("\\s*\\([^\\)]+\\)","",TaxNomVal), # remove words inside parenthesis
    Observer = tolower (Observer),# lower case
    ObserverMod = gsub("\\s*\\([^\\)]+\\)","",Observer),# remove words inside parenthesis) # lower case
    ObserverMod = gsub("\\s*\\(\\)\\s*", "", ObserverMod),# remove parenthesis and space
    ObserverMod = gsub("�", "é", ObserverMod, fixed = TRUE), # remove encoding issue
    ObserverMod = stringi::stri_trans_general(ObserverMod, "Latin-ASCII"),   # removing diacritical marks (accents)
    ObserverMod = gsub("-", " ", ObserverMod ), # remove dash
    ObserverMod = gsub("_", " ", ObserverMod ), # remove underline
    ObserverMod = gsub("[()]", "", ObserverMod), # remove unique parenthesis
    ObserverMod = gsub("/", ",", ObserverMod), # replace /
    ObserverMod = gsub("&", ",", ObserverMod), # replace &
    ObserverMod = gsub("mathilde  lasfargue", "mathilde lasfargue", ObserverMod),
    ObserverMod = gsub("\\?", "e", ObserverMod),# adjust one name
    ObserverMod = gsub("chazaud, p.", "chazaud p.", ObserverMod),# adjust one name
    ObserverMod = gsub("chretien, p.", "chretien p.", ObserverMod),# adjust one name
    ObserverMod = gsub("corradini, p.", "corradini p.", ObserverMod),# adjust one name
    ObserverMod = gsub("davy costaille et michael pasquet", "davy costaille, michael pasquet", ObserverMod),# adjust one name
    ObserverMod = gsub("dufay, c.", "dufay c.", ObserverMod),# adjust one name
    ObserverMod = gsub("durand, g.", "durand g.", ObserverMod),# adjust one name
    ObserverMod = gsub("gaubet, p.", "gaubet p.", ObserverMod),# adjust one name
    ObserverMod = gsub(";", ",", ObserverMod),# adjust one name
    ObserverMod = gsub("s.vrignaud", "s. vrignaud", ObserverMod),# adjust one name
    ObserverMod = gsub("yannig bernard (eliomys), yves suffran (bordeaux métropole)", "yves suffran (bordeaux métropole), yannig bernard (eliomys)", ObserverMod),# adjust one name
    ObserverMod = gsub("pierre grsivard", "pierre grisvard", ObserverMod)) # adjust one name

# cleaning of polygon data
dataPolygonsClean <- dataPolygons %>%
  filter (StatPresen == "présent") %>% # only presence data 
  filter (Maille1 != "") %>% #choose only observations in the scale of 1x1km 
  mutate(#scientific_name = gsub("\\s*\\([^\\)]+\\)","",TaxNomVal), # remove words inside parenthesis
         Observer = tolower (Observer),# lower case
         ObserverMod = gsub("\\s*\\([^\\)]+\\)","",Observer),# remove words inside parenthesis) # lower case
         ObserverMod = gsub("\\s*\\(\\)\\s*", "", ObserverMod),# remove parenthesis and space
         ObserverMod = gsub("�", "é", ObserverMod, fixed = TRUE), # remove encoding issue
         ObserverMod = stringi::stri_trans_general(ObserverMod, "Latin-ASCII"),   # removing diacritical marks (accents)
         ObserverMod = gsub("-", " ", ObserverMod ), # remove dash
         ObserverMod = gsub("_", " ", ObserverMod ), # remove underline
         ObserverMod = gsub("[()]", "", ObserverMod), # remove unique parenthesis
         ObserverMod = gsub("/", ",", ObserverMod), # replace /
         ObserverMod = gsub("&", ",", ObserverMod), # replace &
         ObserverMod = gsub("mathilde  lasfargue", "mathilde lasfargue", ObserverMod),# adjust one name
         ObserverMod = gsub("\\?", "e", ObserverMod),# adjust one name
         ObserverMod = gsub("chazaud, p.", "chazaud p.", ObserverMod),# adjust one name
         ObserverMod = gsub("chretien, p.", "chretien p.", ObserverMod),# adjust one name
         ObserverMod = gsub("corradini, p.", "corradini p.", ObserverMod),# adjust one name
         ObserverMod = gsub("davy costaille et michael pasquet", "davy costaille, michael pasquet", ObserverMod),# adjust one name
         ObserverMod = gsub("dufay, c.", "dufay c.", ObserverMod),# adjust one name
         ObserverMod = gsub("durand, g.", "durand g.", ObserverMod),# adjust one name
         ObserverMod = gsub("gaubet, p.", "gaubet p.", ObserverMod),# adjust one name
         ObserverMod = gsub(";", ",", ObserverMod),# adjust one name
         ObserverMod = gsub("s.vrignaud", "s. vrignaud", ObserverMod),# adjust one name
         ObserverMod = gsub("yannig bernard (eliomys), yves suffran (bordeaux métropole)", "yves suffran (bordeaux métropole), yannig bernard (eliomys)", ObserverMod),# adjust one name
         ObserverMod = gsub("pierre grsivard", "pierre grisvard", ObserverMod)) # adjust one name

# ------------------
# bind polygon (occurrence and atlas), linestring and point data
dataPointsPolygonClean <- rbind(dataPointsClean %>% 
                                  # select the right columns to bind the data
                                  dplyr::select (Year, Month,YearMonth, Maille1, DiscDet, ObserverMod,Observer,TypeGeom,
                                          TaxNomVal,TaxNomVern,DateDebut , dataType),
                                
                                dataLinestringsClean%>% 
                                  # select the right columns to bind the data
                                  dplyr::select (Year, Month,YearMonth, Maille1, DiscDet, ObserverMod,Observer,TypeGeom,
                                                 TaxNomVal,TaxNomVern,DateDebut , dataType),
                                
                                dataPolygonsClean%>% 
                                  # select the right columns to bind the data
                                  dplyr::select (Year, Month,YearMonth, Maille1, DiscDet, ObserverMod,Observer,TypeGeom,
                                                 TaxNomVal,TaxNomVern,DateDebut , dataType)
                                )

# View(data.frame(unique(dataPointsPolygonClean$ObserverMod)[order(unique(dataPointsPolygonClean$ObserverMod))]))
# inconnu here belongs to some institution

# who was doing fieldwork on december and january
#table(dataPointsPolygonClean[which(dataPointsPolygonClean$Month %in% c(1,12)),"ObserverMod"][[1]])[order(table(dataPointsPolygonClean[which(dataPointsPolygonClean$Month %in% c(1,12)),"ObserverMod"][[1]])
#)]

# select the year of interest (2000 onwards) > 2000 observations overall
dataPointsPolygonClean <- dataPointsPolygonClean %>%
  filter (Year >= "2000" & 
            Year <= "2023" &
            Month %in% levels(Month)[2:11] &
            Maille1 != ""
  )

# yearly observations
table(dataPointsPolygonClean$DiscDet,
      dataPointsPolygonClean$Year)

# number of observers associated to each point
# creating a list of unique observers to acknowledge their work
unique_observers <- strsplit (dataPointsPolygonClean$ObserverMod,",") # split names based on the commas
unique_observers <- lapply (unique_observers, function (x) gsub('^ | $', '', x)) # remove the first and the last empty characteres 
length(unique(unlist(unique_observers))) # unlist the IDs

# add the number of observers to the data set
dataPointsPolygonClean$Nobservers <- unlist(lapply (unique_observers,length)) 

#  a tiny number of observations (n=806) had zero observers because the ID was not filled
table( dataPointsPolygonClean$Nobservers )
# explore empty observer names
unique((dataPointsPolygonClean$IdCadreAc[which(dataPointsPolygonClean$Observer =="")]))
# There belong to nine projects - looking metadata we see:
# "756c2c06-7bac-480f-abc1-98b45777e76b" - Inventaire Odonates-Rhopalocères 1999-2002 du PNRLG
# "b85faeaa-ae97-46df-b666-53e3798d0de6" - Inventaire Rhopalocère 2008 du PNRLG
# "12166b2e-5e3c-47aa-a217-b2fecb9ef363" - Inventaires écologiques Grands Projets du Sud-Ouest (GPSO)
# "620dd4ab-9f8d-4186-9fab-baaae127be72" - Edition des bulletins de la Société Linnéenne de Bordeaux
# "a70c4aa3-19b7-4dff-b217-5650f5a6383f" - Étude des invertébrés sur les sites Natura 2000 de l'estuaire de la Gironde
# "7eef371e-9d6a-4fa8-a9bd-0616f8ad1593" - DGA St-Jean-d'Illac
# "bf9630c7-6550-49c2-aa31-3ca277122c45" - Données opportunistes du Parc Naturel Régional Périgord-Limousin (PNRPL)
# "5d1b1f33-0373-4bf0-85ec-ad03773c1fe9" - Saisie de données naturalistes d'observateurs indépendants sur la plateforme de l'Observatoire FAUNA
# "7579f8af-a299-4ef5-a3c1-1910f13135f6" - Elaboration du DOCOB du site Natura 2000 ‘Zones humides de l’arrière dune des pays de Born et de Buch’ (FR7200714)


# these cases occurred for only 10 cells so I (ALLuza) decided to keep these records
length(unique(dataPointsPolygonClean$Maille1 [which(dataPointsPolygonClean$Nobservers>0)])) # N cells with 1+ observers
length(unique(dataPointsPolygonClean$Maille1 [which(dataPointsPolygonClean$Nobservers>=0)])) # N cells with 0+ observers 

#  As such, if the number of observers is zero, because the observer ID/entry was empty, I assumed the existence of at least one observer
dataPointsPolygonClean$Nobservers[dataPointsPolygonClean$Nobservers==0] <- 1

# change yearmonth factor
dataPointsPolygonClean$YearMonth <- gsub ("\\.", "", dataPointsPolygonClean$YearMonth ) # adjust colnames - remove dot

# number of observers per site x year x month ------------------------------------------------------------
obs_table <- (tapply (dataPointsPolygonClean$ObserverMod, # observers per sampling event
                      list (dataPointsPolygonClean$Maille1, # sites
                            dataPointsPolygonClean$YearMonth # survey
                      ),
                      FUN= function(x) length(unique(x)))
)
obs_table[is.na(obs_table)] <- 0 # if sum resulted in NA, replace by zero (ie none species was sampled in this cell and date)
obs_table[which(rownames(obs_table) == "E390N6332"),]

range(obs_table,na.rm=T)
# check this cell
which(apply(obs_table,1,max,na.rm=T)==max(obs_table,na.rm=T))
obs_table [which(apply(obs_table,1,max,na.rm=T)==max(obs_table,na.rm=T)),] # weird

# look the data for that site and year
View(dataPointsPolygonClean[which(dataPointsPolygonClean$Maille1 == "E390N6332"),])
length(unique(dataPointsPolygonClean$ObserverMod[which(dataPointsPolygonClean$Maille1 == "E390N6332")]))

# anonymous observer IDs (simply set numbers)
dataPointsPolygonClean$AnonObserverID<- as.numeric(as.factor(dataPointsPolygonClean$ObserverMod))
range(table((dataPointsPolygonClean$AnonObserverID)))

# ----------------------------------------- 

# cells in the maille
length(unique(cells_NAquitane$CODE_10KM) )

# cells in Bordeaux Metropole
# length(unique(cells_metropole$CODE_10KM))

# cells in the buffer around Bordeaux
length(unique(cells_buffer_bordeaux$CODE_10KM))

# data in each spatial unit
# regionally
nrow(dataPointsPolygonClean[which(unique(dataPointsPolygonClean$Maille1) %in% unique(cells_NAquitane$CODE_10KM) ), ])

# bordeaux metropole
# nrow(dataPointsPolygonClean[which(unique(dataPointsPolygonClean$Maille1) %in% unique(cells_metropole$CODE_10KM) ), ])

# bordeaux buffer
nrow(dataPointsPolygonClean[which(unique(dataPointsPolygonClean$Maille1) %in% unique(cells_buffer_bordeaux$CODE_10KM) ), ])

# of these, the following is observed in the presence-only data of the six species
# table(unique(dataPoints$Maille1) %in% unique(dataPoints$Maille1))
length(unique(dataPointsPolygonClean$Maille1)) # sampled cells

# cells with no sampling
nrow(cells_NAquitane) - length(unique(dataPointsPolygonClean$Maille1))

# cells in the Bordeaux metropole
# table(unique(dataPointsPolygonClean$Maille1) %in% unique(cells_metropole$CODE_10KM))

# cells in the buffer
table(unique(dataPointsPolygonClean$Maille1) %in% unique(cells_buffer_bordeaux$CODE_10KM))

# percentage of the cells with at least one record (of any species / community)
table(unique(dataPointsPolygonClean$Maille1) %in% cells_NAquitane$CODE_10KM)/nrow(cells_NAquitane)

# percentage of cells with zeroes
1-table(unique(dataPointsPolygonClean$Maille1) %in% cells_NAquitane$CODE_10KM)/nrow(cells_NAquitane)

# nobs
nrow (dataPointsPolygonClean)

# number of observations per method
sort(table(dataPointsPolygonClean$TypeGeom), decreasing = T)  # types of data sets

# unique species
sort(unique(dataPointsPolygonClean$TaxNomVal))

# save data to organize occupancy data (array of sites x year x month x species) and observation level covariates
save(dataPointsPolygonClean, 
     cells_NAquitane, 
     cells_buffer_bordeaux,
     file = here ("Processed_data","ProcessedDetectionData.RData"))

# clean workspace
rm(list=ls())
## end
