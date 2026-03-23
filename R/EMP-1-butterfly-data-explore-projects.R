

# ----------------------------------------------------------------

#  Emplore projects
#  s
#  Butterfly data from the whole Nouvelle-Aquitaine region ,at the scale of 1x1 km (systematic grid). 

#  Here we explore the possibility of targeted sampling be more frequent than community-wide sampling (as per Shirey et al. 2018, Methods in Ecol Evol)  

# ------------------------------------------------------------------

# load packages
rm(list=ls())
source ("R/packages.R")

# create a dir to receive the figures
dir.create(here("figures", "empirical"))

# load spatial data (maille) -------------------------------
cells_NAquitane <- st_read(dsn=here ("Data", "SpatialData","Maillage_1x1km"),
                           layer="1x1km_n-a")

# Map Gironde department ----------------------------
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
(a <- ggplot() +
  geom_sf(fill="white")+
  geom_sf(data= cells_NAquitane))

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
    
    # Associate metadata to enable the correct differentiation of data sets
    metadat<-read.csv (here ("Data","SpeciesData",
                                   "Community Data",
                             k,
                             "Metadonnees.csv"), 
                       sep= "\t")
    #unique(metadat$NomJeuDonnees)
    #table(datSpatPts$IdJdd %in% metadat$IdJeuDonnees)
    #table(datSpatPts$IdCadreAc %in% metadat$IdCadre) # they match but the variable names are different
    # change the scd column name to match data and metadata
    
    # match data and metadata
    datSpatPts <- cbind (datSpatPts,
                         metadat [match(datSpatPts$IdJdd, metadat$IdJeuDonnees), c("NomJeuDonnees", "NomCadre")]
    )
    
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
    
    
    # Associate metadata to enable the correct differentiation of data sets
    metadat<-read.csv (here ("Data","SpeciesData",
                             "Community Data",
                             k,
                             "Metadonnees.csv"), 
                       sep= "\t")
    
    #table(datSpatPts$IdCadreAc %in% metadat$IdCadre) # they match but the variable names are different
    # change the scd column name to match data and metadata
    
    # match data and metadata
    datSpatPts <- cbind (datSpatPts,
                         metadat [match(datSpatPts$IdJdd, metadat$IdJeuDonnees), c("NomJeuDonnees", "NomCadre")]
    )
    
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
    
    
    # Associate metadata to enable the correct differentiation of data sets
    metadat<-read.csv (here ("Data","SpeciesData",
                             "Community Data",
                             k,
                             "Metadonnees.csv"), 
                       sep= "\t")
    
    #table(datSpatPts$IdCadreAc %in% metadat$IdCadre) # they match but the variable names are different
    # change the scd column name to match data and metadata
    
    # match data and metadata
    datSpatPol <- cbind (datSpatPol,
                         metadat [match(datSpatPol$IdJdd, metadat$IdJeuDonnees), c("NomJeuDonnees", "NomCadre")]
    )
    
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
                                          TaxNomVal,TaxNomVern,DateDebut , dataType, NomJeuDonnees, NomCadre),
                                
                                dataLinestringsClean%>% 
                                  # select the right columns to bind the data
                                  dplyr::select (Year, Month,YearMonth, Maille1, DiscDet, ObserverMod,Observer,TypeGeom,
                                                 TaxNomVal,TaxNomVern,DateDebut , dataType, NomJeuDonnees, NomCadre),
                                
                                dataPolygonsClean%>% 
                                  # select the right columns to bind the data
                                  dplyr::select (Year, Month,YearMonth, Maille1, DiscDet, ObserverMod,Observer,TypeGeom,
                                                 TaxNomVal,TaxNomVern,DateDebut , dataType, NomJeuDonnees, NomCadre)
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

## standardized and opportunistic data
data_types <- (dataPointsPolygonClean$TypeGeom)
data_types [data_types %in% c("Linéaire", "Transect (Point)", "Transect Sentinelle du Climat (Point)")] <- "Transects"

data.frame (records=sort(table(data_types)),
            Data = names(sort(table(data_types)))) %>%
  ggplot()+
  geom_bar(aes(x=reorder (Data,-records.Freq),y=records.Freq),stat="identity")+
  theme_bw(base_size = 14)+
  theme(axis.text.x = element_text(angle=90,vjust = .65))+
  labs(x="Data type", y="Number of records")

# yearly projects
colSums(table(dataPointsPolygonClean$NomCadre,
      dataPointsPolygonClean$Year)>0)

# species per project
spp_proj_table <- table(dataPointsPolygonClean$NomCadre,
                        dataPointsPolygonClean$TaxNomVal)
spp_proj_table[spp_proj_table>0] <- 1
range(apply (spp_proj_table,1,sum))
dim(spp_proj_table)

data.frame (table(apply (spp_proj_table,1,sum))) %>%
  mutate (Var_cat = cut(as.numeric(Var1), 
                         breaks=c(-Inf, 1, 2, 5, 10, 20, 30, 40, 50, 60, 160), 
                         labels=c("1","2","3-5","6-10","11-20","21-30","31-40","41-50","51-60","61-160"))) %>%
  group_by(Var_cat) %>%
  reframe (Freq = sum(Freq)) %>%
  mutate (Perc = Freq/504) %>%
  #mutate(Freq=sum(Freq))
  ggplot(aes(x=Var_cat, y=Freq))+
  geom_bar(stat="identity")+
  geom_text(aes(x=Var_cat, y=Freq+3,label = round(Perc,3)*100)) +
  labs (x="Number of butterfly taxa (categories)",y="Number of projects")+
  theme_bw()

ggsave(here ("figures", "empirical", "projects_data.png"),width=6,height = 4)
  
# which ones are single species projs
dim(spp_proj_table[which(apply (spp_proj_table,1,sum)==1),]) # 43/504
# two species projects
dim(spp_proj_table[which(apply (spp_proj_table,1,sum)==2),]) # 36/504
# multispecies
dim(spp_proj_table[which(apply (spp_proj_table,1,sum)>2),]) # 425/504

# View
View(spp_proj_table[which(apply (spp_proj_table,1,sum)==1),] %>%
  melt() %>%
  filter (value > 0 )
)

# targeted
dim(dataPointsPolygonClean %>%
  filter (NomCadre %in% "Programme « Papillons menacés des zones humides » en Aquitaine") )
# likely community sampling
dim(dataPointsPolygonClean %>%
      filter (NomCadre %in% "Atlas de la Biodiversité Communal (ABC) de Montmorillon") )

# range
range(rowSums(table(dataPointsPolygonClean$NomCadre,
                    dataPointsPolygonClean$TaxNomVal)>0))
table(rowSums(table(dataPointsPolygonClean$NomCadre,
                    dataPointsPolygonClean$TaxNomVal)>0))==1

# projects with more than 100 species
spp_proj_table[which(apply (spp_proj_table,1,sum)>100),]

# example Données du Parc national des Pyrénées
# years 
dataPointsPolygonClean %>%
  filter (NomCadre == "Données du Parc national des Pyrénées") %>%
  group_by(NomCadre) %>%          # Then, with the filtered data, group it by "bb"
  summarise(Unique_Elements = n_distinct(TaxNomVal)) %>%   # Now summarise with unique elements per group
  ungroup()

# cells
dataPointsPolygonClean %>%
 filter (NomCadre == "Données du Parc national des Pyrénées") %>%
  group_by(NomCadre) %>%          # Then, with the filtered data, group it by "bb"
  summarise(Unique_Elements = n_distinct(Maille1)) %>%   # Now summarise with unique elements per group
  ungroup()

# years 
dataPointsPolygonClean %>%
  filter (NomCadre == "Données du Parc national des Pyrénées") %>%
  group_by(NomCadre) %>%          # Then, with the filtered data, group it by "bb"
  summarise(Unique_Elements = n_distinct(Year),
            min_years = min (Year),
            max_years = max(Year)) %>%   # Now summarise with unique elements per group
  ungroup()



# clean workspace
rm(list=ls())
## end
