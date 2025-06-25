
# --------------------------------------------------------


# Extract spatial covariates used in occupancy models
# Altitude from EU-DEM
# Habitat from CORINE database

# NA = Nouvelle Aquitaine 
# 
# -------------------------------------------------------


# load packages
rm(list=ls())
source ("R/packages.R")


# ggplot theme
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


# -------------------------------------------------------------


# create a dir to receive the figures

dir.create(here("figures", "empirical"))


# ---------------------------------------

# load spatial data
cells_NAquitane <- st_read(dsn=here ("Data", "SpatialData","Maillage_1x1km"),
                           layer="1x1km_n-a")


#load copernicus data
altitude.files <- list.files(path=here ("Data", "Altitude","eudem_buffered", "eudem"),
                             pattern =".tif$", full.names=TRUE)
altitude.europe <- lapply (altitude.files, raster)
altitude.region.interest <- (altitude.europe[[11]])

# crop to decrease the size
# using locator
box <-extent(altitude.region.interest)-400000 # 
b <- as(extent(box), 'SpatialPolygons')
plot(altitude.region.interest)
plot(b,add=T)
plot(cells_NAquitane,add=T)

# crop
altitude.NA <- crop(altitude.region.interest, 
                    b)

# reproject this raster
altitude.NA.Proj <- projectRaster(altitude.NA,
                                  crs = crs(cells_NAquitane))
# check
plot(altitude.NA.Proj)
plot(cells_NAquitane,add=T)

# convert to native format to speed up the extraction
poly_cells<-vect(cells_NAquitane)
rs<-rast(altitude.NA.Proj)

# extract
altitude.NA.cells <-terra::extract(rs,poly_cells)
altitude_stats <- altitude.NA.cells%>%
  group_by(ID) %>%
  summarise(altitude_mean=mean(N2000000E3000000,na.rm=T),
            altitude_sd=sd(N2000000E3000000,na.rm=T),
            #altitude_min=min(N2000000E3000000,na.rm=T),
            #altitude_max=max(N2000000E3000000,na.rm=T)
            ) %>%
  bind_cols(cell_id=poly_cells$ID)

# set values
values(poly_cells) <- cbind (values(poly_cells),
                       alt=altitude_stats$altitude_mean)

# I think it is fine
plot(rs>100)
plot(poly_cells, col = poly_cells$alt>100,add=T)
plot(rs,col=terrain.colors(n=20))

# -----------------------------

# CORINE habitat scheme

# QGIS Generated Color Map Export File
#INTERPOLATION:DISCRETE
#1 230 0 77 255 111 - Continuous urban fabric
#2 255 0 0 255 112 - Discontinuous urban fabric
#3 204 77 242 255 121 - Industrial or commercial units
#4 204 0 0 255 122 - Road and rail networks and associated land
#5 230 204 204 255 123 - Port areas
#6 230 204 230 255 124 - Airports
#7 166 0 204 255 131 - Mineral extraction sites
#8 166 77 0 255 132 - Dump sites
#9 255 77 255 255 133 - Construction sites
#10 255 166 255 255 141 - Green urban areas
#11 255 230 255 255 142 - Sport and leisure facilities
#12 255 255 168 255 211 - Non-irrigated arable land
#13 255 255 0 255 212 - Permanently irrigated land
#14 230 230 0 255 213 - Rice fields
#15 230 128 0 255 221 - Vineyards
#16 242 166 77 255 222 - Fruit trees and berry plantations
#17 230 166 0 255 223 - Olive groves
#18 230 230 77 255 231 - Pastures         <-----------------------
#19 255 230 166 255 241 - Annual crops associated with permanent crops
#20 255 230 77 255 242 - Complex cultivation patterns
#21 230 204 77 255 243 - Land principally occupied by agriculture with significant areas of natural vegetation
#22 242 204 166 255 244 - Agro-forestry areas
#23 128 255 0 255 311 - Broad-leaved forest
#24 0 166 0 255 312 - Coniferous forest
#25 77 255 0 255 313 - Mixed forest
#26 204 242 77 255 321 - Natural grasslands     <-----------------------
#27 166 255 128 255 322 - Moors and heathland   <-----------------------
#28 166 230 77 255 323 - Sclerophyllous vegetation 
#29 166 242 0 255 324 - Transitional woodland-shrub
#30 230 230 230 255 331 - Beaches - dunes - sands
#31 204 204 204 255 332 - Bare rocks
#32 204 255 204 255 333 - Sparsely vegetated areas
#33 0 0 0 255 334 - Burnt areas
#34 166 230 204 255 335 - Glaciers and perpetual snow
#35 166 166 255 255 411 - Inland marshes  <-----------------------
#36 77 77 255 255 412 - Peat bogs <-----------------------
#37 204 204 255 255 421 - Salt marshes  <-----------------------
#38 230 230 255 255 422 - Salines
#39 166 166 230 255 423 - Intertidal flats
#40 0 204 242 255 511 - Water courses
#41 128 242 230 255 512 - Water bodies
#42 0 255 166 255 521 - Coastal lagoons
#43 166 255 230 255 522 - Estuaries
#44 230 242 255 255 523 - Sea and ocean
#48 255 255 255 255 999 - NODATA

#load copernicus data
corine_habitat <- rast(here ("Data", "CorineHabitat","U2018_CLC2018_V2020_20u1", 
                "U2018_CLC2018_V2020_20u1", 
                "U2018_CLC2018_V2020_20u1.tif"),subds=0,lyrs=1)

# coords
transf_cells <- st_transform(cells_NAquitane,
                crs(corine_habitat))

# crop
hab.NA <- crop(corine_habitat, 
               transf_cells)

# find the different types of habitat
# values (hab.NA) <- (ifelse (values(hab.NA) == 15,1,0)) # test grape cover
plot(hab.NA)

# water
water <- hab.NA
# NAs also is water
values(water)[is.na(values(water))] <- 40 # set as water courses

# select water from the list above
values (water) <- (ifelse (values(water) %in% c(39:44),1,0)) # 
plot(water)

# extract water
extract.water <- terra::extract(water, transf_cells,
                                 method='simple',
                                 exact=T,
                                 cellnumbers=T)


# summarize
extract.water.summary <- extract.water %>%
  #filter(is.na(U2018_CLC2018_V2020_20u1)!=T) %>%
  group_by(ID,U2018_CLC2018_V2020_20u1) %>%
  summarise(mean.cover = mean(fraction, na.rm=T),
            max.cover = max(fraction, na.rm=T),
            median.cover = median(fraction, na.rm=T),
            sum.cover = sum(fraction, na.rm=T)
            
  ) %>%
  # table of land use per cell
  reshape::cast (formula=ID~U2018_CLC2018_V2020_20u1,
                 value = "mean.cover",na.rm=T,fill = 0)


# world map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
France <- world %>%
  filter (name == "France")

# crop 
xlim<- c(-2,3) # c(-5,10)
ylim<- c(42.5,48) # c(40,52)

# map
# match and plot
png (here ("figures","empirical","water.png"),width=20,height=25,unit="cm",res=300)

ggplot(data = world %>%
         filter (name == "France")
) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (cells_NAquitane,
                      Cover = extract.water.summary[,"1"]),
          aes (fill=Cover,
               col=Cover),
          alpha=0.5)+
  scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
  scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
  coord_sf(xlim=xlim,ylim=ylim)+
  theme_bw() +
  theme(legend.position = c(0.7,0.2),
        axis.text.y = element_blank(),
        strip.text = element_text(face = "italic"),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.background = element_blank(),
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))

dev.off()


# extract habitat covariates ------------------------------

extract.hab.NA <- terra::extract(hab.NA, transf_cells,
                                 method='simple',
                                 exact=T,
                                 cellnumbers=T)
length(unique(extract.hab.NA$ID))
unique(extract.hab.NA$U2018_CLC2018_V2020_20u1)[order(unique(extract.hab.NA$U2018_CLC2018_V2020_20u1))]

# summarize
extract.hab.NA.summary <- extract.hab.NA %>%
  #filter(is.na(U2018_CLC2018_V2020_20u1)!=T) %>%
  group_by(ID,U2018_CLC2018_V2020_20u1) %>%
  summarise(mean.cover = mean(fraction, na.rm=T),
            max.cover = max(fraction, na.rm=T),
            median.cover = median(fraction, na.rm=T),
            sum.cover = sum(fraction, na.rm=T)
            
  ) %>%
  # table of land use per cell
  reshape::cast (formula=ID~U2018_CLC2018_V2020_20u1,
                 value = "mean.cover",na.rm=T,fill = 0)
  
head(extract.hab.NA.summary)
dim((extract.hab.NA.summary))

# world map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
France <- world %>%
  filter (name == "France")

# crop 
xlim<- c(-2,3) # c(-5,10)
ylim<- c(42.5,48) # c(40,52)

# list of interesing uses
uses <- data.frame (code = c("18","26","27","30","32","35","36","37"),
                    use = c("Pastures", 
                            "Natural grasslands", 
                            "Moors and heathland",
                            "Beaches - dunes - sands",
                            "Sparsely vegetated areas",
                            "Inland marshes",
                            "Peat bogs", 
                            "Salt marshes"))

# land use maps
maps_land_use <- lapply (seq(1,nrow(uses)), function (i)

  # match and plot
  ggplot(data = world %>%
             filter (name == "France")
    ) +
      geom_sf(fill="white")+
      geom_sf(data=cbind (cells_NAquitane,
                          Cover = (extract.hab.NA.summary[,which(colnames(extract.hab.NA.summary) %in% uses[i,1])])),
              aes (fill=Cover,
                   col=Cover),
              alpha=0.5)+
      scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
      scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
      coord_sf(xlim=xlim,ylim=ylim)+
      theme_bw() +
      labs (title = uses[i,2])+
      theme(legend.position = c(0.7,0.2),
            axis.text.y = element_blank(),
            strip.text = element_text(face = "italic"),
            legend.key.size = unit(0.3, 'cm'), #change legend key size
            legend.background = element_blank(),
            legend.key.height = unit(0.3, 'cm'), #change legend key height
            legend.key.width = unit(0.3, 'cm'), #change legend key width
            legend.title = element_text(size=8), #change legend title font size
            legend.text = element_text(size=6))

)

# arrange
grid.arrange(maps_land_use[[1]],maps_land_use[[2]],
             maps_land_use[[3]],maps_land_use[[4]],
             maps_land_use[[5]],maps_land_use[[6]],
             maps_land_use[[7]],maps_land_use[[8]],
             nrow=2,ncol=4
             )


# map of altitude
ggplot(data = world %>%
         filter (name == "France")
) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (cells_NAquitane,
                      Altitude = altitude_stats$altitude_sd),
          aes (fill=Altitude,
               col=Altitude),
          alpha=0.5)+
  scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
  scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
  coord_sf(xlim=xlim,ylim=ylim)+
  theme_bw() +
  theme(legend.position = c(0.7,0.2),
        axis.text.y = element_blank(),
        strip.text = element_text(face = "italic"),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.background = element_blank(),
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=6))

# check observations and covariates
load(here("Processed_data", "ProcessedDetectiondata.RData"))

# species data
species_obs <- (apply (species_table,c(1,3), max,na.rm=T))
species_obs[is.infinite(species_obs)] <- "NS" 

# habitat
grid.arrange (
  data.frame (hab=extract.hab.NA.summary[,"18"],
       species_obs) %>%
  melt(id.vars = "hab") %>%
  mutate(value=as.factor(value)) %>%
  ggplot(aes(x=value,y=hab,fill=value))+
  facet_wrap(~variable,nrow=1)+
  geom_point(position = position_jitterdodge(),
             alpha=0.01)+
  geom_violin(alpha=0.1)+
  theme_minimal()+
    ylab("Cover Pastures")
,
  # habitat
  data.frame (hab=extract.hab.NA.summary[,"35"],
            species_obs) %>%
  melt(id.vars = "hab") %>%
  mutate(value=as.factor(value)) %>%
  ggplot(aes(x=value,y=hab,fill=value))+
  facet_wrap(~variable,nrow=1)+
  geom_point(position = position_jitterdodge(),
             alpha=0.2)+
  geom_violin(alpha=0.1)+
  theme_minimal()+
  ylab("Cover InMarshes")

  ,

  nrow=2)


# elevation
data.frame (hab=altitude_stats$altitude_mean,
            species_obs) %>%
  melt(id.vars = "hab") %>%
  mutate(value=as.factor(value)) %>%
  ggplot(aes(x=value,y=hab,fill=value))+
  facet_wrap(~variable,nrow=1)+
  geom_point(position = position_jitterdodge(),
             alpha=0.05)+
  geom_violin(alpha=0.3)+
  theme_minimal()+
  ylab("Elevation")

# save data -------------------------
save (water,
      extract.water,
      extract.water.summary,
      transf_cells,
      file = here ("Processed_data","Water.RData"))

# save data (slow to obtain)
save (altitude.NA.Proj,
      altitude.NA.cells,
      extract.hab.NA,
      hab.NA,
      transf_cells,
      file = here ("Processed_data","AltitudeHabitatRaw.RData"))

# save for models
save (altitude_stats,
      extract.hab.NA.summary,
      file = here ("Processed_data",
                   "AltitudeHabitatStats.RData"))

