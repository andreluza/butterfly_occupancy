

# ----------------------------------------------------------------

#          Section of plots and maps illustrating the empirical data

# -----------------------------------------------------------------

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
#load(here ("Processed_data",
#           "ProcessedDetectionData.RData"))

# load the objects to be used as input of the multi-season site occupancy models
load (file = here("Processed_data", 
                  "Occupancy_data_spOccupancy.RData"))

# --------------------- Section of plots  ------------------------------

# observations over time
plot_total_records <- data.frame (year  = substr(names(colSums(base_table)),nchar(names(colSums(base_table)))-3,nchar(names(colSums(base_table)))),
                                  records = colSums(base_table)) %>%
  group_by(year) %>%
  summarize(nrecords = sum(records)) %>%
  ggplot(aes(x=year,y=nrecords))+
  geom_point()+
  geom_line(group=1)+
  my_theme+
  labs(x="Year",y="Number of records",
       title="A")+
  theme(axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7,angle=45),
        axis.title = element_text(size=12))

# select data of the "best year" (highest N observations)
surveys_best <- base_table [,grep ("2018",colnames(base_table))]
list_data_surveys <- lapply(seq (1,ncol(surveys_best)), function (i) 
  cbind (cells_NAquitane,
         Surveys =  surveys_best [match(cells_NAquitane$CODE_10KM, rownames(surveys_best)),i],
         Month = colnames(surveys_best)[i]
  )
)

# unlist
list_data_surveys <- do.call(rbind,list_data_surveys)
list_data_surveys$Month <-(factor(list_data_surveys$Month, 
                                             levels=unique(list_data_surveys$Month)))
sum(list_data_surveys$Surveys) # as in the plot_total_records

# map of surveys per month (year 2018)
survey_plot <- ggplot() +
        geom_sf(data= list_data_surveys,
        aes (fill=log(Surveys),
             col=log(Surveys)),
        alpha=0.75)+
        scale_colour_viridis_c(option = "magma",direction=1,na.value = "white")+
        scale_fill_viridis_c(option = "magma",direction=1,na.value = "white")+
        my_theme+
  facet_wrap(~Month,ncol=5)+
        theme(axis.text = element_blank(),
              strip.text.x = element_text(face = "italic",size=8),
              plot.background = element_rect(fill="white"),
              legend.key.size = unit(0.5, 'cm'), #change legend key size
              legend.background = element_blank(),
              legend.key.height = unit(0.5, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.title = element_text(size=10), #change legend title font size
              legend.text = element_text(size=8))

# static maps ---------------------------------
# world map
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
France <- world %>%
  filter (name == "France")

# crop 
xlim<- c(-1.1,-0.2) # c(-5,10)
ylim<- c(44.6,45.1) # c(40,52)

# plot 1
p1a <- ggplot(data = cells_NAquitane) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (cells_NAquitane,
                      det = as.numeric(cells_NAquitane$CODE_10KM %in% dataPointsPolygonClean_subsetSpp$Maille1),
                      Records = (rowSums(base_table)[match(cells_NAquitane$CODE_10KM, rownames(base_table))] 
                      )),
          aes (fill=log(Records),
               col=log(Records)),
          alpha=1)+
  geom_sf(data=dataPointsPolygonClean_subsetSpp %>%
            filter (TaxNomVal %in% sp_list) %>%
            mutate (TaxNomVal=stringr::str_replace(TaxNomVal, " \\s*\\([^\\)]+\\)", ""))
            ,
          fill="orange2",
          colour= "orange2",
          alpha=0.15,size=0.1)+
  facet_wrap(~TaxNomVal,ncol=3)+
  #coord_sf(xlim=xlim,ylim=ylim)+
  #theme_bw() +
  scale_colour_viridis_c(option = "magma",direction=1,end=0.7,na.value = "gray90")+
  scale_fill_viridis_c(option = "magma",direction=1,end=0.7,na.value = "gray90")+
  my_theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
        axis.text.y = element_text(size = 6),
        strip.text.x = element_text(face = "italic",size=9),
        legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.background = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=8), #change legend title font size
        legend.text = element_text(size=8,angle=45))

# number of records per cell (log scale)
p1b <- ggplot(data = cells_NAquitane) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (cells_NAquitane,
                      Records = rowSums(base_table)[match(cells_NAquitane$CODE_10KM, rownames(base_table))]),
          aes (fill=log(Records),
               col=log(Records)),
          alpha=0.75)+
  scale_colour_viridis_c(option = "magma",direction=1,na.value = "gray90")+
  scale_fill_viridis_c(option = "magma",direction=1,na.value = "gray90")+
  #coord_sf(xlim=xlim,ylim=ylim)+
  #theme_bw() +
  #scale_colour_viridis_d(option = "magma",direction=1)+
  theme(#legend.position = c(0.7,0.2),
        axis.text.y = element_blank(),
        strip.text = element_text(face = "italic"),
        legend.key.size = unit(0.8, 'cm'), #change legend key size
        legend.background = element_blank(),
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(0.8, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12,angle=45))+
  labs(title="B")+
  my_theme

# number of surveys months per cell
p1c <- ggplot(data = cells_NAquitane) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (cells_NAquitane,
                      MonthsxYears = rowSums(base_table>0)[match(cells_NAquitane$CODE_10KM, rownames(base_table))]),
          aes (fill=log(MonthsxYears),
               col=log(MonthsxYears)),
          alpha=0.75)+
  scale_colour_viridis_c(option = "magma",direction=-1,na.value = "white")+
  scale_fill_viridis_c(option = "magma",direction=-1,na.value = "white")+
  #coord_sf(xlim=xlim,ylim=ylim)+
  #theme_bw() +
  #scale_colour_viridis_d(option = "magma",direction=1)+
  theme(#legend.position = c(0.7,0.2),
    axis.text.y = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.key.size = unit(0.5, 'cm'), #change legend key size
    legend.background = element_blank(),
    legend.key.height = unit(0.5, 'cm'), #change legend key height
    legend.key.width = unit(0.5, 'cm'), #change legend key width
    legend.title = element_text(size=10), #change legend title font size
    legend.text = element_text(size=8,angle=45))+
  labs(title="B")+
  my_theme

p1c

# distribution of monthly surveys
MontSurv_plot <-table(rowSums(base_table>0)) %>%
  data.frame() %>%
  ggplot (aes(x=(Var1),y=Freq))+
  geom_bar(stat = "identity")+
  labs (y="Number of cells",
        x="Number of Year x Month combinations",
        title="C")+
  theme(axis.text.x = element_text(angle=50,size=6))+
  geom_smooth(aes(x=as.numeric(Var1),y=Freq),method = "glm", se = F, 
              method.args = list(family = "poisson"),
              linetype = "dashed",col="orange")+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=12))+
  scale_y_break(c(200,30000), expand=T,scales=c(1,0.5))+
  #my_theme + 
  theme(axis.text.x = element_text(angle=90,size=8),
                   panel.background =  element_rect(fill="white")
                    )

MontSurv_plot

# plot visits across months
month_plot <- data.frame("NObs"=colSums(base_table)) %>%
  mutate("Survey" = rownames(.),
         "Order" = seq (1,nrow(.)),
         "Month" = substr(colnames(base_table),1,nchar(colnames(base_table))-5),
         "Group" = substr (rownames(.),nchar(rownames(.))-4,nchar(rownames(.)))) %>%
  ggplot (aes(x=reorder(Survey, Order),y=NObs))+
  geom_bar(stat = "identity")+
  #geom_smooth(aes(x=reorder(Survey, Order),y=NObs,group = Group),
  #            method = "glm", se = F, 
  #            formula = 'y ~ poly(x,2)',
  #            method.args = list(family = "poisson"),
  #            linetype = "solid",col="orange",linewidth=0.25)+
  labs (x="",
        y="Number of Observations per Month and Year",
        title="A")+
  #my_theme+
  theme(axis.text.x = element_text(angle=0,size=10,hjust = 1),
        axis.text.y = element_text(size=8)) + 
  coord_polar(start = 0,clip="off",direction = 1,theta="x")+
  scale_x_discrete(breaks = paste0("juin ", seq(2000,2023)), labels = 2000:2023) 
  
month_plot

# distribution of visits across cells
cols_base <- substr(colnames(base_table), 
       nchar(colnames(base_table))-3,
       nchar(colnames(base_table)))

# surveys per year (2018) - the best year
i= "2018" # 
data_visits <- #lapply (unique(cols_base), function (i) # activate if need to evaluate across years
        
    data.frame(
        table(
          rowSums(base_table[,grep(i,colnames(base_table))]>0)
          )
      ) %>%
        mutate(Freq = as.numeric(Freq),
               Var1 = as.numeric(Var1)-1#,
               #Year = as.numeric(i)
               ) 
#)

# melt
#data_visits <- do.call(rbind, data_visits)

#data_visits %>%
#  filter(Year == 2023) %>%
#  summarize(sum(Freq[-1]))

# average number of cells sampled per year
ncells_year <- dataPointsPolygonClean_subsetSpp %>%
  group_by(Year) %>%
  summarise (ncells = length(unique(Maille1))) %>%
  st_drop_geometry()%>%
  as.data.frame()

# calculate the percentages
data_visits$Perc <- data_visits$Freq/length(unique(cells_NAquitane$CODE_10KM))
  
# average across years
data_visits <- data_visits %>%
  group_by(Var1) %>%
  summarise(Freq=(mean(Freq)),
            Perc= mean(Perc))
# check the totals
sum(data_visits$Perc)
sum(data_visits$Freq)

# plot
cells_plot <- data_visits[1:nrow(data_visits),] %>%
  ggplot (aes(x=as.factor(Var1),y=Freq))+
  geom_bar(stat = "identity")+
  labs (y="Number of cells",
        x="Number of Months (2018)",
        title="C")+
  theme(axis.text.x = element_text(angle=50,size=6))+
  geom_smooth(aes(x=Var1+1,y=Freq),method = "glm", se = F, 
              method.args = list(family = "poisson"),
              linetype = "dashed",col="orange")+
  geom_text(aes(x=Var1+1,y=Freq+2000,label = round(Perc,3)*100))+
  theme(axis.text = element_text(size=7),
        axis.title = element_text(size=12),
        panel.background = element_rect(fill="white"))+
  #scale_y_break(c(500,10000), expand=T,scales=c(0.3,8))+
  my_theme

# doser and stoudt distribution of sites x secondary occasions per year
DS_design <- data.frame (Numb = c(1200,1200*0.1),
            SecOcc = (c(1,2)),
            Var1= c(1,0.1)) %>%
  ggplot(aes (x=SecOcc, y=Numb))+
  geom_bar(stat = "identity")+
  labs (#y="Number of cells",
        x=" ",
        title="D")+
  #geom_text(aes(x=Var1+1,y=200,label = round(Var1,3)))+
  #theme(axis.text.x = element_text(angle=50,size=6))+
  geom_smooth(aes(x=SecOcc,y=Numb),method = "glm", se = F, 
              method.args = list(family = "poisson"),
              linetype = "dashed",col="orange")+
  geom_text(aes(x=SecOcc,y=Numb+40,label = round(Var1,2)*100),size=3.5)+
  scale_x_continuous(breaks=c(1,2))+
  my_theme+
  theme(#plot.title = element_text(face="italic"),
        axis.text = element_text(angle=0,size=6),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
# table of surveys
table_cells <- melt(
  
  c ("n=0" = data_visits[1,"Freq"], # single survey cells
     "n=1" = data_visits[2,"Freq"], # single survey cells
     "n=2" =  data_visits[3,"Freq"], # two surveys
     "n=3" = data_visits[4,"Freq"], # three surveys
     "n=4" = data_visits[5,"Freq"], # four surveys
     "n=5" = data_visits[6,"Freq"], # five surveys
     "n=6" = data_visits[7,"Freq"], # six surveys
     "n=7" = data_visits[8,"Freq"], # seven surveys
     "n=8" = data_visits[9,"Freq"], # eight surveys
     "n=9" = data_visits[10,"Freq"], # nine surveys
     "n=10" = data_visits[11,"Freq"] # ten surveys
     
  ) 
) %>%
  dplyr::rename("Cells" = value,
                "Surveys"=L1) %>%
  
  mutate (Cells = round (Cells,0),
          Surveys = gsub (".Freq","",Surveys)
          
  ) %>%
  
  kable(col.names = c("#Cells","Surveyed months"),
        align="c",
        format = "rst") 

# transform into a tableGrob
table_cells<-tableGrob(table_cells)

# arrange
comp_sites <- grid.arrange(
        print(cells_plot),
             table_cells, 
              DS_design,
             nrow=7,ncol=6,
             layout_matrix= rbind(c(1,1,1,1,3,3),
                                  c(1,1,1,1,3,3),
                                  c(1,1,2,2,3,3),
                                  c(1,1,2,2,3,3),
                                  c(1,1,1,1,3,3),
                                  c(1,1,1,1,3,3),
                                  c(1,1,1,1,3,3)))

# plot observations across months
species_plot <- melt (apply(species_table,c(2,3),sum,na.rm=T),as.is=T) %>%
  mutate("Order" = seq (1,nrow(.)),
         "Month" = rep(substr(colnames(species_table),1,nchar(colnames(species_table))-5),6),
         "Group" = substr (as.character(X1),nchar(as.character(X1))-4,
                           nchar(as.character(X1)))) %>%
  mutate (X2=stringr::str_replace(X2, " \\s*\\([^\\)]+\\)", "")) %>%
  ggplot (aes(x=reorder(X1, Order),y=value))+
  facet_wrap(~X2,scales="free_y",nrow = 2)+
  #geom_smooth(aes(x=reorder(X1, Order),y=value,group = Group),
  #            method = "glm", se = F, 
  #            formula = 'y ~ poly(x,2)',
  #            method.args = list(family = "poisson"),
  #            linetype = "solid",col="orange",linewidth=0.25)+
  
  geom_bar(stat = "identity")+
  labs (x="",
        y="Number of Observations per Month and Year",
        title="B")+
  #my_theme+
  theme(axis.text.x = element_text(angle=0,size=4,hjust = 1),
        strip.text = element_text(face = "italic"),
        axis.text.y = element_text(size=5)) + 
  coord_polar(start = 0,clip="off",direction = 1,theta="x")+
  scale_x_discrete(breaks = paste0("juin ", seq(2000,2023)), labels = 2000:2023)   

# plot
png(here ("figures", "empirical", "map_spp_visits.png"),
    width = 30, height = 30, units = "cm",res=400)
  
  grid.arrange(month_plot+labs(title="A"),
               species_plot+labs(title="B"),
               comp_sites,
               p1a+labs(title="E"),
               ncol=6,nrow=7,
               layout_matrix = rbind (c(1,1,1,2,2,2),
                                      c(1,1,1,2,2,2),
                                      c(1,1,1,4,4,4),
                                      c(3,3,3,4,4,4),
                                      c(3,3,3,4,4,4),
                                      c(3,3,3,4,4,4)
               )
)

dev.off()

# plot of detection maps in high res
png(here ("figures", "empirical", "map_spp_highRes.png"),
    width = 30, height = 20, units = "cm",res=600)

             
  p1a
             
             
dev.off()

# histograms
png(here ("figures", "empirical", "records_surveyMonthYears.png"),
    width = 18, height = 18, units = "cm",res=300,bg = "gray")
  
  grid.arrange(plot_total_records,
               p1c,
               print(MontSurv_plot),
               ncol=4,nrow=5,
               layout_matrix = rbind (c(1,1,2,2),
                                      c(1,1,2,2),
                                      c(1,1,2,2),
                                      c(3,3,3,3),
                                      c(3,3,3,3)))

dev.off()

# maps of Bordeaux's buffer ----------------------------

p_BM <- ggplot(data = cells_buffer_bordeaux) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (cells_buffer_bordeaux,
                      MonthsxYears = rowSums(base_table>0)[match(cells_buffer_bordeaux$CODE_10KM, names (rowSums(base_table>0)))]),
          aes (fill=MonthsxYears,
               col=MonthsxYears),
          alpha=0.75)+
  scale_colour_viridis_c(option = "magma",direction=-1,na.value = "gray90")+
  scale_fill_viridis_c(option = "magma",direction=-1,na.value = "gray90")+
  #coord_sf(xlim=xlim,ylim=ylim)+
  #theme_bw() +
  #scale_colour_viridis_d(option = "magma",direction=1)+
  theme(#legend.position = c(0.7,0.2),
    axis.text.y = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.key.size = unit(0.8, 'cm'), #change legend key size
    legend.background = element_blank(),
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(0.8, 'cm'), #change legend key width
    legend.title = element_text(size=14), #change legend title font size
    legend.text = element_text(size=12,angle=45))+
  labs(title="A")+
  my_theme+
  # add Bordeaux
  geom_sf(data=cbind (cells_buffer_bordeaux,
                      MonthsxYears = rowSums(base_table>0)[match(cells_buffer_bordeaux$CODE_10KM, names (rowSums(base_table>0)))]) %>%
            filter (NOM_COMM == "BORDEAUX"),
          aes (fill=MonthsxYears),
          col="black",
            alpha=0.75)

p_BM

# total nomber of records
p_BM_obs <- ggplot(data = cells_buffer_bordeaux) +
  geom_sf(fill="white")+
  geom_sf(data=cbind (cells_buffer_bordeaux,
                      Records = rowSums(base_table)[match(cells_buffer_bordeaux$CODE_10KM, names (rowSums(base_table)))]),
          aes (fill=Records,
               col=Records),
          alpha=0.75)+
  scale_colour_viridis_c(option = "magma",direction=-1,na.value = "gray90")+
  scale_fill_viridis_c(option = "magma",direction=-1,na.value = "gray90")+
  #coord_sf(xlim=xlim,ylim=ylim)+
  #theme_bw() +
  #scale_colour_viridis_d(option = "magma",direction=1)+
  theme(#legend.position = c(0.7,0.2),
    axis.text.y = element_blank(),
    strip.text = element_text(face = "italic"),
    legend.key.size = unit(0.8, 'cm'), #change legend key size
    legend.background = element_blank(),
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(0.8, 'cm'), #change legend key width
    legend.title = element_text(size=14), #change legend title font size
    legend.text = element_text(size=12,angle=45))+
  labs(title="B")+
  my_theme+
  # add Bordeaux
  geom_sf(data=cbind (cells_buffer_bordeaux,
                      Records = rowSums(base_table>0)[match(cells_buffer_bordeaux$CODE_10KM, names (rowSums(base_table>0)))]) %>%
            filter (NOM_COMM == "BORDEAUX"),
          aes (fill=Records),
          col="black",
          alpha=0.75)

p_BM_obs

# save Bordeaux Metropole maps
# histograms
png(here ("figures", "empirical", "Bordeaux_metropole_data.png"),
    width = 30, height = 20, units = "cm", res= 600)

grid.arrange(p_BM,
             p_BM_obs,
             ncol=4,nrow=4,
             layout_matrix = rbind (c(1,1,2,2),
                                    c(1,1,2,2),
                                    c(1,1,2,2),
                                    c(1,1,2,2)))

dev.off()

# clean work space
rm(list=ls())
# end

