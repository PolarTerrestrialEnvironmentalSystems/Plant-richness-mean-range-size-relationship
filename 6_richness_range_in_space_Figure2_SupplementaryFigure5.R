################################################################################

# Author: Ying Liu,Simeon Lisovski
# Email:  ying.liu@awi.de, simeon.lisovski@awi.de
# Date:   May 21 2024

################################################################################
# Description: Here we used the modern plant distribution data downloaded from
# Global Biodiversity Information Facility database (GBIF, https://www.gbif.org/)
# to calculate plant richness to range-size relationship in space, also 
# we used the sedaDNA modern timeslice data to calculate the relationship in space,
# and compare the modern plant taxa range calculated based on 7 lake region and 
# the larger northeast Siberia and Alaska region
################################################################################ 
# load the data and packages

library("dplyr")
library("sf")
library("ggplot2")
library("ggpubr")
library("cowplot")

mydata1<-read.csv("modern_plat_occurance_from_GBIF.csv")
resampleReads_final<-read.csv("resampleReads_final.csv")
lake_coordinate<-read.csv("lake_coordinate.csv")
################################################################################

# dataframe transfer into the sf type, and set a Lambert Azimuthal Equal Area
# (LAEA) projection, the median value is base on the median value of modern 
# plant distribution,and divide the region into 200km*200km grid cells

mydata1<-mydata1[!duplicated(mydata1),]
modern  <- mydata1 %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                crs = 4326) %>% st_shift_longitude()
proj <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                median((modern %>% st_coordinates())[,1]), 
                median((modern %>% st_coordinates())[,2]))

grid   <- modern %>% st_transform(proj) %>% st_bbox() %>% 
  st_as_sfc(crs = proj) %>% 
  st_buffer(150*1000) %>% 
  st_make_grid(200*1000)

df_grid <- data.frame(st_coordinates(grid))

# get the coordination of four vertices

df_grid1 <-data.frame(df_grid$X,df_grid$Y,df_grid$L2)
colnames(df_grid1)<-c("X","Y","L")
df_grid1<-df_grid1[!duplicated(df_grid1),]

df_grid2<-df_grid1%>%group_by(L)%>%
  summarise(X1=min(X),Y1=min(Y),X2=max(X),Y2=max(Y))

# transfer the modern plant coordination into the projection
modern1  <- modern %>% st_transform(proj) %>% st_as_sfc(crs = proj)
modern1<-data.frame(st_coordinates(modern1))


# matching modern species with grid cells
modern1$species<-modern$species

matching_rows_all<-NULL

for (i in 1:nrow(df_grid2)) {
  print(i)
  grid<-df_grid2[i,]
  
  matching_rows <-modern1[(modern1$X)>(grid$X1)&((modern1$X)<(grid$X2))&
                            ((modern1$Y)>(grid$Y1))&((modern1$Y)<(grid$Y2)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$L<-grid$L
    matching_rows_all<-rbind(matching_rows_all,matching_rows)
    
  } else {
    print("fail")
  }
  
}

# claculate the richness and mean range size for every grid

species_range<-matching_rows_all%>%group_by(species)%>%
  summarise(range=n_distinct(L))

matching_rows_all1<-merge(matching_rows_all,species_range,by="species")%>%
  dplyr::select(species,L,range)%>%distinct()


grid_richness_range<-matching_rows_all1%>%group_by(L)%>%
  summarise(richness=n_distinct(species),range=mean(range))
# delete the richness less than 20 because the the plant survey is not so 
# complete, also some grid cover the islands and ocean region

grid_richness_range<-grid_richness_range%>%filter(richness>20)

p1 <- ggplot(data = grid_richness_range, aes(x = richness, y = range)) + 
  geom_point(color="#313695",alpha=0.8) + geom_smooth(method = lm,color="#2166ac") +
  labs(y="mean range size of the 200km*200km grid cell",
       x="richness of the 200km*200km grid cell") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))+
  stat_cor(method = 'spearman', aes(x = richness, y = range))
p1

A<-lm(range~richness,data = grid_richness_range)
summary(A)

# ggsave(p1,filename="richness_range_relationship_in_space_sibala_region.pdf",
#        width = 8,height=8)

# calculate the spatial richness-mean range size relationship in space

modern_sedDNA<-resampleReads_final%>%filter(time_slice==30|time_slice==29)%>%
  group_by(round,lake,NUC_SEQ)%>%summarise(count=sum(N))

NUC_SEQ_range<-modern_sedDNA%>%group_by(round, NUC_SEQ)%>%
  summarise(range=n_distinct(lake))

modern_sedDNA<-merge(modern_sedDNA,NUC_SEQ_range,by=c("round","NUC_SEQ"))


sedDNA_richrange<-modern_sedDNA%>%group_by(round, lake)%>%
  summarise(richness=n_distinct(NUC_SEQ),mean_range=mean(range))


p2 <- ggplot(data = sedDNA_richrange, aes(x = richness, y = mean_range,col=lake)) + 
  geom_point(alpha=1) + geom_smooth(method = lm,color="#2166ac") +
  geom_text(aes(label = lake), hjust = 0, vjust = 0) +
  labs(y="mean range size of 7 lake region based on sedaDNA",
       x="richness of 7 lake region based on sedaDNA") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"),
        legend.position = "none")+
  #stat_cor(method = 'spearman', aes(x = richness, y = mean_range))+
  scale_colour_manual(values = c(Billyakh="#4a6990",Btoko="#8f7700",
      E5="#003c67",Emanda="#a7302f",Ilirney18="#3b3b3b",
      Lele="#7aa6dc",Rauchuagytgyn="#cd534c"))+
      guides(fill = FALSE)
   
p2

B<-lm(mean_range~richness,data = sedDNA_richrange)
summary(B)

# compare the same taxa, get range size from larger northeast Siberia and Alaska
# region, with range size from the seven lake region


# set a Lambert Azimuthal Equal Area(LAEA) projection, the median value is base 
# on lake disreibution,and divide the region into 200km*200km grid cells

lakes  <- lake_coordinate %>% filter(!duplicated(Longitude)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_shift_longitude()

proj_lake <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                     median((lakes %>% st_coordinates())[,1]), 
                     median((lakes %>% st_coordinates())[,2]))


# transfer the modern plant distribution projection
modern1  <- mydata1 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(proj_lake)

modern1<-data.frame(st_coordinates(modern1))

modern1$species<-mydata1$species

# transfer the lake coordination in same projection
lake1  <- lakes %>% st_transform(proj_lake) %>% st_as_sfc(crs =proj_lake)
lake1<-data.frame(st_coordinates(lake1))
lake1$lake<-lakes$lake

lake1$X1<-lake1$X-100*1000
lake1$X2<-lake1$X+100*1000
lake1$Y1<-lake1$Y-100*1000
lake1$Y2<-lake1$Y+100*1000


matching_modern_plant_lake<-NULL
for (i in 1:nrow(lake1)) {
  print(i)
  lake2<-lake1[i,]
  
  matching_rows <-modern1[((modern1$X)>(lake2$X1)) & ((modern1$X)<(lake2$X2))& ((modern1$Y)<(lake2$Y2))& ((modern1$Y)>(lake2$Y1)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$lake<-lake2$lake
    matching_modern_plant_lake<-rbind(matching_modern_plant_lake,matching_rows)
    
  } else {
    print("fail")
  }
  
}


lake_species_range<-matching_modern_plant_lake%>%
  group_by(species)%>%summarise(range_lake=n_distinct(lake))

compare<-merge(lake_species_range,species_range,by="species")


p<-ggplot(compare, aes(x=range_lake, y=range,group=range_lake)) +
  geom_boxplot(color="#2166ac")+
  stat_summary(fun=mean, geom='point',color="#053061", shape=20, size=3)+
  labs(x = "plant taxa range size around seven lakes region",
       y = "plant taxa range size in northeast Siberia and Alaska region")+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))
p

#ggsave(p,filename="compare range in Siberia and lake.pdf",width = 10,height=8)

plot<-plot_grid(p2,p1,p,
                labels=c("b","c","d"),
                label_x = 0.1,label_y = 1,
                ncol = 2,
                align = "v")
plot
#Figure2 b,c,d
ggsave(plot,filename="Figure2_modern_sibala_plant_compare range in Siberia and lake.pdf",
       width = 10,height=8)

################################################################################ 
# For supplementary, try 100km*100km-500km*500km resolution
################################################################################

# try 100km*100km resolution
mydata1<-read.csv("modern_plat_occurance_from_GBIF.csv")
resampleReads_final<-read.csv("resampleReads_final.csv")
lake_coordinate<-read.csv("lake_coordinate.csv")

mydata1<-mydata1[!duplicated(mydata1),]
modern  <- mydata1 %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                crs = 4326) %>% st_shift_longitude()
proj <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                median((modern %>% st_coordinates())[,1]), 
                median((modern %>% st_coordinates())[,2]))

grid   <- modern %>% st_transform(proj) %>% st_bbox() %>% 
  st_as_sfc(crs = proj) %>% 
  st_buffer(50*1000) %>% 
  st_make_grid(100*1000)

df_grid <- data.frame(st_coordinates(grid))

# get the coordination of four vertices

df_grid1 <-data.frame(df_grid$X,df_grid$Y,df_grid$L2)
colnames(df_grid1)<-c("X","Y","L")
df_grid1<-df_grid1[!duplicated(df_grid1),]

df_grid2<-df_grid1%>%group_by(L)%>%
  summarise(X1=min(X),Y1=min(Y),X2=max(X),Y2=max(Y))

# transfer the modern plant coordination into the projection
modern1  <- modern %>% st_transform(proj) %>% st_as_sfc(crs = proj)
modern1<-data.frame(st_coordinates(modern1))


# matching modern species with grid cells
modern1$species<-modern$species

matching_rows_all<-NULL

for (i in 1:nrow(df_grid2)) {
  print(i)
  grid<-df_grid2[i,]
  
  matching_rows <-modern1[(modern1$X)>(grid$X1)&((modern1$X)<(grid$X2))&
                            ((modern1$Y)>(grid$Y1))&((modern1$Y)<(grid$Y2)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$L<-grid$L
    matching_rows_all<-rbind(matching_rows_all,matching_rows)
    
  } else {
    print("fail")
  }
  
}

# claculate the richness and mean range size for every grid

species_range<-matching_rows_all%>%group_by(species)%>%
  summarise(range=n_distinct(L))

matching_rows_all1<-merge(matching_rows_all,species_range,by="species")%>%
  dplyr::select(species,L,range)%>%distinct()


grid_richness_range<-matching_rows_all1%>%group_by(L)%>%
  summarise(richness=n_distinct(species),range=mean(range))
# delete the richness less than 20 because the the plant survey is not so 
# complete, also some grid cover the islands and ocean region

grid_richness_range<-grid_richness_range%>%filter(richness>20)

#Supplementary Figure5a
p11 <- ggplot(data = grid_richness_range, aes(x = richness, y = range)) + 
  geom_point(color="#313695",alpha=0.8) + geom_smooth(method = lm,color="#2166ac") +
  labs(y="modern plant mean range size of the 100km*100km grid cell",
       x="modern plant richness of the 100km*100km grid cell") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))+
  stat_cor(method = 'spearman', aes(x = richness, y = range))
p11

A<-lm(range~richness,data = grid_richness_range)
summary(A)

# compare the same taxa, get range size from larger northeast Siberia and Alaska
# region, with range size from the seven lake region

# set a Lambert Azimuthal Equal Area(LAEA) projection, the median value is base 
# on lake disreibution,and divide the region into 100km*100km grid cells

lakes  <- lake_coordinate %>% filter(!duplicated(Longitude)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_shift_longitude()

proj_lake <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                     median((lakes %>% st_coordinates())[,1]), 
                     median((lakes %>% st_coordinates())[,2]))


# transfer the modern plant distribution projection
modern1  <- mydata1 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(proj_lake)

modern1<-data.frame(st_coordinates(modern1))

modern1$species<-mydata1$species

# transfer the lake coordination in same projection
lake1  <- lakes %>% st_transform(proj_lake) %>% st_as_sfc(crs =proj_lake)
lake1<-data.frame(st_coordinates(lake1))
lake1$lake<-lakes$lake

lake1$X1<-lake1$X-50*1000
lake1$X2<-lake1$X+50*1000
lake1$Y1<-lake1$Y-50*1000
lake1$Y2<-lake1$Y+50*1000


matching_modern_plant_lake<-NULL
for (i in 1:nrow(lake1)) {
  print(i)
  lake2<-lake1[i,]
  
  matching_rows <-modern1[((modern1$X)>(lake2$X1)) & ((modern1$X)<(lake2$X2))& ((modern1$Y)<(lake2$Y2))& ((modern1$Y)>(lake2$Y1)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$lake<-lake2$lake
    matching_modern_plant_lake<-rbind(matching_modern_plant_lake,matching_rows)
    
  } else {
    print("fail")
  }
  
}


lake_species_range<-matching_modern_plant_lake%>%
  group_by(species)%>%summarise(range_lake=n_distinct(lake))

compare<-merge(lake_species_range,species_range,by="species")

#Supplementary Figure5b
p12<-ggplot(compare, aes(x=range_lake, y=range,group=range_lake)) +
  geom_boxplot(color="#2166ac")+
  stat_summary(fun=mean, geom='point',color="#053061", shape=20, size=3)+
  labs(x = "modern plant taxa range size around seven lakes region",
       y = "modern plant taxa range size in northeast Siberia and Alaska region")+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))
#Supplementary Figure5b
p12


plot<-plot_grid(p11,p12,
                labels=c("a","b"),
                label_x = 0.1,label_y = 1,
                ncol = 2,
                align = "v")
plot

ggsave(plot,filename="modern_sibala_plant_compare range in Siberia and lake_100km.pdf",
       width = 13,height=8)

# try 300km*300km resolution

mydata1<-read.csv("modern_plat_occurance_from_GBIF.csv")
resampleReads_final<-read.csv("resampleReads_final.csv")
lake_coordinate<-read.csv("lake_coordinate.csv")

mydata1<-mydata1[!duplicated(mydata1),]
modern  <- mydata1 %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                crs = 4326) %>% st_shift_longitude()
proj <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                median((modern %>% st_coordinates())[,1]), 
                median((modern %>% st_coordinates())[,2]))

grid   <- modern %>% st_transform(proj) %>% st_bbox() %>% 
  st_as_sfc(crs = proj) %>% 
  st_buffer(150*1000) %>% 
  st_make_grid(300*1000)

df_grid <- data.frame(st_coordinates(grid))

# get the coordination of four vertices

df_grid1 <-data.frame(df_grid$X,df_grid$Y,df_grid$L2)
colnames(df_grid1)<-c("X","Y","L")
df_grid1<-df_grid1[!duplicated(df_grid1),]

df_grid2<-df_grid1%>%group_by(L)%>%
  summarise(X1=min(X),Y1=min(Y),X2=max(X),Y2=max(Y))

# transfer the modern plant coordination into the projection
modern1  <- modern %>% st_transform(proj) %>% st_as_sfc(crs = proj)
modern1<-data.frame(st_coordinates(modern1))


# matching modern species with grid cells
modern1$species<-modern$species

matching_rows_all<-NULL

for (i in 1:nrow(df_grid2)) {
  print(i)
  grid<-df_grid2[i,]
  
  matching_rows <-modern1[(modern1$X)>(grid$X1)&((modern1$X)<(grid$X2))&
                            ((modern1$Y)>(grid$Y1))&((modern1$Y)<(grid$Y2)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$L<-grid$L
    matching_rows_all<-rbind(matching_rows_all,matching_rows)
    
  } else {
    print("fail")
  }
  
}

# claculate the richness and mean range size for every grid

species_range<-matching_rows_all%>%group_by(species)%>%
  summarise(range=n_distinct(L))

matching_rows_all1<-merge(matching_rows_all,species_range,by="species")%>%
  dplyr::select(species,L,range)%>%distinct()


grid_richness_range<-matching_rows_all1%>%group_by(L)%>%
  summarise(richness=n_distinct(species),range=mean(range))
# delete the richness less than 20 because the the plant survey is not so 
# complete, also some grid cover the islands and ocean region

grid_richness_range<-grid_richness_range%>%filter(richness>20)

p13 <- ggplot(data = grid_richness_range, aes(x = richness, y = range)) + 
  geom_point(color="#313695",alpha=0.8) + geom_smooth(method = lm,color="#2166ac") +
  labs(y="modern plant mean range size of the 300km*300km grid cell",
       x="modern plant richness of the 300km*300km grid cell") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))+
  stat_cor(method = 'spearman', aes(x = richness, y = range))

#Supplementary Figure5c
p13

A<-lm(range~richness,data = grid_richness_range)
summary(A)

# compare the same taxa, get range size from larger northeast Siberia and Alaska
# region, with range size from the seven lake region


# set a Lambert Azimuthal Equal Area(LAEA) projection, the median value is base 
# on lake disreibution,and divide the region into 300km*300km grid cells

lakes  <- lake_coordinate %>% filter(!duplicated(Longitude)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_shift_longitude()

proj_lake <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                     median((lakes %>% st_coordinates())[,1]), 
                     median((lakes %>% st_coordinates())[,2]))


# transfer the modern plant distribution projection
modern1  <- mydata1 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(proj_lake)

modern1<-data.frame(st_coordinates(modern1))

modern1$species<-mydata1$species

# transfer the lake coordination in same projection
lake1  <- lakes %>% st_transform(proj_lake) %>% st_as_sfc(crs =proj_lake)
lake1<-data.frame(st_coordinates(lake1))
lake1$lake<-lakes$lake

lake1$X1<-lake1$X-150*1000
lake1$X2<-lake1$X+150*1000
lake1$Y1<-lake1$Y-150*1000
lake1$Y2<-lake1$Y+150*1000


matching_modern_plant_lake<-NULL
for (i in 1:nrow(lake1)) {
  print(i)
  lake2<-lake1[i,]
  
  matching_rows <-modern1[((modern1$X)>(lake2$X1)) & ((modern1$X)<(lake2$X2))& ((modern1$Y)<(lake2$Y2))& ((modern1$Y)>(lake2$Y1)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$lake<-lake2$lake
    matching_modern_plant_lake<-rbind(matching_modern_plant_lake,matching_rows)
    
  } else {
    print("fail")
  }
  
}


lake_species_range<-matching_modern_plant_lake%>%
  group_by(species)%>%summarise(range_lake=n_distinct(lake))

compare<-merge(lake_species_range,species_range,by="species")


p14<-ggplot(compare, aes(x=range_lake, y=range,group=range_lake)) +
  geom_boxplot(color="#2166ac")+
  stat_summary(fun=mean, geom='point',color="#053061", shape=20, size=3)+
  labs(x = "modern plant taxa range size around seven lakes region",
       y = "modern plant taxa range size in northeast Siberia and Alaska region")+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))

#Supplementary Figure5d
p14

plot<-plot_grid(p13,p14,
                labels=c("a","b"),
                label_x = 0.1,label_y = 1,
                ncol = 2,
                align = "v")
plot

ggsave(plot,filename="modern_sibala_plant_compare range in Siberia and lake_300km.pdf",
       width = 13,height=8)

# try 400km*400km resolution

mydata1<-read.csv("modern_plat_occurance_from_GBIF.csv")
resampleReads_final<-read.csv("resampleReads_final.csv")
lake_coordinate<-read.csv("lake_coordinate.csv")

mydata1<-mydata1[!duplicated(mydata1),]
modern  <- mydata1 %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                crs = 4326) %>% st_shift_longitude()
proj <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                median((modern %>% st_coordinates())[,1]), 
                median((modern %>% st_coordinates())[,2]))

grid   <- modern %>% st_transform(proj) %>% st_bbox() %>% 
  st_as_sfc(crs = proj) %>% 
  st_buffer(200*1000) %>% 
  st_make_grid(400*1000)

df_grid <- data.frame(st_coordinates(grid))

# get the coordination of four vertices

df_grid1 <-data.frame(df_grid$X,df_grid$Y,df_grid$L2)
colnames(df_grid1)<-c("X","Y","L")
df_grid1<-df_grid1[!duplicated(df_grid1),]

df_grid2<-df_grid1%>%group_by(L)%>%
  summarise(X1=min(X),Y1=min(Y),X2=max(X),Y2=max(Y))

# transfer the modern plant coordination into the projection
modern1  <- modern %>% st_transform(proj) %>% st_as_sfc(crs = proj)
modern1<-data.frame(st_coordinates(modern1))


# matching modern species with grid cells
modern1$species<-modern$species

matching_rows_all<-NULL

for (i in 1:nrow(df_grid2)) {
  print(i)
  grid<-df_grid2[i,]
  
  matching_rows <-modern1[(modern1$X)>(grid$X1)&((modern1$X)<(grid$X2))&
                            ((modern1$Y)>(grid$Y1))&((modern1$Y)<(grid$Y2)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$L<-grid$L
    matching_rows_all<-rbind(matching_rows_all,matching_rows)
    
  } else {
    print("fail")
  }
  
}

# claculate the richness and mean range size for every grid

species_range<-matching_rows_all%>%group_by(species)%>%
  summarise(range=n_distinct(L))

matching_rows_all1<-merge(matching_rows_all,species_range,by="species")%>%
  dplyr::select(species,L,range)%>%distinct()


grid_richness_range<-matching_rows_all1%>%group_by(L)%>%
  summarise(richness=n_distinct(species),range=mean(range))
# delete the richness less than 20 because the the plant survey is not so 
# complete, also some grid cover the islands and ocean region

grid_richness_range<-grid_richness_range%>%filter(richness>20)

p15 <- ggplot(data = grid_richness_range, aes(x = richness, y = range)) + 
  geom_point(color="#313695",alpha=0.8) + geom_smooth(method = lm,color="#2166ac") +
  labs(y="modern plant mean range size of the 400km*400km grid cell",
       x="modern plant richness of the 400km*400km grid cell") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))+
  stat_cor(method = 'spearman', aes(x = richness, y = range))

#Supplementary Figure5e
p15

A<-lm(range~richness,data = grid_richness_range)
summary(A)

# compare the same taxa, get range size from larger northeast Siberia and Alaska
# region, with range size from the seven lake region


# set a Lambert Azimuthal Equal Area(LAEA) projection, the median value is base 
# on lake disreibution,and divide the region into 400km*400km grid cells

lakes  <- lake_coordinate %>% filter(!duplicated(Longitude)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_shift_longitude()

proj_lake <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                     median((lakes %>% st_coordinates())[,1]), 
                     median((lakes %>% st_coordinates())[,2]))


# transfer the modern plant distribution projection
modern1  <- mydata1 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(proj_lake)

modern1<-data.frame(st_coordinates(modern1))

modern1$species<-mydata1$species

# transfer the lake coordination in same projection
lake1  <- lakes %>% st_transform(proj_lake) %>% st_as_sfc(crs =proj_lake)
lake1<-data.frame(st_coordinates(lake1))
lake1$lake<-lakes$lake

lake1$X1<-lake1$X-200*1000
lake1$X2<-lake1$X+200*1000
lake1$Y1<-lake1$Y-200*1000
lake1$Y2<-lake1$Y+200*1000


matching_modern_plant_lake<-NULL
for (i in 1:nrow(lake1)) {
  print(i)
  lake2<-lake1[i,]
  
  matching_rows <-modern1[((modern1$X)>(lake2$X1)) & ((modern1$X)<(lake2$X2))& ((modern1$Y)<(lake2$Y2))& ((modern1$Y)>(lake2$Y1)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$lake<-lake2$lake
    matching_modern_plant_lake<-rbind(matching_modern_plant_lake,matching_rows)
    
  } else {
    print("fail")
  }
  
}


lake_species_range<-matching_modern_plant_lake%>%
  group_by(species)%>%summarise(range_lake=n_distinct(lake))

compare<-merge(lake_species_range,species_range,by="species")


p16<-ggplot(compare, aes(x=range_lake, y=range,group=range_lake)) +
  geom_boxplot(color="#2166ac")+
  stat_summary(fun=mean, geom='point',color="#053061", shape=20, size=3)+
  labs(x = "modern plant taxa range size around seven lakes region",
       y = "modern plant taxa range size in northeast Siberia and Alaska region")+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))

#Supplementary Figure5f
p16


plot<-plot_grid(p15,p16,
                labels=c("a","b"),
                label_x = 0.1,label_y = 1,
                ncol = 2,
                align = "v")
plot

ggsave(plot,filename="modern_sibala_plant_compare range in Siberia and lake_400km.pdf",
       width = 13,height=8)


# try 500km*500km resolution

mydata1<-read.csv("modern_plat_occurance_from_GBIF.csv")
resampleReads_final<-read.csv("resampleReads_final.csv")
lake_coordinate<-read.csv("lake_coordinate.csv")

mydata1<-mydata1[!duplicated(mydata1),]
modern  <- mydata1 %>% st_as_sf(coords = c("Longitude", "Latitude"), 
                                crs = 4326) %>% st_shift_longitude()
proj <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                median((modern %>% st_coordinates())[,1]), 
                median((modern %>% st_coordinates())[,2]))

grid   <- modern %>% st_transform(proj) %>% st_bbox() %>% 
  st_as_sfc(crs = proj) %>% 
  st_buffer(200*1000) %>% 
  st_make_grid(500*1000)

df_grid <- data.frame(st_coordinates(grid))

# get the coordination of four vertices

df_grid1 <-data.frame(df_grid$X,df_grid$Y,df_grid$L2)
colnames(df_grid1)<-c("X","Y","L")
df_grid1<-df_grid1[!duplicated(df_grid1),]

df_grid2<-df_grid1%>%group_by(L)%>%
  summarise(X1=min(X),Y1=min(Y),X2=max(X),Y2=max(Y))

# transfer the modern plant coordination into the projection
modern1  <- modern %>% st_transform(proj) %>% st_as_sfc(crs = proj)
modern1<-data.frame(st_coordinates(modern1))


# matching modern species with grid cells
modern1$species<-modern$species

matching_rows_all<-NULL

for (i in 1:nrow(df_grid2)) {
  print(i)
  grid<-df_grid2[i,]
  
  matching_rows <-modern1[(modern1$X)>(grid$X1)&((modern1$X)<(grid$X2))&
                            ((modern1$Y)>(grid$Y1))&((modern1$Y)<(grid$Y2)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$L<-grid$L
    matching_rows_all<-rbind(matching_rows_all,matching_rows)
    
  } else {
    print("fail")
  }
  
}

# claculate the richness and mean range size for every grid

species_range<-matching_rows_all%>%group_by(species)%>%
  summarise(range=n_distinct(L))

matching_rows_all1<-merge(matching_rows_all,species_range,by="species")%>%
  dplyr::select(species,L,range)%>%distinct()


grid_richness_range<-matching_rows_all1%>%group_by(L)%>%
  summarise(richness=n_distinct(species),range=mean(range))
# delete the richness less than 20 because the the plant survey is not so 
# complete, also some grid cover the islands and ocean region

grid_richness_range<-grid_richness_range%>%filter(richness>20)

p17 <- ggplot(data = grid_richness_range, aes(x = richness, y = range)) + 
  geom_point(color="#313695",alpha=0.8) + geom_smooth(method = lm,color="#2166ac") +
  labs(y="modern plant mean range size of the 500km*500km grid cell",
       x="modern plant richness of the 500km*500km grid cell") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))+
  stat_cor(method = 'spearman', aes(x = richness, y = range))

#Supplementary Figure5g
p17

A<-lm(range~richness,data = grid_richness_range)
summary(A)

# compare the same taxa, get range size from larger northeast Siberia and Alaska
# region, with range size from the seven lake region

# set a Lambert Azimuthal Equal Area(LAEA) projection, the median value is base 
# on lake disreibution,and divide the region into 500km*500km grid cells

lakes  <- lake_coordinate %>% filter(!duplicated(Longitude)) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_shift_longitude()

proj_lake <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                     median((lakes %>% st_coordinates())[,1]), 
                     median((lakes %>% st_coordinates())[,2]))


# transfer the modern plant distribution projection
modern1  <- mydata1 %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(proj_lake)

modern1<-data.frame(st_coordinates(modern1))

modern1$species<-mydata1$species

# transfer the lake coordination in same projection
lake1  <- lakes %>% st_transform(proj_lake) %>% st_as_sfc(crs =proj_lake)
lake1<-data.frame(st_coordinates(lake1))
lake1$lake<-lakes$lake

lake1$X1<-lake1$X-250*1000
lake1$X2<-lake1$X+250*1000
lake1$Y1<-lake1$Y-250*1000
lake1$Y2<-lake1$Y+250*1000


matching_modern_plant_lake<-NULL
for (i in 1:nrow(lake1)) {
  print(i)
  lake2<-lake1[i,]
  
  matching_rows <-modern1[((modern1$X)>(lake2$X1)) & ((modern1$X)<(lake2$X2))& ((modern1$Y)<(lake2$Y2))& ((modern1$Y)>(lake2$Y1)), ]
  
  if (nrow(matching_rows) > 0) {
    matching_rows$lake<-lake2$lake
    matching_modern_plant_lake<-rbind(matching_modern_plant_lake,matching_rows)
    
  } else {
    print("fail")
  }
  
}


lake_species_range<-matching_modern_plant_lake%>%
  group_by(species)%>%summarise(range_lake=n_distinct(lake))

compare<-merge(lake_species_range,species_range,by="species")


p18<-ggplot(compare, aes(x=range_lake, y=range,group=range_lake)) +
  geom_boxplot(color="#2166ac")+
  stat_summary(fun=mean, geom='point',color="#053061", shape=20, size=3)+
  labs(x = "modern plant taxa range size around seven lakes region",
       y = "modern plant taxa range size in northeast Siberia and Alaska region")+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title = element_text(face = "bold", size = 14), 
        axis.text = element_text(size = 12),
        axis.ticks.length = unit(0.2, "cm"))
#Supplementary Figure5h
p18

plot<-plot_grid(p11,p12,p13,p14,p15,p16,p17,p18,
                labels=c("a","b","c","d","e","f","g","h"),
                label_x = 0.1,label_y = 1,
                ncol = 2,
                align = "v")
plot

ggsave(plot,filename="modern_sibala_plant_compare range in Siberia100-500_km.pdf",
       width = 8,height=20)












