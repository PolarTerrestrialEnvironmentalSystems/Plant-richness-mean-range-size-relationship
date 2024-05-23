################################################################################

# Author: Simeon Lisovski,Ying Liu
# Email:  simenon.lisovski@awi.de,ying.liu@awi.de
# Date:   May 21 2024

################################################################################

# Description: Here we calcuted the total plant richness among all the lakes,
# mean range size, and biotic enviromental heterogeneity(beta-diversity)

################################################################################ 
# load the data and the packages

library("dplyr")
library("tidyr")
library("sf")
library("ggplot2")
library("betapart")
coordinate<-read.csv("lake_coordinate.csv")
resampleReads_final<-read.csv("resampleReads_final.csv")[,-c(1)]

################################################################################
# customised a equal area projection (Lambert Azimuthal Equal Area with the 
# center of the Study region),and divided the region into 200 km × 200 km grid

resampleReads_final<-merge(resampleReads_final,coordinate,by="lake")

lakes  <- resampleReads_final %>% filter(!duplicated(Longitude)) %>% 
       st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
       st_shift_longitude()
proj <- sprintf("+proj=laea +lon_0=%f +lat_0=%f", 
                median((lakes %>% st_coordinates())[,1]), 
                median((lakes %>% st_coordinates())[,2]))
grid   <- lakes %>% st_transform(proj) %>% st_bbox() %>% 
                    st_as_sfc(crs = proj) %>% st_buffer(150*1000) %>%
                    st_make_grid(200*1000)

map = suppressWarnings(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    dplyr::select("continent") %>% st_transform(proj)) %>% 
    st_intersection(grid %>% st_union())

ggplot() +
  geom_sf(data = map, fill = "grey90") +
  geom_sf(data = grid, fill = NA) +
  geom_sf(data = lakes, mapping = aes(color = lake)) +
  theme_void()

# calculate the total richness of all the lakes for per timeslice by counting 
# the number of ASVs(amplicon sequence variant) types, for Area of Occupancy (AOO) 
# method,calculating the range size for per ASV by summing lake numbers of the 
# ASV occur; for the Extent of Occurrence (EOO) method,calculating the as the sum of
# the grid cell areas overlapping with the convex hull spanning the lakes in 
# which specific ASVs occurred

RichRange <- resampleReads_final %>% as_tibble() %>%
  group_by(time_slice, round) %>% 
  summarise(richness = length(unique(NUC_SEQ))) %>% 
  ungroup() %>%tibble::rownames_to_column(var = "row_id") %>% 
  group_split(row_id) %>%
  lapply(., function(x) {
    
    dt_sf <- resampleReads_final %>% 
      filter(time_slice==x$time_slice, round == x$round) %>% 
      group_split(NUC_SEQ) %>%
      lapply(., function(s) {
        sfTab <- s %>% filter(!duplicated(Longitude)) %>%
          dplyr::select(lake, Longitude, Latitude) %>% 
          st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
          st_transform(proj)
        
        poly <- grid[(st_intersects(grid, st_convex_hull(sfTab %>% st_union()), 
                                    sparse = FALSE))[,1]]
        
        tibble(round = s$round[1], time_slice = s$time_slice[1], richness = 1,
               range_AOO = nrow(sfTab), 
               range_EOO = as.numeric(sum(poly %>% st_area())/1e6))
      }) %>% Reduce("rbind",.)

    dt_sf %>% group_by(round, time_slice) %>% summarise(richness = sum(richness), 
                      range_AOO = mean(range_AOO), range_EOO = mean(range_EOO))
    
  }) %>% suppressMessages() %>% Reduce("rbind",.)

save(RichRange, file = "RichRange.rda")
load("RichRange.rda")
# compare the mean range size based on the AOO method and EOO method, 
# Supplementary Figure1

ggplot(RichRange, aes(x = range_AOO, y = range_EOO,)) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  xlab("AOO mean range-size")+
  ylab("EOO mean range-size")+
  theme_linedraw()+
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5,size=12),
        axis.text.x = element_text(size = 12),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black"))
ggsave(filename="Supplementary_Figure1.pdf",width = 6,height=6)

# calculate the spatial plant turnover to represent the biotic environmental heterogeneity
resampleList2<- resampleReads_final%>%
  group_by(round,time_slice,lake,NUC_SEQ)%>%
  summarise(sumcount=sum(N))%>%
  pivot_wider(names_from = NUC_SEQ, values_from = sumcount, values_fill = 0)

heterogeneity<-NULL

for (x in unique(resampleList2$round)) {
  beta_all<-NULL
  resampleList3<-resampleList2%>%filter(round==x)%>%ungroup()%>%dplyr::select(-round)
  
  for(i in unique(resampleList3$time_slice)) {
    print(i)
    
    timeslicei<-resampleList3%>%filter(time_slice==i)
    timeslicei<-as.data.frame(timeslicei)
    rownames(timeslicei)<-timeslicei[,c(2)]
    timeslicei<-timeslicei[,-c(1:2)]
    
    timeslicei<-timeslicei[which(colSums(timeslicei)>0)]
    
    timeslicei<- data.frame(ifelse (timeslicei>0, 1, 0))
    beta<-as.data.frame(beta.multi(timeslicei, index.family="jaccard"))
    
    beta$timeslice<-i
    beta$round<-x
    beta_all<-rbind(beta_all,beta)
  }
  heterogeneity<-rbind(heterogeneity,beta_all) 
}

heterogeneity <- heterogeneity %>%
  dplyr::select(heterogeneity = beta.JTU, time_slice = timeslice, round = round)
write.csv(heterogeneity,"heterogeneity.csv")











