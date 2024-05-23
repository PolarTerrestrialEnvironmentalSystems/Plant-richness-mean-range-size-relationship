################################################################################

# Author: Ying Liu,Simeon Lisovski
# Email:  ying.liu@awi.de,simeon.lisovski@awi.de
# Date:   May 21 2024

################################################################################

# Description: Here we calcute the cushion and tree plant abundance

################################################################################ 
# load the data and packages
library(dplyr)
resampleReads_final<-read.csv("resampleReads_final.csv")[,-c(1)]

################################################################################
# select the taxa belong to cushion plant
cushion_count<-resampleReads_final%>%
  filter(scientific_name=="Eritrichium"|scientific_name=="Eremogone capillaris"|
           scientific_name=="Sagina"|scientific_name=="Silene"|
           scientific_name=="Stellaria"|scientific_name=="Potentilla biflora"|
           scientific_name=="Saxifraga"|scientific_name=="Saxifraga oppositifolia"|
           scientific_name=="Diapensia lapponica"|scientific_name=="Draba"|
           scientific_name=="Armeria"|scientific_name=="Phlox hoodii")%>%
  group_by(round,time_slice)%>%summarize(count=sum(N))


total_count<-resampleReads_final%>%group_by(round,time_slice)%>%
  summarise(count=sum(N))

cushion_count$abundance<-(100*cushion_count$count)/(total_count$count)

cushion_count<-cushion_count%>%group_by(round,time_slice)%>%
  summarise(abundance=sum(abundance))

colnames(cushion_count)<-c("round","time_slice","cushion")

# select the taxa belong to cushion plant
tree_count<-resampleReads_final%>%
  filter(scientific_name=="Alnus"|scientific_name=="Alnus alnobetula"|
           scientific_name=="Betula"|scientific_name=="Cornus")%>%
  group_by(round,time_slice)%>%summarize(count=sum(N))


tree_count<-merge(tree_count,total_count,by=c("round","time_slice"))
tree_count$abundance<-(100*tree_count$count.x)/(tree_count$count.y)

tree_count<-tree_count%>%group_by(round,time_slice)%>%
  summarise(abundance=sum(abundance))%>%
  dplyr::select(round = round, time_slice = time_slice, tree = abundance)

#export data

write.csv(cushion_count,"cushion_count.csv")
write.csv(tree_count,"tree_count.csv")

