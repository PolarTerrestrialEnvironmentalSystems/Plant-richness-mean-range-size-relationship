################################################################################

# Author: Ying Liu, Simeon Lisovski
# Email:  ying.liu@awi.de, simeon.lisovski@awi.de
# Date:   May 21 2024

################################################################################

# Description: Here we did the network analysis using the terrestrial plant 
# family based on the pairwise correlation of plant family, only the r-value 
# more than 0.6 (p<0.05) was kept for plot

################################################################################ 
# load the data and packages

library("dplyr")
library("tidyr")
library("Hmisc")
library("igraph")
library("ggplot2")
library("ggraph")
resampleReads_final<-read.csv("resampleReads_final.csv")

################################################################################

data<-resampleReads_final%>%dplyr::select(best_family,scientific_name,N,
                                          time_slice,round)%>%
  group_by(round,time_slice,best_family,scientific_name,)%>%summarise(N=sum(N))

timeslice_age<-data.frame(time_slice=seq(from=1,to=30,by=1),
                          age=seq(from=30000,to=1000,by=-1000))
data<-merge(data,timeslice_age,by="time_slice")%>%dplyr::select(-time_slice)

# remove the aquatic plants,Pteridophyta 

family<-data%>%filter(!best_family=="Acoraceae"&
                        !best_family=="Alismataceae"&
                        !best_family=="Pontederiaceae"&
                        !best_family=="Hydrocharitaceae"&
                        !best_family=="Nymphaeaceae"&!
                        best_family=="Potamogetonaceae"&
                        !best_family=="Typhaceae"&
                        !best_family=="Isoetaceae"&
                        !best_family=="Menyanthaceae"&
                        !best_family=="Ceratophyllaceae"&
                        !best_family=="Haloragaceae"&
                        !best_family=="Lentibulariaceae"&
                        !best_family=="Juncaginaceae"&
                        !best_family=="Andreaeaceae"&
                        !best_family=="Bryaceae"&
                        !best_family=="Dicranaceae"&
                        !best_family=="Grimmiaceae"&
                        !best_family=="Polytrichaceae"&
                        !best_family=="Scapaniaceae"&
                        !best_family=="Ditrichaceae"&
                        !best_family=="Splachnaceae"&
                        !best_family=="Pottiaceae"&
                        !best_family=="Brachytheciaceae"&
                        !best_family=="Harpanthaceae"&
                        !best_family=="Amblystegiaceae"&
                        !best_family=="Brachytheciaceae"&
                        !best_family=="Diphyscium"&
                        !best_family=="Plagiotheciaceae"&
                        !best_family=="Meesiaceae"&
                        !best_family=="Hylocomiaceae"&
                        !best_family=="Rhabdoweisiaceae"&
                        !best_family=="Gymnomitriaceae"&
                        !best_family=="Equisetaceae"&
                        !best_family=="Athyriaceae"&
                        !best_family=="Dryopteridaceae"&
                        !best_family=="Woodsiaceae"&
                        !best_family=="Lycopodiaceae"&
                        !best_family=="Cystopteridaceae"&
                        !best_family=="Thelypteridaceae"&
                        !best_family=="Ophioglossaceae"&
                        !best_family=="Sphagnaceae"&
                        !best_family=="Aulacomniaceae")

# only the plant family more total counts more than 10 was kept 
data2<-family%>%filter(!is.na(best_family))%>%
  group_by(age,best_family)%>%
  summarise(count=sum(N))%>%filter(count>10)%>%
  group_by(age,best_family)%>%
  summarise(count=sum(count))%>%
  pivot_wider(names_from = best_family,values_from = count)%>%
  mutate_all(~ifelse(is.na(.), 0, .))

# only the family occur more than 5 timeslices ware kept
data3 <- data2%>%ungroup%>%dplyr::select(-age)
data31 <- data3
data31[data31>0] <- 1
data3 <- data3[,which(colSums(data31) >= 5)]

# calculate the correlation
corr <- rcorr(as.matrix(data3), type = 'spearman')
r <- corr$r
r[abs(r) < 0.6] <- 0

# p<0.05
p <- corr$P
p <- p.adjust(p, method = 'BH')  
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1

z <- r * p
diag(z) <- 0  
z[z<0]<-0
z <- z[!rowSums(z) == 0, !colSums(z) == 0]

g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

g <- simplify(g)
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
g
plot(g)

# extract the edge
edge <- data.frame(as_edgelist(g))   
edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation,
  type=ifelse(E(g)$correlation > 0,"positive","negative"))

graph <- graph_from_data_frame(edge)

# cluster the taxa
deg <- degree(g, mode="all")
cs<-cluster_optimal(g)
V(g)$community <-cs$membership
plot(cs,g)

ggraph(g, layout_with_graphopt(g, charge=0.03)) +
  geom_edge_link0(edge_linewidth = E(g)$weight*1.5,
                  edge_colour = ifelse(E(g)$correlation > 0.6, "#bababa","#252525")) +
  geom_node_point(aes(fill = degree(g, mode="all"), 
                      size = degree(g, mode="all")), 
                  color = adjustcolor(c("#257ab6","#b2182b"))[cs$membership])+
  geom_node_text(aes(label = name), family = "serif",
                 color= adjustcolor(c("black","black"))[cs$membership])+
  theme(legend.position = "bottom")

ggsave(filename = "Figure6_network.pdf",
       dpi=1000,width = 8,height=7)

cs$membership
cs$names
cushion_group<-data.frame(group=cs$membership,name=cs$names)%>%filter(group==2)
tree_group<-data.frame(group=cs$membership,name=cs$names)%>%filter(group==1)
# cushion group 11,tree group 11
cushion_edge<-edge_list%>%
  filter((edge_list$source%in%cushion_group$name)&(edge_list$target%in%cushion_group$name))%>%
  filter(correlation>0.6)
# cushion edge30

tree_edge<-edge_list%>%
  filter((edge_list$source%in%tree_group$name)&(edge_list$target%in%tree_group$name))%>%
  filter(correlation>0.6)
# tree edge17

# calculate the sum of the possible real link (r>0.6 or r<-0.6)
cushion_link<-r[,colnames(r)%in%cushion_group$name]
cushion_link<-cushion_link[rownames(cushion_link)%in%cushion_group$name,]

cushion_r<-as.data.frame(cushion_link[lower.tri(cushion_link)])
colnames(cushion_r)<-"r"

potential_negative<-cushion_r[cushion_r$r<0]# there is no negative link for the cushion group
# tree

tree_link<-r[,colnames(r)%in%tree_group$name]
tree_link<-tree_link[rownames(tree_link)%in%tree_group$name,]

tree_r<-as.data.frame(tree_link[lower.tri(tree_link)])
colnames(tree_r)<-"r"

potential_negative<-tree_r[tree_r$r<0]# there is no negative link for the tree group

# calculate both the positive and negative interaction, and built the network

# calculate the correlation
corr <- rcorr(as.matrix(data3), type = 'spearman')
r <- corr$r
r[abs(r) < 0.6] <- 0

# p<0.05
p <- corr$P
p <- p.adjust(p, method = 'BH')  

p[p>=0.05] <- (-1)
p[p<0.05 & p>=0] <- 1
p[p<0] <- 0

z <- r * p
diag(z) <- 0  
z <- z[!rowSums(z) == 0, !colSums(z) == 0]

g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g

g <- simplify(g)
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)
g
plot(g)


# extract the edge
edge <- data.frame(as_edgelist(g))   
edge_list <- data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$correlation,
  type=ifelse(E(g)$correlation > 0,"positive","negative"))

graph <- graph_from_data_frame(edge)

# cluster the taxa
deg <- degree(g, mode="all")
cs<-cluster_optimal(g)
V(g)$community <-cs$membership
plot(cs,g)

ggraph(g, layout_with_graphopt(g, charge=0.03)) +
  geom_edge_link0(edge_linewidth = E(g)$weight*1.5,
                  edge_colour = ifelse(E(g)$correlation > 0, "#d6604d","#74add1")) +
  geom_node_point(aes(fill = degree(g, mode="all"), 
                      size = degree(g, mode="all")), 
                  color = adjustcolor(c("#bababa","#bababa","#bababa"))[cs$membership])+
  geom_node_text(aes(label = name), family = "serif",
                 color= adjustcolor(c("black","black","black"))[cs$membership])+
  theme(legend.position = "bottom")+
  theme_void()

ggsave(filename = "network_supple_positive_negaitve_r_more_than0.6.pdf",
       dpi=1000,width = 8,height=5)

# calculate the sum of the possible real link (r>0.6 or r<-0.6)

cushion_link <- data.frame(matrix(NA, nrow = nrow(r), ncol = ncol(r)))

cushion_link[upper.tri(cushion_link, diag = FALSE)] <- r[upper.tri(r, diag = FALSE)]

colnames(cushion_link)<-colnames(r)
rownames(cushion_link)<-rownames(r)

cushion_link<-cushion_link[,colnames(cushion_link)%in%cushion_group$name]

cushion_link[is.na(cushion_link)]<-0

cushion_negative<-sum(cushion_link<0, na.rm = TRUE) #14
cushion_positive<-sum(cushion_link>0, na.rm = TRUE) #30

# tree

tree_link <- data.frame(matrix(NA, nrow = nrow(r), ncol = ncol(r)))

tree_link[upper.tri(tree_link, diag = FALSE)] <- r[upper.tri(r, diag = FALSE)]

colnames(tree_link)<-colnames(r)
rownames(tree_link)<-rownames(r)

tree_link<-tree_link[,colnames(tree_link)%in%tree_group$name]

tree_link[is.na(tree_link)]<-0

#tree_link<-tree_link[!rownames(r)%in%tree_group$name,]

tree_negative<-sum(tree_link<0, na.rm = TRUE) #27
tree_positive<-sum(tree_link>0, na.rm = TRUE) #17

