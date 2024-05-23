################################################################################

# Author: Ying Liu,Simeon Lisovski
# Email:  ying.liu@awi.de, simeon.lisovski@awi.de
# Date:   May 21 2023

################################################################################

# Description: Here we merged all the data for Figure3, and made 5 timeslices
# running time window for the data,and plot the richness to range-size 
# relationship, richness, range-size, heterogeneity, cushion plant
# abundance, tree plant abundance,temperature change along time;also, we test 
# if the different timeslice window influence our result;besides, when load
# RichRange_common.rda, which is the resampling and calculation of only the taxa
# (represented by ASVs) exist in both the glacial period and the Holocene period
# was used for the relationship calculation, aim to remove the effects of extinction.

################################################################################ 
# load data and packages
library("tidyr")
library("dplyr")
library("ggplot2")
library("viridis")
library("arm")
library("cowplot")

load(file="RichRange.rda")
load(file="RichRange_common.rda")

heterogeneity<-read.csv("heterogeneity.csv")
cushion_abundance<-read.csv("cushion_count.csv")
tree_abundance<-read.csv("tree_count.csv")
Tem1<-read.csv("TANN_WAPLS_timeslice.csv")[,-c(1)]
Tem2<-read.csv("NGRIP.csv")
################################################################################
# merge the data togehter
data_resm<-merge(merge(merge(heterogeneity, cushion_abundance, 
                             by = c("round", "time_slice")), 
                       tree_abundance, by = c("round", "time_slice")),
                 RichRange,by=c("round", "time_slice"))

timeslice_age<-data.frame(time_slice=seq(from=1,to=30,by=1),
                          age=seq(from=30000,to=1000,by=-1000))
data<-merge(data_resm,timeslice_age,by="time_slice")

#interpolate the NGRIP value
Tem2<-Tem2%>% 
  group_by(age) %>%
  slice(2) %>%
  ungroup()

interp_Tem2<-as.data.frame(approx(Tem2$age, Tem2$isotope, 
                                        xout = as.numeric(timeslice_age$age)))
colnames(interp_Tem2)<-c("age","isotope")  

data<-merge(data,Tem1,by=c("round","age"))
data<-merge(data,interp_Tem2,by=c("age"))%>%
  dplyr::select(round,age,time_slice,heterogeneity,cushion, tree,richness,
                range_AOO,range_EOO,TANN_WAPLS,TANN_Lower_CI,TANN_Upper_CI,isotope)

{#when load the RichRange_common.rda, select this one
 data<-merge(data,interp_Tem2,by=c("age"))%>%
    dplyr::select(round,age,time_slice,heterogeneity,cushion, tree,richness,
                  range_AOO,TANN_WAPLS,TANN_Lower_CI,TANN_Upper_CI,isotope)
}


# 5 tiemslcie running windows

datWindow <- lapply(seq(30000, 5000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-4000):y))%>% 
    mutate(Window = mean(c((y-4000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+2000)/1000}-{(Window-2000)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

write.csv(datWindow,"datWindow.csv")

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI),
                                isotope=mean(isotope))

# plot the richness and mean range size relationship
p<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 2)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(26))) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.margin = margin(10, 10, 10, 10, "pt"),
    legend.key.size = unit(0.5, "cm")
  )+
  guides(color = guide_legend(ncol = 4,nrow = ))
# Figure4/Supplementary Figure7 when load the RichRange_common.rda
p 

ggsave(p,filename="Figure4_Richness_to_range-size_relationship.pdf",width = 10,height=6)
ggsave(p,filename="Supplementary Figure7.pdf",width = 10,height=6)

# plot the relationship, richness, range , heterogeneity, cushion abundance, 
# tree abundance,temperature together

slope <- lapply(datWindow %>% group_split(Windows), function(x) {
  mod   <- lm(range_AOO~richness, data = x)
  nsim <- nrow(x)
  bsim <- sim(mod, n.sim = nsim) 
  x$slope <- unlist(coef(bsim)[,2])####extract the sigma
  x
}) %>% Reduce("rbind",.)

# Supplementary Table 2
summary <- lapply(datWindow %>% group_split(Windows), function(x) {
  mod   <- lm(range_AOO~richness, data = x)
  summary_mod <- summary(mod)
  adjusted_r_squared <- summary_mod$r.squared
  p_value <- p.adjust(summary_mod$coefficients[2, 4], method = "bonferroni")
  slope <- coef(mod)[2]
  data.frame(window=unique(x$Windows),R2=adjusted_r_squared,p=p_value,slope=slope)
}) %>% Reduce("rbind",.)

slope_AOO<-slope[,c(15,16)]

slope_EOO <- lapply(datWindow %>% group_split(Windows), function(x) {
  mod   <- lm(range_EOO~richness, data = x)
  nsim <- nrow(x)
  bsim <- sim(mod, n.sim = nsim) 
  x$slope <- unlist(coef(bsim)[,2])####extract the sigma
  x
}) %>% Reduce("rbind",.)

# Supplementary Table 3
summary_EOO <- lapply(datWindow %>% group_split(Windows), function(x) {
  mod   <- lm(range_EOO~richness, data = x)
  summary_mod <- summary(mod)
  adjusted_r_squared <- summary_mod$r.squared
  p_value <- p.adjust(summary_mod$coefficients[2, 4], method = "bonferroni")
  slope <- coef(mod)[2]
  data.frame(window=unique(x$Windows),R2=adjusted_r_squared,p=p_value,slope=slope)
}) %>% Reduce("rbind",.)

slope_EOO<-slope_EOO[,c(15,16)]

AOO_EOO<-cbind(slope_AOO,slope_EOO)
colnames(AOO_EOO)<-c("Windows","AOO","Windows1","EOO")
AOO_EOO<-AOO_EOO[,c(1,2,4)]

AOO_EOO1<-AOO_EOO%>%pivot_longer(cols =c("AOO","EOO"),
                                names_to = "slope",
                                values_to = "value")
# Supplementary Figure 2
ggplot(AOO_EOO1) +
  geom_boxplot(aes(x = Windows, y = value, group = Windows),
               outlier.shape = NA,
               color = "#3651a1", fill = "#3651a1", alpha = 1, width = 0)+
    stat_summary(fun = median, geom = "point", 
               aes(x = Windows, y = value, group = Windows),
               color = "#3651a1", size = 1, shape = 16) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_light()+
  theme(
    axis.text.y.left = element_text(),
    axis.text.y.right = element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5, color = "black"), 
    axis.ticks.y = element_line(size = 0.5, color = "black"),
    axis.title = element_text(color = "black", size = 10), 
    panel.border = element_rect(color = "black", linewidth = 0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  facet_wrap(~ slope, scales = "free")

ggsave(filename="Supplementary Figure2.pdf",width = 10,height=6)

slope$age<-as.numeric(as.character(slope$age))
slope1<-slope[order(-slope$age),]
slope1$Windows<-factor(slope1$Windows,levels=unique(slope1$Windows))

slop2<-ggplot(slope1) +
  geom_boxplot(aes(x = Windows, y = slope,color=Windows,fill=Windows), 
               width=0.5,
               alpha = 1,outlier.color = NA)+
  stat_summary(
    data = slope1,
    aes(x = as.factor(Windows), y = slope),
    fun = "mean",
    geom = "point",
    size = 0.5,
  )+
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey",
             size = 0.65) +
  theme_light()+
  scale_color_manual(values = c("#FDE725", "#E5E419", "#C9E120", "#AEDC30",
                                "#93D741", "#7AD151", "#62CB5F", "#4DC36B",
                                "#3ABA76", "#2CB17E", "#22A884", "#1F9F88",
                                "#1F958B", "#228B8D", "#26828E", "#2A788E",
                                "#2E6F8E", "#32658E", "#375A8C", "#3C4F8A",
                                "#414487", "#453882", "#472B7A", "#481E6F",
                                "#471063", "#440154"))+
  scale_fill_manual(values = c("#FDE725", "#E5E419", "#C9E120", "#AEDC30",
                               "#93D741", "#7AD151", "#62CB5F", "#4DC36B",
                               "#3ABA76", "#2CB17E", "#22A884", "#1F9F88",
                               "#1F958B", "#228B8D", "#26828E", "#2A788E",
                               "#2E6F8E", "#32658E", "#375A8C", "#3C4F8A",
                               "#414487", "#453882", "#472B7A", "#481E6F",
                               "#471063", "#440154"))+
  theme(legend.position = "none",
        axis.text.y.left=element_text(),
        axis.text.y.right=element_text(),
        panel.background = element_rect(fill = "white"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.length = unit(0.2, "cm"), 
        axis.ticks.x = element_line(size = 0.5,color = "black"), 
        axis.ticks.y = element_line(size = 0.5,color = "black"),
        axis.title = element_text(color = "black",size = 10), 
        panel.border = element_rect(color = "black",linewidth =0.5),
        strip.text = element_text(),
        axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
  )

slop2

richness<-datWindow%>%group_by(round,Windows)%>%summarise(richness=mean(richness))

richness<-ggplot(richness) +
  geom_boxplot(aes(x = Windows, y = richness, fill = richness), 
               color = "#3651a1", fill = "#3651a1", alpha = 0.9, width = 0.5)+
  theme_light()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
richness

range<-datWindow%>%group_by(round,Windows)%>%summarise(range=mean(range_AOO))

range<-ggplot(range) +
  geom_boxplot(aes(x = Windows, y = range, fill = range), 
               color = "#911934", fill = "#911934", alpha = 0.9, width = 0.5)+
  theme_light()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
range

heterogeneity<-datWindow%>%group_by(round,Windows)%>%summarise(heterogeneity=mean(heterogeneity))
heterogeneity<-ggplot(heterogeneity) +
  geom_boxplot(aes(x = Windows, y = heterogeneity, fill = heterogeneity), 
               color = "#5c2366", fill = "#5c2366", alpha = 0.9, width = 0.5)+
  theme_light()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
heterogeneity

data_summary<-datWindow%>%dplyr::select(round, Windows,cushion,tree)%>%
  group_by(round,Windows)%>%
  summarise(cushion=mean(cushion),tree=mean(tree))%>%group_by(Windows)%>%
  summarise(cushion_sd=sd(cushion),tree_sd=sd(tree),
            cushion=mean(cushion),tree=mean(tree))

cushion<-ggplot(data_summary) +
  geom_bar(aes(x=Windows,y=cushion),fill="#8388ba",stat = "identity",position = "dodge") + 
  scale_y_continuous(name = "cushion plant abundance(%)")+
  theme_bw()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
cushion

tree<-ggplot(data_summary) +
  geom_bar(aes(x=Windows,y=tree),fill="#9e7ca3",stat = "identity",position = "dodge") + 
  scale_y_continuous(name = "tree plant abundance(%)")+
  theme_bw()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
tree

temperature<-ggplot(data=datWindow_mean,aes(x = Windows, y = TANN_WAPLS,group=1)) +
  geom_line()+
  geom_ribbon(aes(ymin = TANN_Lower_CI,
                  ymax = TANN_Upper_CI),
              alpha = 0.2) +
  theme_light()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )+
  labs(y = "temperature")
temperature

isotope<-ggplot(data=datWindow_mean,aes(x = Windows, y = isotope,group=1)) +
  geom_line(size=1)+
  theme_light()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )+
  labs(y = "isotope")
isotope

plot<-plot_grid(richness,range,slop2,heterogeneity,
                cushion,tree,
                temperature,isotope,
                labels=c("a","b","c","d","e","f","g","h"),
                label_x = 0.1,label_y = 1,
                ncol = 1,
                align = "v")
plot

ggsave(plot,filename="Figure3_merged_richness_mean-range_heterogeneity_cushion_tree_temperature_isotope.pdf",width = 6,height=15)

range1<-datWindow%>%group_by(round,Windows)%>%summarise(range=mean(range_EOO))

range1<-ggplot(range1) +
  geom_boxplot(aes(x = Windows, y = range, fill = range), 
               color = "#911934", fill = "#911934", alpha = 0.9, width = 0.5)+
  theme_light()+
  theme(
    axis.text.y.left=element_text(),
    axis.text.y.right=element_text(),
    panel.background = element_rect(fill = "white"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.2, "cm"), 
    axis.ticks.x = element_line(size = 0.5,color = "black"), 
    axis.ticks.y = element_line(size = 0.5,color = "black"),
    axis.title = element_text(color = "black",size = 10), 
    panel.border = element_rect(color = "black",linewidth =0.5),
    strip.text = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
range1
# Supplementary Figure3
# compare the mean range size based on the AOO method and EOO method over time

plot<-plot_grid(richness,range,range1,
                labels=c("a","b","c"),
                label_x = 0.1,label_y = 1,
                ncol = 1,
                align = "v")
plot
ggsave(plot,filename="Supplementary_Figure3.pdf",width = 6,height=8)

# for the richness and mean range size relationship, we test different 
# time window breadth, Supplementary Figure 4

#10 timeslice for 1 time window

heterogeneity<-read.csv("heterogeneity.csv")
cushion_abundance<-read.csv("cushion_count.csv")
tree_abundance<-read.csv("tree_count.csv")
load(file="RichRange.rda")



# merge the data togehter
data_resm<-merge(merge(merge(heterogeneity, cushion_abundance, 
                             by = c("round", "time_slice")), 
                       tree_abundance, by = c("round", "time_slice")),
                 RichRange,by=c("round", "time_slice"))

# 10 timesliced 1 time window

timeslice_age<-data.frame(time_slice=seq(from=1,to=30,by=1),
                          age=seq(from=30000,to=1000,by=-1000))
data<-merge(data_resm,timeslice_age,by="time_slice")

data<-merge(data,Tem1,by=c("round","age"))%>%
  dplyr::select(round,age,time_slice,heterogeneity,cushion, tree,richness,
                range_AOO,range_EOO,TANN_WAPLS,TANN_Lower_CI,TANN_Upper_CI)

datWindow <- lapply(seq(30000, 10000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-9000):y))%>% 
    mutate(Window = mean(c((y-9000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+4500)/1000}-{(Window-4500)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))

# plot the richness and mean range size relationship

p10<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(26))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))
 
p10

###9 timeslice 1 time window

datWindow <- lapply(seq(30000, 9000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-8000):y))%>% 
    mutate(Window = mean(c((y-8000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+4000)/1000}-{(Window-4000)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))


# plot the richness and mean range size relationship

p9<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(26))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))
p9

###8 timeslice 1 time window

datWindow <- lapply(seq(30000, 8000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-7000):y))%>% 
    mutate(Window = mean(c((y-7000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+3500)/1000}-{(Window-3500)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))


# plot the richness and mean range size relationship

p8<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(26))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))

p8

###7 timeslice 1 time window

datWindow <- lapply(seq(30000, 7000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-6000):y))%>% 
    mutate(Window = mean(c((y-6000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+3000)/1000}-{(Window-3000)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))


# plot the richness and mean range size relationship

p7<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(26))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))
  
p7

###6 timeslice 1 time window

datWindow <- lapply(seq(30000, 6000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-5000):y))%>% 
    mutate(Window = mean(c((y-5000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+2500)/1000}-{(Window-2500)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))


# plot the richness and mean range size relationship

p6<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(26))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))
 
p6

###4 timeslice 1 time window

datWindow <- lapply(seq(30000, 4000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-3000):y))%>% 
    mutate(Window = mean(c((y-3000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+1500)/1000}-{(Window-1500)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))


# plot the richness and mean range size relationship

p4<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(27))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))
p4


lm_results <- datWindow %>%
  group_by(Windows) %>%
  summarise(p_value=coef(summary(lm(range_AOO ~ richness)))["richness", "Pr(>|t|)"])
  
  
###3 timeslice 1 time window

datWindow <- lapply(seq(30000, 3000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-2000):y))%>% 
    mutate(Window = mean(c((y-2000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+1000)/1000}-{(Window-1000)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))


# plot the richness and mean range size relationship

p3<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(28))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))
p3

lm_results <- datWindow %>%
  group_by(Windows) %>%
  summarise(p_value=coef(summary(lm(range_AOO ~ richness)))["richness", "Pr(>|t|)"])

###2 timeslice 1 time window

datWindow <- lapply(seq(30000, 2000, by = -1000), function(y) {
  data %>% filter(age %in% c((y-1000):y))%>% 
    mutate(Window = mean(c((y-1000):y)))})%>% 
  Reduce("rbind", .) %>% 
  mutate(Windows = glue::glue("{(Window+500)/1000}-{(Window-500)/1000}ka"))

datWindow <- datWindow %>%
  mutate(age = factor(datWindow$age, levels = unique(datWindow$age)))%>%
  mutate(Windows = factor(datWindow$Windows, levels = unique(datWindow$Windows)))

datWindow_mean<-datWindow[!duplicated(datWindow),]%>%
  group_by(Windows)%>%summarise(TANN_WAPLS=mean(TANN_WAPLS),
                                TANN_Lower_CI=mean(TANN_Lower_CI),                          
                                TANN_Upper_CI=mean(TANN_Upper_CI))


# plot the richness and mean range size relationship

p2<-ggplot(datWindow, aes(x = richness, y = range_AOO, color = as.factor(Windows))) +
  geom_point(size = 1)+
  geom_smooth(method="lm",size = 1,alpha=0.5)+
  scale_color_manual(values = rev(viridis(29))) +
  labs(x="plant taxa richness",y="mean range size")+
  theme(legend.key.size = unit(0.2, "cm"),
        legend.text = element_text(size = 7),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA))+
  guides(color = guide_legend(override.aes = list(size = 0.3),ncol=3))

p2

plot<-plot_grid(p10,p9,p8,p7,
                p6,p4,p3,p2,
                labels=c("10ka time window","9ka time window",
                         "8ka time window","7ka time window",
                         "6ka time window","4ka time window",
                         "3ka time window","2ka time window"),
                label_x = 0.1,label_y = 1,
                ncol = 2,
                align = "v")
plot

ggsave(plot,filename="Supplementary Figure4_merge_different_time_window_relationship.pdf",width = 12,height=15)






