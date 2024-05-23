################################################################################

# Author: Simeon Lisovski,Ying Liu
# Email:  simeon.lisovski@awi.de,ying.liu@awi.de
# Date:   May 25 2024

################################################################################

# Description: Here we did the binomial regression, response variables of 
# richness to range-size slope (R>0.2 for positive relationship and R<-0.2 for 
# negative relationship), cushion plant abundance, and tree plant abundance as 
# the explanatory variables

################################################################################ 
#load the data and R packages
library("dplyr")
library("arm")
library("tidyr")
library("ggplot2")
datWindow<-read.csv("datWindow.csv")

################################################################################
slope <- lapply(datWindow %>% group_split(Windows), function(x) {
  mod   <- lm(range_AOO~richness, data = x)
  nsim <- nrow(x)
  bsim <- sim(mod, n.sim = nsim) 
  x$slope <- unlist(coef(bsim)[,2])####extract the sigma
  # Get summary of the linear model
  mod_summary <- summary(mod)
  
  # Extract R value from summary
  x$R_value <- mod_summary$r.squared
  x
}) %>% Reduce("rbind",.)


slope <- lapply(datWindow %>% group_split(round,Windows), function(x) {
  spearman_r <- cor(x$range_AOO, x$richness, method = "spearman")
  result <- data.frame(round = unique(x$round),
                      window = unique(x$Windows), 
                      spearman_r = spearman_r)
  result$cushion=median(x$cushion)
  result$tree=median(x$tree)
  
  return(result)
  
}) %>% Reduce("rbind",.)

# set the threshold
slope1<-slope%>%filter(spearman_r>0.2|spearman_r<(-0.2))
slope1$group<-ifelse(slope1$spearman_r>0,1,0)

slope_long <- slope1 %>% pivot_longer(cols = c("cushion", "tree"))

ggplot(slope_long, aes(x = value, y = group)) +
  geom_jitter(aes(color = name), width = 0, height = 0.01, alpha = 0.1, size = 1) +
  stat_smooth(aes(color = name), method = "glm", 
              method.args = list(family = binomial), se = TRUE) +
  facet_wrap(~name, scales = "free_x", ncol = 2)+
  scale_color_manual(values = c("#4754b2", "#9244a0"))+
  theme_classic()

model <- glm(group ~ cushion, data =slope1, family = binomial)
summary(model)

model <- glm(group ~ tree, data =slope1, family = binomial)
summary(model)

ggsave(filename = "Figure5_binomial_regression.pdf",
       dpi=1000,width = 6,height=4)

