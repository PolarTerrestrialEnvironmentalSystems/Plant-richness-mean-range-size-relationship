################################################################################

# Author: Thomas Boehmer,Ying Liu
# Email: thomas.boehmer@awi.de, ying.liu@awi.de
# Date: May 21 2024

################################################################################

# Description: Here We reconstructed the mean annual temperature over the past 
# 30,000 years for the study region based on pollen assemblage records from ten lakes

################################################################################ 
#load the data and packages

library("rioja")
library("crayon")
library("palaeoSig")
library("tidyverse")

tann_surrogates <- read.csv("2023-02-01_Tann_Simulated_modern_Matern.csv",  
                            header=TRUE, stringsAsFactors=FALSE)
modern  <- read.csv("2022-11-29_northern_hemispheric_modern_pollenset.csv", 
                    header=TRUE, sep="\t", stringsAsFactors=FALSE)

fossil <- read.csv("2022-01-10_pollen_percentage_hill70_asia.csv", 
                   header=TRUE, sep="\t", stringsAsFactors=FALSE)

# radius for the circle in which modern pollen data is searched; [in km]
distance <- 1000         


wapls_recon_results<- NULL

Dataset_ID<-c("31","788","41138","100273","100327")

runIDs <-Dataset_ID

for(ID in runIDs){
  
# create a subset of one fossil pollen site with the actual id in the loop:
  fossil.subset <- fossil[fossil$Dataset_ID==ID,] 
  
# create a table with all taxa from the subset,use only taxa with countings 
# (exclude Zeros) 
  subset.spp  <- fossil.subset[,-c(1:15)]
  subset.spp1 <- subset.spp[, which(colSums(subset.spp) != 0)] 
  
# extract the coordinates of the actual subset:
  subset.coord <- cbind(fossil.subset$Longitude[1], fossil.subset$Latitude[1])
  
# calculate all geographic distances between the actual fossil site and the 
# modern pollen dataset: 
  geo.dist <- fields::rdist.earth(subset.coord, 
                              modern[,c("Longitude","Latitude")], miles = FALSE) 
  
# write the results in the rdist-column of the modern pollen table:
  r.dist <- geo.dist[1,] 
  
# check which modern pollen datasets are within the geographic range of 1000 km 
# of the fossil pollen subset:
  within.distance <- modern[r.dist < distance, ]
  within.distance <- na.omit(within.distance)
  
# calculating distances among the sampled points in the trainingsdataset:
  geodist.within.distance <- fields::rdist.earth(select(within.distance, 
                                      Longitude, Latitude), miles = FALSE)
  
# reconstruction
# if there are less than 30 samples in range, the reconstruction is skipped    
  if(nrow(within.distance) < 30) {  
    
  }else{
# write all taxa of the modern pollen database in a new table,use only taxa with
# countings (exclude Zeros)  
  spp1 <- within.distance[ ,-c(1:23)]         
  spp2 <- spp1[, which(colSums(spp1) != 0)]   
# calculating the reconstruction for all variables with WA-PLS: 
  climvar <- c("TANN") 
  climate.recon_all <- NULL
  stats.recon_all   <- data.frame(Samples_1000 = nrow(within.distance))
  wapls_components  <- data.frame(TANN=NA)  
    
  for (j in 1:length(climvar)){
      
  stats.recon <- data.frame(min=NA, max=NA, Comp=NA, r2=NA, RMSE=NA, RMSE_training=NA)
  env.var <- within.distance[ ,which(names(within.distance) == climvar[j])] 
# take out the climmate variable that we need to use
      
  stats.recon[1,"min"] <- min(env.var)
  stats.recon[1,"max"] <- max(env.var)
  wapls.var <- suppressWarnings(crossval(WAPLS(sqrt(spp2), env.var)))
# esteblish the relationship of the pollen assemblage and environmental factor,
# and cross-validation
# randomization t-test, test the relationship   
  tt.var <- rand.t.test(wapls.var)
  m.var <- 1
  if(tt.var[2,6] <= -5 & tt.var[2,7] <= 0.05) {m.var <- 2}
      
  wapls_components[1,j] <- m.var
  stats.recon[1,"Comp"] <- m.var
  stats.recon[1,"r2"]   <- tt.var[m.var,2]
  stats.recon[1,"RMSE"] <- tt.var[m.var,1]
      
  pred.var <- suppressWarnings(predict(wapls.var, 
                                sqrt(subset.spp1), npls=m.var, sse=TRUE))
      
# according to the wapls regrassion, predict the environment of the fossil data
      
  if(any(is.na(pred.var$fit[,m.var]))){
        
  del <- which(is.na(pred.var$fit[,m.var]))
        
# if the predict component corresponding to the reconstruction component is na 
# value,then this value will be delete, and do the reconstruction again
        
  subset.spp2 <- subset.spp1[-del, ]
  pred.var <- suppressWarnings(predict(wapls.var, 
                              sqrt(subset.spp2), npls=m.var, sse=TRUE))
        
  climate.recon <- cbind(pred.var$fit[,m.var], 
                         pred.var$SEP.boot[,m.var], 
                        pred.var$v1.boot[,m.var])
        
  colnames(climate.recon) <- c(climvar[j],
                             paste0(climvar[j],"_Error"),
                            paste0(climvar[j],"_se_bootstrap_estimates"))
        
  stats.recon[1,"RMSE_training"] <- pred.var$v2.boot[m.var]
        
  }else{
        
  climate.recon <- cbind(pred.var$fit[,m.var], 
                         pred.var$SEP.boot[,m.var], 
                         pred.var$v1.boot[,m.var])
        
  colnames(climate.recon) <- c(climvar[j],
                             paste0(climvar[j],"_Error"),
                             paste0(climvar[j],"_se_bootstrap_estimates"))
        
  stats.recon[1,"RMSE_training"] <- pred.var$v2.boot[m.var]
        
      }
      
      
  colnames(stats.recon) <- paste0(climvar[j],"_",names(stats.recon))
      
  stats.recon_all   <- as.data.frame(c(stats.recon_all, unlist(stats.recon)))
  climate.recon_all <- cbind(climate.recon_all, climate.recon)
      
  climate.recon_all <- data.frame(climate.recon_all)
      
  if(dim(fossil.subset)[1] != dim(climate.recon_all)[1]
      ){site.info <- fossil.subset[-del,c(1:15)] 
      }else{site.info <- fossil.subset[ ,c(1:15)]}
      
  wapls_recon.table<- cbind(site.info,climate.recon_all) 
      
  wapls_recon.table<- data.frame(wapls_recon.table)
      
  wapls_recon_results<- rbind(wapls_recon_results, wapls_recon.table)
      
    }  
  } 
}   



# america_west
fossil <- read.csv("2022-01-10_pollen_percentage_hill70_namerica_west.csv", 
                   header=TRUE, sep="\t", stringsAsFactors=FALSE)


# radius for the circle in which modern pollen data is searched; [in km]
distance <- 1000         

Dataset_ID <- c(1003,1390,1432,16238,17391)

runIDs <-Dataset_ID

for(ID in runIDs){
  
# create a subset of one fossil pollen site with the actual id in the loop:
  fossil.subset <- fossil[fossil$Dataset_ID==ID,] 
  
# create a table with all taxa from the subset,use only taxa with countings 
# (exclude Zeros) 
  subset.spp  <- fossil.subset[,-c(1:15)]
  subset.spp1 <- subset.spp[, which(colSums(subset.spp) != 0)] 
  
# extract the coordinates of the actual subset:
  subset.coord <- cbind(fossil.subset$Longitude[1], fossil.subset$Latitude[1])
  
# calculate all geographic distances between the actual fossil site and the 
# modern pollen dataset: 
  geo.dist <- fields::rdist.earth(subset.coord, 
                                  modern[,c("Longitude","Latitude")], miles = FALSE) 
  
# write the results in the rdist-column of the modern pollen table:
  r.dist <- geo.dist[1,] 
  
# check which modern pollen datasets are within the geographic range of 1000 km 
# of the fossil pollen subset:
  within.distance <- modern[r.dist < distance, ]
  within.distance <- na.omit(within.distance)
  
# calculating distances among the sampled points in the trainingsdataset:
  geodist.within.distance <- fields::rdist.earth(select(within.distance, 
                                                        Longitude, Latitude), miles = FALSE)
  
# reconstruction
# if there are less than 30 samples in range, the reconstruction is skipped    
  if(nrow(within.distance) < 30) {  
    
  }else{
# write all taxa of the modern pollen database in a new table,use only taxa with
# countings (exclude Zeros)  
    spp1 <- within.distance[ ,-c(1:23)]         
    spp2 <- spp1[, which(colSums(spp1) != 0)]   
# calculating the reconstruction for all variables with WA-PLS: 
    climvar <- c("TANN") 
    climate.recon_all <- NULL
    stats.recon_all   <- data.frame(Samples_1000 = nrow(within.distance))
    wapls_components  <- data.frame(TANN=NA)  
    
    for (j in 1:length(climvar)){
      
      stats.recon <- data.frame(min=NA, max=NA, Comp=NA, r2=NA, RMSE=NA, RMSE_training=NA)
      env.var <- within.distance[ ,which(names(within.distance) == climvar[j])] 
# take out the climmate variable that we need to use
      
      stats.recon[1,"min"] <- min(env.var)
      stats.recon[1,"max"] <- max(env.var)
      wapls.var <- suppressWarnings(crossval(WAPLS(sqrt(spp2), env.var)))
# esteblish the relationship of the pollen assemblage and environmental factor,
# and cross-validation
# randomization t-test, test the relationship   
      tt.var <- rand.t.test(wapls.var)
      m.var <- 1
      if(tt.var[2,6] <= -5 & tt.var[2,7] <= 0.05) {m.var <- 2}
      
      wapls_components[1,j] <- m.var
      stats.recon[1,"Comp"] <- m.var
      stats.recon[1,"r2"]   <- tt.var[m.var,2]
      stats.recon[1,"RMSE"] <- tt.var[m.var,1]
      
      pred.var <- suppressWarnings(predict(wapls.var, 
                                           sqrt(subset.spp1), npls=m.var, sse=TRUE))
      
# according to the wapls regrassion, predict the environment of the fossil data
      
      if(any(is.na(pred.var$fit[,m.var]))){
        
        del <- which(is.na(pred.var$fit[,m.var]))
        
# if the predict component corresponding to the reconstruction component is na 
# value,then this value will be delete, and do the reconstruction again
        
subset.spp2 <- subset.spp1[-del, ]
        pred.var <- suppressWarnings(predict(wapls.var, 
                                   sqrt(subset.spp2), npls=m.var, sse=TRUE))
        
climate.recon <- cbind(pred.var$fit[,m.var], 
                      pred.var$SEP.boot[,m.var], 
                      pred.var$v1.boot[,m.var])
        
colnames(climate.recon) <- c(climvar[j],
                              paste0(climvar[j],"_Error"),
                              paste0(climvar[j],"_se_bootstrap_estimates"))
        
stats.recon[1,"RMSE_training"] <- pred.var$v2.boot[m.var]
        
      }else{
        
      climate.recon <- cbind(pred.var$fit[,m.var], 
                             pred.var$SEP.boot[,m.var], 
                             pred.var$v1.boot[,m.var])
        
      colnames(climate.recon) <- c(climvar[j],
                                   paste0(climvar[j],"_Error"),
                                   paste0(climvar[j],"_se_bootstrap_estimates"))
        
      stats.recon[1,"RMSE_training"] <- pred.var$v2.boot[m.var]
        
      }
      
      
colnames(stats.recon) <- paste0(climvar[j],"_",names(stats.recon))
stats.recon_all   <- as.data.frame(c(stats.recon_all, unlist(stats.recon)))
climate.recon_all <- cbind(climate.recon_all, climate.recon)
      
climate.recon_all <- data.frame(climate.recon_all)
      
if(dim(fossil.subset)[1] != dim(climate.recon_all)[1]
  ){site.info <- fossil.subset[-del,c(1:15)] 
  }else{site.info <- fossil.subset[ ,c(1:15)]}
      
  wapls_recon.table<- cbind(site.info,climate.recon_all) 
      
  wapls_recon.table<- data.frame(wapls_recon.table)
      
  wapls_recon_results<- rbind(wapls_recon_results, wapls_recon.table)
      
    }  
  } 
}   

TANN<-wapls_recon_results%>%dplyr::select(Site_Name,meanAgeBP,TANN)

interpolated_TANN_all<-data.frame(Age=seq(1000,30000,1000))
age_new <-seq(1000,30000,1000)

for (i in unique(TANN$Site_Name)) {
  print(i)
  TANN1<-TANN%>%filter(Site_Name==i) 
  
  interpolated_TANN <- approx(x = TANN1$meanAgeBP, 
                              y = TANN1$TANN, xout = age_new)%>%unlist()%>%
                              as.data.frame()
  
  interpolated_TANN <-as.data.frame(interpolated_TANN[c(31:60),])
  colnames(interpolated_TANN)<-i
  
  interpolated_TANN_all<-cbind(interpolated_TANN_all,interpolated_TANN)
  
}

TANN_WAPLS<-interpolated_TANN_all%>%pivot_longer(2:11)%>%
  rename(Age=Age,lake=name,value=value)%>%
  group_by(Age) %>%
  summarise(
    Mean = mean(value),
    SE = sd(value) / sqrt(n())
  )

summary_data1 <- TANN_WAPLS%>%
  mutate(
    Lower_CI = Mean - qt(0.975, df = n() - 1) * SE,  
    Upper_CI = Mean + qt(0.975, df = n() - 1) * SE  
  )%>%dplyr::select(-SE)

colnames(summary_data1)<-c("age","TANN_WAPLS","TANN_Lower_CI","TANN_Upper_CI")


Tem<-summary_data1
Tem1 <- data.frame()

for (i in 1:100) {
  Tem1 <- rbind(Tem1, Tem) 
}
Tem1 <-as.data.frame(Tem1)

Tem1$round<-rep(1:100,each=30)


write.csv(Tem1,"TANN_WAPLS_timeslice.csv")












