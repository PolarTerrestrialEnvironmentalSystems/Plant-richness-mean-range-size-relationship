################################################################################

# Author: Simeon Lisovski,Ying Liu
# Email: simeon.lisovski@awi.de,ying.liu@awi.de
# Date: May 21 2024

################################################################################

# Description: Here we performed resampling for the sedaDNA dataset as the 
# number of matched ASVs (read counts) differs across lakes and time slices
# after resampling, we make sure same number of ASV (read counts) are selected
# for per lake at every timeslice

################################################################################ 
# load the data and the packages

library("dplyr")
library("data.table")
dataset=read.csv("dataset_for_resampling.csv")[-c(1)]

################################################################################
# only keep the data that toal reads more than 5000 per lake per timeslice

resampleLista<-tidyr::unite(dataset, "lake_time", lake, time_slice, remove = FALSE)

resampleListb<-resampleLista%>%group_by(time_slice,lake)%>%
  summarise(nreads = sum(counts), 
            nseq = n_distinct(NUC_SEQ), 
            nsamples = n_distinct(name_age))

resampleListc<-tidyr::unite(resampleListb, "lake_time", lake, 
                            time_slice, remove = FALSE)%>%
                            filter(nreads>=5000)

resampleLista<-resampleLista%>%
  filter(resampleLista$lake_time %in% resampleListc$lake_time)%>%as.data.table()

#do resampling,repeat 1000 times

resampleReads_final<-NULL
for (r in 1:100) {
  print(r)
  resampleReads_all <-NULL
  for(j in unique(resampleLista$time_slice)) {
    resampleReads<-resampleLista%>%filter(time_slice==j)
    resampleReads1 <- resampleReads[,num:=1:.N
    ][, cbind(.SD, dup = 1:counts), by = "num"][,counts := 1,
    ][,c("NUC_SEQ", "time_slice","lake","scientific_name","best_family")
    ][,.SD[sample(1:.N, 5000)],by = lake
    ][,.(.N), by = c("NUC_SEQ", "time_slice",
    "lake","scientific_name","best_family")]
    resampleReads_all <-rbind(resampleReads_all,resampleReads1,fill=TRUE)
    resampleReads_all$round=r
  }
  resampleReads_final<-rbind(resampleReads_all,resampleReads_final)
}

resampleReads_final<-type.convert(resampleReads_final,as.is=TRUE)

################################################################################
# export the data 
write.csv(resampleReads_final,"resampleReads_final.csv")





