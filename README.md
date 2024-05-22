# Script to Plant interactions shift the relationship of richness to mean range-size relationship
This repository includes the script for "Plant interactions shift the relationship of richness to mean range size from positive to negative at the Glacial–Holocene transition". 

All the scripts for data analysis, as well as for the generation of figures and tables are shown here. The main files included are as follows:

0_resampling

1_calculation: Based on the resampling result, this script calculates plant richness, mean range-size, and heterogeneity. The output includes the results and Supplementary Figure 1.

2_cushion_tree: Cushion plant taxa and tree plant taxa abundance change over time based on the resampling result.

3_merge_plot: Based on the resampling, plant richness to mean range-size relationship was calculated with the 5,000-year running time window, the richness, mean range-size, heterogeneity, and the temperature result, and the isotope result of the NGRIP ice core were all shown in the 5,000-year running time window. The output in this step includes Figure 3, Figure 4, Supplementary Figure 2, Supplementary Figure 4, and Supplementary Figure 7.
              
4_binomial_regression: In this step, we did the binomial regression, response variables of richness to range-size relationship (R>0.2 for positive and R<-0.2 for negative), 
                       cushion plant abundance, and tree plant abundance as the explanatory variables.
5_network_analysis: we did the network analysis using the terrestrial plant family based on the pairwise correlation of plant family, only the r-value more than 0.6 (p<0.05) was kept for plot.
6_richness_range_size_relationship_in_space: Here we used the modern plant distribution data downloaded from the Global Biodiversity Information Facility database (GBIF, https://www.gbif.org/)
to calculate plant richness to range-size relationship in space, in addition, we used the sedaDNA modern timeslice data to calculate the relationship in space, besides, we compared the same modern 
plant taxa range calculated based on two regions: seven lake region and the larger northeast Siberia and Alaska region.
7_temperature_reconstruction: Here We reconstructed the mean annual temperature over the past 30,000 years for the study region based on pollen assemblage records from ten lakes
