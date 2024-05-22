# Script to Plant interactions shift the relationship of richness to mean range-size relationship
This repository includes the script for "Plant interactions shift the relationship of richness to mean range size from positive to negative at the Glacial–Holocene transition". 
All the scripts for data analysis, as well as for the generation of figures and tables are shown here. The main files included are as follows:
0_resampling
1_calculation: Based on the resampling result, this script calculates plant richness, mean range-size, and heterogeneity. The output includes the results and Supplementary Figure 1.
2_cushion_tree: Cushion plant taxa and tree plant taxa abundance change over time based on the resampling result.
3_merge_plot: Based on the resampling, plant richness to mean range-size relationship was calculated with the 5,000-year running time window, the richness, mean range-size, heterogeneity,
              and the temperature result and the isotope result of NGRIP ice core were all shown in the 5,000-year running time window. The output in this step include Figure3, Figure4,
              Supplementary Figure2, Supplementary Figure4, and Supplementary Figure7.
4_binomial_regression: In this step, we did the binomial regression, response variables of 
richness to range-size relationship (R>0.2 for positive and R<-0.2 for negative), 
 cushion plant abundance, and tree plant abundance as the explanatory variables
5_
