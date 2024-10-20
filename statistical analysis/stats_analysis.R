"
Author: Pedro Gonzalez & Xinru Li
Date: Jun. 2022
University of British Columbia
"

library(lattice)
library(onewaytests)
library(ggpubr)
library(dplyr)
library(car)
library(FSA)



################# ================== ###################
################ Levene's Test ##################
################# ================== ###################

# Evaluate the significance of the differences in the proportion of bleaching reports in each bleaching severity category for the reports with priming conditions vs. without priming conditions 
# load data
df <- read.csv("KW&Levene_HSpk1_c749/KW_Levene_test_inputs_DHD2870_c749.csv")

my.data = stack(list(g1=df$bs_dhd2830_ap, g2=df$bs_dhd2830_np))  # replace input testing for other MHW bins
leveneTest(values~ind, data=my.data, center=median)


################# ================== ###################
################ Kruskal-Wallis Test ##################
################# ================== ###################

# Evaluate the significance of the differences in the proportion of bleaching reports in each bleaching severity category for the reports with priming conditions vs. without priming conditions 
kruskal.test(df$bs_dhd2830_ap ~ df$bs_dhd2830_np, data = df)  # replace input testing for other MHW bins

# Evaluate the significance of the differences in the duration and magnitude of the priming periods between different categories of bleaching reports
df <- read.csv("FIGAC5_c749_Dr_DHD3036_ap.csv")

my.data = stack(list(g1=df$Dp_c0, g2=df$Dp_c1, g3=df$Dp_c2, g4=df$Dp_c3))
my.data = stack(list(g1=df$Ap_c0, g2=df$Ap_c1, g3=df$Ap_c2, g4=df$Ap_c3))
kruskal.test(values~ind, data=my.data) 


################# ================== ###################
################ Post Dunn's Test ##################
################# ================== ###################
# Evaluate significance of the differences in the duration and magnitude of the priming periods between different categories of bleaching reports
dunnTest(values~ind, data=my.data)


