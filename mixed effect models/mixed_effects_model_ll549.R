"
Author: Xinru Li
Date: Jul. 2022
University of British Columbia
"

################# ================== ###################
###### Mixed Effects Regression | R Data Analysis ######
################# ================== ###################

require(ggplot2)   # to plot 
require(lme4)      # to perform mixed effect models
library(stargazer) # to make table 
library(MuMIn)     # to compute r2
library(lmerTest)  # to get the p-value
library(sjPlot)    # to plot estimates
library(sjlabelled)
library(sjmisc)
library(ncf)       # to compute MORAN's Index and Mantel test
library(gdata)
library(lmerTest)
library(lme4)
library(nlme)
library(car)
library(glmmTMB)



#############################HSpeak###################################################
dp_MDL3c5_SEL0123_28_nnrm_lhLim49 <- read.csv("MEM_model3_c5_inputs_lhLim49_DHD2856_nnrm_npap_BP.csv")
dp_MDL3c5_SEL0123_56_nnrm_lhLim49 <- read.csv("MEM_model3_c5_inputs_lhLim49_DHD56_nnrm_npap_BP.csv")


########## =============== ##########
### Perform a mixed effects model ###
########## =============== ##########


# perform lmer
m1 <- lmer(BP ~ DHD * Atrc5 +
             (1 | lat), data = dp_MDL3c5_SEL0123_2856_nnrm_lhLim49)
summary(m1)

m2 <- lmer(BP ~ DHD * Dtrc5 + 
             (1 | lat), data = dp_MDL3c5_SEL0123_2856_nnrm_lhLim49)
summary(m2)

m3 <- lmer(BP ~ DHD * Drec5 + 
             (1 | lat), data = dp_MDL3c5_SEL0123_2856_nnrm_lhLim49)
summary(m3)

m4 <- lmer(BP ~ DHD + Dtrc5 + Atrc5 + DHD * Dtrc5 + DHD * Atrc5 +
             (1 | lat), data = dp_MDL3c5_SEL0123_2856_nnrm_lhLim49)
summary(m4)

m5 <- lmer(BP ~ DHD * Atrc5 +
             (1 | lat), data = dp_MDL3c5_SEL0123_56_nnrm_lhLim49)
summary(m1)

m6 <- lmer(BP ~ DHD * Dtrc5 + 
             (1 | lat), data = dp_MDL3c5_SEL0123_56_nnrm_lhLim49)
summary(m2)

m7 <- lmer(BP ~ DHD * Drec5 + 
             (1 | lat), data = dp_MDL3c5_SEL0123_56_nnrm_lhLim49)
summary(m3)

m8 <- lmer(BP ~ DHD + Dtrc5 + Atrc5 + DHD * Dtrc5 + DHD * Atrc5 +
             (1 | lat), data = dp_MDL3c5_SEL0123_56_nnrm_lhLim49)
summary(m4)


