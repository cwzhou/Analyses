# Title: Libraries for Simulation Scripts
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 07.02.2023

# Libraries --------------------------------------------------------------
library(tidyverse)
library(tidyr)
library(survival)
library(dplyr)
# library(randomForestSRC)
library(knitr)
if (local == 1){
  library(itrSurv)
} else{
  # install most recent version of method in the cluster using terminal as needed.
  
  # # Install itrSurv:
  #   # In cluster:
  #   install.packages('~/Dissertation/Paper_1/Paper1_Method/itrSurv_backup/itrSurv_0.1.0.tar.gz', repos = NULL, type='source')
  # library(itrSurv)
  # 
  # devtools::install('/nas/longleaf/home/cwzhou/Dissertation/Paper_1/Paper1_Method/itrSurv_working_package/itrSurv')
  # 
  # R CMD INSTALL -l /nas/longleaf/home/cwzhou/Dissertation/Paper_1/Paper1_Method/itrSurv_working_package/ itrSurv_0.1.0.tar.gz
  # ll -t /nas/longleaf/home/cwzhou/R/x86_64-pc-linux-gnu-library/4.1
  
  # note this is in the 4.1 library, which may be outdated
  library(itrSurv, lib.loc = '/nas/longleaf/home/cwzhou/R/x86_64-pc-linux-gnu-library/4.1')
}
# library(dtrSurv) # we include this within the body scripts directly
library(MASS)
library(reshape2)
library(ggplot2)
library(cowplot)
# for pmcr:
library(rgenoud)
library(rpart)

if (local == 1){
  message("Local Machine")
  library(beepr)
  # beep()
}
# for aipwe:
library(purrr)
library(cmprsk)

message("End of F00.Simulation_Libraries.R")
