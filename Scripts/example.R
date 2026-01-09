
#################################################################################
# Example code to run the Vaganov-Shshkin Model.                                #
# Author: Feng Wang.                                                            #
# Date: January 2026                                                            #
#                                                                               #
# Citation:                                                                     #
# Wang, F., Wise, E.K., Anchukaitis, K.J., Chang, Q., Dannenberg, M.P., 2026.   #
# Evaluation of daily gridded climate products using in situ FLUXNET data and   #
# tree growth modeling. Environmental Research Letters.                         #
#                                                                               #
# We are preparing an R package ("virtualRings") that incorporates this model.  #
# Please cite also the package and the related article when released.           #
#################################################################################

#===========================================
# Preparations
#===========================================

### Clean working environment
rm(list = ls())

### Set up your working directory
setwd("xxx")

#===========================================
# Load functions and packages
#===========================================
source("Scripts/vsmR_v1.0.0.R")
source("Scripts/piece_wise_growth_rate.R")
source("Scripts/day_length.R")

library(lubridate)

#===========================================
# Read in daily climate data for ca691 ITRDB site (climate data from Gridment)
#===========================================
temp <- readRDS("./Data/Gridmet_Tmean_ca691.RDS")
prcp <- readRDS("./Data/Gridmet_Prcp_ca691.RDS")

#===========================================
# Read in VSM parameter list
#===========================================
# Note: User should change parameters in the build_parameters.R file.
# Here, we only use default parameters which may not lead to an optimal simulation result
source("Scripts/build_parameters.R")

#===========================================
# Set up year range and latitude
#===========================================
sy <- 1980      # Start year
ey <- 2011      # End year
phi <- 39.41667 # Latitude of your own site

#===========================================
# Run VSM
#===========================================
example_run <- vsmR(Temp = temp,
                    Prec = prcp,
                    syear = sy,
                    eyear = ey,
                    phi = phi,
                    parameters = parameters) 
                        
