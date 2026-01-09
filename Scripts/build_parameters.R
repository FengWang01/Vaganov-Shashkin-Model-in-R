
#################################################################################
# VSM parameter setup be to used in the main vsmR() function.                   #
# Authors: Feng Wang, Kevin J. Anchukaitis.                                     #
# Date: February 2025                                                           #
# It is mostly 'translated' from the Matlab codes of Anchukaitis et al., 2020.  #
#                                                                               #
# Citation:                                                                     #
# Wang, F., Wise, E.K., Anchukaitis, K.J., Chang, Q., Dannenberg, M.P., 2026.   #
# Evaluation of daily gridded climate products using in situ FLUXNET data and   #
# tree growth modeling. Environmental Research Letters.                         #
#                                                                               #
# We are preparing an R package ("virtualRings") that incorporates this model.  #
# Please cite also the package and the related article when released.           #
#################################################################################

parameters = list(
  # Set to 1 to write headers and program progress to screen, 0 for improved speed and no screen output
  display = 0,  
  ## Main Program and Growth Block Parameters 
  # Growth function parameters for temperature
  Tf = c(3.0000,    # Minimum temperature (C) for growth   
         16.0000,   # Growth rate is max in the range T2-T3 (lower optimal temperature, C)   
         22.0000,   # Growth rate is max in the range T2-T3 (upper optimal temperature, C)
         30.0000),  # Maximum temperature (C) for growth
  
  Wf = c(0.01,      # Minimum soil moisture for growth (v/v)
         0.015,     # The growth rate is max in range W2-W3 (lower optimal soil moisture, v/v)
         0.8,       # The growth rate is max in range W2-W3 (upper optimal soil moisture, v/v)
         0.9),      # Growth is stopped at this soil moisture (v/v)
  
  SNo = 300.0000,   # Initial snowpack (mm)
  Wo = 0.2183,      # Initial soil moisture (v/v)
  Wmax = 0.7000,    # Maximum soil moisture (field capacity, (v/v)
  Wmin = 0.0100,    # Minimum soil moisture (wilting point, v/v)
  rootd = 500.000,  # The root/soil melt depth (mm)
  rated = 0.0051,   # The rate of water drainage from soil
  Pmax  = 30.0000,  # Maximum rate of infiltration water into soil  (mm/day)
  
  # Interception and transpiration parameters
  k = c(0.9000,     # k1 (1-k1) is the interception precipitation by the tree crown
        0.2500,     # 1st coefficient for calculation the transpiration
        0.1800),    # 2nd coefficient for calculation the transpiration
  
  # Soil and snow melting parameters
  Tm = 20.000,      # Sum of temperature for start soil melting (C)
  a = c(10.000,     # 1st coefficient of soil melting
        0.0060),    # 2nd coefficient of soil melting
  Tg = 30.0000,     # Sum of temperature to start growth (C)
  SNr = 3.0000,     # The rate of snow melting (mm/C/day)
  SNmt = 0.0000,    # Minimum temperature for snow melting
  
  # Some switches and variable storage
  # If not specified in the Matlab version (Anchukaitis et al., 2020), a value of zero is given
  K = c(0,          # K[1], Soil melting switch: yes = 1; no = 0
        NA,
        NA,
        1,          # K[4], Snow melting switch: yes = 1; no = 0
        NA,
        NA,
        NA,
        50,         # K[8], Maximum duration (days) of latewood formation
        14,         # K[9] the period (days) over which to sum temperature to calculate start of growth
        10),        # K[10] the period (days) over which to sum temperature to calculate start soil melting
          
  # Parameters for Cambial Block
  ndc = 1,          # Use previous cambium for following year? 1 = yes (dynamic cambium), 0 = no(static);    
  
  # Growth rate parameters
  # If not specified in the Matlab version (Anchukaitis et al., 2020), a value of zero is given
  b = c(0.04,       # b[1], The critical growth rate (Vcr or Vs)
        NA,
        NA,
        0.6,        # b[4], The correction of growth rate (Gr*b[4])
        1.9,        # b[5], The correction of Vo(j*b[5]) and Vmin(j*b[5])
        0.25,       # b[6], Vo[j]= b[6]*j+b[7]
        0.42,       # b[7], b[7] and b[6] determine the function Vo[j]
        2.5,        # b[8], Vmin[j]=(EXP((b[8]+j)*0.4)-5.0)*b[9]
        0.04,       # b[9], b[8] and b[9] determine the function Vmin[j]
        1.00,       # b[10], The growth rate in the S, G2 & M phases 
        8.00,       # b[11], The maximum size of a cell in the G1 phase (SIG1)
        9.00,       # b[12], The maximum size of a cell in the S phase
        9.50,       # b[13], The maximum size of a cell in the G2 phase
        10.00,      # b[14], The maximum size of a cell in the M phase (SIM)
        0.2)        # b[15], The time step in cambium model (1/b(15)) = number of c-cycles/day)
)
