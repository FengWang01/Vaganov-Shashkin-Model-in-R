
#################################################################################
# day_length function to be called in the main vsmR() function.                 #
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

### CALCULATE DAYLENGTH
day_length = function(phi,            # Latitude
                      normalize = T, # If FALSE, daylength in hours will be returned
                      method = 'Friend'){ 
  # Convert latitude to radians
  latr = (pi / 180) * phi  
  
  # Initialize matrices (dataframe) for normalized day length and other variables
  ndl = dtsi = hdl = y = sd = ha = matrix(NA, nrow = 366, ncol = 2)
  wcolumn = 0
  
  for (t in 365:366) {  # Iterate for normal (365) and leap years (366)
    #t=365
    wcolumn = wcolumn + 1
    #lt[wcolumn] = length(1:t) # This is not needed in R
    
    ### NOTE: There is small difference for calculation 
    if(normalize){
      # Solar declination
      sd[1:t, wcolumn] = asin(sin(pi * 23.5 / 180) * sin(pi * (((1:t) - 80) / 180)))
      # Calculate y based on latitude and declination
      y[1:t, wcolumn] = -tan(rep(latr, t)) * tan(sd[1:t, wcolumn])
    }else{
      if(method != 'Friend'){
        # Solar declination
        sd[1:t, wcolumn] = asin(sin(pi * 23.5 / 180) * sin(pi * (((1:t) - 80) / 180)))
        # Calculate y based on latitude and declination
        y[1:t, wcolumn] = -tan(rep(latr, t)) * tan(sd[1:t, wcolumn])
      }else{
        # Equations based on Friend et al., 2022 https://doi.org/10.1038/s41467-022-35451-7
        # Solar declination
        sd[1:t, wcolumn] = (pi / 180) * 23.45 * sin(((pi / 180) * 360 / 365.25) * ((1:t) - 80))
        
        # Calculate y based on latitude and declination
        # This is from http://www.jgiesen.de/astro/suncalc/calculations.htm
        y[1:t, wcolumn] = cos((pi / 180) * 90.833) / (cos(latr) * cos(sd[1:t, wcolumn])) -
          tan(latr) * tan(sd[1:t, wcolumn])
      }
    }
   
    y[y[!is.na(y[, wcolumn]), wcolumn] >= 1, wcolumn] = 1 # NA in Row 366 may lead to errors
    y[y[!is.na(y[, wcolumn]), wcolumn] <= -1, wcolumn] = -1 # NA in Row 366 may lead to errors
    
    # ha is named hdl in the original version of Anchukaitis et al., 2020
    ha[1:t, wcolumn] = acos(y[1:t, wcolumn])
    
    if(normalize){
      # dtsi: Daily total solar insolation = ha*sin(latr)*sin(sd) + cos(latr)*cos(sd)*sin(ha)
      # Note: Solar constant of 1361 W/m2 is not used here, because we calcualte normalized day length
      dtsi[1:t, wcolumn] = (ha[1:t, wcolumn] * sin(latr) * sin(sd[1:t, wcolumn])) +
        (cos(latr) * cos(sd[1:t, wcolumn]) * sin(ha[1:t, wcolumn]))
      
      # ndl: Normalized day length
      # Note that this way to calculate normalized day length also accounts for insolation!
      # Normalization based on day length in hours is not the same!
      ndl[1:t, wcolumn] = dtsi[1:t, wcolumn] / max(dtsi[1:t, wcolumn], na.rm = TRUE)
    }else{
      # hdl: day length in hours
      hdl[1:t, wcolumn] = 2 * (ha[1:t, wcolumn] / (pi / 180)) / 15
    }
  }
  if(normalize){
    return(ndl)
  }else{
    return(hdl)}
}
