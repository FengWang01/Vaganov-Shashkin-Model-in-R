
#################################################################################
# R implementation of VSM Vaganov-Shashkin Cambial Growth Model.                #
# Version: 1.0.0                                                                #                
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

# Inputs:
  # Temp = (366 x number of years) matrix of daily temperatures (C)
  # Prec = (366 x number of years) matrix of daily precipitation (mm)
  # phi = latitude (degrees)
  # syear = first years of Temp and Prec data
  # eyear = last year of Temp and Prec data
  # parameters = structure containing model parameters as fields
  # PAR = our own estimated normalized day length. The model runs well without this.

# Output is a list containing growth and environmental metrics and final simulated ring width

vsmR = function(Temp, 
                    Prec,
                    phi,
                    syear,
                    eyear,
                    parameters
                 ){
  
  #################################
  ## INITIALIZATION AND ERROR CATHCING
  # Unlist parameters to working environment
  list2env(parameters, envir = globalenv())
  if(display){
    print('Vaganov-Shashkin Cambial Growth Model with cell size module [vsmRcell] [version 1.0]')
  }
  
  # Temp and Prec matrices must have the same structure
  if(nrow(Temp) != 366 | nrow(Prec) != 366){
    stop('Climate matrices must be 366 rows x number of years')
  }
  
  if(ncol(Temp) != ncol(Prec)){
    stop('Climate matrices must have the same number of years')
  }
  
  # The number of columns must correspond to the time period defined
  if(ncol(Temp) != (eyear-syear+1)){
    stop('The number of columns in climate matrices must match with the time period')
  }
  
  # Daily growth/environment variables
  # GrT: relative growth rates for temperature
  # GrW: relative growth rates for moisture
  # sm: soil moisture
  # snow: snow depth
  # dep: rooting depth?
  # Gr: overall daily growth rate
  # stm:?
  # GrTrans: growth rate used for transpiration
  # trans: transpiration
  # xmelt: snow melt
  # cellCount
  GrT = GrW = sm = snow = dep = Gr = stm = 
    GrTrans = trans = xmelt = cellCount = matrix(NA, nrow = 366, ncol = (eyear - syear + 1))
  
  # Cell array variables
  # SI:
  # RI:
  # AT:
  # DIV:
  # Allow a maximum number of 200 cells produced per year
  SI = RT = AT = KP = DIV = matrix(NA, nrow = 200, ncol = (eyear - syear + 1))

  # Parse operating variables from parameter and climate data
  # Parameter renaming isn't really necessary but is done to follow the original FORTRAN
  Vm = b[10]     # Fixed growth rate in S, G2, and M phases of mitosis
  Vs = b[1]      # Minimum critical growth, parameterized, otherwise dormancy
  SIG1 = b[11]   # Maximum size of cell in G1 phase
  SIM = b[14]    # Maximum cell size before mitosis occurs
  tc = b[15]     # This is the time step in the cambial model
  tn0 = 1/tc     # ... which determines the cambial model iterations per model day
  
  # set the initial soil depth to the rooting depth from the parameters
  dep[1,1] = rootd
  
  #################################
  ## SOLAR RADIATION
  # Calculate two insolation vectors, one for normal and one for leap year
  # Call the vectors as needed within the year loop depending on value of ydays
  # Preallocate for speed and use the original FORTRAN formulation
  
  # What is PAR? May not need?
  if (!exists('PAR')) {  # Check if 'PAR' exists; if not, calculate based on latitude
    parinput = 0
    # Call day_length function to return normalized daylength: normalize must be TRUE
    ndl = day_length(phi, normalize = TRUE, method = 'VS')
  } else {
    parinput = 1  # Indicate that PAR data is available
  }
  
  #################################
  ## GROWTH MODELING MODULE
  # iyear = is the calendar year the model is currently working on (e.g. 1977, 1978, 1979 ...)
  # cyear = is the integer number of simulation year (e.g. 1st, 2nd, 3rd ...)
  
  iyear = c(syear:eyear) 
  # Begin cycling over years
  
  # Some empty vectors to be filled in the loop
  ndays = ncambium = nring = fday = stday = sday = eday = nxylem = c()
  
  for (cyear in 1:length(iyear)) { 
    #cyear = 7
    if(display){
      print(paste('Modeling', iyear[cyear], '...'))
    }
    if(leap_year(iyear[cyear])){
      ndays[cyear] = 366
    }else{
      ndays[cyear] = 365
    }
    
    #################################
    ## INITIALIZE CAMBIUM
    # Set the initial state of the cambium and a few other things if this the first simulation year
    if(iyear[cyear] == syear || ndc == 0){ 
      # If this is the first year of the simulation, create a generic cambium
      ncambium[cyear] = 1                 # Create a single cambial cell ...
      SI[1,cyear] = SIM/2                 # Start it at half the size it needs to achieve mitosis
      DIV[1,cyear] = 1                    # Make it capable of division
      nring[cyear] = ncambium[cyear]      # Number of cells in the file, calculated as nring = all tracheids + all cambial cells
      RT[1,cyear] = 0                     # Initialize the characteristics of that first cambial cell for the first year
      AT[1,cyear] = 0                     # Initialize the characteristics of that first cambial cell for the first year
      KP[1,cyear] = 0                     # Initialize the characteristics of that first cambial cell for the first year
      sm[1,1] = Wo                        # Initialize soil moisture from the parameter file for the first day of the simulation
      snow[1,1] = SNo                     # Initialize snow depth from the parameter file for the first day of the simulation
      ndiv = 0                            # Number of divisions
    }else{
      # Otherwise, this information is taken from end of previous growth season's cambium
      ncambium[cyear] = ncambium[cyear-1] # New cambium is the cambium from the end of the previous year
      nring[cyear] = ncambium[cyear]      # The number of total cells is equal to the number of cambials cells, to start
      SI[1:ncambium[cyear],cyear] = SI[1:ncambium[cyear],cyear - 1]
      DIV[1:ncambium[cyear],cyear] = DIV[1:ncambium[cyear],cyear - 1]
      RT[1:ncambium[cyear],cyear] = 0
      KP[1:ncambium[cyear],cyear] = 0
      AT[1:ncambium[cyear],cyear] = 0
      ndiv = 0                            # Reset the number of divisions that have occurred
    }  
    
    #################################
    ## INITIATE GROWTH
    # Calculate growing degree days; When it is time to start growth?
    sum_temperature = -Inf # (Re)initialize sum-of-temperature variable
    tt = 1 # Begin on the first day of the year (day counter)

    while(sum_temperature < Tg && tt <= ndays[cyear]){ # Should be <=?
      # As long as sum-of-temperatures has not crossed threshold to begin growth ...
      # Used sum (for Octave) over the period set in the parameter K(9)
      sum_temperature = sum(Temp[tt:(tt+K[9]),cyear]) 
      fday[cyear] = tt+K[9] # Set the first day to begin growth to the current day (in case we cross the threshold)
      tt = tt + 1                      # Increment the day counter
    }
    # Once growing degree days have crossed threshold, cambial activity can begin
    
    #################################
    ## SOIL THAWING
    # Calculate growing degree days; When is it time to start soil thaw?
    sum_st_temperature = -Inf # (Re)initialize sum-of-temperature variable
    st = 1 
    
    while(sum_st_temperature < Tm && st < ndays[cyear]-K[10]){
      # As long as sum-of-temperatures has not crossed threshold to begin thaw ...
      # Sum (for Octave), over the period set in the parameter K(10)
      sum_st_temperature = sum(Temp[st:(st+K[10]),cyear])
      stday[cyear] = st + K[10]
      st = st + 1
      
      if (st == ndays[cyear] - (K[10] + 1)) {
        stday[cyear] = ndays[cyear]
      }  # Condition where soil thaw never begins (uncommon at best)
    }

    # Calculate soil thaw and rooting depth for every day of the year, 
    # which modifies soil moisture availability and other parameters
    if(K[1] == 0){ # If soil thaw is turned 'off' ...
      dep[1:ndays[cyear],cyear] = rootd;  # Rooting depth is just a constant based on the parameters
    }else{
      dep[1:ndays[cyear],cyear] = 0  # The original defaults to refreezing the soil at the beginning of each year
      for (i in (stday[cyear]+1):ndays[cyear]){
        dep[i,cyear] = dep[i-1,cyear] + a[1] * Temp[i-1,cyear] * exp(-a[2] * dep[i-1,cyear])
        if(dep[i,cyear] < 0){ # Error catching, has to be a minimum root depth, of course
          dep[i,cyear] = 0
        }
        if(is.na(dep[i,cyear])){ # Error catching
          dep[i,cyear] = dep[i-1,cyear]
          if(is.na(dep[i,cyear])){
            dep[i,cyear]=rootd
          }
        }
      }
    }
    
    #################################
    ## -- day cycle -- ##
    # Begin cycling over days in a year, starting with the first day (fday) calculated above
    for (t in 1:ndays[cyear]) {
      #t=1
      tn = tn0 # Number of times the cambium simulation will run per day
      
      ## CALCULATE EXTERNAL GROWTH RATES, GrT(t) and GrW(t) based on the piece-wise linear growth function
      ## We'll do this twice, once for "Temp"erature, once for soil moisture ('W'ater)
      
      ## This procedure is simplied because I define a function to do it.
      # First, temperature
      x = Temp[t,cyear]
      GrT[t,cyear] = piece_wise_growth_rate(x = x, para = Tf)
 
      x = sm[t,cyear]
      GrW[t,cyear] = piece_wise_growth_rate(x = x, para = Wf)
      
      # Are we in a leap year?
      if(parinput == 0){ # If we are using normalized daylength ...
        if(ndays[cyear]==366){
          GrE = ndl[,2]
        }else{
          GrE = ndl[,1]
        }
      }else{
        GrE = PAR[,cyear] # If we passed our own estimate (from PAR, etc.) of GrE to the model
      }  
      
      #% modifier for growth rate based on depth, but only if the thaw depth is less than the root depth
      if(dep[t,cyear] < rootd){
        stm[t,cyear] = dep[t,cyear]/rootd
      }else{
        stm[t,cyear] = 1
      }
     
      #################################
      # CAMBIAL MODEL
      # Daily growth rate calculation
      Gr[t,cyear] = GrE[t] * min(GrT[t,cyear],GrW[t,cyear]*stm[t,cyear])
      
      # Start growth?
      # Only enter the growth and cambial blocks if it's time to start growth, otherwise, skip to end of the day cycle
      if(t >= fday[cyear]){
        
        #################################
        ## -- cambial (c) cycle -- ##
        # Begin cambial model cycle (c cycle), at fraction-of-a-day steps (e.g. 5 times per day)
        
        for(c in 1:tn){
          # Begin with the first cell in the file, working from the cambial initial outward
          
          #################################
          ## -- j cycle -- ##
          # While there are still cells in the file to be operated on ...
          # Vmin: Minimum growth rate
          # Vo: Position specific modifier for Gr(t)
          # V: Position specific growth rate
          j = 1
          
          Vmin = Vo = V = c()
          while(j <= nring[cyear]){
            #j=1
            # Minimum growth rate for position in the cell file
            Vmin[j] = b[9] * (exp((b[8] + j*b[5])*0.4) - 5)
            # Position specific modifier for Gr(t)
            Vo[j] = (b[6] * j * b[5]) + b[7]      
            # Position specific growth rate
            V[j] = Vo[j] * Gr[t,cyear] * b[4]
            
            # Proceed through the j-cycle decision tree
            if(DIV[j,cyear] != 0){ # If the cell hasn't differentiated yet ...
              if(SI[j,cyear] <= SIG1){ # If the cell in position j is not yet large enough to divide ...
                if(V[j] < Vs){ # Check if the cell enters dormancy
                  j = j + 1
                }else{
                  if(V[j] < Vmin[j]){ # Check if the cell is growing too slowly, if it is ...
                    DIV[j,cyear]= 0 # The cell can no longer divide
                    RT[j,cyear] = t+(c*tc) # The time it exited the cambium (differentiated) is recorded
                    j = j + 1
                  }else{ # If the cell is growing fast enough ...
                    SI[j,cyear] = SI[j,cyear] + (V[j]*tc) # Increment cell size as a function of growth rate
                    AT[j,cyear] = AT[j,cyear] + tc # Increment time in cambium by cambium time step (fraction of a day)
                    j = j + 1
                  }
                }
              }else{ 
                # If the cell size is larger than the maximum for phase G1, start mitotic cycle (S, G2, and M phases)
                SI[j,cyear] = SI[j,cyear] + (Vm*tc) # Increment the time spent in the cambium
                AT[j,cyear] = AT[j,cyear] + tc # While in the mitotic cycle (~S,G2) grow at constant rate
                
                if(SI[j,cyear] <= SIM){ # Is the cell big enough to divide yet?
                  j = j + 1 # If not, go to the next cell in the file and work on it
                }else{ # If it is, it's time to divide
                  ndiv = ndiv + 1 # Log this division
                  nring[cyear] = nring[cyear] + 1 # Cell Divides! now, increase the number of cells in the ring by one ...
                  SI[(j+2):nrow(SI),cyear] = SI[(j+1):(nrow(SI)-1),cyear]  # and move everything in the file to account for the 'new' cell
                  AT[(j+2):nrow(AT),cyear] = AT[(j+1):(nrow(AT)-1),cyear]  # daughter cells are rows above mother cell
                  RT[(j+2):nrow(RT),cyear]  = RT[(j+1):(nrow(RT)-1),cyear]
                  KP[(j+2):nrow(KP),cyear]  = KP[(j+1):(nrow(KP)-1),cyear]
                  DIV[(j+2):nrow(DIV),cyear] = DIV[(j+1):(nrow(DIV)-1),cyear]
                  
                  # The 'new' cell, in position j
                  RT[j,cyear] = 0 # New cell hasn't yet differentiated ...
                  SI[j,cyear] = SI[j,cyear]/2 # New cell is half the size of the mother cell ...
                  KP[j,cyear] = KP[j,cyear] + 1 # Increment the number of times the cell has divided by one
                  DIV[j,cyear] = DIV[j,cyear] # New cell is capable of division
                  AT[j,cyear] = AT[j,cyear]
                  
                  # The new cell, in position j+1
                  j = j + 1
                  SI[j,cyear] = SI[j - 1,cyear] # New cell is half the size of the mother cell ...
                  AT[j,cyear] = AT[j - 1,cyear]
                  RT[j,cyear] = RT[j - 1,cyear]
                  KP[j,cyear] = KP[j - 1,cyear]
                  DIV[j,cyear] = 1 # New cell is capable of division
                  j = j + 1
                }
              }# finished with current position j, ready to go to the next
            }else{
              j = j + 1
            }
          } # End j cycle
          # Before the beginning of the next cambial (c) cycle, account for state of the cellular file in a multidimensional array
          # not currently operational for speed purposes, since not yet used elsewhere
          # ccells(t,c,cyear) = nring(cyear);                            % number of total cells
          # ccambium(t,c,cyear) = sum(~isnan(DIV(:,cyear)));             % number of those cells which can divide (e.g. are still in the cambium)
          # cxylem(t,c,cyear) = ccells(t,c,cyear) - ccambium(t,c,cyear); % number of differentiated cells
        } # End c-cycle
      } # End of day (t-cycle) here; soil moisture, transpiration, snow melt, soil thaw
      
      #If growth wasn't started yet, but we still need to calculate snow melt, changes in soil moisture, etc., the day cycle jumps to here
      
      #################################
      ## CALCULATE ENVIRONMENTAL VARIABLES
      # Calculate growth rate to use for the transpiration calculation (modified by rooting depth)
      if(t < fday[cyear]){
        GrTrans[t,cyear] = 0 # no transpiration if no growth
      }else{
        GrTrans[t,cyear] = Gr[t,cyear]
      }
      
      # Transpiration
      trt = k[2] * exp(Temp[t,cyear] * k[3])
      trans[t,cyear]  = trt * GrTrans[t,cyear]
      
      # Calculate snow melting and modified precipitation
      # First, did it snow or rain? How much? Is there still snow, and did it melt?
      if(K[4] == 0) {
        snow[t,cyear] = 0 # Snow?
        prsnow = 0 # Precipitation in form of snow
        xmelt[t,cyear] = 0 # Melted snow
        prrain = Prec[t,cyear] # Precipitation in form of rain
      }else if(Temp[t,cyear] <= SNmt){
        prsnow = Prec[t,cyear]
        prrain = 0 # No rain
        xmelt[t,cyear] = 0
      }else{
        prrain = Prec[t,cyear]
        prsnow = 0
        xmelt[t,cyear] = SNr * (Temp[t,cyear] - SNmt) # Snow melts, if there is any.
        if(is.na(xmelt[t,cyear])){
          xmelt[t,cyear] = 0}
        if((snow[t,cyear] - xmelt[t,cyear]) < 0){
          xmelt[t,cyear] = snow[t,cyear]}
        if(xmelt[t,cyear] < 0){
          xmelt[t,cyear] = 0}
      }
      
      # Next, we'll figure out how much snow there will be for the next day
      if(t < ndays[cyear]){
        snow[t+1,cyear] = snow[t,cyear] - xmelt[t,cyear] + prsnow # Snow for the next day is a function of melt and new snowfall
        if(snow[t+1, cyear] < 0){
          snow[t+1,cyear] = 0 # Make sure we don't have negative snow depths
        }
      }else if(t == ndays[cyear] && iyear[cyear] != eyear){
        # The last day of year
        snow[1,cyear+1] = snow[t,cyear] - xmelt[t,cyear] + prsnow # Snow for the next day is a function of melt and new snowfall
        if(snow[1,cyear+1] < 0){
          snow[1,cyear+1] = 0 # Make sure we don't have negative snow depths
        }
      }
        
      # Calculate soil moisture for the next day
      K[6] = min(k[1] * prrain, Pmax) # Soil moisture recharge cannot exceed maximum daily infiltration
      xm = min(xmelt[t,cyear],Pmax) # Soil moisture recharge from snow melt cannot exceed maximum daily infiltration
      w = sm[t,cyear] * dep[t,cyear] * (1 - rated) + K[6] - trans[t,cyear] + xm # Available water used for soil moisture calculation
        
      if(is.na(w) | w <=0){w=0} # Error catching, available snow can't be less than zero
        
      if(t < ndays[cyear]){
        if(dep[t+1,cyear] > 0){#... as long as soil isn't completely frozen ... 
          sm[t+1,cyear] = w/dep[t+1,cyear]
          if(sm[t+1,cyear] <= Wmin){ # Error catching
            sm[t+1,cyear] = Wmin}
          if(sm[t+1,cyear] >= Wmax){ # Error catching
            sm[t+1,cyear] = Wmax}
          if(is.na(sm[t+1,cyear])){
            sm[t+1,cyear] = Wmin}
          }else{
            sm[t+1,cyear] = sm[t,cyear]}
        # Check for soil moisture above Wmax or below Wmin
        # This is the last day of the year, so sm data for next day needs to be saved in sm[1,cyear+1]
        }else if(t == ndays[cyear] && iyear[cyear] != eyear){ 
          if(dep[t,cyear] > 0){
            sm[1,cyear+1] = w/dep[t,cyear]
          }else{
            sm[1,cyear+1] = Wmin}
          if(sm[1,cyear+1] <= Wmin){ # Error catching
            sm[1,cyear+1] = Wmin}
          if(sm[1,cyear+1] >= Wmax){ # Error catching
            sm[1,cyear+1] = Wmax}
          if(is.na(sm[1,cyear+1])){
            sm[1,cyear+1] = Wmin}
        } 
      cellCount[t,cyear] = nring[cyear]
    } # End day (t) cycle
  
  # Cambial state
  # Accounting for the state of the cambium and other annually resolved variables
  i = which(!is.na(DIV[,cyear]))
  ncambium[cyear] = sum(DIV[i,cyear]) # Number of cells in the cambium for the start of the next year
    
  ii = intersect(which(RT[, cyear] > 0), which(!is.na(RT[, cyear]))) # Look for non-NaN and non-zero values within RT
  if(length(ii) > 0){ # Not empty
    # Record start of growth (sday) as time of first differentiation
    sday[cyear] = RT[ii[length(ii)],cyear] 
    # Record end of growth (eday) as time of last differentiation
    eday[cyear] = RT[ii[1],cyear] 
  }else{
    sday[cyear] = eday[cyear] = NA
  }
    
  # Number of xylem (non cambium) cells
  nxylem[cyear] = nring[cyear] - ncambium[cyear] # Number of tracheid cells in the ring for the year
  } # End of year cycle
  
 ## CELL SIZE MODULE
 #if(cellsize){}
  
 ## FINAL ACCOUNTING AND OUTPUTS
 # Create two kinds of normalized ring width chronology using number of rings
 trw = nxylem/mean(nxylem, na.rm = T) # calculate tree-ring width by cell number (good approximation)
 trws = colSums(SI, na.rm = T)/mean(colSums(SI, na.rm = T)) # calculate tree-ring width by cell 'size' (ring width as an of integration of Gr as filtered through the cambial model)

 output = list()
 # Write output data to a single structure
 output$startYear = syear
 output$endYear = eyear
 output$simulationLength = (eyear - syear + 1)
 output$latitude = phi
 output$trw = trw
 output$trws = trws
 output$nXylem = nxylem
 output$nCambium = ncambium
 output$Gr = Gr
 output$GrT = GrT
 output$GrW = GrW
 output$GrE = GrE
 output$transpiration = trans
 output$sm = sm
 output$snowdepth = snow
 output$snowmelt = xmelt
 output$soilDepth = dep
 output$fday = fday
 output$startDay = sday
 output$endDay = eday
 output$parameters = parameters
 output$difftime = RT
 output$cellCount = cellCount
  
 return(output)
}
