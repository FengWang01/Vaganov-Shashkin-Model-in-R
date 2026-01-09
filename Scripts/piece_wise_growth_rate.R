
#################################################################################
# Piece_wise_growth_rate function to be called in the main vsmR() function.     #
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

### CALCULATE EXTERNAL GROWTH RATES, GrT(t) and GrW(t) based on the piece-wise linear growth function
# See inset of Figure 1 in Anchukaitis et al., 2020 for details
piece_wise_growth_rate = function(x,          # Input daily temperature or soil moisture
                                  para        # Defined parameters of min, optimal lower, optimal upper, and max 
                                  ){
  #para = parameters$Tf
  if(is.na(x)){
    gr = NA
  }else{
    if(x < para[1]){
      # If temp/soil moisture  is smaller than Tmin, relative growth rate is zero
      gr = 0
    }else if((x >= para[1]) && (x <= para[2])){
      # If temp/soil moisture is between min and optimal lower, relative growth rate is the ratio
      gr = (x - para[1])/(para[2] - para[1])
    }else if((x >= para[2]) && (x <= para[3])){
      # If temp/soil moisture is between optimal lower and optimal upper, relative growth rate is 1
      gr = 1
    }else if((x >= para[3]) && (x <= para[4])){
      gr = (para[4] - x)/(para[4] - para[3])
    }else{
      gr = 0
    }
  }
}

