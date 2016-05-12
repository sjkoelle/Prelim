#this function computes the node conditional partition integral from an eta1, an eta2, 
#the sufficient statistics function, and the base measure function
Anode = function(eta1Anode = eta1[i,s], eta2Anode = eta2[i,s], 
                 sufficientstatistic = sufficientstatistic, 
                 basemeasure = basemeasure, 
                 measure = "lebesgue",
                 family = "exponential"){
  
  integrand = function(x){
    return(exp(eta1Anode * sqrt(sufficientstatistic(x, family = family)) +
                 eta2Anode * sufficientstatistic(x, family = family) + 
                 basemeasure(x, family = family)))
  }
  
  #count data sum until high values
  if(measure == "counting"){
    partitionintegral = mean(integrand(0:1000))
  }
  
  #lebesgue
  if(measure == "lebesgue"){
    partitionintegral = integral(integrand, 0, xmax = Inf)
  }
  
  return(partitionintegral)
  #return(integrand)
}