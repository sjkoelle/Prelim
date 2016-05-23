

#P proportional evaluates the part of the probability that depends on parameters and data
#base measure cancels out in the ratio
Pproportional = function(x, Phi, Theta, family, sufficientstatistic){
  output = exp(Theta %*% sqrt(sufficientstatistic(x, family)) +  
                 sqrt(sufficientstatistic(x, family = family)) %*% Phi %*% sqrt(sufficientstatistic(x, family = family)))
  return(output)
}

# PproportionalBM = function(x, Phi, Theta, family, sufficientstatistic, basemeasure){
#   output = exp(Theta %*% (sign(x) * sqrt(sufficientstatistic(x, family))) +  
#                  (sign(x) * sqrt(sufficientstatistic(x, family))) %*% Phi %*% (sign(x) * sqrt(sufficientstatistic(x, family)))
#                + sum(basemeasure(x, family)))
#   return(output)
# }

PproportionalBM = function(x, Phi, Theta, family, sufficientstatistic, basemeasure){
  output = exp(Theta %*% (sign(x) * sqrt(sufficientstatistic(x, family))) +  
                 (sign(x) * sqrt(sufficientstatistic(x, family))) %*% Phi %*% (sign(x) * sqrt(sufficientstatistic(x, family)))
               + sum(mapply(FUN = basemeasure, x = as.list(x), MoreArgs = list(family = family))))
  return(output)
}


sufficientstatistic = function(x, family){
  
  
  
  if(family == "Poisson"){
    # ss = function(x){
    return(x)
    #}
  }
  if(family == "Gaussian"){
    #ss = function(x){
    return(x^2)
    # }
  }
  if(family == "exponential"){
    #ss = function(x){
    return(x)
    #}
  }
  
  # return(ss)
}


basemeasure = function(x, family){
  if(family == "exponential"){
    return(0)
  }
  
  if(family == "Poisson"){
    return(log(1/factorial(x)))
  }
  
  
  
  if(family == "Gaussian"){
    return(log(2*pi) * (-1/2))
    #return((2*pi) ^ (-1/2))
  }
  # if(family = "Ising"){
  
  
}

