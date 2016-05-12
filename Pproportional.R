#P proportional evaluates the part of the probability that depends on parameters and data
#base measure cancels out in the ratio
Pproportional = function(x, Phi, Theta, family, sufficientstatistic){
  output = exp(Theta %*% sqrt(sufficientstatistic(x, family)) +  
                 sqrt(sufficientstatistic(x, family = family)) %*% Phi %*% sqrt(sufficientstatistic(x, family = family)))
  return(output)
}

PproportionalBM = function(x, Phi, Theta, family, sufficientstatistic, basemeasure){
  output = exp(Theta %*% (sign(x) * sqrt(sufficientstatistic(x, family))) +  
                 (sign(x) * sqrt(sufficientstatistic(x, family))) %*% Phi %*% (sign(x) * sqrt(sufficientstatistic(x, family)))
               + sum(basemeasure(x, family)))
  return(output)
}
