basemeasure = function(x, family){
  if(family == "exponential"){
      return(1)
  }
  
  if(family == "Poisson"){
    return(log(1/factorial(x)))
  }
  
 
  
  if(family == "Gaussian"){
    return((2*pi)^(-1/2))
  }
 # if(family = "Ising"){
    
    
}
