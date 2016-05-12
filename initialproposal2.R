#generates proposed particle from the independent exponential family distribution given by diagonal entries of Phi
initialproposal = function(initializingparameters,
                                    family = "Poisson", 
                                    nparticles = nparticles){
  if(family == "Poisson"){
    output = mcmapply(FUN = rpois, lambda = initializingparameters, MoreArgs  = list(n = nparticles))
  }
  
  if(family == "Gaussian"){
    output = mcmapply(FUN = rnorm, sd)
  }
  return(output)
}

#each particle generates proposed particle from the independent exponential family distribution centered at the current iterate
newproposal = function(family= "Poisson", initializingparameters = (particlehistory[i-1,1,])){
  
  if(family == "Poisson"){
    output = mapply(FUN = rpois, lambda = initializingparameters, MoreArgs = list(n=1))
    }
  return(output)
}
