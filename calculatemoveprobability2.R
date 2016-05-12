calculatemoveprobability = function(problaw = PproportionalBM,
                                    proposedmove = proposedmoves[1,],
                                    currentposition = particlehistory[i-1,1,],
                                    Phii = Phihistory[i,,], 
                                    Thetai = Thetahistory[i,],
                                    sufficientstatistic = sufficientstatistic,
                                    family = "Poisson", 
                                    basemeasure = basemeasure){
  
  
  probratio = problaw(x = proposedmove, 
                       Phi = Phii, 
                       Theta = Thetai,
                       sufficientstatistic = sufficientstatistic,
                       family = family,
                       basemeasure = basemeasure) /
    problaw(x = currentposition, 
            Phi = Phii, 
            Theta = Thetai,
            sufficientstatistic = sufficientstatistic,
            family = family,
            basemeasure = basemeasure)
  
  #if currentposition is zero, sometimes this will be NaN
  propratio = prod(mapply(FUN = dpois, x = currentposition, lambda = proposedmove)) / 
    prod(mapply(FUN = dpois, x = proposedmove, lambda = currentposition))
    
  output = propratio * probratio
  return(output)
}

calculateweights = function(problaw = Pproportional, 
                            particlehistory = particlehistory,
                            sufficientstatistic = sufficientstatistic, 
                            Phihistory = Phihistory, 
                            Thetahistory = Thetahistory){
  
  weights = array(dim = c((nsteps),nparticles))
  for(i in 2:dim(particlehistory)[1]){
    weights[i-1,] = mcmapply(FUN = Pproportional, 
                            x = split(particlehistory[i,,],
                            row(particlehistory[i,,])),
                            MoreArgs = list(Phi = Phihistory[i,,], 
                            sufficientstatistic = sufficientstatistic, 
                            Theta = Thetahistory[i,], 
                            family = "Poisson")) / 
                    mcmapply(FUN = Pproportional, 
                             x = split(particlehistory[i,,], row(particlehistory[i,,])),
                              MoreArgs = list(Phi = Phihistory[i-1,,], 
                            sufficientstatistic = sufficientstatistic, 
                            Theta = Thetahistory[i-1,],
                            family = "Poisson"))
  }
  output = apply(FUN = prod, MARGIN = 2, X = weights)
  return(output)
}