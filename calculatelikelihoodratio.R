calculatemoveprobability = function(problaw = PproportionalBM, 
                                    proposedmoves = split(proposedmove, row(proposedmove)),
                                    currentposition = split(particlehistory[i-1,,], row(particlehistory[i-1,,])),
                                    Phii = Phihistory[i,,], 
                                    Thetai = Thetahistory[i,],
                                    sufficientstatistic = sufficientstatistic,
                                    family = "Poisson", 
                                    basemeasure = basemeasure){
  
  output = mcmapply(FUN = problaw, 
                             x = proposedmoves,
                             MoreArgs = list(Phi = Phii,
                                             sufficientstatistic = sufficientstatistic,
                                             Theta = Thetai, 
                                             family = "Poisson",
                                             basemeasure = basemeasure)) /
    mcmapply(FUN = problaw,
             x = currentposition,
             MoreArgs = list(Phi = Phii, 
                             sufficientstatistic = sufficientstatistic, 
                             Theta = Thetai, 
                             family = "Poisson",
                             basemeasure = basemeasure))
  
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