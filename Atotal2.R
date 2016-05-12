#use annealed importance sampling to compute normalizing constant 
Atotal = function(Theta, Phi, nsteps, family){
  
  p = length(Theta)
  Phidiag = Phi
  Phioffdiag = Phi
  for(s in 1:p){
    for(t in 1:p){
      if(s==t){
      Phidiag[s,t] = Phi[s,t]
      Phioffdiag[s,t] = 0
      }
      else{
        Phidiag[s,t] = 0
        Phioffdiag[s,t] = Phi[s,t]
      }
    }
  }
  
  
  
  #Annealed Importance Sampling
  nsteps = 100
  nparticles = 500
  Phihistory = array(dim = c(nsteps+1, p, p))
  Thetahistory = array(0,dim = c(nsteps+1,p))
  Phihistory[1,,] = Phidiag
  for(i in 1:nsteps){
    Phihistory[i+1,,] = Phidiag + (i / nsteps) * Phioffdiag
    Thetahistory[i+1,] = (i / nsteps) * Theta
  }
  
  particlehistory = array(dim = c((nsteps+1),nparticles, p))
  if(family == "Poisson"){
    #start by generating n particles of dimension p from p independent exponential families
    particlehistory[1,,] = initialproposal(initializingparameters = as.list(diag(Phi)), family = family, nparticles = nparticles)
    for(i in 2:(nsteps+1)){
      print(i)
      #propose n moves of dimension p from p independent exponential families
     proposedmoves = t(mcmapply(FUN = newproposal,
                              initializingparameters = split(particlehistory[i-1,,], row(particlehistory[i-1,,])),
                              MoreArgs = list(family= family)))
     #Calculate the probabilities of these n moves
     moveprobabilities        = mcmapply(FUN = calculatemoveprobability,
                                 proposedmove = split(proposedmoves, row(proposedmoves)),
                                 currentposition = split(particlehistory[i-1,,], row(particlehistory[i-1,,])),
                                 MoreArgs = list(Phii = Phihistory[i,,], 
                                                 Thetai = Thetahistory[i,],
                                                 sufficientstatistic = sufficientstatistic,
                                                 family = "Poisson", 
                                                 basemeasure = basemeasure,
                                                 problaw = PproportionalBM))

      #update particle paths
      particlehistory[i,,] = particlehistory[i-1,,]
      acceptancethreshold = runif(n = nparticles, min = 0, max = 1)
      particlehistory[i,moveprobabilities > acceptancethreshold,] = proposedmoves[moveprobabilities > acceptancethreshold,]
    }
    
    #compute the weights
    #only need P proportional since base measures cancel out
    #what about A
    
    #paths that are correlated have high weight
    particlechainweight = calculateweights(problaw = Pproportional,
                                           particlehistory = particlehistory,
                                          sufficientstatistic = sufficientstatistic,
                                           Phihistory = Phihistory,
                                           Thetahistory = Thetahistory)
    }
    
  }



# P = function(x, Phi, Theta, sufficientstatistic, basemeasure, family){
#   sum(Theta * sqrt(sufficientstatistic(x))) + sqrt(sufficientstatistic(x)) * Phi * sqrt(sufficientstatistic(x)) 
#   + sum(basemeasure(x)) - Atotal(Phi, Theta, sufficientstatistic, basemeasure, family)
# }


