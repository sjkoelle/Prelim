#eta is n x p
#this function takes in eta1, eta1, data, and basemeasure, and returns the gradient of Anode wrt eta1, eta2
#for all observations and all nodes
getEtaGradient = function(etas, sufficientstatistics, sufficientstatistic, basemeasure, e = 0.0001){
  eta1 = etas[[1]]
  eta2 = etas[[2]]
  n = dim(eta1)[1]
  p = dim(eta1)[2]
  An = array(dim = c(n,p))
  Aneta1 = array(dim = c(n,p))
  Aneta2 = array(dim = c(n,p))
  for(i in 1:n){
    for(s in 1:p){
      An[i,s] = Anode(eta1Anode = eta1[i,s],eta2Anode =  eta2[i,s], sufficientstatistic, basemeasure, measure = "lebesgue",
                      family = "exponential")
      Aneta1[i,s] = Anode(eta1Anode = (eta1[i,s] + e),eta2Anode =   eta2[i,s], sufficientstatistic, basemeasure, measure = "lebesgue",
      family = "exponential")
      Aneta2[i,s] = Anode(eta1Anode = eta1[i,s],eta2Anode =  (eta2[i,s] + e), sufficientstatistic, basemeasure, measure = "lebesgue",
                          family = "exponential")
    }
  }
  
  dAndeta1 = (Aneta1 - An) / e
  dAndeta2 = (Aneta2 - An) / e
  
  dPdeta1 = 1/2 * sufficientstatistics^(-1/2) - dAndeta1
  dPdeta2 = sufficientstatistics - dAndeta2
  
  #return gradient of negative log likelihood
  result = list(-dPdeta1, -dPdeta2)
  return(result)
}

#this function takes in our dAdetas and returns dAdtheta, which is the mean
#this outputs some direction vector of length p
getThetaGradient = function(etagradient, sufficientstatistics){
  dtheta = apply(FUN = mean, X = etagradient[[1]], MARGIN = 2)
  return(dtheta)
}

#phi gradient is more complicated
#we apply the chain rule
#this outputs a direction matrix of size pxp
#we take the mean of the chain ruled directions
getPhiGradient = function(etagradient, sufficientstatistics){
  n = dim(etagradient[[1]])[1]
  p = dim(etagradient[[1]])[2]
  output = matrix(0, nrow= p, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      output[j,] = output[j,] + (2*sqrt(sufficientstatistics[i,])*etagradient[[1]][i,])/n
      output[j,j] = output[j,j] - (2*sqrt(sufficientstatistics[i,j])*etagradient[[1]][i,j])/n
      output[j,j] = output[j,j] + etagradient[[2]][i,j]/n
    }
  }
  return(output)
}

getEtas = function(Phi, Theta, sufficientstatistics){
  n = dim(sufficientstatistics)[1]
  p = dim(sufficientstatistics)[2]
  eta1 = array(dim = c(n,p))
  eta2 = array(dim = c(n,p))
  for(i in 1:n){
    for(s in 1:p){
      eta2[i,s] = Phi[s,s]
      eta1[i,s] = Theta[s] + 2*Phi[s,-s] %*% sqrt(sufficientstatistics[i,-s])
    }
  }
  return(abind(eta1,eta2, along = 3))
}