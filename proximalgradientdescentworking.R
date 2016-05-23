#load data
setwd("/Users/samsonkoelle/Desktop/Prelim/New Code 051216")
rawdata = read.csv("airlinedata.csv", header = FALSE)
rawdatanames = read.csv("airport-geocodes.csv")
rawdata = as.matrix(rawdata)
Phi = matrix(0, nrow = 30, ncol = 30)
initPhi = -1 / apply(FUN = mean, X = rawdata, MARGIN = 2) 
diag(Phi) <- initPhi
Theta = array(0, dim = 30)

sortedcors = cbind(order(-Phi) %%30, ((order(-Phi) - (order(-Phi) %%30)) / 30) + 1)
sortedpcor = cbind(order(-cor(rawdata))%%30, ((order(-cor(rawdata))- (order(-cor(rawdata))%%30)) / 30)+1)
rawdatanames[[3]][]

rawdatanames[[3]][as.vector(sortedcors[1,])]
nodewiseProximalGradientDescent = function(X, Phi, Theta, index, maxIter = 1000, maximalPhi = -0.01, lambda= 1, maxBeta = 20){
  
  Ptestsave = array(dim = c(maxIter,maxBeta))
  Phicol = Phi[,index]
  eta1 = getEta1(X, theta = Theta, phi = Phicol, index)
  eta2 = getEta2(X, phi= Phicol, index=index)
  Pnodecurrent = Pnode(X, eta1 = eta1, eta2 = eta2, index = index) + lambda * sum(abs(Phicol[-index]))
  for(i in 1:maxIter){
    print(i)
    g = array(dim = 30)
    g[index] = dPnodeindexphi(X,eta1,eta2,index);
    #temp = mapply(FUN = dAnode2, eta1pt = eta2, eta2pt = eta1, MoreArgs = list(eps = 0.0001))
    temp = mapply(FUN = dAnode2fixed, eta1 = eta1, eta2 = eta2, MoreArgs = list(eps = 0.0001))
    g[-index] = dPnodeotherphi(X = X, dA2 = temp, index = index)
    dtheta = dPnodetheta(X, dA2 = temp, index = index)
    g = append(g,dtheta)
    
    for(beta in 0:maxBeta){
      print(beta)
      step = 0.25 ^ beta
      testPhi = Phicol - step*g[1:30]
      testTheta = Theta[index] - step*dtheta
      
      #Proximal projection
      testPhi[index] = min(testPhi[index],maximalPhi)
      testPhi[-index] = mapply(FUN = softThreshold, value = testPhi[-index], MoreArgs = list(lam = lambda))
      testEta1 = getEta1single(X, theta = testTheta,phi = testPhi,index = index)
      testEta2 = getEta2(X,testPhi,index)
      PnodeTest = Pnode(X,testEta1,testEta2,index) + lambda*sum(abs(testPhi[-index]))
      print(PnodeTest)
      Ptestsave[i,beta] = PnodeTest
      gDev = 1/step * append(Phicol - testPhi, Theta[index] - testTheta)
      if(!is.na(PnodeTest)){
       if(PnodeTest < Pnodecurrent && PnodeTest <= Pnodecurrent - step* g %*% gDev + step/2 * (gDev %*% gDev)){
        relDiff = (PnodeTest - Pnodecurrent) / abs(Pnodecurrent)
        #print(
        #fprintf('s = %d, iter = %d, lineIter = %d, cur f = %g, test f = %g, relDiff = %g in %g s\n', ...
        #        s, iter, beta+1, fCur, fTest, relDiff, toc(tStart));
        Pnodecurrent = PnodeTest
        Phicol = testPhi;
        Theta[index] = testTheta;
        eta1 = testEta1;
        eta2 = testEta2;
        break;
       }
      }
    }
    if(beta == maxBeta){
      break;
    }
  }
  
  return(list(Phicol, testTheta))
}

bbb = mcmapply(FUN = nodewiseProximalGradientDescent, index = 1:30, mc.cores = 4, MoreArgs = list(X = rawdata, Phi = Phi, 
                                                                                                Theta = Theta, 
                                                                                                maxBeta = 20, 
                                                                                                maxIter = 1000, 
                                                                                                lambda = .0000464))


Thetafinal = array(dim = 30)
Phifinal = array(dim = c(30,30))
for(i in 1:30){
  Phifinal[,i] = bbb[,i][[1]]
  Thetafinal[i] = bbb[,i][[2]]
}


diagPhifinal = diag(Phifinal)
nodiagPhifinal = Phifinal - diag(diagPhifinal)
nodiagPhifinalaverage = (nodiagPhifinal + t(nodiagPhifinal)) / 2
Phifinal = nodiagPhifinalaverage  + diag(diagPhifinal)



asdf = nodewiseProximalGradientDescent(X = rawdata, Phi = Phi, 
                                       Theta = Theta, 
                                       maxBeta = 20, 
                                       maxIter = 1000, lambda = .0000464, index = 1)       

asdf = nodewiseProximalGradientDescent(X = rawdata, Phi = Phi, 
                                       Theta = Theta, 
                                       maxBeta = 20, 
                                       maxIter = 1000, lambda = 1, index = 1) 

asdf = nodewiseProximalGradientDescent(X = rawdata, Phi = Phi, 
                                       Theta = Theta, 
                                       maxBeta = 20, 
                                       maxIter = 1000, lambda = .0001, index = 1) 

softThreshold = function(value=0, lam = 0){
  output = sign(value) * max((abs(value) - lam),0)
  return(output)
}

dPnodetheta = function(X, dA2 = temp, index = index){
  output = (-1/n)* (sum(sqrt(X[,index])) - sum(dA2))
  return(output)
}


#dfo
dPnodeotherphi = function(X, index, dA2 = temp){
  output = -1/n * (2 * t(sqrt(X[,-index])) %*% sqrt(X[,index]) - 2*t(sqrt(X[,-index])) %*% dA2)
  return(output)
}

dPnodethetafixed = function(X, dA2 = temp, index = index){
  output = (-1/n)* (sum(sqrt(X[,index])) - sum(dA1))
  return(output)
}


#dfo
dPnodeotherphifixed = function(X, index, dA2 = temp){
  output = -1/n * (2 * t(sqrt(X[,-index])) %*% sqrt(X[,index]) - 2*t(sqrt(X[,-index])) %*% dA1)
  return(output)
}

dAnode1fixed = function(eta1, eta2, eps){
  output = (Anode(eta1 = (eta1 + eps), eta2 = eta2) - Anode(eta1 = eta1, eta2 = eta2) )/ eps
  return(output)
}
dAnode2fixed = function(eta1, eta2, eps){
  output = (Anode(eta1 = (eta1), eta2 = (eta2+eps)) - Anode(eta1 = eta1, eta2 = eta2) )/ eps
  return(output)
}
Anodefixed = function(eta1, eta2){
  output = sqrt(pi) * eta1 * exp((eta1^2) / (-4*eta2)) / (2*(-eta2)^1.5) * (1 - erf(-eta1 /
                                                                                              (2*sqrt(-eta2)))) + 1/-eta2
  return(log(output))
}
dAnode1 = function(eta1pt, eta2pt, eps){
  output = (Anode(eta1pt = (eta1pt + eps), eta2pt = eta2pt) - Anode(eta1pt = eta1pt, eta2pt = eta2pt) )/ eps
  return(output)
}
dAnode2 = function(eta1pt, eta2pt, eps){
  output = (Anode(eta1pt = (eta1pt), eta2pt = (eta2pt+eps)) - Anode(eta1pt = eta1pt, eta2pt = eta2pt) )/ eps
  return(output)
}
Anode = function(eta1pt, eta2pt){
  output = sqrt(pi) * eta2pt * exp((eta2pt^2) / (-4*eta1pt)) / (2*(-eta1pt)^1.5) * (1 - erf(-eta2pt /
                                                                                              (2*sqrt(-eta1pt)))) + 1/-eta1pt
  return(log(output))
}
#this is "dfs"
#eps = 0.001
dPnodeindexphi = function(X, eta1, eta2, index){
  #eta1 is a list
  output =  -1/n * (sum(X[,index]) - sum(mapply(FUN = dAnode1, eta1pt = eta2, eta2pt = eta1, MoreArgs = list(eps = 0.001))))
return(output)
}
dPnodeindexphifixed = function(X, eta1, eta2, index){
  #eta1 is a list
  output =  -1/n * (sum(X[,index]) - sum(mapply(FUN = dAnode2fixed, eta1 = eta1, eta2 = eta2, MoreArgs = list(eps = 0.001))))
  return(output)
}



#why do we switch eta1 and eta2?
#this is "f"
Pnode = function(X,eta1,eta2,index){
  n = dim(X)[1]
  output = (-1/n)*(eta1 %*% (sqrt(X[,index])) + eta2 %*% X[,index] - sum(mapply(FUN = Anode,
                                                                                eta1pt = eta2, eta2pt = eta1)))
  return(output)
}

getEta2 = function(X, phi = Phicol, index = index){
  output = phi[index]
  return(rep(output,365))
}

getEta1 = function(X, theta = Theta, phi = Phicol, index){
  output = theta[index] + 2*t(phi[-index]) %*% t(sqrt(X[,-index]))
  return(output)
}

getEta1single = function(X, theta = Theta, phi = Phicol, index){
  output = theta + 2*t(phi[-index]) %*% t(sqrt(X[,-index]))
  return(output)
}

logPpt = function(Xpt, Theta = Theta, Phi = Phi, A){
  Xpt = as.numeric(Xpt)
  print(Xpt)
  output = (Theta %*% (sqrt(Xpt)) + (sqrt(Xpt)) %*% Phi %*% (sqrt(Xpt)) + 
            sum(mapply(FUN = basemeasure, Xpt, MoreArgs = list(family = "exponential"))) - A)
  return(output)
}

logP = function(X, Theta = Theta, Phi = Phi, A = A){
  output = sum(mapply(FUN = logPpt, Xpt = split(X, rownames(X)), MoreArgs = list(Phi = Phi, Theta = Theta, A = A)))
  return(output)
}




