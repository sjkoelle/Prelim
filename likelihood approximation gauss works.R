approximatePartition = function(theta = theta, Phi, nSamplesPerIter, nAnnealing, nSamples){

  p = dim(Phi)[1]
  diagPhi = diag(diag(Phi))
  offPhi = Phi - diagPhi
  #Y = array(0, nSamples)
  #logwVec = zeros(nSamples,1);
  asdf = mcmapply(FUN = annealingParticle, mc.cores = 4, index = 1:nSamples, MoreArgs = list(theta = theta, 
                                                                                             diagPhi = diagPhi, 
                                                                                             offPhi = offPhi,
                                                                                             nSamplesPerIter= nSamplesPerIter,
                                                                                             nAnnealing = nAnnealing))
  
  Aexp = sum(-log(-diag(Phi)))
  logwvec = array(dim = nSamples)
  for(i in 1:nSamples){
     logwvec[i] =asdf[,i][[2]]
  }
  A = log(mean(exp(logwvec - max(logwvec)))) + max(logwvec) + Aexp
  return(A)
}

annealingParticle = function(theta = theta, diagPhi = diagPhi, offPhi = offPhi, index, nAnnealing = nAnneadling,
                               nSamplesPerIter = nSamplesPerIter){
     logw = 0
      cur = array(0, p)
      distSeq = seq(from = 0, to = 1, by = 1/nAnnealing)
      for(i in 1:(nAnnealing-1)){
        curTheta = distSeq[i] * theta
        curPhi = distSeq[i] * offPhi + diagPhi
        nextTheta = distSeq[i+1] * theta
        nextPhi = distSeq[i+1] * offPhi + diagPhi
        if(i == 1){
          cur = mapply(rexp, rate = 1/-diag(curPhi), MoreArgs = list(n= 1))
        } else {
          cur = mcmcsample(curTheta, curPhi, nSamplesPerIter = nSamplesPerIter, cur = cur);
        }
        
        sqrCur = sqrt(cur);
        logw = logw + (nextTheta-curTheta) %*% sqrCur + sqrCur %*% (nextPhi - curPhi) %*% sqrCur
      
    }
    return(list(sqrCur, logw))
  }


#mcmc sample is for use with annealed importance sampling
#curTheta and curPhi refer to annealing intermediates
#cur is the initialization point of the sampling, which is relevant since we typically don't run many steps
mcmcsample = function(curTheta, curPhi, nSamplesPerIter, cur){
  
  sqrCur = sqrt(cur)
  for(j in 1:nSamplesPerIter){
    print(j)
    proposal = -cur * log(runif(min = 0, max = 1, n = length(cur)))
  
    sqrp = sqrt(proposal)
    Pprop = exp(curTheta %*% sqrp + (sqrp) %*% curPhi %*% sqrp)
    Pcur = exp(curTheta %*% sqrCur + (sqrCur) %*% curPhi %*% sqrCur)
    
    Qcur = prod(1/proposal * exp(-cur/proposal))
    Qprop = prod(1/cur * exp(-proposal/cur))
    
    acceptProb = Pprop/Pcur * Qcur/Qprop
    print(acceptProb)
    if(!is.na(acceptProb)){
    if(acceptProb > runif(1)){
      cur = proposal
      sqrCur = sqrp
    }
    }
    #Y[,j] = cur
  }
  
  return(cur)
}

#switch in symmetric Gaussian proposal density with var = 1
mcmcsampleGaussian = function(curTheta, curPhi, nSamplesPerIter, cur){
  
  sqrCur = cur
  recordpos = array(dim = c(length(curTheta), nSamplesPerIter))
  if(nSamplesPerIter > 0){
    for(j in 1:nSamplesPerIter){
      print(j)
      proposal = mapply(FUN = rnorm, mean = sqrCur, MoreArgs = list(n = 1, sd = 1))
      #-cur * log(runif(min = 0, max = 1, n = length(cur)))
      
      #sqrp = sqrt(proposal)
      sqrp = proposal
      Pprop = exp(curTheta %*% sqrp + (sqrp) %*% curPhi %*% sqrp)
      Pcur = exp(curTheta %*% sqrCur + (sqrCur) %*% curPhi %*% sqrCur)
      
  #     Qcur = prod(1/proposal * exp(-cur/proposal))
  #     Qprop = prod(1/cur * exp(-proposal/cur))
  #     
      #acceptProb = Pprop/Pcur * Qcur/Qprop
      acceptProb = Pprop/Pcur
      
      print(acceptProb)
      if(!is.na(acceptProb)){
        if(acceptProb > runif(1)){
          #ur = proposal
          sqrCur = sqrp
        }
      }
      #Y[,j] = cur
      #recordpos[,j] = sqrCur
    }
    return(sqrCur)
  } else {
    return(sqrCur)
  }
}

#(curTheta %*% recordpos) * recordpos %*% (curPhi %*% recordpos)

#smoothScatter(recordpos[1,1000:10000], recordpos[2,1000:10000])



annealingParticleGaussian = function(theta = Theta, diagPhi = diag(diag(Phi)), offPhi = Phi - (diag(diag(Phi))),
                                     index, nAnnealing = nAnneadling,
                             nSamplesPerIter = nSamplesPerIter){
  #initialize weight
  logw = 0
  #initialize particle to be annealed
  cur = array(0, length(theta))
  #initialize mixing factors
  distSeq = seq(from = 0, to = 1, by = 1/nAnnealing)
  #at each annealing step, move the particle, and compute weight
  for(i in 1:(nAnnealing)){
    curTheta = distSeq[i] %*% theta
    curPhi = distSeq[i] * offPhi + diagPhi
    nextTheta = distSeq[i+1] * theta
    nextPhi = distSeq[i+1] * offPhi + diagPhi
    if(i == 1){
      #cur = mapply(rnorm, mean = solve((-2*(diagPhi+offPhi))) %*% curTheta , MoreArgs = list(sd = 1, n= 1))
      cur = mapply(rnorm, mean = solve((-2*(diagPhi+offPhi))) %*% t(curTheta) , 
                   sd = sqrt(diag(solve(-2*curPhi))), MoreArgs = list( n= 1))
    } else {
      cur = mcmcsampleGaussian(curTheta, curPhi, nSamplesPerIter = nSamplesPerIter, cur = cur)
    }
    
    sqrCur = cur
    logw = logw + (nextTheta-curTheta) %*% sqrCur + sqrCur %*% (nextPhi - curPhi) %*% sqrCur
    
  }
  return(list(sqrCur, logw))
}

approximatePartitionGaussian = function(theta = theta, Phi, nSamplesPerIter, nAnnealing, nSamples){
  
  p = dim(Phi)[1]
  diagPhi = diag(diag(Phi))
  offPhi = Phi - diagPhi
  #Y = array(0, nSamples)
  #logwVec = zeros(nSamples,1);
  asdf = mcmapply(FUN = annealingParticleGaussian, mc.cores = 4, index = 1:nSamples, MoreArgs = list(theta = theta, 
                                                                                             diagPhi = diagPhi, 
                                                                                             offPhi = offPhi,
                                                                                             nSamplesPerIter= nSamplesPerIter,
                                                                                             nAnnealing = nAnnealing))
  
  Aexp = sum(log(sqrt(diag(solve(-2*diagPhi)))))
  logwvec = array(dim = nSamples)
  for(i in 1:nSamples){
    logwvec[i] = asdf[,i][[2]]
  }
 # A = log(mean(exp(logwvec - max(logwvec)))) + max(logwvec) + Aexp
  A = log(mean(exp(logwvec - max(logwvec)))) + max(logwvec) + Aexp
  return(A)
}
