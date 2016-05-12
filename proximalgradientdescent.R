#take gradient step, then soft threshold (phi off-diagonal)
prox = function(phiinit = Phi, thetainit = Theta, phigradient = phigradient, thetagradient = thetagradient, stepsize = 0.001,
                lambda = lambda){
  Phinew = Phi - stepsize * phigradient
  Thetanew = Theta - stepsize * thetagradient
  Phioffdiag = Phi
  for(s in 1:p){
    for(t in 1:p){
      if(s == t){
        Phinew[s,t] = Phinew[s,t]
      } else {
        Phinew[s,t] = soft.threshold(Phinew[s,t], lambda)
      }
    }
  }
  output = list(Thetanew, Phinew)
}