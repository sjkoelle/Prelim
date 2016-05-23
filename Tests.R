#simulate multivariate Gaussian random variables using SQGM
grid = seq(from = -4, to = 4, by = 1)

#Recall, Theta = sigma^-1 * mu, Phi = -1/2 sigma^-1
#Phi must be negative definite!
Phi = -1/2
Theta = 0 
prob = array(dim = c(length(grid)))
for(i in 1:length(grid)){
  prob[i] = PproportionalBM(x= grid[i], Theta = Theta, Phi = Phi, family = "Gaussian",
                            sufficientstatistic = sufficientstatistic, basemeasure = basemeasure)
}

#bivariate
#sigma (1,.5,.5,1) gives
Phi = matrix(c(-2/3,1/3,1/3,-2/3), nrow = 2)
Theta = c(0,0)
prob = array(dim = c(length(grid),length(grid)))
prob2 = prob
for(i in 1:length(grid)){
  for(j in 1:length(grid)){
    prob[i,j] = PproportionalBM(x= c(grid[i], grid[j]), Theta = Theta, Phi = Phi, family = "Gaussian",
                                sufficientstatistic = sufficientstatistic, basemeasure = basemeasure)
    
    prob2[i,j] = dmvnorm(c(grid[i], grid[j]), mean = rep(0, 2), sigma = matrix(c(1,1/2,1/2,1), nrow = 2), log = FALSE)
  }
}

prob/prob2
plot_ly(z = prob2, type = "surface")

#bivariate
sigma = matrix(c(10,-2,-2,10), nrow = 2)
Phi = -.5*solve(sigma)
Theta = c(0,0)
prob = array(dim = c(length(grid),length(grid)))
prob2 = prob
for(i in 1:length(grid)){
  for(j in 1:length(grid)){
    prob[i,j] = PproportionalBM(x= c(grid[i], grid[j]), Theta = Theta, Phi = Phi, family = "Gaussian",
                                sufficientstatistic = sufficientstatistic, basemeasure = basemeasure)
    
    prob2[i,j] = dmvnorm(c(grid[i], grid[j]), mean = rep(0, 2), sigma = sigma, log = FALSE)
  }
}

prob/prob2

##########





prob = prob * Atotal(Theta = Theta, Phi = Phi, family = "Poisson")



#simulate correlated exponential random variables using SQGM
expgrid = seq(from = 0, to = 0.4, by = 0.01)

prob = array(dim = c(length(expgrid),length(expgrid)))
for(i in 1:length(expgrid)){
  for(j in 1:length(expgrid)){
    prob[i,j] = Pprop(c(i,j), Theta, Phi, sufficientstatistic, basemeasure, family = "exponential")
    prob[i,j] = prob * exp(basemeasure(c(i,j), family = "exponential"))
  }
}
prob = prob * Atotal(Theta = Theta, Phi = Phi, basemeasure = basemeasure, family = "exponential")

#simulate correlated poisson random variables using SQGM
poisgrid = seq(from = 0, to = 9, by = 1)
bm = basemeasure(family = "Poisson")

#univariate
Phi = 1
Theta = 0 
prob = array(dim = c(length(poisgrid)))
for(i in 1:length(poisgrid)){
  prob[i] = Pproportional(x= i, Theta = Theta, Phi = Phi, sufficientstatistic = sufficientstatistic)
  prob[i] = prob[i] * exp(sum(bm(i)))
}

#bivariate
Phi = matrix(c(1,1,1,1), nrow = 2)
Theta = c(0,0)
prob = array(dim = c(length(poisgrid),length(poisgrid)))
for(i in 1:length(poisgrid)){
  for(j in 1:length(poisgrid)){
    prob[i,j] = Pproportional(x= c(i,j), Theta = Theta, Phi = Phi, sufficientstatistic = sufficientstatistic)
    prob[i,j] = prob[i,j] * exp(sum(bm(c(i,j))))
  }
}
prob = prob * Atotal(Theta = Theta, Phi = Phi, family = "Poisson")


