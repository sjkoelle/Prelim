getdobjectivedeta = function(eta1, eta2, ss, sufficientstatistic, basemeasure, e = 0.0001){

  An = Anode(eta1Anode = eta1,eta2Anode =  eta2, sufficientstatistic, basemeasure, measure = "lebesgue",
                     family = "exponential")
  Aneta1 = Anode(eta1Anode = (eta1 + e),eta2Anode =   eta2, sufficientstatistic, basemeasure, measure = "lebesgue",
                         family = "exponential")
  Aneta2 = Anode(eta1Anode = eta1,eta2Anode =  (eta2 + e), sufficientstatistic, basemeasure, measure = "lebesgue",
                         family = "exponential")
  dAndeta1 = (Aneta1 - An) / e
  dAndeta2 = (Aneta2 - An) / e
  dObjectivedeta1 = -(sqrt(ss) - dAndeta1)
  dObjectivedeta2 = -(ss - dAndeta2)
  dObjectivedeta = c(dObjectivedeta1,  dObjectivedeta2)
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

#phi gradient is more complicated
#we apply the chain rule
#this outputs a direction matrix of size pxp
#we take the mean of the chain ruled directions
getdobjectivedphi = function(etagradient, sufficientstatistics){
  n = dim(etagradient[,,1])[1]
  p = dim(etagradient[,,1])[2]
  output = matrix(0, nrow= p, ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      output[j,] = output[j,] + (2*sqrt(sufficientstatistics[i,])*etagradient[i,,1])
      output[j,j] = output[j,j] - (2*sqrt(sufficientstatistics[i,j])*etagradient[i,j,1])
      output[j,j] = output[j,j] + etagradient[i,j,2]
    }
  }
  return(output)
}

# #phi gradient is more complicated
# #we apply the chain rule
# #this outputs a direction matrix of size pxp
# #we take the mean of the chain ruled directions
# getdobjectivedphi = function(etagradient, sufficientstatistics){
#   n = dim(etagradient[,,1])[1]
#   p = dim(etagradient[,,1])[2]
#   output = matrix(0, nrow= p, ncol = p)
#   
#     for(s in 1:p){
#       
#       for(t in 1:p){
#         if(s == t){
#           output[s,t] = mean(etagradient[,s,2])
#         } else{
#           output[s,t] = 0
#           for(i in 1:n){
#           for(j in 1:n){
#             for(v in 1:p){
#           output[s,t] = output[s,t] + etagradient[i,v,1] * deta1dphi[i,v,j,t]
#             }
#           }
#         }
#       }
#     }
#   }
#   
#         output[j,] = output[j,] + (2*sqrt(sufficientstatistics[i,])*etagradient[i,,1])
#         output[j,j] = output[j,j] - (2*sqrt(sufficientstatistics[i,j])*etagradient[[1]][i,j])
#         output[j,j] = output[j,j] + etagradient[[2]][i,j]
#         
#             
#     }
#   }
#   return(output)
# }

