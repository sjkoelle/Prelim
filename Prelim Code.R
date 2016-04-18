####Stat 572 Prelim Exam
#Sam Koelle
#UW Stat Department

setwd("/Users/samsonkoelle/Desktop/Prelim")
rawdata = read.csv("On_Time_On_Time_Performance_2014_1.csv")

averagedelaytimes = array(dim = c(30, 365))

get30busiestairports = function(rawdata){
  flightsperairport = table(rawdata$OriginAirportID)
  top30 = order(-table(rawdata$OriginAirportID))[1:30]
  top30IDs = names(flightsperairport)[top30]
#   flightsfromtop30 = which(rawdata$OriginAirportID %in% top30IDs)
#   top30data = rawdata[flightsfromtop30,]
}

for(i in 1:30){
  print(i)
  airporti = rawdata[which(rawdata$OriginAirportID == top30IDs[i]),]
  datevec = as.numeric(strftime(as.character(airporti$FlightDate), format = "%j"))
  for(k in 1:365){
    print(k)
    airportidayk = which(datevec == k)
    averagedelaytimes[i,k] = mean(airporti[airportidayk,]$DepDelay, na.rm = TRUE)
  }
}



#p l1-regularized node-wise regressions
# the purpose of these regressions is to minimize the penalized square root convex objective function
# the partition function is part of the objective function, and we are going to approximate it
#they use the quadgk matlab function, which hopefully is somewhat similar to the quadgk R function
# Anode  has two parameters, nu1 and nu2, so by slightly varying nu1 and nu2, we can get the gradient