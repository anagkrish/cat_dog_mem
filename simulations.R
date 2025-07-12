#--------------------------------------------------------------------
# Simulations exploring the behaviour of ridge estimation
# script by Dr. Michael Noonan
#--------------------------------------------------------------------


#Load in the requisite packages into the root environment
library(ctmm)
library(raster)
library(terra)
library(sp)
library(sf)
library(grDevices)
library(foreach)
library(doParallel)

#Source in functions from file
source("simulation_functions.R")

#--------------------------------------------------------------------
# Set up the cores for the parallelisation
#--------------------------------------------------------------------

#Reg. multiple cores for doParallel
nCores <- 4

registerDoParallel(nCores)

#Check that it's setup correctly (should match what you've asked for)
getDoParWorkers()


#--------------------------------------------------------------------
# Simulations varying the sampling duration
#--------------------------------------------------------------------

#Specify an OUF model for simulation
#One day in seconds
ds <- 86400

#Spatial variance in m^2
sig <- 100000

#Positional autocorrelation timescale
tau_p <- ds

#Velocity autocorrelation timescale
tau_v <- ds/2

#sampling duration in days
nd <- c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)

#The number of locations per day
pd <- 24

#Specify the number of replicates per sampling frequency
nReps <- 100

#Create an empty list to store the results
res <- list()

#Loop over the sampling frequencies
for(i in 1:length(pd)){
  
  #For each sampling frequency, repeat  simulation n times
  x <- foreach(j=1:nReps,
               .combine='rbind',
               .packages = c('ctmm'),
               .export = c("ds",
                           "sig",
                           "tau_p",
                           "tau_v",
                           "nd",
                           "pd")) %dopar% {
                             
                             #Specify the movement model
                             model <- ctmm(tau=c(tau_p,tau_v),
                                           isotropic=TRUE,
                                           sigma=sig,
                                           mu=c(0,0))
                             
                             #sampling times
                             st <- 1:(nd[i]*pd)*(ds/pd)
                             
                             #Simulate some data from the pre-defined model
                             sim <- simulate(model,t=st)
                             
                             #Estimate ridges
                             results <- ridges(sim)
                             
                             #Output results
                             x <- c(nd[i], results)
                           } #Closes the parallelisation
  
  
  #Store results in a list
  res[[i]] <- x
  
  #Save the results as they come out
  save(res, file = "duration_sims.rda")
  
  #Print out an indicator of progress
  print(pd[i])
  
}

#Some data carpentry to get the results into a more usable format
res <- do.call(rbind, res)

res <- data.frame("duration" = res[,1],
                  "akde_area_low" = res[,2],
                  "akde_area_ML" = res[,3],
                  "akde_area_high" = res[,4],
                  "akde_ridge_length_low" = res[,5],
                  "akde_ridge_length_ML" = res[,6],
                  "akde_ridge_length_high" = res[,7],
                  "akde_ridge_density_low" = res[,8],
                  "akde_ridge_density_ML" = res[,9],
                  "akde_ridge_density_high" = res[,10],
                  "bb_area" = res[,11],
                  "bb_ridge_length" = res[,12],
                  "ridge_density" = res[,13])

write.csv(res, file = "duration_sims.csv")


#--------------------------------------------------------------------
# Simulations varying the sampling frequency
#--------------------------------------------------------------------

#Specify an OUF model for simulation
#One day in seconds
ds <- 86400

#Spatial variance in m^2
sig <- 100000

#Positional autocorrelation timescale
tau_p <- ds

#Velocity autocorrelation timescale
tau_v <- ds/2

#sampling duration in days
nd <- 100

#Create a vector specifying the number of locations per day
pd <- c(2, 4, 8, 16, 32, 64, 128, 256)

#Specify the number of replicates per sampling frequency
nReps <- 100

#Create an empty list to store the results
res <- list()

#Loop over the sampling frequencies
for(i in 1:length(pd)){
  
  #For each sampling frequency, repeat  simulation n times
  x <- foreach(j=1:nReps,
               .combine='rbind',
               .packages = c('ctmm'),
               .export = c("ds",
                           "sig",
                           "tau_p",
                           "tau_v",
                           "nd",
                           "pd")) %dopar% {
                             
                             #Specify the movement model
                             model <- ctmm(tau=c(tau_p,tau_v),
                                           isotropic=TRUE,
                                           sigma=sig,
                                           mu=c(0,0))
                             
                             #sampling times
                             st <- 1:(nd*pd[i])*(ds/pd[i])
                             
                             #Simulate some data from the pre-defined model
                             sim <- simulate(model,t=st)
                             
                             #Estimate ridges
                             results <- ridges(sim, model)
                             
                             #Output results
                             x <- c(pd[i], results)
                           } #Closes the parallelisation
  
  #Store results in a list
  res[[i]] <- x
  
  #Save the results as they come out
  save(res, file = "frequency_sims.rda")
  
  #Print out an indicator of progress
  print(pd[i])
  
}

#Some data carpentry to get the results into a more usable format
res <- do.call(rbind, res)

res <- data.frame("frequency" = res[,1],
                  "akde_area_low" = res[,2],
                  "akde_area_ML" = res[,3],
                  "akde_area_high" = res[,4],
                  "akde_ridge_length_low" = res[,5],
                  "akde_ridge_length_ML" = res[,6],
                  "akde_ridge_length_high" = res[,7],
                  "akde_ridge_density_low" = res[,8],
                  "akde_ridge_density_ML" = res[,9],
                  "akde_ridge_density_high" = res[,10],
                  "bb_area" = res[,11],
                  "bb_ridge_length" = res[,12],
                  "ridge_density" = res[,13])

write.csv(res, file = "frequency_sims.csv")


#--------------------------------------------------------------------
# Simulations varying the movement tortuosity
#--------------------------------------------------------------------

#Specify an OUF model for simulation
#One day in seconds
ds <- 86400

#Spatial variance in m^2
sig <- 100000

#Positional autocorrelation timescale
tau_p <- ds

#Velocity autocorrelation timescale
tau_v <- c(ds/2, ds/4, ds/8, ds/16, ds/32, ds/64, ds/128, ds/256)

#sampling duration in days
nd <- 100

#The number of locations per day
pd <- 24

#The number of replicates per sampling frequency
nReps <- 100

#Create an empty list to store the results
res <- list()

#Loop over the sampling frequencies
for(i in 1:length(tau_v)){
  
  #For each sampling frequency, repeat  simulation n times
  x <- foreach(j=1:nReps,
               .combine='rbind',
               .packages = c('ctmm'),
               .export = c("ds",
                           "sig",
                           "tau_p",
                           "tau_v",
                           "nd",
                           "pd")) %dopar% {
                             
                             #Specify the movement model
                             model <- ctmm(tau=c(tau_p,tau_v[i]),
                                           isotropic=TRUE,
                                           sigma=sig,
                                           mu=c(0,0))
                             
                             #sampling times
                             st <- 1:(nd*pd)*(ds/pd)
                             
                             #Simulate some data from the pre-defined model
                             sim <- simulate(model,t=st)
                             
                             #Estimate ridges
                             results <- ridges(sim)
                             
                             #Output results
                             x <- c(tau_v[i], results)
                           } #Closes the parallelisation
  
  
  #Store results in a list
  res[[i]] <- x
  
  #Save the results as they come out
  save(res, file = "tortuosity_sims.rda")
  
  #Print out an indicator of progress
  print(pd[i])
  
}

#Some data carpentry to get the results into a more usable format
res <- do.call(rbind, res)

res <- data.frame("tau_v" = res[,1],
                  "akde_area_low" = res[,2],
                  "akde_area_ML" = res[,3],
                  "akde_area_high" = res[,4],
                  "akde_ridge_length_low" = res[,5],
                  "akde_ridge_length_ML" = res[,6],
                  "akde_ridge_length_high" = res[,7],
                  "akde_ridge_density_low" = res[,8],
                  "akde_ridge_density_ML" = res[,9],
                  "akde_ridge_density_high" = res[,10],
                  "bb_area" = res[,11],
                  "bb_ridge_length" = res[,12],
                  "ridge_density" = res[,13])

write.csv(res, file = "tortuosity_sims.csv")

