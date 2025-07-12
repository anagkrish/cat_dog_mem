#Calculate ridges on an individual movement track using CTMM package

# load libraries
library(ctmm)
library(sp)
library(sf)
library(grDevices)

track <- "load/movement/data/here"

telemetry <- as.telemetry(track)
SVF <- variogram(track)
GUESS <- ctmm.guess(track,variogram=SVF,interactive=F)

#fit model
FIT <- ctmm.select(telemetry, GUESS, 
                        level = 0.95) #non-default argument used to speed up model fitting for 1500+ individuals (made some fits unstable, see methods for more info))
UD <- akde(telemetry, CTMM = FIT)

#Estimate ridges
RIDGE <- ctmm:::ridges.UD(UD)
ridge_indicator <- RIDGE$Indicator
