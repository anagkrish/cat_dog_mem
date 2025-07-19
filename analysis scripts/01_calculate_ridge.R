#Calculate ridges on an individual movement track using CTMM package

# load libraries
library(ctmm)
library(sp)
library(sf)
library(grDevices)

#movement data is any csv file
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

#estimate ridge density
#ridge object is a list with 1) average width and 2) indicator matrix
ridge_indicator <- RIDGE$Indicator #original ridge matrix 

SP.UD <- SpatialPolygonsDataFrame.UD(UD)
threshold <- 0.5 #we tried out different thresholds and they did not influence the clade difference estimate

#get UDs at diff thresholds
UD_High <- SpatialPolygonsDataFrame.UD(UD)
UD_High <- tryCatch( { sf::st_as_sf(UD_High) },
                     error=function(e) { UD_High <- 0 })
UD_High <- tryCatch( { sf::st_make_valid(UD_High$geometry[3]) },
                     error=function(e) { UD_High <- 0 })

UD_Mean <- SpatialPolygonsDataFrame.UD(UD)
UD_Mean <- tryCatch( { sf::st_as_sf(UD_Mean) },
                     error=function(e) { UD_Mean <- 0 })
UD_Mean <- tryCatch( { sf::st_make_valid(UD_Mean$geometry[2]) },
                     error=function(e) { UD_Mean <- 0 })

UD_Low <- SpatialPolygonsDataFrame.UD(UD)
UD_Low <- tryCatch( { sf::st_as_sf(UD_Low) },
                    error=function(e) { UD_Low <- 0 })
UD_Low <- tryCatch( { sf::st_make_valid(UD_Low$geometry[1]) },
                    error=function(e) { UD_Low <- 0 })

#get ridges as contourLines
ridges <- grDevices::contourLines(UD$r$x, UD$r$y, ridge_indicator, level=threshold)
lines <- list()
for (i in seq_along(ridges)) {
  lines[[i]] <- Lines(list(Line(cbind(ridges[[i]]$x, ridges[[i]]$y))), 
                      ID = as.character(i))
}

#convert to spatiallines
all_lines <- SpatialLines(lines, proj = SP.UD@proj4string)
all_lines <- sf::st_as_sf(all_lines)

#get all ridge lengths and densities
ridges_high <- tryCatch({sf::st_intersection(all_lines, UD_High) },
                        error= function(e) { ridges_high <- 0 })
#convert to sf bc rgeos is being retired
ridges_high <- tryCatch( { sf::st_as_sf(ridges_high) },
                         error=function(e) { ridges_high <- 0 })
length_h <-  tryCatch({sum(sf::st_length(ridges_high)) },
                      error= function(e) { length_h <- 0 })
dens_h <- length_h^2/area_ud[[3]]

ridges_mean <- tryCatch({sf::st_intersection(all_lines, UD_Mean) },
                        error= function(e) { ridges_high <- 0 })
#convert to sf bc rgeos is being retired
ridges_mean <- tryCatch( { sf::st_as_sf(ridges_mean) },
                         error=function(e) { ridges_mean <- 0 })
length_m <-  tryCatch({sum(sf::st_length(ridges_mean)) },
                      error= function(e) { length_m <- 0 })
dens_m <- length_m^2/area_ud[[2]]

ridges_low <- tryCatch({sf::st_intersection(all_lines, UD_Low) },
                       error= function(e) { ridges_low <- 0 })
#convert to sf bc rgeos is being retired
ridges_low <- tryCatch({sf::st_as_sf(ridges_low) },
                       error= function(e) { ridges_low <- 0 })
length_l <-  tryCatch({sum(sf::st_length(ridges_low)) },
                      error= function(e) { length_l <- 0 })
dens_l <- length_l^2/area_ud[[1]]

