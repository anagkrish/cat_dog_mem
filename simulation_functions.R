library(ctmm)
library(adehabitatHR)
library(raster)
library(terra)
library(sp)
library(sf)
library(grDevices)

#--------------------------------------------------------------------
# Function for estimating ridges from telemetry object
#--------------------------------------------------------------------

ridges <- function(telemetry, CTMM = model, threshold = 0.5){
  
  #-------------------------------
  #Estimate the home range
  #-------------------------------

  FIT <- ctmm.fit(telemetry, model)
  UD <- akde(telemetry, CTMM = FIT)
  area_ud <- summary(UD, units=F)$CI[1,]
  
  #-------------------------------
  #Estimate ridge density on HR
  #-------------------------------
  
  #Estimate ridges
  RIDGE <- ctmm:::ridges.UD(UD)
  ridge_indicator <- RIDGE$Indicator
  
  threshold <- 0.5
  
  #get UDs at diff thresholds
  UD_High <- SpatialPolygonsDataFrame.UD(UD)
  UD_High <- sf::st_make_valid(st_as_sf(UD_High)[3,])
  
  UD_Mean <- SpatialPolygonsDataFrame.UD(UD)
  UD_Mean <- sf::st_make_valid(st_as_sf(UD_Mean)[2,])
  
  UD_Low <- SpatialPolygonsDataFrame.UD(UD)
  UD_Low <- sf::st_make_valid(st_as_sf(UD_Low)[1,])
  
  threshold=0.5
  #get ridges as contourLines
  ridges <- grDevices::contourLines(UD$r$x, UD$r$y, ridge_indicator, level=threshold)
  lines <- list()
  for (i in seq_along(ridges)) {
    lines[[i]] <- Lines(list(Line(cbind(ridges[[i]]$x, ridges[[i]]$y))), 
                        ID = as.character(i))
  }
  
  #convert to spatiallines
  all_lines <- SpatialLines(lines)
  all_lines <- sf::st_make_valid(st_as_sf(all_lines))
  
  #get all ridge lengths and densities
  ridges_high <- tryCatch({ terra::intersect(terra::vect(all_lines), terra::vect(UD_High)) },
                          error= function(e) { ridges_high <- 0 })
  #convert to sf bc rgeos is being retired
  ridges_high <- tryCatch( { sf::st_as_sf(ridges_high) },
                           error=function(e) { ridges_high <- 0 })
  length_h <-  tryCatch({sum(sf::st_length(ridges_high)) },
                        error= function(e) { length_h <- 0 })
  dens_h <- length_h^2/area_ud[[3]]
  
  ridges_mean <- tryCatch({ terra::intersect(terra::vect(all_lines), terra::vect(UD_Mean)) },
                          error= function(e) { ridges_mean <- 0 })
  #convert to sf bc rgeos is being retired
  ridges_mean <- tryCatch({ sf::st_as_sf(ridges_mean) },
                          error= function(e) { ridges_mean <- 0 })
  length_m <-  tryCatch({sum(sf::st_length(ridges_mean)) },
                        error= function(e) { length_m <- 0 })
  dens_m <- length_m^2/area_ud[[2]]
  
  ridges_low <- tryCatch({ terra::intersect(terra::vect(all_lines), terra::vect(UD_Low)) },
                         error= function(e) { ridges_low <- 0 })
  #convert to sf bc rgeos is being retired
  ridges_low <- tryCatch({ sf::st_as_sf(ridges_low) },
                         error= function(e) { ridges_low <- 0 })
  length_l <-  tryCatch({sum(sf::st_length(ridges_low)) },
                        error= function(e) { length_l <- 0 })
  dens_l <- length_l^2/area_ud[[1]]
  
  #-------------------------------
  #process simulated tracks with adehabitatHR
  #-------------------------------
  
  #needed to allow adehabitatHR to read simulated tracks as ltraj obj
  t_date <- seq(from=as.POSIXct("2000-1-1 0:00"), by="hour", length.out = length(telemetry$t))
  
  dat <- telemetry %>%
    data.frame() %>%
    mutate(t=t_date)
  
  #convert to ltraj
  DATA <- as.ltraj(cbind(dat$x, dat$y),
                   dat$t, id="sim")
  
  t=720
  l=0
  h=10 #units are m
  
  #estimate Diffusion coefficient
  vv <- BRB.D(DATA,
              Tmax=t*60, #12 hour upper threshold (depends on sampling frequency and should be biologically accurate)
              Lmin=l) #unsure of biological grounds of lmin with simulated tracks
  
  #estimate brb UD
  UD.BB <- BRB(DATA, D = vv, Tmax = 720*60, Lmin = 0, hmin=h)
  #use acast function to turn dataframe into x by y grid similar to PDF from ctmm
  grid <- acast(cbind(UD.BB@coords, UD.BB@data), x~y, value.var='dens') 
  
  #-------------------------------
  #create an empty UD object that can function as a ctmm ud to run through the ridge function
  #-------------------------------

  UD.BB.RIDGE <- list()
  UD.BB.RIDGE$PDF <- grid #grid of densities is analogous to PDF in ctmm
  #hmin is analogous to bandwith in ctmm
  UD.BB.RIDGE$dr <- UD.BB@grid@cellsize
  #print("got bb") 
  
  #-------------------------------
  #Estimate ridge density on HR
  #------------------------------
  
  RIDGE_BB <- ctmm:::ridges.UD(UD.BB.RIDGE)
  ridge_indicator_bb <- RIDGE_BB$Indicator
  
  threshold <- 0.5
  
  #with added error messages to catch misbehaving cases
  ridges <- grDevices::contourLines(as.numeric(dimnames(grid)[[1]]),
                                    as.numeric(dimnames(grid)[[2]]), ridge_indicator_bb, level=threshold)
  lines <- list()
  for (i in seq_along(ridges)) {
    lines[[i]] <- Lines(list(Line(cbind(ridges[[i]]$x, ridges[[i]]$y))),
                        ID = as.character(i))
  }
  
  #convert ridges to to spatiallines
  crs <- telemetry@info$projection
  all_lines <- SpatialLines(lines, proj = CRS(crs))
  all_lines <- sf::st_as_sf(all_lines)

  #get BRB UD as sp and convert to sf
  SP.BB <- getverticeshr(UD.BB, percent = 95)
  SF.BB <- sf::st_as_sf(SP.BB)
  area_bb <- st_area(SF.BB)
  
  ridges_mean <- tryCatch({sf::st_intersection(all_lines, SF.BB) },
                          error= function(e) { ridges_mean <- 0 })
  ##convert to sf bc rgeos is being retired
  ridges_mean <- tryCatch( { sf::st_as_sf(ridges_mean) },
                           error=function(e) { ridges_mean <- 0 })
  bb_length_m <-  tryCatch({sum(sf::st_length(ridges_mean)) },
                           error= function(e) { bb_length_m <- 0 })
  bb_dens_m <- bb_length_m^2/area_bb
  
  
  #-------------------------------
  #Compile and return results
  #-------------------------------
  
  results <- c(area_ud[[1]],
               area_ud[[2]],
               area_ud[[3]],
               length_l,
               length_m,
               length_h,
               dens_l,
               dens_m,
               dens_h,
               area_bb,
               bb_length_m,
               bb_dens_m)
  
  return(results)
  
  
} #Closes the function

#test <- ridges(telemetry)
