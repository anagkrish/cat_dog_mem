#FIGURE 1 MAIN

library(tidyverse)
library(ctmm)
library(lubridate)
library(parsedate)
library(dplyr)
library(geodist)
library(raster)
library(ks)
library(geosphere)
library(foreach)
library(snow)
library(doSNOW)

load("movement/data/here")

#only read in individuals used for this figure
ridge <- read_csv("ridge.csv") %>%
  filter(`kept (y/n)` == "y") %>%
  filter(updatedstudy %in% c("Vanak", "Young", "Sekercioglu"))

#modified 7/12/25 to plot ridges w ks algorithm

#For panels A-D
################## example workflow for ridge density calculations
#puma concolor F70 from Young study
#ridge plots are now displayed with the slower, more accurate ks algorithm
#this algorithm was not feasible to run on all 1528 individuals, but can be used to more accurately visualize the ridges

#set projection because this is what we use when we display all individuals on terrain - this way orientations
#line up between example and terrain visualizations

track <- movementdata %>% filter(individual.local.identifier=="F70") %>% as.telemetry
GUESS <- ctmm.guess(track,interactive=F)
FIT <- ctmm.select(track, GUESS, trace = TRUE, level = 0.95, cores = 2)
UD <- akde(track,FIT,debias=FALSE,weights=FALSE,grad=TRUE) # used for ks defaults

#set up for ks algorithm
x <- data.frame(track[,c('x','y')])
H <- UD$H
w <- UD$weights * nrow(x)
gridsize <- dim(UD$PDF)
fhat <- ks::kde(x,H=H,w=w,gridsize=gridsize)
density.cutoff <- ks::contourLevels(fhat,cont=95)
threshold <- ks::contourProbs(fhat, cont=65)

xmin <- c(UD$r$x[1],UD$r$y[1])
xmax <- c(UD$r$x[gridsize[1]],UD$r$y[gridsize[2]])

# this takes forever but looks good
RIDGE_KS <- ks::kdr(x,H=H,w=w,fhat=fhat,density.cutoff=density.cutoff,
                    gridsize=gridsize,xmin=xmin,xmax=xmax,pre=FALSE)

segments <- RIDGE_KS[["end.points"]] %>%
  group_split(segment)
lines_ks <- list()

for (i in seq_along(segments)) {
  lines_ks[i] <- Lines(list(Line(cbind(segments[[i]]$x, segments[[i]]$y))),
                       ID = as.character(i))
}

#convert to spatiallines
all_lines_ks <- SpatialLines(lines_ks, proj = SpatialPolygonsDataFrame.UD(UD)@proj4string)
all_lines_ks <- sf::st_as_sf(all_lines_ks)

UD_Mean <- SpatialPolygonsDataFrame.UD(UD)
UD_Mean <- tryCatch( { sf::st_as_sf(UD_Mean) },
                     error=function(e) { UD_Mean <- 0 })
UD_Mean <- tryCatch( { sf::st_make_valid(UD_Mean$geometry[2]) },
                     error=function(e) { UD_Mean <- 0 })

ridges_mean_ks <- tryCatch({sf::st_intersection(all_lines_ks, UD_Mean) },
                           error= function(e) { ridges_mean_ks <- 0 })
#convert to sf bc rgeos is being retired
ridges_mean_ks <- tryCatch( { sf::st_as_sf(ridges_mean_ks) },
                            error=function(e) { ridges_mean_ks <- 0 })

proj <- str_replace(UD@info$projection, "units=m", "units=km")
ridges_km <- spTransform(as_Spatial(ridges_mean_ks),proj)

#plot ks calculated ridges over track
ctmm::plot(track, col="black", pch=19, cex=0.25, error=F,
           UD=UD, col.UD="#FF6969", cex.axis=2.5, cex.lab=2.5) #D: ridges w/ UD and GPS tracks as point
plot(ridges_km, col="#960000", lwd=5, add=T)


#For panels E-G
################# calculate and plot multiple ridges in same extent

inds <- (movementdata %>%
           filter(updatedstudy=="Vanak") #%>% (studies are Drouilly, Young and Vanak)
         #unite(id, sp, study, col="id2", sep="") %>% 
          #drouilly inds
         #filter(id2 %in% c("BBJackal_Forest", "BBJackal_Rain","Lucky"))
          #vanak inds
         #filter(id2 %in% c("Jackal 03 (Nanda)", #"Jackal 07 (Marker)", "Jackal 09 (Roma)", 
         # "Jungle cat 09 (Dagdu)", "Jungle cat 10 (Bheegi Billi)", "Jungle cat 12 (Cavity)", "Jungle cat 13 (Kamini)"))
          #young inds
         #filter(id2 %in% c("C023", "C026", "C028", "F108", "F157","F159", "F70"))
)$id


#modifications to tracks as necessary
mod_to_tracks <- function(track) {

  id <- unique(track$"individual.local.identifier")[[1]]
  species <- unique(track$"individual.taxon.canonical.name")[[1]]
  study <- unique(track$"study.id")[[1]]

  updatedstudy <- filter(updatedstudynames, ID==id, Species==species, Study==study)$`Updated Study`
  # print(paste(id, species, study, updatedstudy))

  ############# below are lines to edit problematic individuals as needed

  if (updatedstudy=="Drouilly") {
    if (id%in%c("BBJackal_Rain")) {
      track <- crop_range_res(track)
    }
    if (id=="Lucky") {
      track <- crop_range_res(track) %>%
        filter(location.lat > -32.845)
    }
  }

  if (updatedstudy=="Vanak") {
    if (id=="Jackal 02 (Zoom)") {
      track <- crop_range_res(track)
    }
    if(id%in%c("Jungle cat 08 (Sultan)", "Jungle Cat 17 (Momo)")) {
      track <- track %>%
        drop_na(location.lat) %>%
        crop_range_res()
    }
  }

  if (updatedstudy=="Young") {
    if (id%in%c("C037","B004","B009")) {
      track <- crop_range_res(track)
    }
  }

  return(track)
}

#minsampling
get_minsampling <- function(individual) {
  minsampling <- min(diff(individual$t))
  return(minsampling) #pretty sure units are in seconds
}

#resample inds with minsampling <60s
resample <- function(track) {
  
  if (get_minsampling(as.telemetry(track)) < 60) {
    track_rs <- track %>%
      mutate(timestamp_rs = sapply(timestamp, cut,
                                   breaks = "1 min")) %>% #special feature in lubrdiate pkg!
      #manually search for any cases in which the 00:00:00 has been dropped
      mutate(timestamp_rs = ifelse(grepl("[0-9]{4}-[0-9]{2}-[0-9]{2}$",timestamp_rs),
                                   paste(as.character(timestamp_rs), "00:00:00"),
                                   as.character(timestamp_rs)),
             #timestamp_rs = ymd_hms(timestamp_rs)
      ) %>%
      #collapse to individual by resampled time points but save og time points
      distinct(timestamp_rs, .keep_all = "TRUE") %>%
      dplyr::select(-c("timestamp")) %>%
      rename("timestamp" = timestamp_rs) #added to resample sketchy individuals
    
    return(track_rs)
  }
  
}

#crop non-range resident inds
crop_range_res <- function(track) {

  median <- ctmm::median(as.telemetry(track))
  dists <- track %>%
    mutate(dists=NA)

  #calculate distances from each point to median
  for (i in seq_along(dists$location.long)) {

    dists$dists[i] <- distm(c(dists$location.long[[i]], dists$location.lat[[i]]),
                            c(median$longitude[1], median$latitude[1]), fun = distHaversine)
  }


  #calculate median absolute deviation
  mad <- mad(dists$dists)

  #pull out all points 5 mad thresholds away from centroid; idk if this is the best threshold, can change if need be
  subset <- subset(dists, dists < (mad*5))
  return(subset)

}

mod_to_tracks <- function(track) { #edit non-range res individuals
  
  id <- unique(track$"individual.local.identifier")[[1]]
  species <- unique(track$"individual.taxon.canonical.name")[[1]]
  study <- unique(track$"study.id")[[1]]
  
  updatedstudy <- filter(updatedstudynames, ID==id, Species==species, Study==study)$`Updated Study`
  # print(paste(id, species, study, updatedstudy))
  
  ############# below are lines to edit problematic individuals as needed
  
  if (updatedstudy=="Drouilly") {
    if (id%in%c("BBJackal_Rain")) {
      track <- crop_range_res(track)
    }
    if (id=="Lucky") {
      track <- crop_range_res(track) %>%
        filter(location.lat > -32.845)
    }
  }
  
  if (updatedstudy=="Vanak") {
    if (id=="Jackal 02 (Zoom)") {
      track <- crop_range_res(track)
    }
    if(id%in%c("Fox 25 (Chandra)","Jungle cat 08 (Sultan)", "Jungle Cat 17 (Momo)")) {
      track <- track %>%
        drop_na(location.lat) %>%
        crop_range_res()
    }
  }
  
  if (updatedstudy=="Young") {
    if (id%in%c("C037","B004","B009")) {
      track <- crop_range_res(track)
    }
  }
  
  return(track)
}

#####

data <- alldata %>%
  filter(individual.local.identifier%in%c("selected/individuals/here"))

data <- split(data, data$individual.local.identifier)
DATA <- list()
DATA <- list()

for (i in seq_along(data)) {
  
  DATA[[i]] <- data[[i]] %>%
    mod_to_tracks() %>%
    resample() %>%
    as.telemetry()
  
}

FITS <- list()
UDS <- list()

for (i in seq_along(DATA)) {
  
  ctmm::projection(DATA[[i]]) <- median(DATA)
  GUESS <- ctmm.guess(DATA[[i]],interactive=F)
  FITS[[i]] <- ctmm.select(DATA[[i]], GUESS, trace = TRUE, level = 0.95, cores = 2)
  UDS[[i]] <- akde(DATA[[i]],FITS[[i]],debias=FALSE,weights=FALSE,grad=TRUE) # used for ks defaults
  
}

#get all UDS to plot together
SP.UDS <- SpatialPolygonsDataFrame.UD(UDS[[1]])
for (i in seq_along(UDS)) {
  i=i+1
  SP <- SpatialPolygonsDataFrame.UD(UDS[[i]])
  SP.UDS <- rbind(SP.UDS, SP)
}

threshold <- 0.5

CTMM_RIDGE <- list()
KS_RIDGE <- list()

for (k in seq_along(UDS)) {
  
  UD <- UDS[[k]]
  
  #calculate ridge density!
  RIDGE <- ctmm:::ridges.UD(UD)
  ridge_indicator <- RIDGE$Indicator
  
  SP.UD <- SpatialPolygonsDataFrame.UD(UD)
  area_ud <- summary(UD, units=F)$CI[1,]

  #get uds at different thresholds
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
  
  #get all ridge lengths and densities
  
  #with added error messages to catch misbehaving cases
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
  
  CTMM_RIDGE[[k]] <- ridges_high
  print(paste("ctmm", k))
  
  ##KS method
  
  # ks doesn't well automate anything
  x <- data.frame(DATA[[k]][,c('x','y')])
  H <- UD$H
  w <- UD$weights * nrow(x)
  gridsize <- dim(UD$PDF)
  fhat <- ks::kde(x,H=H,w=w,gridsize=gridsize)
  density.cutoff <- ks::contourLevels(fhat,cont=95)
  threshold <- ks::contourProbs(fhat, cont=65)
  
  xmin <- c(UD$r$x[1],UD$r$y[1])
  xmax <- c(UD$r$x[gridsize[1]],UD$r$y[gridsize[2]])
  
  # this takes forever but looks good
  RIDGE_KS <- ks::kdr(x,H=H,w=w,fhat=fhat,density.cutoff=density.cutoff,
                      gridsize=gridsize,xmin=xmin,xmax=xmax,pre=FALSE)
  
  KS_RIDGE[[k]] <-  RIDGE_KS
  
}

SP.UDS <- SpatialPolygonsDataFrame.UD(UDS[[1]]) #SpatialPolygonsDataFrame.UD(UDS[[1]])
for (i in seq_along(UDS)) {
  i=i+1
  SP <- SpatialPolygonsDataFrame.UD(UDS[[i]])
  #crs(SP) <- crs(SpatialPolygonsDataFrame.UD(UDS[[1]]))
  SP.UDS <- rbind(SP.UDS, SP)
}

ext = extent(SP.UDS)

utmcoor<-SpatialPoints(rbind(cbind(ext[1], ext[3]),
                             cbind(ext[1], ext[4]),
                             cbind(ext[2], ext[3]),
                             cbind(ext[2], ext[4])), 
                       proj4string=crs(SP.UDS[1]))

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

##non-satellite imagery downloaded using extent from Stadia/Stamenmaps and saved as raster
background <- raster::stack("stamenmap.tif")
ras <- terra::project(terra::rast(background), crs(utmcoor, asText=T))
ras_crop <- terra::crop(ras, ext)

#example for vanak landscape
for(i in seq_along(UDS)) { print(UDS[[i]]@info$identity) }

raster::plotRGB(ras, axes=T, mar=1.5, buffer=F)
plot(SP.UDS, border=c("#D4B59E", "black","#818181"), add=T)
plot(KS_RIDGE[[1]], col="blue", type="l", lwd=4, add=T)
plot(KS_RIDGE[[2]], col="blue", type="l", lwd=4, add=T)
plot(KS_RIDGE[[3]], col="red", type="l", lwd=4, add=T)
plot(KS_RIDGE[[4]], col="red", type="l", lwd=4, add=T)
plot(KS_RIDGE[[5]], col="red", type="l", lwd=4, add=T)


