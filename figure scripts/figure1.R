library(tidyverse)
library(ctmm)
library(lubridate)
library(parsedate)
library(dplyr)
library(geodist)
library(raster)
library(geosphere)
library(foreach)
library(snow)
library(doSNOW)

#load full csv
ridge <- read_csv("ridge.csv") %>%
  filter(`kept (y/n)` == "y") %>% #drop all dropped individuals
  mutate(pursuit = as.factor(pursuit),
         disruptfast = as.factor(disruptfast),
         slowwalking = as.factor(slowwalking), 
         #convert to factors bc they're read in as numerical
         log_mass = log(`mass (g)`),
         log_hr = log(`area_ud_est (m^2)`),
         inv_ess = 1/ess,
         log_roughness = log(mean_roughness),
         log_hfi = log(mean_hfi),
         log_dhi_gpp = log(mean_dhi_gpp),
         point = as.factor(1:length(id)),
         log_ridge = log(`ridge_dens_est (1/m)`)) %>%
  #standardize vars
  mutate(log_mass_st=mosaic::zscore(log_mass),
         log_hr_st=mosaic::zscore(log_hr),
         log_roughness_st=mosaic::zscore(log_roughness),
         seasonality_dhi_gpp_st=mosaic::zscore(seasonality_dhi_gpp),
         speed_est_st=mosaic::zscore(`speed_est (m/s)`, na.rm=T),
         mean_treecover=mean_treecover/100)

load("movement/data/here")

#For panels A-D
################## example workflow for ridge density calculations
#puma concolor F70 from Young study

#set projection because this is what we use when we display all individuals on terrain - this way orientations
#line up between example and terrain visualizations

proj <- as.telemetry(alldata %>%
                       unite(individual.local.identifier, individual.taxon.canonical.name,
                             study.id, col="id", sep="", remove=F) %>%
                       filter(id==str_replace("C023Canis latransMahoney", ".rda", "")))

tracks <- alldata %>% filter(individual.local.identifier == "F70", study.id=="Mahoney")

dat <- as.telemetry(tracks)
ctmm:::projection(dat) <- median(proj) #align projections

#workflow to calculate ridges
SVF <- variogram(dat, CI="Gauss") #gauss CI for more accurate error bars
GUESS <- ctmm.guess(dat,variogram=SVF,interactive=F)
BEST_FIT <- ctmm.select(dat, GUESS2, level = 0.95, cores = 2, trace=2)
UD <- akde(data=dat, CTMM=BEST_FIT, grid=list(extent(proj)), 
           debias=FALSE, weights=FALSE)
RIDGE <- ctmm:::ridges.UD(UD2)

SP.UD <- SpatialPolygonsDataFrame.UD(UD2)
threshold <- 0.5

#get UD at low medium and high 95% core area
UD_High <- SpatialPolygonsDataFrame.UD(UD)
UD_High@polygons[[1]] <- NULL
UD_High@polygons[[1]] <- NULL
UD_High@plotOrder <- c(1L)
UD_High <- rgeos::gMakeValid(UD_High)

UD_Mean <- SpatialPolygonsDataFrame.UD(UD)
UD_Mean@polygons[[3]] <- NULL
UD_Mean@polygons[[1]] <- NULL
UD_Mean@plotOrder <- c(1L)
UD_Mean <- rgeos::gMakeValid(UD_Mean)

UD_Low <- SpatialPolygonsDataFrame.UD(UD)
UD_Low@polygons[[3]] <- NULL
UD_Low@polygons[[2]] <- NULL
UD_Low@plotOrder <- c(1L)
UD_Low <- rgeos::gMakeValid(UD_Low) 

#get ridges as contourLines
ridges <- grDevices::contourLines(UD$r$x, UD$r$y, RIDGE, level=threshold)

lines <- list()
for (i in seq_along(ridges)) {
  lines[[i]] <- Lines(list(Line(cbind(ridges[[i]]$x, ridges[[i]]$y))), 
                      ID = as.character(i))
}

#convert to spatiallines
all_lines <- SpatialLines(lines, proj = SP.UD@proj4string)

ridges_mean <- tryCatch({ raster::intersect(all_lines, UD_Mean) },
                         error= function(e) { ridges_mean <- 0 })

#transform axes into km to plot
ridges_km <- spTransform(ridges_mean,"+proj=aeqd +lat_0=38.6480207305513 +lon_0=-112.074775967984 +x_0=0 +y_0=0 +datum=WGS84 +units=km")

#plots
plot(dat, col="black", error=F, type="l",cex.axis=1.15, cex.lab=1.15) #A: GPS tracks as line

xlim=c(0, tail(SVF2$lag,1))
plot(SVF,BEST_FIT, xlim=xlim, cex.axis=1.15, cex.lab=1.15) #B: variogram

plot(dat, col="black", UD=UD2, col.DF="#FF6969", cex.axis=1.15, cex.lab=1.15) #C: UD with GPS tracks as point

plot(dat, col="black", error=F, UD=UD, col.DF="#FF6969", cex.axis=1.15, cex.lab=1.15) #D: ridges w/ UD and GPS tracks as point
plot(ridges_km, col="#960000",lwd=3,add=T)


#For panels E-G
################# calculate ridges in same extent

inds <- (ridge %>%
           filter(updatedstudy=="Drouilly") #%>% (studies are Drouilly, Young and Vanak)
         #unite(id, sp, study, col="id2", sep="") %>% 
          #drouilly inds
         #filter(id2 %in% c("BBJackal_Forest", "BBJackal_Rain","Lucky"))
          #vanak inds
         #filter(id2 %in% c("Fox 13 (Vasu)", "Fox 22 (Ding Dong)", "Fox 23 (300)", "Jackal 03 (Nanda)",
         #"Jackal 07 (Marker)", "Jackal 09 (Roma)", "Jungle cat 09 (Dagdu)", "Jungle cat 10 (Bheegi Billi)",
         #"Jungle cat 12 (Cavity)", "Jungle cat 13 (Kamini)"))
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

#minsampling
get_minsampling <- function(individual) {
  minsampling <- min(diff(individual$t))
  return(minsampling) #pretty sure units are in seconds
}

#resample inds with minsampling <60s
resample <- function(track) {

  if (get_minsampling(as.telemetry(track)) < 60) {
    track_rs <- track %>%
      mutate(timestamp_rs = sapply(ymd_hms(timestamp), cut,
                                   breaks = "1 min")) %>% #special feature in lubrdiate pkg!
      distinct(timestamp_rs, .keep_all = "TRUE") %>% #collapse to individual resampled time points
      dplyr::select(-c("timestamp")) %>%
      rename("timestamp" = timestamp_rs) #added to resample sketchy individuals

    return(track_rs)
  }

  else {
    return(track)
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

#meant to be run in parallel, returns list of UD, ridge, UD, ridge etc for each ind

#set projection to first inds, not sure how else to go about this :(
#allows for multiple UDs to be plotted in same extent

proj <- as.telemetry(alldata %>%
                       unite(individual.local.identifier, individual.taxon.canonical.name,
                             study.id, col="id", sep="", remove=F) %>%
                       filter(id==str_replace(inds[1], ".rda", "")))

landscape_ridges <- list()

for (i in inds) {

  track <- alldata %>%
    unite(individual.local.identifier, individual.taxon.canonical.name,
          study.id, col="id", sep="", remove=F) %>%
    mutate(id=str_replace(id, "/","")) %>%
    filter(id==str_replace(inds[i], ".rda$", "")) %>%
    dplyr::select(-c("id")) %>%
    #above two lines for ramesh inds bc the naming system sucks
    resample() %>% #resample if less than a minute minsampling
    mod_to_tracks() %>%
    unite(individual.local.identifier, individual.taxon.canonical.name,
          col="individual.local.identifier", sep=" ", remove=F)

  track <- track%>%
    as.telemetry()

  print(head(track))
  ctmm:::projection(track) <- ctmm::median(proj) #set projections to be the same for all

  print(head(track))

  SVF <- variogram(track)

  GUESS <- ctmm.guess(track,variogram=SVF,interactive=F)
  BEST_FIT <- ctmm.select(track, GUESS, trace = TRUE, level = 0.95, cores = 2)
  UD <- akde(track, BEST_FIT, debias=FALSE, weights=FALSE)

  #calculate ridge density!
  RIDGE <- ctmm:::ridges.UD(UD)
  SP.UD <- SpatialPolygonsDataFrame.UD(UD)
  threshold <- 0.5

  #get UDs at diff thresholds
  UD_High <- SpatialPolygonsDataFrame.UD(UD)
  UD_High@polygons[[1]] <- NULL
  UD_High@polygons[[1]] <- NULL
  UD_High@plotOrder <- c(1L)
  UD_High <- rgeos::gMakeValid(UD_High)

  UD_Mean <- SpatialPolygonsDataFrame.UD(UD)
  UD_Mean@polygons[[3]] <- NULL
  UD_Mean@polygons[[1]] <- NULL
  UD_Mean@plotOrder <- c(1L)
  UD_Mean <- rgeos::gMakeValid(UD_Mean)

  UD_Low <- SpatialPolygonsDataFrame.UD(UD)
  UD_Low@polygons[[3]] <- NULL
  UD_Low@polygons[[2]] <- NULL
  UD_Low@plotOrder <- c(1L)
  UD_Low <- rgeos::gMakeValid(UD_Low)

  #get ridges as contourLines
  ridges <- grDevices::contourLines(UD$r$x, UD$r$y, RIDGE, level=threshold)

  lines <- list()
  for (i in seq_along(ridges)) {
    lines[[i]] <- Lines(list(Line(cbind(ridges[[i]]$x, ridges[[i]]$y))),
                        ID = as.character(i))
  }

  #convert to spatiallines
  all_lines <- SpatialLines(lines, proj = SP.UD@proj4string)

  ridges_mean <- tryCatch({ raster::intersect(all_lines, UD_Mean) },
                          error= function(e) { ridges_mean <- 0 })

  #add to list
  landscape_ridges <- c(landscape_ridges, list(UD), list(ridges_mean))

}

################# for each landscape, get uds/sp.uds and ridges

###### drouilly 
#species: Canis aureus, Caracal caracal
#individuals: BBJackal_Forest, BBJackal_Rain, Lucky

AKDES <- readRDS("drouilly/akdes/here")
#naming is slightly confusing; this is just the spatialpolygons version of the akde output above
UDS <- readRDS("drouilly/uds/here") 
ridges <- readRDS("drouilly/ridges/here")

ext = extent(SpatialPolygons(lapply(UDS, function(x){x@polygons[[3]]})))
utmcoor<-SpatialPoints(rbind(cbind(ext[1]-1000, ext[3]-1000),
                             cbind(ext[1]-1000, ext[4]+1000),
                             cbind(ext[2]+1000, ext[3]-1000),
                             cbind(ext[2]+1000, ext[4]+1000)), 
                       proj4string=crs(UDS[[1]]))

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

##satellite imagery downloaded via google earth engine sentinel-2
#drouilly_ras <- terra::rast("satellite/raster/here") 
#drouilly_ras <- terra::project(drouilly_ras, crs(utmcoor, asText=T))
#drouilly_ras <- terra::crop(drouilly_ras, extent(utmcoor))

##non-satellite imagery downloaded from Stadia/Stamenmaps and saved as raster
drouilly_ras <- terra::rast("terrain/raster/here") 
plotRGB(drouilly_ras, axes=T, mar=1.5, buffer=F) #just background

#plot ridges on top of background
par(bg=NA, cex=0.68)

#change names to match dataset (sp canis mesomelas --> lupulella mesomelas)
AKDES[[1]]@info$identity <- "BBJackal_Forest Lupulella mesomelas"
AKDES[[2]]@info$identity <- "BBJackal_Rain Lupulella mesomelas"

plotRGB(drouilly_ras, axes=T, mar=1.5, buffer=F) #for background terrain

clade <- (ridge %>% 
            rowwise() %>% 
            mutate(id=paste(id,sp, sep=" ")) %>% 
            filter(updatedstudy=="Drouilly", id==AKDES[[1]]@info$identity))$clade

if(is_empty(clade)==TRUE) {
  col=NA } else if(clade=="canidae") { 
    col="blue"} else if(clade=="felidae") { 
      col="red"}

if(is_empty(clade)==FALSE) {
  
  plot(SpatialPolygons(UDS[[1]]@polygons[1]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border="#D4B59E")
  plot(SpatialPolygons(UDS[[1]]@polygons[2]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border="black")
  plot(SpatialPolygons(UDS[[1]]@polygons[3]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border = "#818181")
  plot(drouillyridge[[1]], col=col, xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000), 
       axes=T, add=T)
  
}

print(paste(AKDES[[1]]@info$identity, clade))

for (i in 2:length(ridges)) {
  
  clade <- (ridge %>%
              rowwise() %>% 
              mutate(id=paste(id,sp, sep=" ")) %>% 
              filter(updatedstudy=="Drouilly", 
                     id==AKDES[[i]]@info$identity))$clade
  
  if(is_empty(clade)==TRUE) {
    col=NA } else if(clade=="canidae") { 
      col="blue"} else if(clade=="felidae") { 
        col="red"}
  
  if(is_empty(clade)==FALSE) {
    plot(SpatialPolygons(UDS[[i]]@polygons[1]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, axes=T, add=T, border="#D4B59E")
    plot(SpatialPolygons(UDS[[i]]@polygons[2]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, axes=T, add=T, border="black")
    plot(SpatialPolygons(UDS[[i]]@polygons[3]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, axes=T, add=T, border="#818181")
    plot(ridges[[i]], col=col, xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000), axes=T, add=TRUE)
    
  }
  
  print(paste(i, AKDES[[i]]@info$identity, clade))
  
}

dev.off()

#####
#vanak
#species: Vulpes bengalensis, Canis aureus, Felis Chaus
#individuals: Fox 13 (Vasu), Fox 22 (Ding Dong), Fox 23 (300), Jackal 03 (Nanda), Jackal 07 (Marker),
  #Jackal 09 (Roma), Jungle cat 09 (Dagdu), Jungle cat 10 (Bheegi Billi), Jungle cat 12 (Cavity),
  #Jungle cat 13 (Kamini)

AKDES <- readRDS("vanak/akdes/here")
UDS <- readRDS("vanak/uds/here")
ridges <- readRDS("vanak/ridges/here")

ext = extent(SpatialPolygons(lapply(UDS, function(x){x@polygons[[3]]})))
utmcoor<-SpatialPoints(rbind(cbind(ext[1]-1000, ext[3]-1000),
                             cbind(ext[1]-1000, ext[4]+1000),
                             cbind(ext[2]+1000, ext[3]-1000),
                             cbind(ext[2]+1000, ext[4]+1000)), 
                       proj4string=crs(UDS[[1]]))

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

##satellite imagery downloaded via google earth engine sentinel-2
#vanak_ras <- terra::rast("satellite/raster/here") 
#vanak_ras <- terra::project(vanak_ras, crs(utmcoor, asText=T))
#vanak_ras <- terra::crop(vanak_ras, extent(utmcoor))

##non-satellite imagery downloaded from Stadia/Stamenmaps and saved as raster
vanak_ras <- terra::rast("terrain/raster/here") 
plotRGB(vanak_ras, axes=T, mar=1.5, buffer=F) #just background

par(bg=NA) #comment out to get background back

plotRGB(vanak_ras, axes=T, mar=1.5, buffer=F) #for background terrain

#plot(ind_ras) #for all other rasts
clade <- (ridge %>% 
            rowwise() %>% #for ramesh
            mutate(id=paste(id,sp, sep=" ")) %>% #for ramesh
            filter(updatedstudy=="Vanak", id==AKDES[[1]]@info$identity))$clade
sp <- (ridge %>% 
         rowwise() %>% #for ramesh
         mutate(id=paste(id,sp, sep=" ")) %>% #for ramesh
         filter(updatedstudy=="Vanak", id==AKDES[[1]]@info$identity))$sp

if(is_empty(clade)==TRUE) {
  col=NA } else if(clade=="canidae") { 
    col="blue"} else if(clade=="felidae") { 
      col="red"}

if(is_empty(clade)==FALSE) {
  
  plot(SpatialPolygons(UDS[[1]]@polygons[1]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border="#D4B59E")
  plot(SpatialPolygons(UDS[[1]]@polygons[2]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border="black")
  plot(SpatialPolygons(UDS[[1]]@polygons[3]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border = "#818181")
  plot(ridges[[1]], col=col, xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000), axes=T, add=T)
  
  if (sp=="Vulpes bengalensis") { text(UDS[[1]]@bbox[[1]], UDS[[1]]@bbox[[2]], "2", col="blue") }
  if (sp=="Canis aureus") { text(UDS[[1]]@bbox[[1]], UDS[[1]]@bbox[[2]], "1", col="blue") }
  
}

print(paste(AKDES[[1]]@info$identity, clade))

for (i in 2:length(ridges)) {
  
  clade <- (ridge %>%
              mutate(id=paste(id,sp, sep=" ")) %>% #for ramesh
              filter(updatedstudy=="Vanak", id==AKDES[[i]]@info$identity))$clade
  sp <- (ridge %>% 
           rowwise() %>% #for ramesh
           mutate(id=paste(id,sp, sep=" ")) %>% #for ramesh
           filter(updatedstudy=="Vanak", id==AKDES[[i]]@info$identity))$sp
  
  if(is_empty(clade)==TRUE) {
    col=NA } else if(clade=="canidae") { 
      col="blue"} else if(clade=="felidae") { 
        col="red"}
  
  if(is_empty(clade)==FALSE) {
    plot(SpatialPolygons(UDS[[i]]@polygons[1]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, add=T, border="#D4B59E")
    plot(SpatialPolygons(UDS[[i]]@polygons[2]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, add=T, border="black")
    plot(SpatialPolygons(UDS[[i]]@polygons[3]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, add=T, border="#818181")
    plot(ridges[[i]], col=col, xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000), add=TRUE)
    
    if (sp=="Vulpes bengalensis") { text(UDS[[i]]@bbox[[3]]-350, UDS[[i]]@bbox[[4]]-500, "2", col="blue") }
    if (sp=="Canis aureus") { text(UDS[[i]]@bbox[[3]], UDS[[i]]@bbox[[4]]-1200, "1", col="blue") }
    
    
  }
  
  print(paste(i, AKDES[[i]]@info$identity, clade))
  
}

#####
#young
#species: Canis latrans, Puma concolor
#individuals: C023, C026, C028, F108, F157, F159, F70
AKDES <- readRDS("young/akdes/here")
UDS <- readRDS("young/uds/here")
ridges <- readRDS("young/ridges/here")

ext = extent(SpatialPolygons(lapply(UDS, function(x){x@polygons[[3]]})))
utmcoor<-SpatialPoints(rbind(cbind(ext[1]-1000, ext[3]-1000),
                             cbind(ext[1]-1000, ext[4]+1000),
                             cbind(ext[2]+1000, ext[3]-1000),
                             cbind(ext[2]+1000, ext[4]+1000)), 
                       proj4string=crs(UDS[[1]]))

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

##satellite imagery downloaded via google earth engine sentinel-2
#young_ras <- terra::rast("satellite/raster/here") 
#young_ras <- terra::project(young_ras, crs(utmcoor, asText=T))
#young_ras <- terra::crop(young_ras, extent(utmcoor))

##non-satellite imagery downloaded from Stadia/Stamenmaps and saved as raster
young_ras <- terra::rast("terrain/raster/here") 
plotRGB(young_ras, axes=T, mar=1.5, buffer=F) #just background


#add ridges
par(bg=NA)

plotRGB(young_ras, axes=T, mar=1.5, buffer=F) #for background terrain

clade <- (ridge %>% 
            rowwise() %>% 
            mutate(id=paste(id,sp, sep=" ")) %>% 
            filter(updatedstudy=="Young", id==AKDES[[1]]@info$identity))$clade

if(is_empty(clade)==TRUE) {
  col=NA } else if(clade=="canidae") { 
    col="blue"} else if(clade=="felidae") { 
      col="red"}

if(is_empty(clade)==FALSE) {
  
  plot(SpatialPolygons(UDS[[1]]@polygons[1]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border="#D4B59E")
  plot(SpatialPolygons(UDS[[1]]@polygons[2]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border="black")
  plot(SpatialPolygons(UDS[[1]]@polygons[3]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
       lwd = 1.5, axes=T, add=T, border = "#818181")
  plot(ridges[[1]], col=col, xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000), 
       axes=T, add=T)
  
}

print(paste(AKDES[[1]]@info$identity, clade))

par(xpd=TRUE)
legend(-26592, 16500, "Washington, USA\n\n\n", box.col = "black", bg = "white", cex=0.9, adj = c(0.12,0.1)) #vanak
text(-24000, 12500, substitute(paste(italic("Canis latrans"))), 
     col="blue", cex=0.9, adj = c(0.12)) #vanak
text(-24000, 10000, substitute(paste(italic("Puma concolor"))),
     col="red", cex=0.9, adj = c(0.12)) #vanak

for (i in 2:length(youngridge)) {
  
  clade <- (ridge %>%
              rowwise() %>% #for ramesh
              mutate(id=paste(id,sp, sep=" ")) %>% #for ramesh
              filter(updatedstudy=="Young", id==youngAKDES[[i]]@info$identity))$clade
  
  if(is_empty(clade)==TRUE) {
    col=NA } else if(clade=="canidae") { 
      col="blue"} else if(clade=="felidae") { 
        col="red"}
  
  if(is_empty(clade)==FALSE) {
    plot(SpatialPolygons(UDS[[i]]@polygons[1]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, axes=T, add=T, border="#D4B59E")
    plot(SpatialPolygons(UDS[[i]]@polygons[2]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, axes=T, add=T, border="black")
    plot(SpatialPolygons(UDS[[i]]@polygons[3]), xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000),
         lwd = 1.5, axes=T, add=T, border="#818181")
    plot(ridges[[i]], col=col, xlim=c(ext[1]-1000, ext[2]+1000), ylim=c(ext[3]-1000,ext[4]+1000), axes=T, add=TRUE)
    
  }
  
  print(paste(i, AKDES[[i]]@info$identity, clade))
  
}

#figures 

