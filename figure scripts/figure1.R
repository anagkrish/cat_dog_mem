library(tidyverse)
library(ctmm)
library(lubridate)
library(parsedate)
library(dplyr)
library(geodist)
library(raster)
library(geosphere)
library(ggmap)
library(poisspatial)
library(foreach)
library(snow)
library(doSNOW)

#---CLEAN CSV SO ALL YOU HAVE TO DO IS LOAD IT
ridge <- read_csv("allridge.csv") %>%
  filter(sp!="Canis dingo") %>% #we're dropping dingos now
  mutate(sp = ifelse(sp=="Canis lupus x lycaon", "Canis lupus", sp), #and merging lupus x lycaon with lupus
         sp = ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp),
         sp = ifelse(sp=="Pseudalopex vetulus", "Lycalopex vetula", sp)) %>% #rename species
  mutate(phylo = ifelse(phylo=="Canis_mesomelas", "Lupulella_mesomelas", phylo),
         phylo = ifelse(phylo=="Pseudalopex_vetulus", "Lycalopex_vetula", phylo)) %>%
  mutate(pack_hunting = as.factor(pack_hunting),
         log_mass = log10(mass),
         log_hr = log(area_ud_est),
         inv_ess = 1/ess,
         log_roughness = log(mean_roughness),
         log_hfi = log(mean_hfi),
         log_dhi_gpp = log(mean_dhi_gpp),
         point = as.factor(1:length(id)),
         log_ridge = log(ridge_dens_est),
         hunting_movement = as.factor(hunting_movement),
         hunting_movement = fct_collapse(hunting_movement,
                                         `Mixed Strategies` = c("Mixed Strategies","Mixed Strategies\r\n")),
         hunting_cooperativity = as.factor(hunting_cooperativity),
         pursuit = ifelse(hunting_movement=="Pursuit", 1, 0),
         pursuit = as.factor(pursuit),
         disruptfast = ifelse(hunting_movement %in% c("Disruptive Fast hunting", "Mixed Strategies"), 
                              1, 0),
         disruptfast = as.factor(disruptfast),
         slowwalking = ifelse(hunting_movement %in% c("Slow walking", "Mixed Strategies"), 
                              1, 0),
         slowwalking = as.factor(slowwalking)) %>%
  rename("mean_road_cover"=mean_road_dens)

################# calculate ridges in same extent
load("movement/data/here")

inds <- (ridge %>%
  filter(updatedstudy=="Abrahms") #%>%
  #filter(sp %in% c("Chrysocyon brachyurus", "Puma concolor")) #oliveira-santos
  #filter(sp %in% c("Canis latrans", "Lynx rufus")) #for clark
  #filter(study=="Mahoney") #drops bobcats; for young
  #filter(grepl("MV", id)) #for prugh MV
  #filter(grepl("NE", id)) #for prugh NE
  #unite(id, sp, study, col="id2", sep="") %>% #for ramesh & drouilly
    #ramesh inds
  #filter(id2 %in% c("2Canis mesomelasRamesh", "2Caracal caracalRamesh", "4Caracal caracalRamesh",
    #"6Canis mesomelasRamesh", "ElevenLeptailurus servalRamesh", "FiveFiveLeptailurus servalRamesh",
    #"FiveLeptailurus servalRamesh", "FourteenLeptailurus servalRamesh", "NineLeptailurus servalRamesh",
    #"OneLeptailurus servalRamesh", "SevenLeptailurus servalRamesh", "SevensevenLeptailurus servalRamesh",
    #"SixLeptailurus servalRamesh", "TenLeptailurus servalRamesh", "TentenLeptailurus servalRamesh",
    #"ThreeLeptailurus servalRamesh", "TwelveLeptailurus servalRamesh", "TwoLeptailurus servalRamesh"))
    #drouilly inds
  #filter(id2 %in% c("BBJackal_EskimoCanis mesomelasDrouilly", "BBJackal_ForestCanis mesomelasDrouilly",
    #"BBJackal_LunaCanis mesomelasDrouilly", "BBJackal_RainCanis mesomelasDrouilly",
    #"LuckyCaracal caracalDrouilly", "MoonshineCaracal caracalDrouilly"))
    )$id

#modifications to tracks as necessary
mod_to_tracks <- function(track) {

  id <- unique(track$"individual.local.identifier")[[1]]
  species <- unique(track$"individual.taxon.canonical.name")[[1]]
  study <- unique(track$"study.id")[[1]]

  updatedstudy <- filter(updatedstudynames, ID==id, Species==species, Study==study)$`Updated Study`
  # print(paste(id, species, study, updatedstudy))

  ############# below are lines to edit problematic individuals as needed

  if (updatedstudy=="Clark") {

    # if (id%in%c("L19","L19_a","L19_b")) {
    # break
    #}
    if (id%in%c("L12","L27","R3","C263")) {
      track <- crop_range_res(track)
    }
    if (id=="C248"){
      track <- track %>%
        distinct(timestamp, .keep_all=TRUE)
    }
    if (id=="L30"){
      track <- track %>%
        filter(timestamp > ymd_hms("2022-01-01 07:00:39"))
    }
  }

  if (updatedstudy=="Drouilly") {
    if (id%in%c("BBJackal_Rain")) {
      track <- crop_range_res(track)
    }
    if (id=="Lucky") {
      track <- crop_range_res(track) %>%
        filter(location.lat > -32.845)
    }
  }

  if (updatedstudy=="Oliveira-Santos.Dataset1") {
    if (id%in%c("150011","150041","150102","150312","150402","150462",
                "150552","150681","163181","164820","164886","164900",
                "164957","1649671","165164","165194","1651941","165224","165252","1652521")) {
      break
      #skip Cerdocyon thous duplicates in OS dataset1
    }
    if (id%in%c("CAN13","CAN49", "Grupo005_Id001","Grupo005_Id032")) {
      track <- crop_range_res(track)
    }
  }

  if (updatedstudy=="Prugh") {
    if (id%in%c("C_NECOY20F","B_MVBOB54F","B_MVBOB66M","B_NEBOB33M","B_NEBOB35M")) {
      track <- crop_range_res(track)
    }
    if (id=="C_MVCOY98M") {
      track <- crop_range_res(track) %>%
        filter(timestamp < ymd_hms("2020-09-01 00:00:00"))
    }
    if (id=="B_NEBOB23M") {
      track <- crop_range_res(track) %>%
        filter(location.long < -118)
    }
  }

  if (updatedstudy=="Sekercioglu") {
    if (id%in%c("Bilge")) {
      track <- crop_range_res(track)
    }
  }

  if (updatedstudy=="Ramesh") {

    if (id%in%c("4")&species%in%c("Canis mesomelas")) {
      track <- crop_range_res(track)
    }
    if (id=="4" & species=="Caracal caracal") {
      track <- crop_range_res(track)
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
utmcoor<-SpatialPoints(rbind(cbind(ext[1], ext[3]),
                             cbind(ext[1], ext[4]),
                             cbind(ext[2], ext[3]),
                             cbind(ext[2], ext[4])), 
                       proj4string=crs(UDS[[1]]))

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

background <- get_stamenmap(c(extent(longlatcoor)[1]-0.01, extent(longlatcoor)[3]-0.02,
                              extent(longlatcoor)[2]+0.01, extent(longlatcoor)[4]+0.01), zoom=12, maptype="terrain")

background <- poisspatial::ps_ggmap_to_raster(background)
crs(background) <- longlatcoor

drouilly_ras <- terra::project(terra::rast(background), crs(utmcoor, asText=T))
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

par(xpd=TRUE)
legend(c(-136093.1, -118000), c(-52679.85, -47479.85), " Western Cape, South Africa\n\n", 
       box.col = "black", bg = "white", cex=0.85, adj = c(0.12,0.1)) #vanak
text(-134093.1, -49679.85, substitute(paste(italic("Lupulella mesomelas"))), 
     col="blue", cex=0.85, adj = c(0.12)) #vanak
text(-134290.1, -51079.85, substitute(paste(italic("Caracal caracal"))),
     col="red", cex=0.85, adj = c(0.12)) #vanak

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
utmcoor<-SpatialPoints(rbind(cbind(ext[1], ext[3]),
                             cbind(ext[1], ext[4]),
                             cbind(ext[2], ext[3]),
                             cbind(ext[2], ext[4])), 
                       proj4string=crs(UDS[[1]]))

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

background <- get_stamenmap(c(extent(longlatcoor)[1]-0.01, extent(longlatcoor)[3]-0.02,
                              extent(longlatcoor)[2]+0.01, extent(longlatcoor)[4]+0.01), zoom=12, maptype="terrain")

background <- poisspatial::ps_ggmap_to_raster(background)
crs(background) <- longlatcoor

vanak_ras <- terra::project(terra::rast(background), crs(utmcoor, asText=T))
plotRGB(vanak_ras, axes=T, mar=1.5, buffer=F)

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

par(xpd=TRUE)
legend(10000, 300, title ="Clade", c("Canidae", "Felidae"), fill=c("blue","red"), box.col = "black", bg = "white",
       cex=0.9, text.col="black") #vanak
legend(-8700, 950, "Maharashtra, India   \n\n\n\n", box.col = "black", bg = "white", cex=0.9, adj = c(0.12,0.1)) #vanak
text(-8080, -200, substitute(paste(italic("Canis aureus(1)\nVulpes bengalensis(2)"))), 
     col="blue", cex=0.9, adj = c(0.12)) #vanak
text(-8180, -700, substitute(paste(italic("Felis chaus"))),
     col="red", cex=0.9, adj = c(0.12)) #vanak

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
utmcoor<-SpatialPoints(rbind(cbind(ext[1], ext[3]),
                             cbind(ext[1], ext[4]),
                             cbind(ext[2], ext[3]),
                             cbind(ext[2], ext[4])), 
                       proj4string=crs(UDS[[1]]))

longlatcoor<-spTransform(utmcoor,CRS("+proj=longlat"))

background <- get_stamenmap(c(extent(longlatcoor)[1]-0.01, extent(longlatcoor)[3]-0.01,
                              extent(longlatcoor)[2]+0.01, extent(longlatcoor)[4]+0.05), zoom=12, maptype="terrain")

background <- poisspatial::ps_ggmap_to_raster(background)
crs(background) <- longlatcoor

young_ras <- terra::project(terra::rast(background), crs(utmcoor, asText=T))
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

#put them together
library(cowplot)

vanak <- ggdraw() +  draw_image("/vanak/plot/here", clip="on")

young <- ggdraw() + 
  theme(plot.background = element_rect(fill="transparent", color = NA)) + 
  draw_image("young/plot/here", clip="on")

drouilly <- ggdraw() + 
  theme(plot.background = element_rect(fill="transparent", color = NA)) + 
  draw_image("drouilly/plot/here", clip="on")

#FINAL figure 1!
plot_grid(vanak, NULL,
          plot_grid(young, drouilly,
                    nrow=1, ncol=2, 
                    rel_widths=c(0.85,1),
                    labels=c("B)","C)"),
                    label_size = 30,
                    label_x = 0.065, label_y = 0.86),
          nrow=3, ncol=1,
          rel_heights=c(1, -0.3, 1.3),
          labels=c("A)","",""),
          label_size = 30,
          label_x = 0.05, label_y = 0.92)


