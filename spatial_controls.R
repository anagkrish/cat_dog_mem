library(tidyverse)
library(ctmm)
library(sp)
library(sf)
library(raster)
library(rgdal)
library(rgeos)
library(maptools)
library(geosphere)
library(lubridate)
library(ncdf4)

"%ni%" <- Negate("%in%")

#load movement object alldata for all 1528 individuals
load("movement/data/here")

#copied from getridgefromfits, edit tracks so you extract correct extent for cropped individuals
mod_to_tracks <- function(track) {
  
  id <- track$individual.local.identifier[[1]]
  updatedstudy <- filter(updatedstudynames, 
                         ID==id, 
                         Species==track$individual.taxon.canonical.name[[1]],
                         Study==track$study.id[[1]])$`Updated Study`
  
  ############# below are lines to edit problematic individuals as needed
  if (updatedstudy=="Belant-Beyer") {
    if (id%in%c("W07")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Berteaux") {
    if (id%in%c("BORR","BVOB","JVOJ","OBBB")) {
      track <- crop_range_res(track)
    }
    if (id=="ORRR") {
      track <- track %>%
        filter(location.long > NA, #coords withheld 
               location.lat < NA) #coords withheld
    }
  }
  
  if(updatedstudy=="Clark") {
    if (id%in%c("L12","L27","R3","C263")) {
      track <- crop_range_res(track)
    }
    
    if (id=="C248"){
      track <- track %>%
        distinct(timestamp, .keep_all=TRUE)
    }
    
  }
  
  if (updatedstudy=="Conner") {
    if (id%in%c("F46","M21", "F30")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Cristescu") {
    if (id%in%c("NCM1","NCM13","NCM8","NCM9")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Darlington") {
    if (id%in%c("C30","C8")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Drouilly") {
    if (id%in%c("BBJackal_Rain")) {
      track <- crop_range_res(track)
    }
    
    if (id=="Lucky") {
      track <- crop_range_res(track) %>%
        filter(location.lat > NA) #coords withheld
    }
  }
  
  if (updatedstudy=="Ferrell") {
    if (id%in%c("B21")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Frair") {
    if (id=="M5") {
      track <- track %>%
        filter(timestamp > ymd_hms("2006-01-01 18:28:59")) %>%
        crop_range_res()
    }
  }
  
  if (updatedstudy=="Fryxell") {
    
    if (id=="148.66") {
      track <- crop_range_res(track)
    }
    
    if (id=="148.81") {
      track <- track %>%
        mutate(timestamp=ymd_hms(timestamp)) %>%
        filter(location.long > NA, #coords withheld 
               location.lat < NA) #coords withheld
    }
    
  }
  
  if (updatedstudy=="Garthe") {
    if (id%in%c("REDFOX-2015-01-BHKoog",
                "REDFOX-2016-03-BHKoog",
                "REDFOX-2018-03-Sylt")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Getz-Bellan") {
    if (id%in%c("CM18","CM26","CM36","CM95", "CM83")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Hebblewhite") {
    if (id%in%c("B065","J031")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Hernandez-Blanco") {
    if (id=="Tsetseg") {
      track <- crop_range_res(track)
    }
    if(id=="Safar") {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Jachowski") {
    if (id%in%c("F33","F39")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Jackson") {
    if (id%in%c("NAN00177")) {
      track <- crop_range_res(track)
    }
    
    if (id=="NAN00178"){
      track <- track %>%
        mutate(timestamp = ymd_hms(timestamp)) %>%
        filter(timestamp > ymd_hms("2016-06-14 20:45:32"))
    }
  }
  
  if (updatedstudy=="Kays") {
    
    if (id%in%c("30822")) {
      track <- crop_range_res(track)
    }
    
    if (id=="30839") {
      track <- track %>%
        filter(timestamp > ymd_hms("2011-08-01 00:03:07"))
    }
  }
  
  if (updatedstudy=="VanDerWeyde-Kral") {
    if (id%in%c("Kealiboka","Matsoshetsi")) {
      track <- crop_range_res(track)
    }
  }
  
  
  if (updatedstudy=="Morato") {
    if (id%in%c("Panthera_11", "Panthera_114", "Panthera_13", "Panthera_27","Panthera_34", "Panthera_36",
                "Panthera_44")) {
      track <- crop_range_res(track)
    }
    
    if (id=="Panthera_106"){
      track <- track %>%
        filter(timestamp<=ymd_hms("2010-01-01 10:00:00"))
    }
    
    if (id=="Panthera_111"){
      track <- track %>%
        filter(timestamp<=ymd_hms("2009-11-18 04:00:00"))
    }
    
    if (id=="Panthera_43"){
      track <- track %>%
        filter(timestamp<=ymd_hms("2013-02-28 11:00:00"))
    }
    
    if (id=="Panthera_36"){
      track <- track %>%
        filter(timestamp<=ymd_hms("2000-10-18 15:50:00"))
    }
    
  }
  
  if (updatedstudy=="Azevedo.Lemos") {
    if (id%in%c("150011","150041","150102","150312","150402","150462",
                "150552","150681","163181","164820","164886","164900",
                "164957","1649671","165164","165194","1651941","165224",
                "165252","1652521")) {
      
      #return empty df for tracks to skip (can't use next since not formally in loop)
      track <- data.frame(location.lat=0, location.long=0, 
                          timestamp=ymd_hms("1970-01-01 00:00:00"),
                          individual.local.identifier=track$individual.local.identifier[[1]], 
                          individual.taxon.canonical.name=track$individual.taxon.canonical.name[[1]],
                          study.id=track$study.id[[1]])
      
      #skip Cerdocyon thous duplicates in OS dataset1
    }
    
    if (id%in%c("CAN13","CAN49", "Grupo005_Id001","Grupo005_Id032")) {
      track <- crop_range_res(track)
    }
    
    if (id=="shack") {
      track <- crop_range_res(track)
    }
    
  }
  
  if (updatedstudy=="Oliveira-Santos") {
    if (id%in%c("150041_GustavoCT")) {
      track <- crop_range_res(track)
    }
    if (id=="kayapo") {
      track <- track %>%
        filter(location.lat > NA) #coords withheld
    }
  }
  
  if (updatedstudy=="Palomares") {
    if (id%in%c("Garfio", "Patsuezo")) {
      
      track <- crop_range_res(track)
      
    } 
  }
    
    if (updatedstudy=="Palacios Gonzalez") {
      if (id%in%c("Lynx_Ketamina", "Lynx_Llerena", 
                "Lynx_Miera", "Lynx_Negral", "Lynx_Nitrógeno")) {
      track <- crop_range_res(track)
    }
    
    if (id=="Lynx_Neruda") {
      track <- crop_range_res(track)
    }
    
  }
  
  if (updatedstudy=="Patterson.B") {
    
    if (id%in%c("W97143","W97155","W97177","W97185","W97188","W97189","W97311","W97359",
                "H01", "T21", "W215", "W1605")) {
      track <- crop_range_res(track)
    }
    
    if(id=="T50") { 
      track <- track %>%
        filter(location.lat > NA) #coords withheld
    }
    
    if (track$individual.taxon.canonical.name[[1]]=="Canis lupus x lycaon") {
      if (id %in% c("W97104","W97107","W97110","W97141","W97142","W97147",
                    "W97162","W97303","W97305","W97306","W97356","W97362")) {
        
        track <- data.frame(location.lat=0, location.long=0, 
                            timestamp=ymd_hms("1970-01-01 00:00:00"),
                            individual.local.identifier=track$individual.local.identifier[[1]], 
                            individual.taxon.canonical.name=track$individual.taxon.canonical.name[[1]],
                            study.id=track$study.id[[1]])
      }
    }
    
  }
  
  if (updatedstudy=="Prugh") {
    if (id%in%c("C_MVBOB91M", "C_NEBOB13F", "C_MVBOB80M", "C_MVBOB66M", "C_NEBOB35M", "C_NEBOB5M",
                "C_MVBOB67M", "C_NEBOB11M", "C_MVBOB87M", "C_MVBOB99F", "C_NEBOB6F", "C_MVBOB88M",
                "C_MVBOB85F", "C_MVBOB51M", "C_MVBOB62M", "C_MVBOB69F", "C_MVBOB90M", "C_NEBOB7F",
                "C_MVBOB76M", "C_NEBOB33M", "C_NEBOB8M", "C_MVBOB77M", "C_NEBOB38M", "C_MVBOB83M",
                "C_NEBOB45M", "C_MVBOB52M", "C_NEBOB25F", "C_NEBOB37M", "C_NEBOB41F", "C_NEBOB16M",
                "C_NEBOB32F", "C_NEBOB10F")) {
      
      track <- data.frame(location.lat=0, location.long=0, 
                          timestamp=ymd_hms("1970-01-01 00:00:00"),
                          individual.local.identifier=track$individual.local.identifier[[1]], 
                          individual.taxon.canonical.name=track$individual.taxon.canonical.name[[1]],
                          study.id=track$study.id[[1]]) #skip duplicated coyotes
    }
    
    if (id%in%c("B_MVBOB54F","B_MVBOB66M","B_NEBOB33M","B_NEBOB35M", "C_NECOY20F")) {
      track <- crop_range_res(track)
    }
    
    if (id=="C_MVCOY98M") { 
      track <- crop_range_res(track) %>%
        filter(timestamp < ymd_hms("2020-09-01 00:00:00"))
    }
    
    if (id=="B_NEBOB23M") { 
      track <- crop_range_res(track) %>%
        filter(location.long < NA) #coords withheld
    }
    
  }
  
  
  if (updatedstudy=="Ramesh") {
    if (id=="4") { #in this case both 4s need cropping which is lucky
      track <- crop_range_res(track)
    }
  }
  
  
  if (updatedstudy=="Roshier") {
    if (id%in%c("Kev","Mem")) {
      track <- crop_range_res(track)
    }
    if (id=="Ada"){ 
      track <- crop_range_res(track) %>%
        filter(location.long < NA) #coords withheld
    }
  }
  
  if (updatedstudy=="Sekercioglu") {
    if (id%in%c("Bilge")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Serieys") {
    if (id=="Titan") {
      track <- crop_range_res(track)
    }
    
    if(id=="Xolani") {
      track <- track %>%
        filter(location.lat < NA) #coords withheld
    }
  }
  
  
  if (updatedstudy=="Serieys.AromasHills") {
    if (id%in%c("Bobcat_B31M","Bobcat_B35F","Bobcat_B38M")) {
      track <- crop_range_res(track)
    }
  }
  
  
  if (updatedstudy=="Serieys.CoyoteValley") {
    if (id=="B23M_5624") {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="SimonTrinzenHerrmannGotz") {
    if (id%in%c("79","130","127","92")) {
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
    
    if (id=="Jungle cat 06 (Bubbly)") {
      track <- track %>%
        filter(location.long > NA) #coords withheld
    }
    
    if (id=="Fox 09 (Broken Tail)") {
      track <- crop_range_res(track)
    }
    
    if (id=="Fox 20 (Bijli)") {
      track <- track %>%
        filter(location.lat > NA) #coords withheld
    }
    
  }
  
  if (updatedstudy=="Walton") {
    if (id%in%c("Mattias","Viktor","Oskar_1")) {
      track <- crop_range_res(track)
    }
    
    if (id=="Bengt") {
      track <- track %>%
        filter(timestamp<ymd_hms("2015-02-02 07:39:00"))
    }
    if (id=="Spank") {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Wheeldon") {
    if (id%in%c("PEC053","PEC059","PEC070","PEC082","PEC105",
                "PEC112","PEC114")) {
      track <- crop_range_res(track)
    }
    
    if(id=="PEC141") {
      track <- crop_range_res(track)
    }
    
  }
  
  if (updatedstudy=="Wilmers") {
    if (id%in%c("Astrid","Charlotte","Hedley","19F","35M","43F","54M","66M","95F","Chunga")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Young") {
    if (id%in%c("C037","B004","B009")) {
      track <- crop_range_res(track)
    }
  }
  
  return(track)
}

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

### LINKS TO SPATIAL DATA

#rougness: 
#link to dataset: http://www.earthenv.org/topography
#source: https://doi.org/10.1038/sdata.2018.40
#how to get: download from link with specifications:
#Dataset: "Roughness", Aggregation: "Median", Sources: "GMTED 2010", Resolution: "1 km"

#tree cover:
#link to dataset: http://www.earthenv.org/landcover
#source: https://doi.org/10.1111/geb.12182
#how to get: download cover classes 1-4 (needleleaf trees, evergreen broadleaf, deciduous broadleaf, others) and
#sum all all rasters to get total tree cover.

#human footprint index
#link to dataset: https://mountainscholar.org/handle/10217/216207
#source: https://doi.org/10.1088/1748-9326/abe00a
#how to get: go to source and download 2019 map

#road cover from merged osm/groads
#link to dataset: https://github.com/scabecks/humanfootprint_2000-2013/blob/master/spatial_data/pressure_layers/roads_osm_groads_union.tif
#source: https://doi.org/10.1016/j.oneear.2020.08.009
#how to get: download from link and reproject into longlat
#raster(terra::project(terra::rast("roads_osm_groads_union.tif"), crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")))

#dynamic habitat index (gross primary productivity)
#link to dataset: https://silvis.forest.wisc.edu/data/dhis/
#source: https://doi.org/10.1016/j.rse.2018.12.009
#how to get: read in as stack and get third layer (this one corresponds to variation/seasonality)
#dhi_gpp <- stack("gpp_dhi_combined-v5/dhi_gppqa_f.tif") 
#dhi_gpp <- dhi_gpp[[3]]

data_list <- alldata %>%
  group_split(id) #split into lists for each individual

vals <- list()
names <- c()

for (i in seq_along(data_list)) {
  
  ind <- data_list[[i]] %>%
    mod_to_tracks()
  
  id <- ind$id[[1]]
  
  if (length(ind$timestamp) < 2) { next } #skip objects with too few data
  
  ind_tel <- ind %>%
    as.telemetry() 
  
  coordinates(ind)<- ~location.long+location.lat
  crs(ind) <- ind_tel@info$projection #set projection from telemetry object
  
  #plot(ind)
  #ind@bbox #automatically has bbox? that's useful
  
  #get vals
  lat <- c(ind@bbox[[2]], ind@bbox[[4]])
  lon <- c(ind@bbox[[1]], ind@bbox[[3]])
  
  #make shapefile of bounding box
  minbox <- data.frame(lon, lat) %>%
    st_as_sf(coords = c("lon", "lat")) %>% 
    st_bbox() %>% 
    st_as_sfc() %>%
    as(Class="Spatial") #convert to spatialpolygons df
  
  ind_extent <- raster::extent(minbox)
  print(ind_extent)
  
  #if the individual tracks do not line up, return null for raster
  if(is.null(intersect(ind_extent, extent(curr_rast)))==TRUE) {
    
    next
    
  }
  
  ind_raster <- raster::crop(curr_rast, ind_extent)
  indvals <- getValues(ind_raster)
  
  #update df
  vals = c(vals, list(indvals))
  names = c(names, id)
  
}

names(vals) <- names

#get mean value over each individual extent
means <- lapply(vals, mean, na.rm=TRUE) %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("id")
