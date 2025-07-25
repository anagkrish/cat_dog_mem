#part 2 of calculating ridge densities; uses saved model fits from calculate_fits.R to get ridge densities

#load libraries
library(tidyverse)
library(ctmm)
library(lubridate)
library(stringr)
library(geosphere)
library(raster)
library(sp)
library(sf)
library(grDevices)
library(snow)
library(doSNOW)
#library(Rmpi)

"%ni%" <- Negate("%in%")

#load movement object alldata for all 1528 individuals
#note: movement data used in this manuscript is archived in Movebank. see manuscript for more details.
load("movement/data/here")

#list of individuals/species/study names -- used to filter individuals
updatedstudynames <- read_csv("updated/study/here")

#########################################

fitsfolder <- "/all/fits/final" #name of folder where all model fits are stored (saved in calculate_fits.R)

ALLRIDGE <- data.frame("id"=NA, "sp"=NA, "study"=NA,
                       "updatedstudy"=NA,
                       "mod_name"=NA, "n_fixes"=NA,
                       "duration"=NA, "year"=NA,
                       "resampled?"=NA, "minsampling"=NA, "ess"=NA,
                       "speed_low"=NA,"speed_est"=NA,"speed_high"=NA,
                       "mean_dist_day_low"=NA, "mean_dist_day_est"=NA, "mean_dist_day_high"=NA,
                       "area_mod_low"=NA, "area_mod_est"=NA, "area_mod_high"=NA,
                       "area_ud_low"=NA, "area_ud_est"=NA, "area_ud_high"=NA,
                       "tau_p_low"=NA, "tau_p_est"=NA, "tau_p_high"=NA,
                       "tau_v_low"=NA, "tau_v_est"=NA, "tau_v_high"=NA,
                       "sigma_p"=NA, "b_l_s"=NA,
                       "ridge_length_low"=NA, "ridge_length_est"=NA, "ridge_length_high"=NA,
                       "ridge_dens_low"=NA, "ridge_dens_est"=NA, "ridge_dens_high"=NA) %>% drop_na()


#FUNCTIONS ####################################################################

calc_from_fits <- function(file) {

mod_to_tracks <- function(track) {
  
  id <- unique(track$"individual.local.identifier")[[1]]
  species <- unique(track$"individual.taxon.canonical.name")[[1]]
  study <- unique(track$"study.id")[[1]]
  
  updatedstudy <- filter(updatedstudynames, ID==id, Species==species, Study==study)$`Updated Study`
  print(paste(id, species, study, updatedstudy))
  
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
          filter(location.long > NA, 
                 location.lat < NA) #coords withheld
      }
      
    }
    
    if(updatedstudy=="Clark") {
      if (id%in%c("L12","L27","R3","C263")) {
        track <- crop_range_res(track)
      }
      
      if (id=="L19") { 
        track <- track %>%
          split_for_mean(dt="2018-07-10 22:00:53")
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
    
    if (updatedstudy=="Conner") {
      if (id%in%c("F46","M21")) {
        track <- crop_range_res(track)
      }
      if (id=="M7") {
        track <- track %>%  split_for_mean(dt="2014-11-30 00:15:00")
      }
      if (id=="F30") {
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
      if (id=="BBJackal_Rooky") {
        track <- crop_range_res(track) %>% 
          split_for_mean(dt="2014-12-20 00:00:00")
      }
    }
  
  if (updatedstudy=="Frair") {
    if (id=="M5") {
      track <- track %>%
        filter(timestamp > ymd_hms("2006-01-01 18:28:59")) %>%
        crop_range_res()
    }
    if (id=="Gray1") {
      track <- track %>% slice(-c(210)) #edits for speed
    }
  }
  
  if (updatedstudy=="Ferrell") {
    if (id%in%c("B21")) {
      track <- crop_range_res(track)
    }
    
    if (id=="B42") {
      track <- track %>%
        slice(-c(53, 179)) #edits to get rid of unrealistic speed calculations
    }
    
  }
  
  if (updatedstudy=="Fryxell") {
    
    if (id=="148.71") {
        track <- track %>% 
          split_for_mean(dt="2010-07-22 22:00:34")
    }
    
    if (id=="148.66") {
      track <- crop_range_res(track)
    }
    
    if (id=="148.81") {
      track <- track %>%
        mutate(timestamp=ymd_hms(timestamp)) %>%
        filter(location.long > NA, 
               location.lat < NA) #coords withheld
        split_for_mean(dt="2010-10-15 18:00:00")
    }
    
    }
    
    if (updatedstudy=="Garthe") {
      if (id%in%c("REDFOX-2015-01-BHKoog",
                  "REDFOX-2016-03-BHKoog",
                  "REDFOX-2018-03-Sylt")) {
        track <- crop_range_res(track)
      }
      if(id=="RACDOG-2019-01-Sylt") {
        track <- track %>%
          split_for_mean(dt="2019-12-16 02:00:00")
      }
    }
    
    if (updatedstudy=="Getz-Bellan") {
      if (id%in%c("CM18","CM26","CM36","CM95")) {
        track <- crop_range_res(track)
      }
      
      if(id=="CM83") {
        track <- crop_range_res(track) %>%
          split_for_mean(dt="2010-05-01 00:00:45")
      }
    }
    
    if (updatedstudy=="Hebblewhite") {
      if (id%in%c("B065","J031")) {
        track <- crop_range_res(track)
      }
      if(id=="JW01"){ 
        track <- track %>%
          split_for_mean(dt="2010-12-01 11:00:55")
      }
    }
    
    if (updatedstudy=="Hernandez-Blanco") {
      if (id=="Tsetseg") {
        track <- crop_range_res(track)
      }
      if(id=="Safar") {
        track <- crop_range_res(track) %>%
          split_for_mean(dt="2012-01-31 22:02:00")
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
          filter(timestamp >= ymd_hms("2016-06-14 20:45:32"))
      }
      
    }
    
    if (updatedstudy=="Kays") {
      if (id=="30822") {
        track <- crop_range_res(track)
      }
      if (id=="30839") {
        track <- track %>%
	 mutate(timestamp = ymd_hms(timestamp)) %>%
         filter(timestamp > ymd_hms("2011-08-01 00:03:07"))
      }
      if (id=="31754") {
        track <- track %>%
          split_for_mean(dt="2012-01-15 05:02:20")
      }
    }
    
    if (updatedstudy=="VanDerWeyde-Kral") {
      if (id%in%c("Kealiboka","Matsoshetsi")) {
        track <- crop_range_res(track)
      }
    }
   
   if (updatedstudy=="Lantham") {
     if (id=="3") {
       track <- track %>%
         split_for_mean(dt="2006-11-15 01:01:47")
     }
   }

   if (updatedstudy=="Lang") {
      if (id=="Felis_200") {
        track <- track %>%
          split_for_mean(dt="2021-04-01 00:01:00")
      }
   }

  if (updatedstudy=="Mannil-Kont") {
    if (id=="450") {
      track <- track %>%
        split_for_mean(dt="2017-09-01 09:01:11")
    }
  }
  
    if (updatedstudy=="Morato") {
      if (id%in%c("Panthera_11", "Panthera_114", "Panthera_13", "Panthera_27",
		  "Panthera_34","Panthera_44")) {
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
          split_for_mean(dt="2013-02-28 11:00:00")
      }
      
      if (id=="Panthera_5") {
        track <- track %>%
          split_for_mean(dt="2009-08-29 02:09:00")
      }
      
      if (id=="Panthera_53") {
        track <- track %>%
          split_for_mean(dt="2013-06-10 04:00:00")
      }
      
      if (id=="Panthera_59") {
        track <- track %>%
          split_for_mean(dt="2013-06-30 07:01:00")
      }
      
      if (id=="Panthera_82") {
        track <- track %>%
          split_for_mean(dt="2014-10-01 04:29:00")
      }
      
      if (id=="Panthera_36"){
        track <- track %>%
          filter(timestamp<=ymd_hms("2000-10-18 15:50:00")) %>%
          crop_range_res()
      }
      
    }
    
    if (updatedstudy=="Azevedo.Lemos") { #originally Oliveira-Santos.Dataset1
      if (id%in%c("150011","150041","150102","150312","150402","150462",
                  "150552","150681","163181","164820","164886","164900",
                  "164957","1649671","165164","165194","1651941","165224","165252","1652521")) {
        break
        #skip Cerdocyon thous duplicates in OS dataset1
      }
      
      if (id%in%c("CAN13","CAN49", "Grupo005_Id001","Grupo005_Id032")) {
        track <- crop_range_res(track)
      }
      
      if (id=="shack") {
        track <- crop_range_res(track)
      }
      
    }
    
    if (updatedstudy=="Oliveira-Santos") { #originally Oliveira-Santos.Dataset2
      if (id%in%c("150041_GustavoCT")) {
        track <- crop_range_res(track)
      }
      if (id=="kayapo") {
        track <- track %>%
          filter(location.lat > NA) #coords withheld
      }
    }
    
    if (updatedstudy=="Palomares") { #originally grouped with palacios gonzales no idea why
      if (id%in%c("Garfio", "Patsuezo")) {
        track <- crop_range_res(track) 
      }
    }
         
  if (updatedstudy=="Palacios Gonzalez") { #originally grouped w palomares, no idea why
    if (id%in%c("Lynx_Ketamina", "Lynx_Llerena", 
                  "Lynx_Miera", "Lynx_Negral", "Lynx_Nitrógeno")) {
        track <- crop_range_res(track)
      }
      
      if (id=="Lynx_Mayo") {
        track <- track %>%
          split_for_mean("2016-06-01 05:00:00 ")
      }
      
      if (id=="Lynx_Neruda") {
        track <- crop_range_res(track) %>%
          split_for_mean("2017-08-31 23:25:00")
      }
      
    }
  
  if (updatedstudy=="Patterson.BD") {
    if (id=="Diana") {
      track <- track %>%
        slice(-c(66,69,70,72)) #edits to get rid of unrealistic speed calculations
    }
  }
  
  if (updatedstudy=="Patterson.B") {
    
    if (id%in%c("W97143","W97155","W97177","W97185","W97189","W97311","W97359",
                "H01", "T21", "W215", "W1605")) {
      track <- crop_range_res(track)
    }
    
    if(id=="W97188") {
      track <- track %>% crop_range_res() %>% slice(-c(103))
    }
    
    if(id=="W97363") {
      track <- track %>% slice(-c(70, 2421, 3211)) #edits to get rid of unrealistic speed calculations
    }
    
    if(id=="W97381") {
      track <- track %>% slice(-c(2493)) #edits to get rid of unrealistic speed calculations
    }
    
    if(id=="W97303") {
      track <- track %>%
        split_for_mean(dt="2012-04-07 04:59:00")
    }
    
    if(id=="T02") {
      track <- track %>%
        split_for_mean(dt="2005-12-01 23:33:00")
    }
    
    if(id=="T50") { 
      track <- track %>%
        filter(location.lat > NA) #coords withheld 
    }
   
    if (species=="Canis lupus x lycaon") {
      if (id %in% c("W97104","W97107","W97110","W97141","W97142","W97147",
                    "W97162","W97303","W97305","W97306","W97356","W97362")) {
        break #skip old lupus x lycaon individuals (but keep the updated lupus inds)
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
        break #skip duplicated coyotes
      }
    
      if (id%in%c("C_NECOY20F","B_MVBOB54F","B_MVBOB66M","B_NEBOB33M","B_NEBOB35M")) {
        track <- crop_range_res(track)
      }
      
      if (id=="C_MVCOY98M") { 
        track <- crop_range_res(track) %>%
          filter(timestamp < ymd_hms("2020-09-01 00:00:00"))
      }
      
      if (id=="C_NECOYaF") { 
        track <- track %>%
          split_for_mean(dt="2019-01-31 20:00:00")
      }
      
      if (id=="B_NEBOB23M") { 
        track <- crop_range_res(track) %>%
          filter(location.long < NA) #coords withheld
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
      if (id=="Azure") {
        track <- track %>%
          split_for_mean(dt="2016-06-05 00:00:00")
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
      if (id=="106") {
        track <- track %>%
          split_for_mean(dt="2017-04-30 01:00:49")
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
      
      if (id=="PEC009") {
        track <- track %>%
          split_for_mean(dt="2010-12-01 15:02:45")
      }
      
      if(id=="PEC141") {
        track <- crop_range_res(track) %>%
          split_for_mean(dt="2012-12-10 00:01:01")
      }
      
    }
    
    if (updatedstudy=="Wilmers") {
      if (id%in%c("Astrid","Charlotte","Hedley",
		  "19F","35M","43F","54M","66M","95F")) {
        track <- crop_range_res(track)
      }
      
      if(id=="Chunga") {
        track <- crop_range_res(track) %>%
          split_for_mean(dt="2016-02-28 02:30:00")
      }
      if (id=="Merimela") {
        track <- track %>%
          split_for_mean(dt="2015-01-28 02:30:00")
      }
      if (id=="16M") {
        track <- track %>%
          split_for_mean(dt="2011-05-30 22:00:53")
      }
      if (id=="31M") {
        track <- track %>%
          split_for_mean(dt="2012-09-10 22:01:38")
      }
      if (id=="42M") {
        track <- track %>%
          split_for_mean(dt="2015-03-01 22:00:42")
      }
      if (id=="56M") {
        track <- track %>%
          split_for_mean(dt="2015-09-01 22:00:42")
      }
    }
    
  if (updatedstudy=="Young") {
    if (id%in%c("C037","B004","B009")) {
      track <- crop_range_res(track)
    }
    
    if(id=="B007") {
      track <- track %>%
        split_for_mean(dt="2014-09-20 02:01:20")
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

#get duration
get_duration <- function(individual) { #individual name
  duration <- as.duration(interval(as.Date(individual$timestamp[1], #change to 2 for inds where first val is 00:00:00
                                           format = "%Y-%m-%d %H:%M:%S"),
                                   as.Date((tail(individual,1)$timestamp),
                                           #change to tail(individual,2)$timestamp[1] for
                                           #inds where LAST val is 00:00:00
                                           format='%Y-%m-%d %H:%M:%S')))@.Data/604800
  
  #gets duration in secs and converts to weeks (604800 seconds in a WEEK bitch)
  
  #duration <- parse_number(duration)
  year <- as.character(as.Date(individual$timestamp[1], format = "%Y-%m-%d %H:%M:%S"))
  
  return(list(duration,year))
  
}

#split track into two "ranges" for individuals w seasonal/spatiotemporal shifts that are basically two home ranges
split_for_mean <- function(track, dt) {
  
  #modified from fit split for mean function to just crop tracks (since we already have fits)
  if(grepl("_a", file)==TRUE) {
    
    track <- track %>%
      mutate(timestamp = ymd_hms(timestamp)) %>%
      filter(timestamp < ymd_hms(dt))
  }

  if(grepl("_b", file)==TRUE) {
    
    track <- track %>%
      mutate(timestamp = ymd_hms(timestamp)) %>%
      filter(timestamp > ymd_hms(dt))
    
  } else {
     
    track <- track

  }

  return(track)
  
}

#crop out extraneous points to home range
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

get_stats <- function(track, resamp, nametag) { #nametag is just a way to differentiate split tracks
  
  id <- unique(track$"individual.local.identifier")[[1]]
  species <- unique(track$"individual.taxon.canonical.name")[[1]]
  study <- unique(track$"study.id")[[1]]
  updatedstudy <- filter(updatedstudynames, ID==id, Species==species, Study==study)$`Updated Study`
  
  track <- track%>%
    as.telemetry()
  
  #BEST_FIT loaded as file
  dur <- get_duration(track)
  minsampling <- get_minsampling(track)/60 #convert to mins
  mod_name <- summary(BEST_FIT)$name[[1]]
  ess <- summary(BEST_FIT)$DOF["area"][[1]]
  n_fixes <- length(track$timestamp)
  
  print(paste("ok",id, species, study, updatedstudy))
  
  #akde and hr area (debias and weights off)
  UD <- akde(track, BEST_FIT, debias=FALSE, weights=FALSE)
  
  #general summary info (hr area, tau p, v and sigma)
  area_model <- summary(BEST_FIT, units=F)$CI[1,]
  area_ud <- summary(UD, units=F)$CI[1,]
  
  if ((mod_name == "OU anisotropic") || (mod_name == "OU") 
      || (mod_name == "OU anisotropic error") || (mod_name == "OU error")) {
    
    tau_p <- summary(BEST_FIT, units=F)$CI[2,]
    
    #so the summary function is returning values for OU/OU anisotropic tau_v and sigma 
    #but since they were excluded before I'm going to keep doing that unless told otherwise
    tau_v <- c(NA,NA,NA) 
    sigma_p <- NA
    b_l_s <- NA
    
    #speed is inf if not ouf/ouf anisotropic so manually insert na and ask Chris
    speed <- c(NA,NA,NA)
    dists <- c(NA,NA,NA)
  }
  
  else if ((mod_name == "OUF anisotropic") || (mod_name == "OUF")
           || (mod_name == "OUF anisotropic error") || (mod_name == "OUF error")
           || (mod_name == "OUf anisotropic") || (mod_name == "OUf")
           || (mod_name == "OUf anisotropic error") || (mod_name == "OUf error")) {
    
    tau_p <- summary(BEST_FIT, units=F)$CI[2,]
    tau_v <- summary(BEST_FIT, units=F)$CI[3,]
    sigma_p <- ctmm:::area.covm(BEST_FIT$sigma)
    
    b_l_s <- sqrt((tau_v[[2]]/tau_p[[2]]) *  sigma_p)
    
    speed <- speed(BEST_FIT, units=F)$CI
    
    #dist/day sloppy version
    #mod from https://doi.org/10.1186/s40462-019-0177-1 s2
    track$day <-cut(track$timestamp, breaks="day") 
    days <- unique(track$day)
    
    #An empty list to fill with the results
    dists <- list()
    
    #Loop over the number of days
    for(i in 1:length(days)){
      
      #select data for the day in question
      day.dat <- track[which(track$day == days[i]),]
      samp.time <- diff(c(day.dat$t[1], day.dat$t[nrow(day.dat)])) #get duration of sampling time
      ctmm_dist <- speed*samp.time # get the estimated distance travelled (in m) 
      
      x <- c(ctmm_dist[1], #min
             ctmm_dist[2], #est dist
             ctmm_dist[3]) #max
      names(x) <- c("min.dist", "est.dist", "max.dist")
      dists[[i]] <- x
      
    }
    
    dists <- as.data.frame(do.call(rbind , dists)) %>%
      colMeans #get mean dist travelled per day 
  }
  
  else {
    
    print("ok")
    #IID/IID error 
    
    #tau_p not in IID
    tau_p <- c(NA,NA,NA)
    
    #so the summary function is returning values for OU/OU anisotropic tau_v and sigma 
    #but since they were excluded before I'm going to keep doing that unless told otherwise
    tau_v <- c(NA,NA,NA)
    sigma_p <- NA
    b_l_s <- NA
    
    #speed is inf if not ouf/ouf anisotropic so manually insert na and ask Chris
    speed <- c(NA,NA,NA)
    dists <- c(NA,NA,NA)
    
  }
  
  #calculate ridge density!
  RIDGE <- ctmm:::ridges.UD(UD)
  
  ridge_indicator <- RIDGE$Indicator #original ridge matrix
  avg_width <- RIDGE$Ave.Width #avg width used to calculate percent ridge use probability
  
  bandwidth <- sqrt(eigen(UD$H)$values) #chris is suspicious that avg width is related to bandwidth 
  
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
  
  hr_width <- sqrt(BEST_FIT$sigma@par['minor']) #hr width for buffer
  
  #creates buffer around ridges and calculates how much the individual actually uses those tracks
  
  buffer <- sf::st_buffer(ridges_mean, dist=avg_width) #buffer using avg width
  buffer_reduced10 <- sf::st_buffer(ridges_mean, dist=avg_width/10) #buffer using avg width
  buffer_reduced100 <- sf::st_buffer(ridges_mean, dist=avg_width/100) #buffer using avg width
  
  buffer_frac5 <- sf::st_buffer(ridges_mean, dist=(hr_width*0.05)) #buffer using fraction of hr width
  buffer_frac10 <- sf::st_buffer(ridges_mean, dist=(hr_width*0.1)) #buffer using fraction of hr width
  buffer_frac25 <- sf::st_buffer(ridges_mean, dist=(hr_width*0.25)) #buffer using fraction of hr width
  
  buffer_fixed10 <- sf::st_buffer(ridges_mean, dist=10) #buffer using fixed width
  buffer_fixed50 <- sf::st_buffer(ridges_mean, dist=50) #buffer using fixed width
  buffer_fixed100 <- sf::st_buffer(ridges_mean, dist=100) #buffer using fixed width
  buffer_fixed200 <- sf::st_buffer(ridges_mean, dist=200) #buffer using fixed width
  
  pmf <- ctmm::raster(UD, DF="PMF") #export probability mass function
  pmf_mask <- raster::mask(pmf, buffer)
  pmf_mask_frac5 <- raster::mask(pmf, buffer_frac5)
  pmf_mask_frac10 <- raster::mask(pmf, buffer_frac10)
  pmf_mask_frac25 <- raster::mask(pmf, buffer_frac25)
  pmf_mask_reduced10 <- raster::mask(pmf, buffer_reduced10)
  pmf_mask_reduced100 <- raster::mask(pmf, buffer_reduced100)
  pmf_mask_fixed10 <- raster::mask(pmf, buffer_fixed10)
  pmf_mask_fixed50 <- raster::mask(pmf, buffer_fixed50)
  pmf_mask_fixed100 <- raster::mask(pmf, buffer_fixed100)
  pmf_mask_fixed200 <- raster::mask(pmf, buffer_fixed200)
  
  ridge_probability <- raster::cellStats(pmf_mask, "sum")
  ridge_prob_frac5 <- raster::cellStats(pmf_mask_frac5, "sum")
  ridge_prob_frac10 <- raster::cellStats(pmf_mask_frac10, "sum")
  ridge_prob_frac25 <- raster::cellStats(pmf_mask_frac25, "sum")
  ridge_prob_reduced10 <- raster::cellStats(pmf_mask_reduced10, "sum")
  ridge_prob_reduced100 <- raster::cellStats(pmf_mask_reduced100, "sum")
  ridge_prob_fixed10 <- raster::cellStats(pmf_mask_fixed10, "sum")
  ridge_prob_fixed50 <- raster::cellStats(pmf_mask_fixed50, "sum")
  ridge_prob_fixed100 <- raster::cellStats(pmf_mask_fixed100, "sum")
  ridge_prob_fixed200 <- raster::cellStats(pmf_mask_fixed200, "sum")
  
  returned <- c(paste(id, nametag, sep=""), 
                species, study, updatedstudy, mod_name, n_fixes,
                dur[1][[1]], dur[2][[1]], resamp, minsampling, ess,
                speed[[1]], speed[[2]], speed[[3]],
                dists[[1]], dists[[2]], dists[[3]],
                area_model[[1]], area_model[[2]], area_model[[3]],
                area_ud[[1]], area_ud[[2]], area_ud[[3]],
                tau_p[[1]], tau_p[[2]], tau_p[[3]],
                tau_v[[1]], tau_v[[2]], tau_v[[3]],
                sigma_p, b_l_s,
                length_l, length_m, length_h,
                dens_l, dens_m, dens_h, avg_width, bandwidth[[1]], bandwidth[[2]], 
                ridge_probability, ridge_prob_frac5, ridge_prob_frac10, ridge_prob_frac25,
                ridge_prob_reduced10, ridge_prob_reduced100,
                ridge_prob_fixed10, ridge_prob_fixed50, ridge_prob_fixed100, ridge_prob_fixed200)
  
  print(returned)
  return(returned)
  
}
# actual loop to run ###################################################################

    print(file)
    load(paste(folder, "/", file, sep = ""))
    resamp <- "n" #keep track of which individuals were resampled

#  for (file in list.files(folder)) {
    if(grepl("_a", file)==TRUE) {
  
      #note if resampled
      if (get_minsampling(as.telemetry(alldata %>%
                                       filter(id==str_replace(str_replace(file, "_a", ""), ".rda$", "")) %>%
                                       dplyr::select(-c("id")))) < 60) { resamp <- "y" }
      
      track <- alldata %>%
        filter(id==str_replace(str_replace(file, "_a", ""), ".rda$", "")) %>%
        dplyr::select(-c("id")) %>%
        resample() %>% #resample if less than a minute minsampling
        mod_to_tracks()
      
      print(track)	
      returned <- get_stats(track, resamp, "_a")
      
      return(returned)
      #next
    }
    
    if (grepl("_b", file)==TRUE) {
      
      #note if resampled
      if (get_minsampling(as.telemetry(alldata %>%
                                       filter(id==str_replace(str_replace(file, "_a", ""), ".rda$", "")) %>%
                                       dplyr::select(-c("id")))) < 60) { resamp <- "y" }
      
      
      track <- alldata %>%
        filter(id==str_replace(str_replace(file, "_b", ""), ".rda$", "")) %>%
        dplyr::select(-c("id")) %>%
        resample() %>% #resample if less than a minute minsampling
        mod_to_tracks()
      
      print(track)
      return(get_stats(track, resamp, "_b"))

    }

    else {
      
      #note if resampled
      if (get_minsampling(as.telemetry(alldata %>%
                                       filter(id==str_replace(str_replace(file, "_a", ""), ".rda$", "")) %>%
                                       dplyr::select(-c("id")))) < 60) { resamp <- "y" }
      
      #select id by recreating og id and then dropping that column
      track <- alldata %>%
        mutate(id=str_replace(id, "/","")) %>% #specifically for serieys.coyotevalley B03F_5622/5974 and B05F_56285620
        filter(id==str_replace(file, ".rda$", "")) %>%
        dplyr::select(-c("id")) %>%
        resample() %>% #resample if less than a minute minsampling
        mod_to_tracks() #make necessary mods
      
      print(track)
      return(get_stats(track, resamp, ""))
      
    }

#}    

}

#runs everything but with a special case to print actual returned error (works better for parallel runs)
workerfun <- function(i) {
  tryCatch({
    calc_from_fits(i)
  },
  error=function(e) {
    print(e)
    stop(e)
  })
}

#####################################################################
#code to make it parallel
#actually run the code

final_list <- list()

#make it parallel!!!
np <- mpi.universe.size() - 1
cl <- snow::makeCluster(np, type = "MPI", outfile ="")

clusterEvalQ(cl, library("tidyverse"))
clusterEvalQ(cl, library("ctmm"))
clusterEvalQ(cl, library("lubridate"))
clusterEvalQ(cl, library("dplyr"))
clusterEvalQ(cl, library("geosphere"))
clusterEvalQ(cl, library("raster"))
clusterEvalQ(cl, library("sp"))
clusterEvalQ(cl, library("rgeos"))
clusterEvalQ(cl, library("grDevices"))
clusterEvalQ(cl, library("stringr"))

snow::clusterExport(cl, list= c("updatedstudynames"), envir=environment())
snow::clusterExport(cl, list= c("alldata"), envir=environment())
snow::clusterExport(cl, list= c("calc_from_fits"), envir=environment())

ptm <- proc.time()

#for 0.5 fitsfolder[c(5,6,8,9,17,19,25,26,33,36,39,40,46,48,50,52,54,56,57)]
#for 0.75 fitsfolder)[c(17,25,26,36,39,52)
#for OUf/OUf anisotropic inds[c(5,7,19,21,33,48,53)]

for (folder in list.dirs(fitsfolder)) { #each study is saved in its own folder; loop through each study/ individual

    print(folder)
    snow::clusterExport(cl, list= c("folder"), envir=environment())
    result <- snow::clusterApply(cl = cl, list.files(folder), workerfun) #calc_from_fits)
    final_list <- append(final_list, result)

}

for(i in 1:length(final_list)) {
  ALLRIDGE <- structure(rbind(ALLRIDGE, final_list[[i]]), .Names = names(ALLRIDGE))
}

proc.time()-ptm

#Rmpi::mpi.close.Rslaves()

print("cluster done")

#proc.time() - ptm
print("writing to rda")

write_csv(ALLRIDGE, file="ridgedatstats.csv")

print("done")

#output of this file is ridge.csv (see 'data' folder)

