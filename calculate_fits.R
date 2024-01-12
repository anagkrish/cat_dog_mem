#part 1 of calculating ridge densities; calculates model fits on ALL individuals (this is the longest part)
#meant to be run in parallel for full effect

#load all libraries
library(tidyverse)
library(ctmm)
library(parsedate)
library(lubridate)
library(snow)
library(doSNOW)
library(geosphere)
#library(Rmpi)

"%ni%" <- Negate("%in%")

#load object alldata (all movement data for 1528 individuals)
load("movement/data/here")

#list of individuals/species/study names -- used to filter individuals
updatedstudy <- read_csv("updated/study/here")

#########################################

#set up dataframe
df <- data.frame(matrix(ncol = 24, nrow = 0))

cols <- c("ID", "Species", "Study", "Updated Study", "Model", "Min Sampling", "Duration", "Year",
	  "Effective_Sample_Size_Model",  "Effective_Sample_Size_AKDE",
	  "Home_Range_Area_Model (lower)", "Home_Range_Area_Model (mean)", "Home_Range_Area_Model (upper)",
	  "Home_Range_Area_AKDE (lower)", "Home_Range_Area_AKDE (mean)", "Home_Range_Area_AKDE (upper)", 
	  "Tau_p (lower)", "Tau_p (mean)", "Tau_p (upper)","Tau_v (lower)", "Tau_v (mean)", "Tau_v (upper)", 
          "Sigma_P", "Ballistic_Length_Scale")

colnames(df) <- cols

#########################################
#space to clean up select individuals and clean up tracks and whatnot
study <- c() #list studies of interest here
inds <- c() #list individuals of interest here

#select inds to fit (can go by study?)
ids <- updatedstudy %>% #filter by whatever criteria
	#filter(`Updated Study`%in%study) %>%
	#filter(ID%in%inds)
	unite("id", c(ID, Species), sep="_") %>%
	dplyr::select("id")

print(ids)

#########################################
#get normal tracks
tracks <- alldata %>%
  dplyr::select(timestamp, `location.long`,`location.lat`, `individual.local.identifier`,
                `individual.taxon.canonical.name`, `study.id`) %>%
  unite("id", c(`individual.local.identifier`, `individual.taxon.canonical.name`), sep="_", remove=F) %>%
  filter(id%in%ids$id) %>%
  group_by(paste(`individual.local.identifier`, `individual.taxon.canonical.name`)) %>%
  group_split

##########################################

print(length(tracks))

n <- length(tracks)
k <- 50 ## your LEN
sublists <- split(tracks, rep(1:ceiling(n/k), each=k)[1:n])

##########################################
#MAIN FUNCTION (to be run with clusterApply)

get_stats <- function(track) {  

  id <- unique(track$"individual.local.identifier")[[1]]
  species <- unique(track$"individual.taxon.canonical.name")[[1]]
  study <- unique(track$"study.id")[[1]]
  
  write(id, stdout())
  write(species, stdout()) #get species as well
  write(study, stdout())

  updatedstudy <- filter(updatedstudy, ID==id, Species==species, Study==study)$`Updated Study`
  print(paste(id, species, study, updatedstudy))

  write(updatedstudy, stdout())

  #sub-functions 
  
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
               timestamp_rs = ymd_hms(timestamp_rs)) %>%
        #collapse to individual by resampled time points but save og time points
        distinct(timestamp_rs, .keep_all = "TRUE") %>%
        dplyr::select(-c("timestamp")) %>%
        rename("timestamp" = timestamp_rs) #added to resample sketchy individuals
      
      return(track_rs)
    }
    
    else {
      return(track)
    }
    
  }

  #split track into two "ranges" for individuals w seasonal/spatiotemporal shifts that are basically two home ranges
  split_for_mean <- function(track, dt) {
 
     track1 <- track %>%
        filter(timestamp < ymd_hms(dt)) %>%
        mutate(individual.local.identifier=paste(individual.local.identifier, "a", sep="_"))

     track2 <- track %>%
        filter(timestamp > ymd_hms(dt)) %>%
        mutate(individual.local.identifier=paste(individual.local.identifier, "b", sep="_"))

      return(rbind(track1, track2))

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

  #function to get all mod fits so that you can rerun as many times as need be
  mod_stats <- function(tel,id,species,study) { #input telemetry object

    min_sampling <- get_minsampling(dat)/60 #outputs in seconds, divide by 60 to get minutes
  
    #save image of raw tracks
    png(file=paste("rawtracks/", str_replace(id, "/", ""),
                 species,study,".png",sep=""))
    plot(dat, col=color(dat, by="time"))
    dev.off()
    
    get_duration <- function(individual) { #individual name
	    
	    duration <- as.duration(interval(as.Date(individual$timestamp[1], #change to 2 for inds where first val is 00:00:00
                                            format = "%Y-%m-%d %H:%M:%S"),
                                   as.Date((tail(individual,1)$timestamp),
                                            #change to tail(individual,2)$timestamp[1] for
                                            #inds where LAST val is 00:00:00
                                            format='%Y-%m-%d %H:%M:%S')))@.Data/604800

          #gets duration in secs and converts to weeks (604800 seconds in a week)
          year <- as.Date(individual$timestamp[1], format = "%Y-%m-%d %H:%M:%S")

          return(list(duration,year))
    }

    times <- get_duration(dat)
    duration <- times[[1]]
    year <- as.character(times[[2]])
    
    print(paste(id,species,"got durations/minsampling"))

  SVF <- variogram(dat)
  GUESS <- ctmm.guess(dat,variogram=SVF,interactive=F)
  
  ## to fit error model for resampled individuals and IID fits, uncomment below
  #uere(dat) <- 10
  #GUESS <- ctmm.guess(dat,CTMM=ctmm(error=TRUE),interactive=F)
  
  #fit model with default CTMM parameters (phREML)
  ALL_FITS <-
    ctmm.select(
      dat,
      GUESS,
      verbose = TRUE,
      level = 0.95, #non-default argument used to speed up model fitting for 1500+ individuals (made some fits unstable, see methods for more info)
      cores = 2
    )

  DATA <- dat
  # Finds best model and gets its name
  BEST_FIT <- ALL_FITS[[1]]

  #save fit from each individuals!! important for later
  save(BEST_FIT, file=paste("allfits/study/", str_replace(id, "/", ""),
                            species,study,".rda",sep=""))

  print(paste(id,species,"got fits"))

  best_model_name <- summary(BEST_FIT)$name[[1]]

  #units = FALSE ensures area is m^2, and tau measurements are in seconds
  BEST_FIT_DATA <- summary(BEST_FIT, units = FALSE)

  if (best_model_name == "inactive") {
    return(NULL)
  }

  #save variogram plot
  xlim=c(0, tail(SVF$lag,1))

  print(xlim)

  png(file=paste("variograms/", str_replace(id, "/", ""),
                 species,study,".png",sep=""))
  plot(SVF,BEST_FIT, xlim=xlim)
  dev.off()

  UD <- NA
  
  #AKDEc (default) with debias OFF and weights OFF
  UD <- akde(DATA,BEST_FIT,debias=FALSE,weights=FALSE,trace=TRUE)
  AKDE_DATA <- summary(UD, units = FALSE)

  # Sets up variables to be used in the dataframe
  # (vals will depend on which best model was selected)
  low_area_model <- NA
  est_area_model <- NA
  high_area_model <- NA

  low_area_akde <- NA
  est_area_akde <- NA
  high_area_akde <- NA

  low_tau_p <- NA
  est_tau_p <- NA
  high_tau_p <- NA

  low_tau_v <- NA
  est_tau_v <- NA
  high_tau_v <- NA
  sigma_p <- NA

  # ballistic length scale: sqrt(est_tau_v/est_tau_p * sigm_p)
  b_l_s <- NA

  # Set values
  effec_samp_size_model <- BEST_FIT_DATA$DOF["area"][[1]]
  effec_samp_size_akde <- AKDE_DATA$DOF["area"][[1]]

  low_area_model <- BEST_FIT_DATA$CI[1, ][[1]]
  est_area_model <- BEST_FIT_DATA$CI[1, ][[2]]
  high_area_model <- BEST_FIT_DATA$CI[1, ][[3]]

  low_area_akde <- AKDE_DATA$CI[1, ][[1]]
  est_area_akde <- AKDE_DATA$CI[1, ][[2]]
  high_area_akde <- AKDE_DATA$CI[1, ][[3]]

  if ((best_model_name == "OU anisotropic") || (best_model_name == "OU")) {
    low_tau_p <- BEST_FIT_DATA$CI[2,][[1]]
    est_tau_p <- BEST_FIT_DATA$CI[2,][[2]]
    high_tau_p <- BEST_FIT_DATA$CI[2,][[3]]
  } else if ((best_model_name == "OUF anisotropic") || (best_model_name == "OUF")) {

    low_tau_p <- BEST_FIT_DATA$CI[2,][[1]]
    est_tau_p <- BEST_FIT_DATA$CI[2,][[2]]
    high_tau_p <- BEST_FIT_DATA$CI[2,][[3]]

    low_tau_v <- BEST_FIT_DATA$CI[3,][[1]]
    est_tau_v <- BEST_FIT_DATA$CI[3,][[2]]
    high_tau_v <- BEST_FIT_DATA$CI[3,][[3]]
    sigma_p <- ctmm:::area.covm(BEST_FIT$sigma)

    # ballistic length scale: sqrt(tau_v/tau_p * sigm_p)
    b_l_s <- sqrt((est_tau_v/est_tau_p) *  sigma_p)
  }

  returned <- c(id, species, study, updatedstudy, best_model_name, min_sampling, duration, year,
                effec_samp_size_model,  effec_samp_size_akde,
                low_area_model, est_area_model, high_area_model,
                low_area_akde, est_area_akde, high_area_akde,
                low_tau_p, est_tau_p, high_tau_p, low_tau_v, est_tau_v,  high_tau_v, sigma_p, b_l_s)

  return(returned)
  
 }

#THESE LINES EDIT INDIVIDUALS THAT WERE FLAGGED AS PROBLEMATIC BUT FIXEABLE
############# below are lines to edit problematic individuals as needed
if (updatedstudy=="Belant-Beyer") {
if (id%in%c("W07")) {
    track <- crop_range_res(track)
}
}

if (updatedstudy=="Berteaux") {
if (id%in%c("BORR","BVOB","JVOJ","OBBB","JVMJ")) {
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
    track <- split_for_mean(dt="2014-11-30 00:15:00")
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
}

if (updatedstudy=="Ferrell") {
  if (id%in%c("B21")) {
      track <- crop_range_res(track)
  }
  if (id=="B42") {
    track <- track %>%
      slice(-c(53, 179)) #edits for speed calc
  }
}

if (updatedstudy=="Fryxell") {
  if (id=="148.71") {
    track <- track %>%
      split_for_mean(dt="2010-07-22 22:00:34")
  }
  
  #iid inds added nov 27 23
  if (id=="148.66") {
    track <- crop_range_res(track)
  }
  
  if (id=="148.81") {
    track <- track %>%
      mutate(timestamp=ymd_hms(timestamp)) %>%
      filter(location.long > NA, 
             location.long < NA) %>% #coords withheld
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
if (id%in%c("Tsetseg")) {
    track <- crop_range_res(track)
}
if(id=="Killi") {
   track <- track %>%
	   split_for_mean(dt="2016-07-15 08:01:00")
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
if (id%in%c("30822")) {
    track <- crop_range_res(track)
}
if (id=="30839") {
  track <- track %>%
    mutate(timestamp = ymd_hms(timestamp)) %>%
    filter(timestamp > ymd_hms("2011-08-01 00:03:07"))
}
if (id=="31754") {
  track <- track %>%
	  split_for_mean(dt="2012-02-01 05:02:20")
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

if (updatedstudy=="Young") {
if (id%in%c("C037","B004","B009")) {
    track <- crop_range_res(track)
}

if(id=="B007") {
   track <- track %>%
           split_for_mean(dt="2014-09-20 02:01:20")
}
}
  
if (updatedstudy=="Mannil-Kont") {
  if (id=="450") {
    track <- track %>%
      split_for_mean(dt="2017-09-01 09:01:11")
  }
}

if (updatedstudy=="Morato") {
if (id%in%c("Panthera_114","Panthera_11","Panthera_13","Panthera_27","Panthera_34",
	    "Panthera_44")) {
    track <- crop_range_res(track)
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
	   split_for_mean(dt="2013-07-10 04:00:00")
}

if (id=="Panthera_82") {
   track <- track %>%
	   split_for_mean(dt="2015-05-01 04:29:00")
}

if (id=="Panthera_106"){
  track <- track %>%
          filter(timestamp<=ymd_hms("2009-11-21 10:00:00"))
}

if (id=="Panthera_111"){
  track <- track %>%
          filter(timestamp<=ymd_hms("2009-11-18 04:00:00"))
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
  return(c(id, species, study, updatedstudy, "duplicate!", NA, NA, NA, 
	 NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
  break
#skip Cerdocyon thous duplicates in OS dataset1
}

if (id%in%c("CAN13","CAN49","Grupo005_Id001","Grupo005_Id032", "shack")) {
    track <- crop_range_res(track)
  }
}

if (updatedstudy=="Oliveira-Santos") { #originally Oliveira-Santos.Dataset1
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

if (updatedstudy=="Prugh") {
if (id%in%c("C_MVBOB91M", "C_NEBOB13F", "C_MVBOB80M", "C_MVBOB66M", "C_NEBOB35M", "C_NEBOB5M",
	    "C_MVBOB67M", "C_NEBOB11M", "C_MVBOB87M", "C_MVBOB99F", "C_NEBOB6F", "C_MVBOB88M",
	    "C_MVBOB85F", "C_MVBOB51M", "C_MVBOB62M", "C_MVBOB69F", "C_MVBOB90M", "C_NEBOB7F",
	    "C_MVBOB76M", "C_NEBOB33M", "C_NEBOB8M", "C_MVBOB77M", "C_NEBOB38M", "C_MVBOB83M",
	    "C_NEBOB45M", "C_MVBOB52M", "C_NEBOB25F", "C_NEBOB37M", "C_NEBOB41F", "C_NEBOB16M",
	    "C_NEBOB32F", "C_NEBOB10F")) {
  return(c(id, species, study, updatedstudy, "duplicate!", NA, NA, NA,
         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
  break #skip duplicated coyotes
}

if (id%in%c("B_MVBOB66M","B_MVBOB54F","B_NEBOB33M","B_NEBOB35M","C_NECOY20F")) {
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

if (id=="C_NECOYaF") {
    track <- track %>%
      split_for_mean(dt="2019-01-31 20:00:00")
}

}

if (updatedstudy=="Patterson.BD") {
  if (id=="Diana") {
    track <- track %>%
      slice(-c(66,69,70,72))
  }
}
  
if (updatedstudy=="Patterson.B") {

if (id%in%c("W97155","W97143","W97177", "W97189","W97311","W97359",
	    "H01","T21","W1605","W215")) {
    track <- crop_range_res(track)
}
  
if(id=="W97188") {
  track <- track %>% crop_range_res() %>% slice(-c(103)) #edits to speed
}
  
  if(id=="W97363") {
    track <- track %>% slice(-c(70, 2421, 3211)) #edits to speed
  }
  
if(id=="W97303") {
   track <- track %>%
	   split_for_mean(dt="2012-04-07 04:59:00")
}
if(id=="T50") { 
   track <- track %>%
	   filter(location.lat > NA) #coords withheld 
}
if(id=="T02") {
   track <- track %>%
	   split_for_mean(dt="2005-12-01 23:33:00")
}
  
if (species=="Canis lupus x lycaon") {
  if (id %in% c("W97104","W97107","W97110","W97141","W97142","W97147",
		"W97162","W97303","W97305","W97306","W97356","W97362")) {
	  return(c(id, species, study, updatedstudy, "duplicate!", NA, NA, NA,
         NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
	  break #skip old lupus x lycaon individuals (but keep the updated lupus inds)
	}
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
if (id%in%c("Titan")) {
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
if (id%in%c("Bobcat_B31M","Bobcat_B38M","Bobcat_B35F")) {
    track <- crop_range_res(track)
}
}


if (updatedstudy=="Serieys.CoyoteValley") {
if (id%in%c("B23M_5624")) {
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
if (id%in%c("Jackal 02 (Zoom)")) {
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
}

if (updatedstudy=="Wheeldon") {
if (id%in%c("PEC059","PEC053","PEC070","PEC082","PEC105",
	    "PEC112","PEC114")) {
    track <- crop_range_res(track)
}
if(id=="PEC141") {
    track <- crop_range_res(track) %>%
	    split_for_mean(dt="2012-12-10 00:01:01")
}
if (id=="PEC009") {
    track <- track %>%
	    split_for_mean(dt="2010-12-01 15:02:45")
}
}

if (updatedstudy=="Wilmers") {
if (id%in%c("Astrid","Charlotte","Hedley","95F","35M","19F","43F","54M","66M")) {
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
if (id=="56M") {
    track <- track %>%
	    split_for_mean(dt="2015-09-01 22:00:42")
}
if (id=="31M") {
    track <- track %>%
	    split_for_mean(dt="2012-09-10 22:01:38")
}
if (id=="42M") {
    track <- track %>%
	    split_for_mean(dt="2015-03-01 22:00:42")
}
if (id=="16M") {
    track <- track %>%
	    split_for_mean(dt="2011-05-30 22:00:53")
}
}

#############
#for inds where you need to temporally split
print(paste(id,length(as.telemetry(track)),sep=" "))

if (length(as.telemetry(track)) == 2) { #very specific but allows for selection of individuals that have split tracks (hopefully we won't have to deal with cases upwards of 2)
  returned <- list() #initialize list to output both halves of model fit
  
  #merge to get one track; allows for calculation of og stats
  #and can get projection of total data to put individual splits in same projection (easier to visualize when plotting)
  tot_track <- track %>%
	 mutate(individual.local.identifier=str_replace(individual.local.identifier, "_a", "")) %>%
	 mutate(individual.local.identifier=str_replace(individual.local.identifier, "_b", ""))

  #save variogram and tracks for total ind before splits (for comparison purposes)
  #so each ind in this category should have three tracks/vgs/fits associated with
  if(get_minsampling(as.telemetry(tot_track)) < 60) {

     print(paste(id,species,"refitting"))

    dat <- as.telemetry(resample(tot_track), timeformat = "auto", timezone = "UTC", projection = NULL, datum = NULL,
                        timeout = Inf, na.rm = "row", mark.rm = FALSE, keep = FALSE, drop = TRUE)

    print("resampled!")
    print(head(dat))

  } else {
          print(paste(id, species,"not refitting"))

           dat <- as.telemetry(tot_track, timeformat = "auto", timezone = "UTC", projection = NULL, datum = NULL,
                        timeout = Inf, na.rm = "row", mark.rm = FALSE, keep = FALSE, drop = TRUE)

          print(head(dat))

  }
  returned <- append(returned, list(mod_stats(dat, id, species, study)))
  print(paste("split tot track:", returned, sep="/n"))

  #then split ind and continue
  track <- track %>%
    group_by(`individual.local.identifier`) %>%
    group_split

  for (i in seq_along(track)) {
    split <- track[[i]]
    id <- unique(split$"individual.local.identifier")[[1]]

    if(get_minsampling(as.telemetry(split)) < 60) {

       print(paste(id,species,"refitting"))

       dat <- as.telemetry(resample(split), timeformat = "auto", timezone = "UTC", projection = NULL, datum = NULL,
                        timeout = Inf, na.rm = "row", mark.rm = FALSE, keep = FALSE, drop = TRUE)
    
      print("resampled!")
      print(head(dat))

    } else {
          print(paste(id, species,"not refitting"))

	   dat <- as.telemetry(split, timeformat = "auto", timezone = "UTC", projection = NULL, datum = NULL,
                        timeout = Inf, na.rm = "row", mark.rm = FALSE, keep = FALSE, drop = TRUE)

          print(head(dat))

  }

  ctmm::projection(dat) <- ctmm::projection(as.telemetry(tot_track)) #set split projection to projection calculated above 
  returned <- append(returned, list(mod_stats(dat, id, species, study)))
}

print(paste("split track:", returned, sep="/n"))
return(returned)

} else {

#for all inds that are not being split
#resample if minsampling less than a minute
if(get_minsampling(as.telemetry(track)) < 60) {

    print(paste(id,species,"refitting"))
    
    dat <- as.telemetry(resample(track), timeformat = "auto", timezone = "UTC", projection = NULL, datum = NULL,
			timeout = Inf, na.rm = "row", mark.rm = FALSE, keep = FALSE, drop = TRUE)
    
    print("resampled!")
    print(head(dat))

  } else {
      
      print(paste(id, species,"not refitting"))
      
      dat <- as.telemetry(track, timeformat = "auto", timezone = "UTC", projection = NULL, datum = NULL,
                        timeout = Inf, na.rm = "row", mark.rm = FALSE, keep = FALSE, drop = TRUE)
      
      print(head(dat))

  }

  returned <- list(mod_stats(dat, id, species, study))

  print(paste("regular track:", returned, sep="/n"))
  return(returned)

}
}

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
snow::clusterExport(cl, list= c("updatedstudy"), envir=environment())

ptm <- proc.time()

print(sublists)
for (curr_list in sublists) {
 
  result <- snow::clusterApply(cl = cl, curr_list, get_stats)
  print(length(result))
  print(length(result[[1]]))
  
  #edits, allows u to deal with split inds with list() results outputs
  for (i in seq_along(result)) {
    if (length(result[[i]])>1) {
      for (j in seq_along(result[[i]])) {
      	      print(result[[i]][[j]])
	final_list <- append(final_list,list(result[[i]][[j]]))
    }
  }
  else { 
    final_list <- append(final_list, list(result[[i]][[1]])) 
  }
  }

}

print(final_list)

proc.time()-ptm

print("cluster done")

print("writing to rda")

for(i in 1:length(final_list)) {
  df <- structure(rbind(df, final_list[[i]]), .Names = names(df))
}

write.csv(df, "model_fit_output.csv")

print("done!")

Rmpi::mpi.close.Rslaves()

