#phenology plots

#load libraries
library(tidyverse)
library(ctmm)
library(lubridate)
library(parsedate)
library(dplyr)

alldata <- load("movement/data/here")

#set up dataframe -- count days for each individual
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
  
  if (updatedstudy=="Bertreaux") {
    if (id%in%c("BORR","BVOB","JVOJ","OBBB","ORRR")) {
      track <- crop_range_res(track)
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
    if (id%in%c("F46","M21")) {
      track <- crop_range_res(track)
    }
  }
  
  if (updatedstudy=="Critescu") {
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
        filter(location.lat > -32.845)
    }
  }
  
  if (updatedstudy=="Ferrell") {
    if (id%in%c("B21")) {
      track <- crop_range_res(track)
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
    if (id%in%c("CM18","CM26","CM36","CM95")) {
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
    if (id=="30839") {
      track <- track %>%
        filter(timestamp > ymd_hms("2011-08-01 00:03:07"))
    }
  }
  
  if (updatedstudy=="Kral") {
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
  
  if (updatedstudy=="Oliveira-Santos.Dataset1") {
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
  }
  
  if (updatedstudy=="Oliveira-Santos.Dataset2") {
    if (id%in%c("150041_GustavoCT")) {
      track <- crop_range_res(track)
    }
    if (id=="kayapo") {
      track <- track %>%
        filter(location.lat > -18.34)
    }
  }
  
  if (updatedstudy=="Palomares") {
    if (id%in%c("Garfio", "Patsuezo", "Lynx_Ketamina", "Lynx_Llerena", 
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
        filter(location.lat > 48) 
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
    
    if (id%in%c("B_MVBOB54F","B_MVBOB66M","B_NEBOB33M","B_NEBOB35M")) {
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
        filter(location.long < 141.1)
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
        filter(location.lat < -34)
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
  
  if (updatedstudy=="Tatler Kalamurina") {
    
    if(id=="JT38") { 
      track <- crop_range_res(track) %>%
        filter(location.long < 137.9828, location.lat < -27.85)
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
  
  if (updatedstudy=="Walton") {
    if (id%in%c("Mattias","Viktor","Oskar_1")) {
      track <- crop_range_res(track)
    }
    
    if (id=="Bengt") {
      track <- track %>%
        filter(timestamp<ymd_hms("2015-02-02 07:39:00"))
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
    if (id%in%c("Astrid","Charlotte","Hedley","19F","35M","43F","54M","66M","95F")) {
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

get_duration <- function(individual) { #individual name
  
  duration <- as.duration(interval(as.Date(individual$timestamp[1], #change to 2 for inds where first val is 00:00:00
                                           format = "%Y-%m-%d %H:%M:%S"),
                                   as.Date((tail(individual,1)$timestamp),
                                           #change to tail(individual,2)$timestamp[1] for
                                           #inds where LAST val is 00:00:00
                                           format='%Y-%m-%d %H:%M:%S')))@.Data/604800
  
  #gets duration in secs and converts to weeks (604800 seconds in a WEEK bitch)
  
  #duration <- parse_number(duration)
  year <- as.Date(individual$timestamp[1], format = "%Y-%m-%d %H:%M:%S")
  
  return(list(duration,year))
}

data_list <- alldata %>%
  group_split(id) #split into lists for each individual

#from https://stackoverflow.com/questions/7176870/create-a-vector-of-all-dates-in-a-given-year
get_days <- function(year){
  seq(as.Date(paste(year, "-01-01", sep="")), as.Date(paste(year, "-12-31", sep="")), by="+1 day")
}

names <- c()
days <- list(get_days(year("2000-01-01")))

for (i in seq_along(data_list)) {
  
  ind <- data_list[[i]] %>%
    mod_to_tracks()
  
  dur <- get_duration(ind)
  
  #split up individuals with a duration >1 year to count their years seperately
  if (dur[[1]] > 52) {
    
    ind <- ind %>% 
      mutate(year=year(timestamp)) %>%
      unite(c(id, year), col="id", remove=F)
    
    ind_years <- group_split(ind, by=ind$year)
    
    for (j in seq_along(ind_years)) {
      
      ind_days <- ind_years[[j]] %>%
        #mutate to random year because year doesn't matter
        mutate(timestamp = gsub("^[^-]*", "2000", timestamp)) %>%
        mutate(timestamp = as_date(timestamp)) %>%
        distinct(timestamp, .keep_all=T) %>%
        dplyr::select(timestamp)
      
      if(length(ind_days$timestamp)==1) {next}
      
      names <- c(names, unique(ind_years[[j]]$id))
      days <- c(days, list(ind_days$timestamp))
      
    }
    
  } else {
    
    ind_days <- ind %>%
      #mutate to random year because year doesn't matter
      mutate(timestamp = gsub("^[^-]*", "2000", timestamp)) %>%
      mutate(timestamp = as_date(timestamp)) %>%
      distinct(timestamp, .keep_all=T) %>%
      dplyr::select(timestamp)
    
    if(length(ind_days$timestamp)==1) {next}
    
    names <- c(names, unique(ind$id))
    days <- c(days, list(ind_days$timestamp))
    
  }
  
}

names(days) <- c("year", names)

#align all days 
days_dat <- data.frame("name"=NA, "timeseries"=NA) %>%
  drop_na()

for (i in seq_along(days)) {
  
  ind <- data.frame("name"=as.character(names(days)[[i]]), "timeseries"=days[[i]])
  
  days_dat <- days_dat %>%
    rbind(ind)
  
}

#get all clean ridge individuals
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
         mean_treecover=mean_treecover/100) %>%
  unite(c("id", "sp", "updatedstudy"), col="tojoin", sep="", remove=F) %>%
  dplyr::select("tojoin", "id", "sp", "updatedstudy", "clade")

#join with full dataframe for species data
days_full <- days_dat %>%
  separate(col=name, into = c("name","yeartracked"), sep="_(?=[0-9]{4}$)", remove=T) %>% #regex is to get last undersscore
  left_join(ridge, by=c("name"="tojoin")) %>%
  mutate(Species = replace(Species, Species == "Panthera pardus tulliana", "Panthera pardus"),
         Species = replace(Species, Species == "Panthera tulliana", "Panthera pardus"),
         Species = replace(Species, Species == "Panthera pardus orientalis", "Panthera pardus"),
         Species = replace(Species, Species == "Panthera pardus ciscaucasica", "Panthera pardus"),
         Species = replace(Species, Species == "Lycalopex vetulus", "Pseudalopex vetulus"),
         Species = replace(Species, Species == "Panthera tigris tigris", "Panthera tigris"),
         Species = replace(Species, Species == "Canis lupus dingo", "Canis dingo")) %>%
  relocate(c(ID, Species, `Updated Study`, Study), .after=yeartracked) %>%
  right_join(dplyr::select(ridge, c("id", "sp", "updatedstudy", "clade")),
             by=c("ID"="id", "Species"="sp", "Updated Study"="updatedstudy"))

#check individual counts (sum$n should be 1219)
sum <- days_full %>% group_by(`Updated Study`, Species) %>% summarize(n=length(unique(ID)))
sum(sum$n) #1219

#by clade
days_full %>%
  group_by(clade, timeseries) %>%
  mutate(dens = length(timeseries)) %>%
  ungroup() %>%
  group_by(clade) %>%
  mutate(dens = dens/length(clade)) %>%
  ggplot() +
  stat_smooth(aes(x = timeseries, y = dens, color=clade), method = "gam", formula = y ~ s(x, bs = "cc")) +
  coord_polar() +
  scale_x_date(date_breaks="month", date_labels="%b") +
  scale_color_manual(name="Clade", values=c("canidae"="blue","felidae"="red"), labels=c("Canidae", "Felidae")) +
  scale_y_continuous(limits = c(0,0.003)) +
  labs(x=NULL, y="Density of fixes on day of year") +
  theme_bw() 

#by most divergent species
days_full %>%
  mutate(clade=ifelse(clade=="canidae", "Canidae", "Felidae")) %>%
  filter(Species %in% c("Cuon alpinus", "Cerdocyon thous", "Leopardus wiedii", "Canis lupus",
                        "Puma concolor", "Leopardus pardalis", "Vulpes lagopus")) %>% #select most divergent species
  group_by(Species, timeseries) %>%
  mutate(dens = length(timeseries)) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(dens = dens/length(Species)) %>%
  ggplot() +
  stat_smooth(aes(x = timeseries, y = dens, color=Species), method = "gam", formula = y ~ s(x, bs = "cc")) +
  coord_polar() +
  scale_x_date(date_breaks="month", date_labels="%b") +
  scale_fill_manual(values=as.vector(pals::polychrome(7))) +
  labs(x=NULL, y="Density of fixes on day of year") +
  theme_bw() +
  facet_wrap(~clade,ncol=1) +
  theme(text=element_text(size=14),legend.position="bottom")
