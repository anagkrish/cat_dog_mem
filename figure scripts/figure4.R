#code for figure 4 (shared landscapes)
#########

#load libraries
library(tidyverse)
library(ctmm)
library(lubridate)
library(parsedate)
library(dplyr)
library(ape)
library(phytools)
library(MASS)
library(metafor)

source("ranef.rma.mv.R") #hacked metafor function from Chris Fleming
"%ni%" <- Negate("%in%")

#########################
#cleaned final dataset (n=1219 without dingos) - cleaned in file get_new_ridges.Rmd
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


#download tree from: https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre
tree_all <- read.nexus("/tree/file/here")
#recommended to use the dna only one but it doens't have two species that we need!!
#https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre

tree_all$tip.label <- sapply(tree_all$tip.label, function(x) str_extract(x, "[^_]*_[^_]*"))
#rename canis mesomelas and pseudalopex vetulus in phylogeny
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Canis_mesomelas"), "Lupulella_mesomelas")
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Pseudalopex_vetulus"), "Lycalopex_vetulus")

#fit landscape level mods + clade diff tests

#filter data to studies and add within study indicator variable
study <- ridge %>%
  mutate(updatedstudygroups = updatedstudy) %>%
  mutate(updatedstudygroups = replace(updatedstudygroups, 
                                      updatedstudy %ni% c("Abrahms", "Clark", "Drouilly", "Young", 
                                                          "Azevedo.Lemos", "Prugh", "Ramesh",
                                                          "Sekercioglu", "Vanak"),
                                      NA),
         updatedstudygroups = replace(updatedstudygroups, study=="Young", NA),
         updatedstudygroups = replace(updatedstudygroups, 
                                      updatedstudygroups=="Clark"&grepl("C2",id), NA),
         updatedstudygroups = replace(updatedstudygroups, 
                                      updatedstudygroups=="Clark"&grepl("OR",id), NA),
         updatedstudygroups = replace(updatedstudygroups, 
                                      updatedstudygroups=="Prugh"&grepl("MV",id), "Prugh MV"),
         updatedstudygroups = replace(updatedstudygroups, 
                                      updatedstudygroups=="Prugh"&grepl("NE",id), "Prugh NE")) %>%
  relocate(updatedstudygroups, .after=updatedstudy) %>%
  filter(updatedstudy%in%c("Abrahms", "Clark", "Drouilly", "Mahoney", "Azevedo.Lemos",
                           "Prugh", "Ramesh", "Vanak", "Sekercioglu","Young")) %>%
  filter(sp%in%c("Canis latrans", "Lynx rufus", "Lycaon pictus", "Panthera leo", "Lupulella mesomelas", "Caracal caracal",
                 "Chrysocyon brachyurus", "Puma concolor", "Leptailurus serval", "Canis lupus", "Lynx lynx",
                 "Vulpes bengalensis", "Canis aureus", "Felis chaus")) %>% #get all sp 
  filter(!grepl("C2", id), #drop puma concolor from clark bc it's not same landscape
         study != "Young") #drop young bobcats bc its not same landscape

#set up phylo correlation matrix for this study
phylo <- drop.tip(tree_all, c(which(str_extract(tree_all$tip.label, "[^_]*_[^_]*") %ni%
                                      c(unique(ridge$phylo),
                                        "Phataginus_tetradactyla",
                                        "Hyaena_hyaena",
                                        "Herpestes_sanguineus",
                                        "Enhydra_lutris",
                                        "Halichoerus_grypus",
                                        "Ursus_americanus")))) #keep all species in studies + 6 outgroups

ph_corr = vcv(corPagel(1, phylo), corr = TRUE) #correlation matrix
study <- study[order(match(study$phylo, phylo$tip.label)),] #reorder data to match correlation matrix

# marginalize covariance matrix down to the species we have data on (not sure if necessary)
ph_corr = vcv(corPagel(1, phylo), corr = TRUE) #CORRELATION MATRIX for brms
SPECIES <- unique(study$phylo)
PHS <- rownames(ph_corr)
IN <- PHS[PHS %in% SPECIES]
ph_corr <- ph_corr[IN,IN] #get full corrmat and then trim (not totally sure if this is relevant)

#random effects for each species (phylo) and each individual (point)
random <- list(~ 1|phylo, ~ 1|point) # these will be the only fitted variances in metafor, the latter is a normal regression error that is not included by default
R <- list(phylo=ph_corr)

#final model fit with all vars BUT spatial and species level!!
mods <- ~ log_mass_st + pursuit + disruptfast + slowwalking + log_hr_st + inv_ess + log_roughness_st + mean_treecover + mean_hfi + seasonality_dhi_gpp_st

#fit model
FIT <- metafor::rma.mv(log_ridge,V=0,mods=mods,random=random,R=R,data=study)
summary(FIT) 

#clade diff estimates
EFF <- ranef2(FIT)$phylo
REST <- EFF$est # random effects
RCOV <- EFF$COV # error covariance

CANINE <- study$clade=="canidae"
CANINE <- unique(study$phylo[CANINE])
FELINE <- study$clade=="felidae"
FELINE <- unique(study$phylo[FELINE])

iCOV <- ctmm:::PDsolve(RCOV)
dW.dl <- iCOV %*% cbind( names(REST) %in% CANINE , names(REST) %in% FELINE )
M <- rbind( colSums( dW.dl[CANINE,,drop=FALSE] ) , colSums( dW.dl[FELINE,,drop=FALSE] ) )
lambda <- c( solve(M) %*% c(1,-1) )
W <- c(dW.dl %*% lambda)
names(W) <- rownames(dW.dl)

# W[CANINE,] # check +1
sum(W[CANINE])
#
# W[FELINE,] # check -1
sum(W[FELINE])

# point estimate
DIFF <- c(W %*% REST)
VAR.DIFF <- c(W %*% RCOV %*% W)

# relative impact on ridge density (back transform)
ctmm:::norm.ci(DIFF, VAR.DIFF)
exp(ctmm:::norm.ci(DIFF, VAR.DIFF))

#predict vals for landscape (aggregate and get empirical mean and confidence intervals)
study_ests <- do.call(data.frame,
                      aggregate(cbind(inv_ess, `ridge_dens_est (1/m)`) ~ sp + updatedstudy + updatedstudygroups, 
                                data=study,
                                FUN = function(x) c(mn = mean(x), se = (sd(x)/sqrt(length((x))))))) %>%
  left_join(dplyr::select(ridge, c("sp","updatedstudy","clade")),
            by=c("sp"="sp","updatedstudy"="updatedstudy")) %>%
  distinct(sp, updatedstudygroups, .keep_all=T) %>%
  mutate(ridge_dens_est_cil = `ridge_dens_est (1/m).mn` - (`ridge_dens_est (1/m).se`*qnorm(0.975)),
         ridge_dens_est_ciu = `ridge_dens_est (1/m).mn` + (`ridge_dens_est (1/m).se`*qnorm(0.975)))

#aggregate variables to predict ridge dens/species/site
sp_level_ests <- aggregate(cbind(log_mass_st, log_roughness_st, mean_treecover,
                                 mean_hfi, log_hr_st, inv_ess, seasonality_dhi_gpp_st, `ridge_dens_est (1/m`) ~ sp + updatedstudygroups, #mean_dhi_ndvi,  log_dhi_gpp,
                           data=study,
                           FUN = mean) %>%
  left_join(dplyr::select(ridge, c("sp", "clade", "pursuit", "disruptfast", "slowwalking", "hunting_movement")),
            by=c("sp"="sp")) %>%
  distinct()

#use model to predict ridge density (mean) for each species
for (i in seq_along(sp_level_ests$sp)) {
  
  #get avg vars
  species <- sp_level_ests$sp[i]
  study <- sp_level_ests$updatedstudygroups[i]

  SLOPE.EST <- REST[[str_replace(species, " ", "_")]]
  SLOPE.VAR <- RCOV[str_replace(species, " ", "_"), str_replace(species, " ", "_")]
  
  #get avg vars
  log_mass <- (sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$log_mass_st
  log_roughness <- (sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$log_roughness_st
  mean_treecover <- (sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$mean_treecover
  mean_hfi <- (sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$mean_hfi
  seasonality_dhi_gpp <- (sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$seasonality_dhi_gpp_st
  log_hr <- (sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$log_hr_st
  inv_ess <- (sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$inv_ess
  #indicator vars for movement (need to set vars manually >:())
  pursuit <- ifelse((sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$pursuit==1, 1, 0)
  disruptfast <- ifelse((sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$disruptfast==1, 1, 0)
  slowwalking <- ifelse((sp_level_ests %>% filter(sp==species, updatedstudygroups==study))$slowwalking==1, 1, 0)
  
  GRAD <- c(1, log_mass, pursuit, disruptfast, slowwalking,
            log_hr, inv_ess, log_roughness, mean_treecover, mean_hfi, seasonality_dhi_gpp)
  
  EST <- GRAD %*% FIT$beta + SLOPE.EST
  
  #species var
  VAR <- GRAD %*% FIT$vb %*% GRAD + SLOPE.VAR
  SE <- sqrt(VAR)
  
  #then you can calculate normal CIs and back-transform to the ridge density.
  # ctmm:::norm.ci(EST[,1], SE[,1])
  ctmm:::norm.ci(EST[,1], SE[,1])
  exp(ctmm:::norm.ci(EST[,1], SE[,1]))
  
  sp_level_ests$pred_ridge[i] <- exp(ctmm:::norm.ci(EST[,1], SE[,1]))[2]
  
}

#join dataframes; less hassle this way
study_ests <- study_ests %>%
  left_join(dplyr::select(sp_level_ests, -c("clade")), #sp, updatedstudygroups, pred_ridge), 
            by=c("sp"="sp", "updatedstudygroups"="updatedstudygroups"))

######
#create figure 4

#manually labeller (study biome/country)
study_labs = c(Abrahms = "Tropical Savanna\nBotswana",
               Clark = "Temperate Conifer Forest\nOregon, United States",
               Drouilly = "Desert & Xeric Shrubland\nSouth Africa",
               `Azevedo.Lemos` = "Tropical Savanna/\nTropical Moist Broadleaf Forest\nBrazil",
               `Prugh MV` = "Temperate Conifer Forest/\nTemperate Savanna\nWashington, United States",
               `Prugh NE` = "Temperate Conifer Forest\nWashington, United States",
               Ramesh = "Tropical Savanna\nSouth Africa",
               Vanak = "Desert & Xeric Shrubland\nIndia",
               Sekercioglu = "Temperate Broadleaf & Mixed Forest\nTürkiye",
               Young = "Desert & Xeric Shrubland/\nTemperate Conifer Forest\nUtah, United States")

#figure 4 version 1 (with confidence intervals)
#actual version was put together by Dr. Anshuman Swain

study_ests %>% 
  mutate(sp=ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp),
         sp=str_replace(sp, " ", "\n"),
         sp=factor(sp, levels=c("Lycaon\npictus", "Panthera\nleo", "Canis\nlatrans", "Lynx\nrufus",
                                "Lupulella\nmesomelas", "Caracal\ncaracal", "Chrysocyon\nbrachyurus", 
                                "Leptailurus\nserval", "Canis\nlupus", "Lynx\nlynx", "Canis\naureus", 
                                "Vulpes\nbengalensis","Felis\nchaus", "Puma\nconcolor"))) %>%
  ggplot(mapping=aes(x=sp, y=ridge_dens_est.mn, color=clade)) +
  geom_point() +
  geom_errorbar(mapping=aes(ymin=ridge_dens_est_cil, ymax=ridge_dens_est_ciu, color=clade), width=0.5) +
  geom_point(mapping=aes(x=sp, y=pred_ridge, color=clade), 
             size=4, shape = 24) + #get mean pred_ridge, they're actually very close
  scale_color_manual(name="Clade", values=c("canidae"="blue","felidae"="red"), labels=c("Canidae", "Felidae")) +
  labs(x="Species", y="Ridge Density", color="Clade", 
       caption = "Ridge density threshold: 0.5 \n Triangles represent predicted ridge values") +
  theme_bw() +
  facet_wrap(~updatedstudygroups, ncol=2, scales="free", labeller = as_labeller(study_labs)) + 
  theme(text=element_text(size=15),
        plot.margin = unit(c(1.3, 0.2, 0.2, 0.2), "cm"),
        legend.direction = "horizontal",
        legend.position = c(0.15, 1.06), # c(0,0) bottom left, c(1,1) top-right.
        legend.background = element_rect(fill = "white", colour = NA))
