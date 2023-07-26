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

#biome info extracted from lat/long from here:
#https://www.arcgis.com/apps/View/index.html?appid=144b1d74a5964d728b25aeb0542de485
studylocations <- data.frame(study=c("Abrahms", "Clark", "Drouilly", "Mahoney", "Oliveira-Santos", "Prugh", "Ramesh", 
                                     "Vanak", "Sekercioglu","Young"),
                             country=c("Botswana", "United States", "South Africa", "United States", "Brazil", 
                                       "United States", "South Africa", "India", "Turkey", "United States"),
                             biome=c("Savanna", "Temperate Forest", "Desert", "Desert/Temperate Forest",
                                     "Savanna", "Temperate Grasslands/Forest", "Grasslands",
                                     "Desert", "Temperate Forest", "Grasslands")) %>%
  unite(c(country, study), col="label", sep="\n", remove=F)

#########################
#---CLEAN CSV SO ALL YOU HAVE TO DO IS LOAD IT
ridge <- read_csv("allridge.csv") %>%
  filter(sp!="Canis dingo") %>% #we're dropping dingos now
  #update species names and phylo indicators (set this up earlier)
  mutate(sp = ifelse(sp=="Canis lupus x lycaon", "Canis lupus", sp),
         sp = ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp),
         sp = ifelse(sp=="Pseudalopex vetulus", "Lycalopex vetula", sp)) %>% 
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
         hunting_cooperativity = as.factor(hunting_cooperativity))

#download tree from: https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre
tree_all <- read.nexus("/tree/file/here")
#recommended to use the dna only one but it doens't have two species that we need!!
#https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre

tree_all$tip.label <- sapply(tree_all$tip.label, function(x) str_extract(x, "[^_]*_[^_]*"))
#rename canis mesomelas and pseudalopex vetulus in phylogeny
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Canis_mesomelas"), "Lupulella_mesomelas")
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Pseudalopex_vetulus"), "Lycalopex_vetula")

#fit landscape level mods + clade diff tests

#filter data to studies and add within study indicator variable
study <- ridge %>%
  mutate(updatedstudygroups = updatedstudy) %>%
  mutate(updatedstudygroups = replace(updatedstudygroups, 
                                      updatedstudy %ni% c("Abrahms", "Clark", "Drouilly", "Young", 
                                                          "Oliveira-Santos.Dataset1", 
                                                          "Prugh", "Ramesh","Sekercioglu", "Vanak"),
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
  filter(updatedstudy%in%c("Abrahms", "Clark", "Drouilly", "Mahoney", "Oliveira-Santos.Dataset1",
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
mods <- ~ updatedstudygroups + inv_ess

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

#predict vals for landscape (aggregate and get raw data mean/se)
study_ests <- do.call(data.frame,
                      aggregate(cbind(inv_ess, ridge_dens_est) ~ sp + updatedstudy + updatedstudygroups, 
                                data=study,
                                FUN = function(x) c(mn = mean(x), se = (sd(x)/sqrt(length((x))))))) %>%
  left_join(dplyr::select(ridge, c("sp","updatedstudy","clade")),
            by=c("sp"="sp","updatedstudy"="updatedstudy")) %>%
  mutate(updatedstudy=ifelse(updatedstudy=="Oliveira-Santos.Dataset1", "Oliveira-Santos", updatedstudy)) %>%
  distinct(sp, updatedstudygroups, .keep_all=T) %>% #there are join warnings but I checked and they're all just duplicates, we're not losing any data here! 
  mutate(ridge_dens_est_cil = ridge_dens_est.mn - (ridge_dens_est.se*qnorm(0.975)),
         ridge_dens_est_ciu = ridge_dens_est.mn + (ridge_dens_est.se*qnorm(0.975)))

#use model to predict ridge density (mean) for each species
for (i in seq_along(study_ests$sp)) {
  
  species <- study_ests$sp[i]
  st <- study_ests$updatedstudygroups[i]
  
  #get species-level vars (kinda messy)
  abrahms <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Abrahms", 1, 0)
  clark <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Clark", 1, 0)
  drouilly <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Drouilly", 1, 0)
  oliveira <- ifelse((study_ests %>% filter(updatedstudygroups==st,
                                            sp==species))$updatedstudygroups=="Oliveira-Santos.Dataset1", 1, 0)
  prughmv <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Prugh MV", 1, 0)
  prughne <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Prugh NE", 1, 0)
  ramesh <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Ramesh", 1, 0)
  sekercioglu <- ifelse((study_ests %>% filter(updatedstudygroups==st, 
                                               sp==species))$updatedstudygroups=="Sekercioglu", 1, 0)
  vanak <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Vanak", 1, 0)
  young <- ifelse((study_ests %>% filter(updatedstudygroups==st, sp==species))$updatedstudygroups=="Young", 1, 0)
  inv_ess <- (study_ests %>% filter(updatedstudygroups==st, sp==species))$inv_ess.mn
  
  #code to estimate species, from Chris
  SLOPE.EST <- REST[[str_replace(species, " ", "_")]]
  SLOPE.VAR <- RCOV[str_replace(species, " ", "_"), str_replace(species, " ", "_")]
  
  GRAD <- c(1, clark, drouilly, oliveira, prughmv, prughne, ramesh, sekercioglu, vanak, young, inv_ess)
  
  EST <- GRAD %*% FIT$beta + SLOPE.EST
  
  #species var
  VAR <- GRAD %*% FIT$vb %*% GRAD + SLOPE.VAR
  SE <- sqrt(VAR)
  
  #then you can calculate normal CIs and back-transform to the ridge density.
  # ctmm:::norm.ci(EST[,1], SE[,1])
  ctmm:::norm.ci(EST[,1], SE[,1])
  exp(ctmm:::norm.ci(EST[,1], SE[,1]))
  
  study_ests$pred_ridge[i] <- exp(ctmm:::norm.ci(EST[,1], SE[,1]))[2]

}

######
#create figure 4

#manually labeller (study biome/country)
study_labs = c(Abrahms = "Savanna\nBotswana",
               Clark = "Temperate Evergreen Forest\nUnited States",
               Drouilly = "Desert\nSouth Africa",
               `Oliveira-Santos` = "Savanna\nBrazil",
               `Prugh MV` = "Temperate Grasslands/Forest\nUnited States",
               `Prugh NE` = "Temperate Evergreen Forest\nUnited States",
               Ramesh = "Mountain Grassland\nSouth Africa",
               Vanak = "Desert\nIndia",
               Sekercioglu = "Temperate Evergreen Forest\nTurkey",
               Young = "Desert/Temperate Evergreen Forest\nUnited States")

#figure 4 version 1 (with confidence intervals)
study_ests %>% 
  mutate(updatedstudygroups=ifelse(updatedstudygroups=="Oliveira-Santos.Dataset1",
                                   "Oliveira-Santos",updatedstudygroups),
         sp=str_replace(sp, " ", "\n"), #mutate for formatting (sorry it's so ugly)
         sp=factor(sp, levels=c("Lycaon\npictus", "Panthera\nleo", "Canis\nlatrans", "Lynx\nrufus",
                                "Lupulella\nmesomelas", "Caracal\ncaracal", "Chrysocyon\nbrachyurus", 
                                "Leptailurus\nserval", "Canis\nlupus", "Lynx\nlynx", 
                                "Canis\naureus", "Vulpes\nbengalensis","Felis\nchaus",
                                "Puma\nconcolor"))) %>%
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

#figure 4 version 2 (with violin plots)
study %>%
  mutate(updatedstudygroups=ifelse(updatedstudygroups=="Oliveira-Santos.Dataset1", 
                                   "Oliveira-Santos", updatedstudygroups)) %>%
  mutate(sp=str_replace(sp, " ", "\n")) %>%
  mutate(sp=factor(sp, levels=c("Lycaon\npictus", "Panthera\nleo", "Canis\nlatrans", "Lynx\nrufus",
                                "Lupulella\nmesomelas", "Caracal\ncaracal", "Chrysocyon\nbrachyurus", 
                                "Leptailurus\nserval", "Canis\nlupus", "Lynx\nlynx", 
                                "Canis\naureus", "Vulpes\nbengalensis", "Felis\nchaus","Puma\nconcolor"))) %>%
  ggplot(mapping=aes(x=sp,
                     y=ridge_dens_est, fill=clade)) +
  geom_violin(alpha = 0.5) +
  geom_point(data = mutate(study_ests, 
                           updatedstudygroups=ifelse(updatedstudygroups=="Oliveira-Santos.Dataset1", 
                                                     "Oliveira-Santos", updatedstudygroups),
                           sp=str_replace(sp, " ", "\n"),
                           sp=factor(sp, levels=c("Lycaon\npictus", "Panthera\nleo", "Canis\nlatrans", "Lynx\nrufus",
                                                  "Lupulella\nmesomelas", "Caracal\ncaracal", "Chrysocyon\nbrachyurus", 
                                                  "Leptailurus\nserval", "Canis\nlupus", "Lynx\nlynx", 
                                                  "Canis\naureus", "Vulpes\nbengalensis", "Felis\nchaus","Puma\nconcolor"))), 
             mapping=aes(x=sp, 
                         y=pred_ridge, fill=clade), 
             size=4, shape = 24) + #get mean pred_ridge, they're actually very close
  scale_fill_manual(name="Clade", values=c("canidae"="blue","felidae"="red"), labels=c("Canidae", "Felidae")) +
  labs(x="Species", y="Ridge Density", color="Clade", 
       caption = "Ridge density threshold: 0.5 \n Triangles represent predicted ridge values") +
  theme_bw() +
  facet_wrap(~updatedstudygroups, ncol=2, scales="free_x", labeller = as_labeller(study_labs)) + 
  theme(text=element_text(size=15),
        legend.position = c(0.1, 0.96), # c(0,0) bottom left, c(1,1) top-right.
        legend.background = element_rect(fill = "white", colour = NA))
