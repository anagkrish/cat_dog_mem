library(tidyverse)
library(ctmm)
library(lubridate)
library(parsedate)
library(dplyr)
library(ape)
library(phytools)
library(ggtree)
library(MASS)
library(metafor)
library(geodist)
library(raster)
library(geosphere)
library(ggmap)
library(rosm)
library(poisspatial)

source("~ranef.rma.mv.R")
"%ni%" <- Negate("%in%")

#########################
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

###E1. Effective Sizes
#to get ess, run best fit model on each subset of data and then run below code:
# subs = "shared landscape" #change to name of data subset
# # effsizes <- data.frame(subset=NA, var=NA, low=NA, est=NA, high=NA, type=NA, color=NA) %>%
# #   drop_na()
# 
# effsizes <- effsizes %>%
#   rbind(c(subs, "Phylogenetic", append((exp(qnorm(1-0.05/2)*c(-1,1)*sqrt(FIT$sigma2[[1]])))-1, NA, after=1), 
#           "Biological Variance", "Biological Variance")) %>%
#   rbind(c(subs, "Individual", append((exp(qnorm(1-0.05/2)*c(-1,1)*sqrt(FIT$sigma2[[2]])))-1, NA, after=1), 
#           "Biological Variance", "Biological Variance")) %>%
#   rbind(c(subs, "Phylogenetic\n(Clade)", exp(ctmm:::norm.ci(DIFF, VAR.DIFF))-1, 
#           "Model Coefficient", "Biological Variance"))
# 
# tidy <- broom::tidy(FIT) %>%
#   mutate(ub = estimate + (std.error*qnorm(0.975)),
#          lb = estimate - (std.error*qnorm(0.975)))
# 
# for (i in seq_along(tidy$term)) {
# if (tidy$term[[i]] %in% c("speed_est_st")) {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]], 
#                       get_coefs(filter(tidy, term==tidy$term[[i]])), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# if (tidy$term[[i]] %in% c("log_mass_st")) {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]], 
#                       get_coefs(filter(tidy, term==tidy$term[[i]])), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# if (tidy$term[[i]] == "pursuit1") {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]],
#                       get_coefs(filter(tidy, term==tidy$term[[i]])), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# if (tidy$term[[i]] == "disruptfast1") {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]], 
#                       get_coefs(filter(tidy, term==tidy$term[[i]])), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# if (tidy$term[[i]] == "slowwalking1") {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]],
#                       get_coefs(filter(tidy, term==tidy$term[[i]])), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# 
# if (tidy$term[[i]] %in% c("log_hr_st")) {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]], 
#                       (get_coefs(filter(tidy, term==tidy$term[[i]]))), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# 
# if (tidy$term[[i]] %in% c("log_roughness_st")) {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]], 
#                       (get_coefs(filter(tidy, term==tidy$term[[i]]))), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# if (tidy$term[[i]] == "mean_hfi") {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]], 
#                       (get_coefs(filter(tidy, term==tidy$term[[i]]))), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# if (tidy$term[[i]] == "mean_treecover") {
#   effsizes <- rbind(effsizes, 
#                     c(subs, tidy$term[[i]], 
#                       get_coefs(filter(tidy, term==tidy$term[[i]])), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# if (tidy$term[[i]] %in% c("seasonality_dhi_gpp_st")) {
#   effsizes <- rbind(effsizes, 
#                     #c(subs, "Dynamic Habitat Index\n(Gross Primary Productivity)\n(kg C/m²))", 
#                     c(subs, tidy$term[[i]], 
#                       get_coefs(filter(tidy, term==tidy$term[[i]])), 
#                       "Model Coefficient", "Model Coefficient"))
# }
# 
# }
# 
# colnames(effsizes) = c("subset", "var", "low", "est", "high", "type", "color")

#alternatively, save all effect sizes to a csv and load it here
effsizes <- read_csv("effect/size/here") %>%
  mutate(type=ifelse(type=="Model Coefficient", "Model Coefficients and Effects", type),
         color=ifelse(color=="Model Coefficient", "Model Coefficients and Effects", color)) %>%
  mutate(var=ifelse(var=="Human Footprint Index (%)", "Human Footprint\nIndex (%)", var))

#modify for display and naming purposes
effsizes <- effsizes %>% 
  mutate(type=ifelse(type=="Model Coefficient", "Model Coefficients and Effects", type),
         color=ifelse(color=="Model Coefficient", "Model Coefficients and Effects", color),
         name=var) %>%
  mutate(name=ifelse(var=="log_mass_st", "Mass (g)\\*",
                     ifelse(var=="pursuit1", "Hunting Movement<br>(Pursuit)",
                     ifelse(var=="disruptfast1", "Hunting Movement<br>(Disruptive Fast)",
                     ifelse(var=="slowwalking1", "Hunting Movement<br>(Slow Walking)",
                     ifelse(var=="log_hr_st", "Home Range<br>Area (m²)\\*",
                     ifelse(var=="log_roughness_st", "Terrain<br>Roughness (m)\\*",
                     ifelse(var=="mean_hfi", "Human Footprint Index (%)",
                     ifelse(var=="mean_treecover", "Tree Cover (%)",
                     ifelse(var=="seasonality_dhi_gpp_st", 
                            "Dynamic Habitat Index<br>(Gross Primary Productivity)<br>(kg C/m²))\\*",
                     ifelse(var=="speed_est", "Speed<br>(m/s)\\*", name)))))))))))


`SUBSET NAME` <- effsizes %>% #plot object name corresponds to subset: full, dna, speed, clean, year
  filter(subset%in%c("subset name here")) %>%
  mutate(subset = factor(subset, levels = c("full", "dna only", "speed inc", 
                                            "clean", "duration > 1yr", "same landscape"),
                         labels = c("Full", "DNA Only Tree", "Speed Included", 
                                    "No Preprocessing", "Year or Longer", "Shared Landscapes"))) %>%
  mutate(textcol = ifelse(var %in% c("Phylogenetic", "Individual", "Phylogenetic<br>(Clade)"), "**", 
                          #phylo/rand effects
                          ifelse(var %in% c("log_mass_st", "log_hr_st", "log_roughness_st"), "***",
                          #log transformed vars
                          ifelse(var %in% c("pursuit1", "disruptfast1", "slowwalking1"), "_", 
                          #indicator vars
                          ifelse(var %in% c("mean_hfi", "mean_treecover", "seasonality_dhi_gpp_st", "speed_est"), "", NA))))) %>%
                          #no transformation (not including standardization)
  mutate(label = paste(textcol, name, textcol, sep = ""),
         label = fct_reorder(label, rev(sort(as.character(label))))) %>%
  ggplot(mapping=aes(x=est, y = label, color = color)) +
  geom_point() +
  geom_errorbar(mapping=aes(xmin=low, xmax=high)) +
  geom_vline(mapping=aes(xintercept=0), linetype="dashed") +
  labs(x="Effect Size (Percent Increase in Ridge Density)", y = NULL, color = NULL) +
  scale_x_continuous(labels = scales::percent) +
  scale_color_manual(values=c("Model Coefficients and Effects"="#5ac85a","Biological Variance"="#6d2f9c")) +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.position = "none",
        text=element_text(size=15),
        axis.text.y=element_markdown()) +
  facet_grid(type~., scales="free", space="free") +
  ggforce::facet_col(type~., scales = 'free', space = 'free')

p <- cowplot::plot_grid(dna + labs(subtitle="\nDNA Only Tree (N=1183, C=15, F=17)") + theme(axis.title.x=element_blank()), 
                        speed + labs(subtitle="\nSpeed Included (N=1011, C=16, F=18)") + theme(axis.title.x=element_blank()), 
                        clean + labs(subtitle="\nNo Preprocessing (N=1064, C=16, F=18)") + theme(axis.title.x=element_blank()), 
                        year + labs(subtitle="\nYear or Longer (N=338, C=13, F=16)") + theme(axis.title.x=element_blank()),
                        shared + labs(subtitle="\nShared Landscapes (N=216, C=7, F=7)") + theme(axis.title.x=element_blank()),
                        ncol=2,
                        rel_heights=(c(1,1,0.4)),
                        labels=c("A)", "B)", "C)", "D)", "E)"),
                        label_size=20,
                        hjust = -2,
                        vjust = 2)  +
  theme(panel.background = element_rect(fill="white"), text=element_text(size=35))

cowplot::ggdraw(cowplot::add_sub(p, label="Effect Size (Percent Increase in Ridge Density)", hjust = 0.5, color="black"))

###E2. Same Species Different Landscapes
speciesdifflandscape <- ridge %>% filter(sp%in%c("Canis latrans", "Canis lupus", "Vulpes vulpes", "Lynx rufus",
                                                 "Puma concolor", "Caracal caracal", "Lupulella mesomelas", "Felis silvestris",
                                                 "Lycaon pictus", "Lynx lynx", "Panthera leo"))

##fit mod
phylo <- drop.tip(tree_all, c(which(str_extract(tree_all$tip.label, "[^_]*_[^_]*") %ni%
                                      c(unique(ridge$phylo),
                                        "Phataginus_tetradactyla",
                                        "Hyaena_hyaena",
                                        "Herpestes_sanguineus",
                                        "Enhydra_lutris",
                                        "Halichoerus_grypus",
                                        "Ursus_americanus"))))


ph_corr = vcv(corPagel(1, phylo), corr = TRUE) #CORRELATION MATRIX for brms
speciesdifflandscape <- speciesdifflandscape[order(match(speciesdifflandscape$phylo, 
                                                         phylo$tip.label)),] #reorder data to match phylogeny

# marginalize covariance matrix down to the species we have data on (not sure if necessary)
######################################################
# marginalize covariance matrix down to the species we have data on (not sure if necessary)
# except on steroids bc now it's per landscape

SPECIES <- unique(speciesdifflandscape$phylo)
PHS <- rownames(ph_corr)
IN <- PHS[PHS %in% SPECIES]
ph_corr <- ph_corr[IN,IN]

#random effects for each species (phylo) and each individual (point)
random <- list(~ 1|phylo, ~ 1|point) # these will be the only fitted variances in metafor, the latter is a normal regression error that is not included by default
R <- list(phylo=ph_corr)

mods <- ~ log_hr + inv_ess

FIT <- metafor::rma.mv(log_ridge,V=0,mods=mods,random=random,R=R,data=speciesdifflandscape)
summary(FIT)

EFF <- ranef2(FIT)$phylo
REST <- EFF$est # random effects
RCOV <- EFF$COV # error covariance

CANINE <- speciesdifflandscape$clade=="canidae"
CANINE <- unique(speciesdifflandscape$phylo[CANINE])
FELINE <- speciesdifflandscape$clade=="felidae"
FELINE <- unique(speciesdifflandscape$phylo[FELINE])

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

#get phylo model predictions  for landscape
study_ests <- do.call(data.frame,
                      aggregate(cbind(log_mass_st,
                                      log_hr_st, inv_ess, log_roughness_st,
                                      mean_treecover,  mean_hfi, 
                                      seasonality_dhi_gpp_st,
                                      ridge_dens_est) ~ sp + updatedstudy, 
                                data=speciesdifflandscape,
                                FUN = mean)) %>%
  left_join(dplyr::select(speciesdifflandscape, 
                          c("sp", "clade", "pursuit", 
                            "disruptfast", "slowwalking", "hunting_movement"))) %>%
  distinct(sp, updatedstudy, .keep_all=T)


for (i in seq_along(study_ests$sp)) {
  
  species <- study_ests$sp[i]
  st <- study_ests$updatedstudy[i]
  
  SLOPE.EST <- REST[[str_replace(species, " ", "_")]]
  SLOPE.VAR <- RCOV[str_replace(species, " ", "_"), str_replace(species, " ", "_")]
  
  #get avg vars
  log_mass <- (study_ests %>% filter(sp==species, updatedstudy==st))$log_mass_st
  log_hr <- (study_ests %>% filter(sp==species, updatedstudy==st))$log_hr_st
  inv_ess <- (study_ests %>% filter(sp==species, updatedstudy==st))$inv_ess
  log_roughness <- (study_ests %>% filter(sp==species, updatedstudy==st))$log_roughness_st
  mean_treecover <- (study_ests %>% filter(sp==species, updatedstudy==st))$mean_treecover
  mean_hfi <- (study_ests %>% filter(sp==species, updatedstudy==st))$mean_hfi
  seasonality_dhi_gpp <- (study_ests %>% filter(sp==species, updatedstudy==st))$seasonality_dhi_gpp_st
  #indicator vars for movement (need to set vars manually >:())
  pursuit <- ifelse((study_ests %>% filter(sp==species, updatedstudy==st))$pursuit==1, 1, 0)
  disruptfast <- ifelse((study_ests %>% filter(sp==species, updatedstudy==st))$disruptfast==1, 1, 0)
  #slowwalking dropped due to redundancy (probably due to small subset of species)
  # slowwalking <- ifelse((study_ests %>% filter(sp==species, updatedstudy==st))$slowwalking==1, 1, 0)
  
  GRAD <- c(1, log_mass, pursuit, disruptfast, 
            log_hr, inv_ess, log_roughness, mean_treecover, mean_hfi,seasonality_dhi_gpp)
  
  EST <- GRAD %*% FIT$beta + SLOPE.EST
  
  #species var
  VAR <- GRAD %*% FIT$vb %*% GRAD + SLOPE.VAR
  SE <- sqrt(VAR)
  
  ctmm:::norm.ci(EST[,1], SE[,1])
  exp(ctmm:::norm.ci(EST[,1], SE[,1]))
  
  study_ests$pred_ridge[i] <- exp(ctmm:::norm.ci(EST[,1], SE[,1]))[2]

}


biomes <- readxl::read_excel("biomes.csv") %>% #biomes extracted from RESOLVE ecoregions map in arcGIS online
  mutate(identifier = replace_na(as.character(identifier), ""),
         sp=ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp))

study_ests_plot <- study_ests %>%
  left_join(biomes, by = c("sp"="sp", "updatedstudy"="updatedstudy")) %>%
  mutate(sp=ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp)) %>%
  mutate(updatedstudy=ifelse(updatedstudy=="Oliveira-Santos.Dataset1", "Oliveira-Santos", updatedstudy)) %>%
  unite("breaks", sp, updatedstudy, sep=" ", remove=F) %>%
  unite("label", `alpha code`, biome, sep="\n", remove=T) %>%
  unite("label", label, identifier, sep=" ", remove=T)

levels <- c(rbind(FELINE, CANINE))
levels <- levels %>%
  str_replace_all("_"," ") %>%
  replace(which(levels=="Canis mesomelas"), "Lupulella mesomelas") %>%
  unique()
levels <- append(levels, " ", after=10)

#get landscape-level predictions for each mod
species_summary_all <- data.frame("mod"=NA, "species"=NA, "term"=NA, 
                                  "estimate"=NA, "std.error"=NA, "statistic"=NA, "p.value"=NA, "aic"=NA, "bic"=NA) %>%
  drop_na()

species_ests_all <- data.frame("mod"=NA, "species"=NA, "term"=NA, 
                               "estimate"=NA, "std.error"=NA, "statistic"=NA, "p.value"=NA) %>%
  drop_na()

species_dat <- speciesdifflandscape %>%
  left_join(biomes, by = c("sp"="sp", "updatedstudy"="updatedstudy")) %>%
  mutate(sp=ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp))  %>%
  unite("landscape", country, biome, sep=" ", remove=T) %>%
  split(speciesdifflandscape$sp)

for (species in species_dat) {
  
  mod_sparse <- lm(log_ridge ~ updatedstudy, data=species)
  
  species_summary <- broom::tidy(mod_sparse) %>%
    cbind("mod"="indicator only", "species"=unique(species$sp), "aic"=AIC(mod_sparse), "bic"=BIC(mod_sparse)) %>%
    relocate(c(mod, species), .before = term)
  
  species_summary_all <- rbind(species_summary_all, species_summary)
  
  species_ests <- aggregate(cbind(log_hr, inv_ess, `ridge_dens_est (1/m)`) ~ sp + updatedstudy + landscape, 
                            data=species,
                            FUN = mean)
  
  species_ests <- species_ests %>%
    cbind(predict.lm(mod_sparse, newdata=species_ests, interval="prediction")) %>%
    rename("pred_ridge_est" = "fit",
           "pred_ridge_lower" = "lwr",
           "pred_ridge_upper" = "upr") %>%
    mutate(pred_ridge_lower = exp(pred_ridge_lower),
           pred_ridge_est = exp(pred_ridge_est),
           pred_ridge_upper = exp(pred_ridge_upper)) %>%
    cbind("mod"="indicator only") %>%
    relocate(mod, .before = sp)
  
  species_ests_all <- rbind(species_ests_all, species_ests)
  
  mod_full <- lm(log_ridge ~ updatedstudy + log_hr + inv_ess, data=species)
  
  species_summary <- broom::tidy(mod_full) %>%
    cbind("mod"="indicator + model", "species"=unique(species$sp), "aic"=AIC(mod_full), "bic"=BIC(mod_full)) %>%
    relocate(c(mod, species), .before = term)
  
  species_summary_all <- rbind(species_summary_all, species_summary)
  
  species_ests <- aggregate(cbind(log_hr, inv_ess, ridge_dens_est) ~ sp + updatedstudy + landscape, 
                            data=species,
                            FUN = mean)
  
  species_ests <- species_ests %>%
    cbind(predict.lm(mod_full, newdata=species_ests, interval="prediction")) %>%
    rename("pred_ridge_est" = "fit",
           "pred_ridge_lower" = "lwr",
           "pred_ridge_upper" = "upr") %>%
    mutate(pred_ridge_lower = exp(pred_ridge_lower),
           pred_ridge_est = exp(pred_ridge_est),
           pred_ridge_upper = exp(pred_ridge_upper)) %>%
    cbind("mod"="indicator + model") %>%
    relocate(mod, .before = sp)
  
  species_ests_all <- rbind(species_ests_all, species_ests)
  
}

# species_summary_all_table <- species_summary_all %>%
#   mutate(term=str_replace(term, "updatedstudy","")) %>%
#   left_join(dplyr::distinct(do.call(rbind, species_dat), sp, updatedstudy, landscape), 
#             by = c("species"="sp",
#                    "term"="updatedstudy")) %>%
#   relocate(landscape, .after=term)

species_ests_all_toplot <- species_ests_all %>%
  left_join(dplyr::distinct(speciesdifflandscape, sp, clade), by = c("sp"="sp")) %>%
  mutate(clade=ifelse(sp=="Lupulella mesomelas", "canidae", clade)) %>%
  unite("breaks", sp, updatedstudy, sep=" ", remove=F)

#make figure E2
speciesdifflandscape %>%
  mutate(sp=ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp),
         clade=paste(clade, "full", sep=" ")) %>%
  unite("breaks", sp, updatedstudy, sep=" ", remove=F) %>%
  mutate(sp=factor(sp, levels=levels)) %>% #make sure names are standardized to this subset of data
  ggplot(mapping=aes(x=breaks,
                     y=ridge_dens_est, fill=clade)) +
  #distribution of points per species per landscape
  geom_violin(alpha = 0.5, show.legend=FALSE) +
  #est from simple (non phylo) model
  geom_point(data= mutate(filter(species_ests_all_toplot, mod == "indicator + model"),
                          sp=factor(sp, levels=levels),
                          clade=paste(clade, "landscape", sep=" ")), 
             mapping=aes(x=breaks, 
                         y=pred_ridge_est, color=clade), 
             size=4, shape=5, stroke=2) +
  #error from simple (non phylo) model
  geom_errorbar(data= mutate(filter(species_ests_all_toplot, mod == "indicator + model"),
                             sp=factor(sp, levels=levels),
                             clade=paste(clade, "landscape", sep=" ")), 
                mapping=aes(ymin=pred_ridge_lower, ymax=pred_ridge_upper, color=clade), width=0.35) +
  #est from phylo model
  geom_point(mapping=aes(x=breaks, 
                         y=pred_ridge, fill=clade), 
             size=4, shape = 24) + #get mean pred_ridge, they're actually very close
  scale_color_manual(name=NULL, values=c("canidae full"=NULL,
                                         "canidae landscape"="blue",
                                         "felidae full"=NULL,
                                         "felidae landscape"="red"),
                     labels=c(NULL, "Canidae (Species-specific model\nestimates/confidence intervals)", 
                              NULL, "Felidae (Species-specific model\nestimates/confidence intervals)")) +
  geom_point(data= mutate(study_ests_plot,
                          sp=ifelse(sp=="Canis mesomelas", "Lupulella mesomelas", sp),
                          sp=factor(sp, levels=levels),
                          clade=paste(clade, "full", sep=" ")), 
             mapping=aes(x=breaks, 
                         y=pred_ridge, color=NULL, fill=clade), 
             size=4, shape = 24) + #get mean pred_ridge, they're actually very close
  scale_fill_manual(name="Clade", values=c("canidae full"="blue",
                                           "canidae landscape"=NULL,
                                           "felidae full"="red",
                                           "felidae landscape"=NULL),
                    labels=c("Canidae (Phylogenetic model estimates)", NULL,
                             "Felidae (Phylogenetic model estimates)", NULL)) +
  scale_shape_manual(name=NULL, values=c("canidae full"=24,
                                         "canidae landscape"=4,
                                         "felidae full"=24,
                                         "felidae landscape"=4),
                     labels=c(NULL, NULL, NULL, NULL)) +
  scale_x_discrete(breaks = c(study_ests_plot$breaks), labels = c(study_ests_plot$label)) +
  geom_text(data=mutate(study_ests_plot,
                        sp=factor(sp, levels=levels)), 
            aes(label = c(rep("",41), "***", rep("", 5), "*", ""), 
                y = c(rep(1,41), 40, rep(1,5), 190, 1)), size=10) +
  labs(x=NULL, y="Ridge Density",fill=NULL,color=NULL) +
  theme_bw() +
  facet_wrap(~sp, ncol=2, nrow=6, scales="free", drop=F) +
  theme(text=element_text(size=19),
        legend.position = c(0.22, 0.065), #c(0,0) bottom left, c(1,1) top-right.
        strip.background = element_blank(),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.direction="vertical") +
  guides(fill = guide_legend(order=1, override.aes = list(shape = 24, stroke=1)),
         color = guide_legend(order=2, override.aes = list(shape = 5)))


###E3. Species level Predicted Means (full dataset, best model)

#aggregate predictors (will lead to some inaccuracy bc some species are across mult landscape but its the best we can do)
sp_level_ests <- aggregate(cbind(log_mass_st, log_roughness_st, mean_treecover, 
                                 mean_hfi, log_hr_st, inv_ess, 
                                 seasonality_dhi_gpp_st, ridge_dens_est) ~ sp, #mean_dhi_ndvi,  log_dhi_gpp,
                           data=ridge,
                           FUN = mean) %>%
  left_join(dplyr::select(ridge, c("sp", "clade", "pursuit", "disruptfast", "slowwalking", "hunting_movement")), 
            by=c("sp"="sp")) %>%
  distinct(sp, .keep_all=T)

#make sure best fit model for whole dataset has been loaded and run
#as well as corresponding clade difference estimator
for (i in seq_along(sp_level_ests$sp)) {
  
  species <- sp_level_ests$sp[i]
  
  SLOPE.EST <- REST[[str_replace(species, " ", "_")]]
  SLOPE.VAR <- RCOV[str_replace(species, " ", "_"), str_replace(species, " ", "_")]
  
  #get avg vars
  log_mass_st <- (sp_level_ests %>% filter(sp==species))$log_mass_st
  log_roughness_st <- (sp_level_ests %>% filter(sp==species))$log_roughness_st
  mean_treecover <- (sp_level_ests %>% filter(sp==species))$mean_treecover
  mean_hfi <- (sp_level_ests %>% filter(sp==species))$mean_hfi
  seasonality_dhi_gpp_st <- (sp_level_ests %>% filter(sp==species))$seasonality_dhi_gpp_st
  log_hr_st <- (sp_level_ests %>% filter(sp==species))$log_hr_st
  inv_ess <- (sp_level_ests %>% filter(sp==species))$inv_ess
  #indicator vars for movement (need to set vars manually >:())
  pursuit <- ifelse((sp_level_ests %>% filter(sp==species))$pursuit==1, 1, 0)
  disruptfast <- ifelse((sp_level_ests %>% filter(sp==species))$disruptfast==1, 1, 0)
  slowwalking <- ifelse((sp_level_ests %>% filter(sp==species))$slowwalking==1, 1, 0)
  
  GRAD <- c(1, log_mass_st, pursuit, disruptfast, slowwalking,
            log_hr_st, inv_ess, log_roughness_st, mean_treecover, mean_hfi, seasonality_dhi_gpp_st)
  
  
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

#make violin plot with actual vs predicted values
ggplot() +
  geom_violin(data = ridge, mapping=aes(x=sp,y=`ridge_dens_est (1/m)`, fill=clade), alpha=0.5, width = 1.5) +
  geom_point(data = sp_level_ests, mapping=aes(x=sp, y=pred_ridge, color=clade, shape=hunting_movement), size=5) +
  scale_fill_manual(name="Clade", values=c("canidae"="blue","felidae"="red"), labels=c("Canidae", "Felidae")) +
  scale_color_manual(name="Clade", values=c("canidae"="blue","felidae"="red"), labels=c("Canidae", "Felidae")) +
  scale_shape_manual(name="Hunting Strategy",
                     values=c("Disruptive Fast hunting"=15,
                              "Mixed Strategies"=16,
                              "Pursuit"=17,
                              "Slow walking"=18),
                     labels=c("Disruptive Fast hunting"="Disruptive Fast Hunting",
                              "Mixed Strategies"="Mixed Strategies",
                              "Pursuit"="Pursuit",
                              "Slow walking"="Slow Walking")) +  
  labs(x="Species", y="Ridge Density", caption = "Points represent predicted species-level ridge values") +
  theme_bw() +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=90)) +
  guides(fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
         color = guide_legend(order=1),
         shape = guide_legend(order = 2))

