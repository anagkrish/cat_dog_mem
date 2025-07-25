#FIGURE 2 MAIN

library(tidyverse) 
library(rsample)
library(lubridate)
library(emmeans)
library(effects)
library(ape)
library(phytools)
library(ggtree)
library(brms)
library(MASS)
library(nlme)
library(metafor)
library(distributions3)

"%ni%" <- Negate("%in%")
source("ranef.rma.mv.R") #hacked metafor function from Chris Fleming

#load full csv
ridge <- read_csv("ridge.csv") %>%
  filter(`kept (y/n)` == "y") %>% #drop all dropped individuals
  mutate(pursuit = as.factor(pursuit),
         disruptfast = as.factor(disruptfast),
         slowwalking = as.factor(slowwalking), #convert to factors bc they're read in as numerical
         log_mass = log(mass),
         log_hr = log(area_ud_est),
         inv_ess = 1/ess,
         log_roughness = log(mean_roughness),
         log_hfi = log(mean_hfi),
         log_dhi_gpp = log(mean_dhi_gpp),
         point = as.factor(1:length(id)),
         log_ridge = log(ridge_dens_est),
         log_ridge_length = log(ridge_length_est),
         logit_ridge_prob = log(ridge_use_prob/(1-ridge_use_prob)),
         logit_ridge_prob_reduced10 = log(ridge_prob_reduced10/(1-ridge_prob_reduced10)),
         logit_ridge_prob_reduced100 = log(ridge_prob_reduced100/(1-ridge_prob_reduced100)),
         logit_ridge_prob_fixed10 = log(ridge_prob_fixed10/(1-ridge_prob_fixed10)),
         logit_ridge_prob_fixed50 = log(ridge_prob_fixed50/(1-ridge_prob_fixed50)),
         logit_ridge_prob_fixed100 = log(ridge_prob_fixed10/(1-ridge_prob_fixed100))) %>%
  mutate(log_mass_st=mosaic::zscore(log_mass),
         log_hr_st=mosaic::zscore(log_hr), 
         log_roughness_st=mosaic::zscore(log_roughness),
         seasonality_dhi_gpp_st=mosaic::zscore(seasonality_dhi_gpp),
         speed_est_st=mosaic::zscore(speed_est, na.rm=T),
         mean_treecover=mean_treecover/100)

#download tree from: https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre
tree_all <- read.nexus("/tree/file/here")
#recommended to use the dna only one but it doens't have two species that we need!!
#https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre

tree_all$tip.label <- sapply(tree_all$tip.label, function(x) str_extract(x, "[^_]*_[^_]*"))
#rename canis mesomelas and pseudalopex vetulus in phylogeny
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Canis_mesomelas"), "Lupulella_mesomelas")
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Pseudalopex_vetulus"), "Lycalopex_vetulus")

#drop all tips except species in study and some other carnivora ingroups
#phataginus: pangolin (closest outgroup relative to cats and dogs)
#hyaena: hyena (ingroup for cats)
#herpestes: mongoose (ingroup for cats)
#ursus: bear (further along carnivora branch)
#halichoerus: grey seal (further along carnivora branch)
#enhydra: sea otter (furthest along carnivora branch)
# phylo_dna <- drop.tip(tree_dna, c(which(str_extract(tree_dna$tip.label, "[^_]*_[^_]*") %ni%
#                                       c(unique(ridge$phylo),
#                                                               "Phataginus_tetradactyla",
#                                                               "Hyaena_hyaena",
#                                                               "Herpestes_sanguineus",
#                                                               "Enhydra_lutris",
#                                                               "Halichoerus_grypus",
#                                                               "Ursus_americanus"))))

phylo <- drop.tip(tree_all, c(which(str_extract(tree_all$tip.label, "[^_]*_[^_]*") %ni%
                                      c(unique(ridge$phylo),
                                        "Phataginus_tetradactyla",
                                        "Hyaena_hyaena",
                                        "Herpestes_sanguineus",
                                        "Enhydra_lutris",
                                        "Halichoerus_grypus",
                                        "Ursus_americanus"))))

plotTree(phylo)
ph_corr = vcv(corPagel(1, phylo), corr = TRUE) #CORRELATION MATRIX for brms
ridge <- ridge[order(match(ridge$phylo, phylo$tip.label)),] #reorder data to match phylogeny

######################################################
#set up covariance matrix for metafor analysis

# marginalize covariance matrix down to the species we have data on (not sure if necessary)
SPECIES <- unique(ridge$phylo)
PHS <- rownames(ph_corr)
IN <- PHS[PHS %in% SPECIES]
ph_corr <- ph_corr[IN,IN] #get full corrmat and then trim (not totally sure if this is relevant)

#random effects for each species (phylo) and each individual (point)
random <- list(~ 1|phylo, ~ 1|point) # these will be the only fitted variances in metafor, the latter is a normal regression error that is not included by default
R <- list(phylo=ph_corr)

#### single mod fit
#best fit all
mods <- mods <- ~ log_mass_st + pursuit + disruptfast + slowwalking + log_hr_st + inv_ess + log_roughness_st + mean_treecover + mean_hfi + seasonality_dhi_gpp_st
FIT <- metafor::rma.mv(log_ridge,V=0,mods=mods,random=random,R=R,data=ridge)
summary(FIT)

# hacked metafor function to give covariance matrix (source from file)
EFF <- ranef2(FIT)$phylo
REST <- EFF$est # random effects
RCOV <- EFF$COV # error covariance

CANINE <- ridge$clade=="canidae"
CANINE <- unique(ridge$phylo[CANINE])
FELINE <- ridge$clade=="felidae"
FELINE <- unique(ridge$phylo[FELINE])

# minimum-error difference effect (more complicated than the usual relations)
# W == c(w.canine,w.feline) # weights to optimize
# sum(w.canine) == +1 && sum(w.feline) == -1
# so that E[W %*% EST] == mean(canine) - mean(feline)
# Lagrangian for optimal weights W
# W %*% RCOV %*% W + (lambda.canine*CANINE + lambda.feline*FELINE) %*% W - lambda.canine + lambda.feline
# solutions
# RCOV[CANINE,] %*% W == lambda.canine * CANINE
# RCOV[FELINE,] %*% W == lambda.feline * FELINE
# W == iCOV %*% c(lambda.c,lambda.f)

iCOV <- ctmm::pd.solve(RCOV)
dW.dl <- iCOV %*% cbind( names(REST) %in% CANINE , names(REST) %in% FELINE )
M <- rbind( colSums( dW.dl[CANINE,] ) , colSums( dW.dl[FELINE,] ) )
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

sp_level_ests <- aggregate(cbind(log_mass_st, log_roughness_st, mean_treecover, 
                                 mean_hfi, log_hr_st, inv_ess, 
                                 seasonality_dhi_gpp_st, `ridge_dens_est (1/m)`) ~ sp, 
                           data=ridge,
                           FUN = mean) %>%
  left_join(dplyr::select(ridge, c("sp", "clade", "pursuit", "disruptfast", "slowwalking", "hunting_movement")), 
            by=c("sp"="sp")) %>%
  distinct(sp, .keep_all=T)

for (i in seq_along(sp_level_ests$sp)) {
  
  species <- sp_level_ests$sp[i]
  
  SLOPE.EST <- REST[[str_replace(species, " ", "_")]]
  SLOPE.VAR <- RCOV[str_replace(species, " ", "_"), str_replace(species, " ", "_")]
  
  #get avg vars
  log_mass_st <- (sp_level_ests %>% filter(sp==species))$log_mass_st
  log_roughness_st <- (sp_level_ests %>% filter(sp==species))$log_roughness_st
  mean_treecover <- (sp_level_ests %>% filter(sp==species))$mean_treecover
  mean_hfi <- (sp_level_ests %>% filter(sp==species))$mean_hfi
  log_hr_st <- (sp_level_ests %>% filter(sp==species))$log_hr_st
  inv_ess <- (sp_level_ests %>% filter(sp==species))$inv_ess
  seasonality_dhi_gpp <- (sp_level_ests %>% filter(sp==species))$seasonality_dhi_gpp
  pursuit <- ifelse((sp_level_ests %>% filter(sp==species))$pursuit==1, 1, 0)
  disruptfast <- ifelse((sp_level_ests %>% filter(sp==species))$disruptfast==1, 1, 0)
  slowwalking <- ifelse((sp_level_ests %>% filter(sp==species))$slowwalking==1, 1, 0)
  
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

phylo$tip.label <- replace(phylo$tip.label, which(phylo$tip.label=="Canis_mesomelas"), "Lupulella_mesomelas")
phylo$tip.label <- replace(phylo$tip.label, which(phylo$tip.label=="Pseudalopex_vetulus"), "Lycalopex_vetulus")
names(phylo$tip.label) <- NULL

#avg ridge vals
d_avg <- sp_level_ests %>%
  #filter(sp %ni% c("Canis dingo", "Canis lupus x lycaon")) %>%
  mutate(sp=str_replace(sp, " ", "_")) %>%
  dplyr::select(c("sp", "ridge_dens_est")) %>%
  rbind(cbind(sp="Phataginus_tetradactyla",`ridge_dens_est (1/m)`=NA),
        cbind(sp="Halichoerus_grypus",`ridge_dens_est (1/m)`=NA),
        cbind(sp="Enhydra_lutris",`ridge_dens_est (1/m)`=NA),
        cbind(sp="Herpestes_sanguineus",`ridge_dens_est (1/m)`=NA),
        cbind(sp="Ursus_americanus",`ridge_dens_est (1/m)`=NA),
        cbind(sp="Hyaena_hyaena",`ridge_dens_est (1/m)`=NA)) #this paired with color scheme manually forces them to be grey

svl_avg <- parse_number(d_avg$`ridge_dens_est (1/m)`)
names(svl_avg) <- d_avg$sp

#predicted ridge vals
d_pred <- sp_level_ests %>%
  filter(sp %ni% c("Canis dingo", "Canis lupus x lycaon")) %>%
  mutate(sp=str_replace(sp, " ", "_")) %>%
  dplyr::select(c("sp", "pred_ridge")) %>%
  rbind(cbind(sp="Phataginus_tetradactyla",pred_ridge=NA),
        cbind(sp="Halichoerus_grypus",pred_ridge=NA),
        cbind(sp="Enhydra_lutris",pred_ridge=NA),
        cbind(sp="Herpestes_sanguineus",pred_ridge=NA),
        cbind(sp="Ursus_americanus",pred_ridge=NA),
        cbind(sp="Hyaena_hyaena",pred_ridge=NA)) #this paired with color scheme manually forces them to be grey

svl_pred <- parse_number(d_pred$pred_ridge)
names(svl_pred) <- d_pred$sp

phylo_dat <- as.matrix(cbind("Average"=svl_avg, "Predicted"=svl_pred), row.names=names(svl_pred))

#heatmap with avg vs predicted!
#version with ggtree (more flexible visualization)

names <- data.frame(label = phylo$tip.label, label2 = str_replace(phylo$tip.label, "_", " "))
phylo_plot <- full_join(phylo, names, by = "label")

#get individual plots of dog/cat ridges for fig. 2
load("vanak/fox13 (vasu) data")

#calculate akde and ks ridges with workflow from figure 1 or calculatefits/getridgefromfits
AKDES <- load("precalc/akde")
allridge <- load("precalc/ridges")

projection(allridge@proj4string)

ridges_km <- spTransform(as_Spatial(ridges_mean_ks),"+proj=aeqd +lat_0=18.3247561459615 +lon_0=74.6171183357356 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")
AKDES@info$identity

par(mar = c(5, 6, 1, 1))
plot(UD, col.UD="grey50", col.grid=NA, cex.lab=3, cex.axis=3, yaxt="n",xaxt = "n")
plot(ridges_km, col="blue", lwd=2, add=T)
text(-15.7, 1.5, "Vulpes bengalensis", cex=2.5)
axis(1, at = c(-16, -15, -14),
     labels = c(-1, 0, 1), cex.axis=3)
axis(2, at = c(-1, -0, 1),
     labels = c(-1, 0, 1), cex.axis=3)

#felid ridge 
load("sekercioglu/ani/lynx/data")

#calculate akde and ridges with workflow from figure 1 or calculatefits/getridgefromfits
AKDES <- load("precalc/akde")
allridge <- load("precalc/ridges")

projection(allridge@proj4string)

ridges_km <- spTransform(allridge,"+proj=aeqd +lat_0=40.2058246972533 +lon_0=42.6315888788603 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")
AKDES@info$identity

par(mar = c(6, 6, 1, 1))
plot(UD, col.UD="grey50", col.grid=NA, cex.lab=3, cex.axis=3)
plot(ridges_km, col="red", lwd=2, add=T)
text(-8, 9.3, "Lynx lynx", cex=2.5)

#put it all together
p <- ggtree(phylo_plot) +
  scale_y_reverse() +
  geom_tiplab(aes(label=label2)) + 
  geom_cladelab(node=63, label="Canidae", align=TRUE, offset = -38.5, angle=270, 
                offset.text = 0.5, textcolor="blue") +
  geom_cladelab(node=45, label="Felidae", align=TRUE, offset = -38.5, angle=270, 
                offset.text = 0.5, textcolor="red") +
  theme_tree2() 

p <- revts(p) +
  scale_x_continuous(breaks=c(-60,-50,-40,-30,-20,-10,0), labels = abs) +
  labs(x="MYA") +
  theme(text=element_text(size=18))

gheatmap(p, phylo_dat,
         offset = 20,
         width=0.3,
         legend_title="Ridge Density") +
  colorspace::scale_fill_continuous_sequential(palette="Greens3", na.value="white") +
  theme(plot.margin = unit(c(0.1, 8, 0, 0.5), "cm"),
        legend.position = c(1.03,0.45), legend.direction="vertical",
        legend.key = element_rect(color="black", fill = NA),
        legend.background=element_blank(),
        legend.key.height = unit(2, 'cm'), text=element_text(size=16),
        axis.title.x = element_text(hjust=0.03, vjust=15)) +
  coord_cartesian(clip="off") +
  guides(fill = guide_colorbar(override.aes=list(fill=NA),
                               title.position = "right", title.theme = element_text(angle=270),
                               title.hjust = 0.5, title.vjust = -2.5)) +
  ggimage::geom_image(x=65, y=-6.5, 
             image="dog.png", 
             size=0.4) +
  ggimage::geom_image(x=65, y=-35.5, 
                      image="cat.png", 
                      size=0.4) +
  labs(fill="Ridge Density")


