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

#download tree from: https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre
tree_all <- read.nexus("/tree/file/here")
#recommended to use the dna only one but it doens't have two species that we need!!
#https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre

tree_all$tip.label <- sapply(tree_all$tip.label, function(x) str_extract(x, "[^_]*_[^_]*"))
#rename canis mesomelas and pseudalopex vetulus in phylogeny
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Canis_mesomelas"), "Lupulella_mesomelas")
tree_all$tip.label <- replace(tree_all$tip.label, which(tree_all$tip.label=="Pseudalopex_vetulus"), "Lycalopex_vetula")

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
mods <- ~ log_mass + pursuit + disruptfast + slowwalking + log_hr + inv_ess + log_roughness + mean_treecover + mean_road_cover

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

iCOV <- ctmm:::PDsolve(RCOV)
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

sp_level_ests <- aggregate(cbind(log_mass, log_roughness, mean_treecover, 
                                 mean_road_cover, log_hr, inv_ess, log_dhi_gpp, ridge_dens_est) ~ sp, #mean_dhi_ndvi,  log_dhi_gpp,
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
  log_mass <- (sp_level_ests %>% filter(sp==species))$log_mass
  log_roughness <- (sp_level_ests %>% filter(sp==species))$log_roughness
  mean_treecover <- (sp_level_ests %>% filter(sp==species))$mean_treecover
  mean_road_cover <- (sp_level_ests %>% filter(sp==species))$mean_road_cover
  log_hr <- (sp_level_ests %>% filter(sp==species))$log_hr
  inv_ess <- (sp_level_ests %>% filter(sp==species))$inv_ess
  pursuit <- ifelse((sp_level_ests %>% filter(sp==species))$pursuit==1, 1, 0)
  disruptfast <- ifelse((sp_level_ests %>% filter(sp==species))$disruptfast==1, 1, 0)
  slowwalking <- ifelse((sp_level_ests %>% filter(sp==species))$slowwalking==1, 1, 0)
  
  GRAD <- c(1, log_mass, pursuit, disruptfast, slowwalking,
            log_hr, inv_ess, log_roughness, mean_treecover, mean_road_cover)
  
  EST <- GRAD %*% FIT$beta + SLOPE.EST
  
  #species var
  VAR <- GRAD %*% FIT$vb %*% GRAD + SLOPE.VAR
  SE <- sqrt(VAR)
  
  #then you can calculate normal CIs and back-transform to the ridge density.
  # ctmm:::norm.ci(EST[,1], SE[,1])
  ctmm:::norm.ci(EST[,1], SE[,1])
  exp(ctmm:::norm.ci(EST[,1], SE[,1]))
  
  sp_level_ests$pred_ridge[i] <- exp(ctmm:::norm.ci(EST[,1], SE[,1]))[2]
  #print(paste(sp, sp_level_ests$mean_ridge[i], sp_level_ests$pred[i]))
}

phylo$tip.label <- replace(phylo$tip.label, which(phylo$tip.label=="Canis_mesomelas"), "Lupulella_mesomelas")
phylo$tip.label <- replace(phylo$tip.label, which(phylo$tip.label=="Pseudalopex_vetulus"), "Lycalopex_vetula")
names(phylo$tip.label) <- NULL

#avg ridge vals
d_avg <- sp_level_ests %>%
  #filter(sp %ni% c("Canis dingo", "Canis lupus x lycaon")) %>%
  mutate(sp=str_replace(sp, " ", "_")) %>%
  dplyr::select(c("sp", "ridge_dens_est")) %>%
  rbind(cbind(sp="Phataginus_tetradactyla",ridge_dens_est=NA),
        cbind(sp="Halichoerus_grypus",ridge_dens_est=NA),
        cbind(sp="Enhydra_lutris",ridge_dens_est=NA),
        cbind(sp="Herpestes_sanguineus",ridge_dens_est=NA),
        cbind(sp="Ursus_americanus",ridge_dens_est=NA),
        cbind(sp="Hyaena_hyaena",ridge_dens_est=NA)) #this paired with color scheme manually forces them to be grey

svl_avg <- parse_number(d_avg$ridge_dens_est)
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

#get individual ridge/UD plots

for (i in seq_along(AKDES)) {
  
  length_m <-  tryCatch({rgeos::gLength(allridge[[i]]) },
                        error= function(e) { length_m <- 0 })
  dens_m <- length_m^2/summary(AKDES[[i]], units=F)$CI[[2]]
  print(paste(AKDES[[i]]@info$identity, summary(AKDES[[i]], units=F)$CI[[2]], length_m, dens_m))
  
  
}

#load akdes/allridge from pre-calculated files (will add og code soon)

png("dog.png")
plot(dogAKDE, units=F, col.DF="grey50")
plot(dogridge, add=T)
dev.off()

png("cat.png")
plot(catAKDE, units=F, col.DF="grey50")
plot(catridge, add=T)
dev.off()

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
  theme(text=element_text(size=14))

gheatmap(p, phylo_dat,
         offset = 20,
         width=0.3,
         legend_title="Ridge Density") +
  #scale_fill_brewer(type="qual", palette=1, na.value="white") +
  scale_fill_gradientn(colors=viridis::plasma(n=1000), 
                       na.value = "white") +
  theme(plot.margin = unit(c(0.1, 8, 0, 0.5), "cm"),
        legend.position = c(1.03,0.45), legend.direction="vertical",
        legend.key = element_rect(color="black", fill = NA),
        legend.key.height = unit(2, 'cm'), text=element_text(size=10),
        axis.title.x = element_text(hjust=0.03, vjust=15)) +
  coord_cartesian(clip="off") +
  guides(fill = guide_colorbar(override.aes=list(fill=NA),
                               title.position = "right", title.theme = element_text(angle=270),
                               title.hjust = 0.5, title.vjust = -7.3)) +
  ggimage::geom_image(x=65, y=-5, 
             image="dog.png", 
             size=0.3) +
  ggimage::geom_image(x=65, y=-38, 
                      image="cat.png", 
                      size=0.3) +
  labs(fill="Ridge Density")


