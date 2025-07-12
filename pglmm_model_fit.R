#phylogenetic GLMM model fit workflow (see figure scripts for more examples)
#code last modified jan 15 2024

#load necessary libraries
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
  filter(`kept (y/n)` == "y") %>% #keep only individuals acceptable for analysis
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

#phylogeny from https://vertlife.org/data/mammals/
#download tree from: https://github.com/n8upham/MamPhy_v1/blob/master/_DATA/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre
tree_all <- read.nexus("/tree/file/here")
#recommended to use the dna only one but we had to use the full tree
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
