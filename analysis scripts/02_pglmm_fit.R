#phylogenetic GLMM model fit workflow (see figure scripts for more examples)
#fits both ridge density and ridge probability (with different back transformations)
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

######################################################
#### RIDGE DENSITY MODEL FIT
mods <- mods <- ~ log_mass_st + pursuit + disruptfast + slowwalking + log_hr_st + inv_ess + log_roughness_st + mean_treecover + mean_hfi + seasonality_dhi_gpp_st

#ridge density
FIT <- metafor::rma.mv(log_ridge,V=0,mods=mods,random=random,R=R,data=ridge)
#ridge probability
#FIT <- metafor::rma.mv(logit_ridge_prob_fixed50,V=0,mods=mods,random=random,R=R,data=ridge)

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

# relative impact on ridge probability (back transform)
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

#########################
# relative impact on ridge density (back transform)

ctmm:::norm.ci(DIFF, VAR.DIFF)
exp(ctmm:::norm.ci(DIFF, VAR.DIFF))

# relative impact on ridge probability (back transform)

# inv.logit <- function(x) {
#   exp(x)/(1+exp(x))
# }
# inv.logit(ctmm:::norm.ci(DIFF, VAR.DIFF))

############################
# save model fit coefficients to a spreadsheet (can be later used to make figures easily)
# 
# mod_comps <- data.frame(n=length(FIT[["data"]]$id),
#                          response_var=as.character(FIT[["call"]][["yi"]]), version=NA,
#                         terms=NA, aic=NA, bic=NA,
#                         clade_reldiff_low = NA, clade_reldiff_est = NA, clade_reldiff_high = NA,
#                         clade_reldiff_zval = NA, clade_reldiff_pval = NA,
#                         phylo_coeff=NA, point_coeff=NA,
#                         intrcpt_coeff=NA, intrcpt_se=NA, intrcpt_p=NA,
#                         log_mass_coeff=NA, log_mass_se=NA, log_mass_p=NA,
#                         log_homerange_coeff=NA, log_homerange_se=NA, log_homerange_p=NA,
#                         inv_ess_coeff=NA, inv_ess_se=NA, inv_ess_p=NA,
#                         log_roughness_coeff=NA, log_roughness_se=NA, log_roughness_p=NA,
#                         treecover_coeff=NA, treecover_se=NA, treecover_p=NA,
#                         mean_hfi_coeff=NA, mean_hfi_se=NA, mean_hfi_p=NA,
#                         road_cover_coeff=NA, road_cover_se=NA, road_cover_p=NA,
#                         seasonality_gpp_coeff = NA, seasonality_gpp_se = NA, seasonality_gpp_p = NA,
#                         speed_est_coeff=NA, speed_est_se=NA, speed_est_p=NA,
#                         movement_disruptfast_coeff=NA, movement_disruptfast_se=NA, movement_disruptfast_p=NA,
#                         movement_pursuit_coeff=NA, movement_pursuit_se=NA, movement_pursuit_p=NA,
#                         movement_slow_coeff=NA, movement_slow_se=NA, movement_slow_p=NA,
#                         cooperativity_solitary_coeff=NA, cooperativity_solitary_se=NA, cooperativity_solitary_p=NA)
# 
# r=length(mod_comps$version)
# 
# dat=FIT$data
# 
# mod_comps$version[r]="full"
# 
# tidy <- broom::tidy(FIT)
# 
# mod_comps$phylo_coeff[r]=FIT$sigma2[1] #corresponds to random coeff 1 (phylo), can't find se??
# mod_comps$point_coeff[r]=FIT$sigma2[2] #corresponds to random coeff 2 (point)
# 
# mod_comps$intrcpt_coeff[r]=(tidy%>% filter(term %in% c("intercept", "overall")))$estimate
# mod_comps$intrcpt_se[r]=(tidy%>% filter(term %in% c("intercept", "overall")))$std.error
# mod_comps$intrcpt_p[r]=(tidy%>% filter(term %in% c("intercept", "overall")))$p.value
# 
# if (grepl("log_mass", as.character(terms)) ==TRUE) {
#     mod_comps$log_mass_coeff[r]=(tidy%>% filter(term%in%c("log_mass","log_mass_st")))$estimate
#     mod_comps$log_mass_se[r]=(tidy%>% filter(term%in%c("log_mass","log_mass_st")))$std.error
#     mod_comps$log_mass_p[r]=(tidy%>% filter(term%in%c("log_mass","log_mass_st")))$p.value
# }
# 
# if (grepl("log_roughness", as.character(terms)) ==TRUE) {
#     mod_comps$log_roughness_coeff[r]=(tidy%>% filter(term%in%c("log_roughness", "log_roughness_st")))$estimate
#     mod_comps$log_roughness_se[r]=(tidy%>% filter(term%in%c("log_roughness", "log_roughness_st")))$std.error
#     mod_comps$log_roughness_p[r]=(tidy%>% filter(term%in%c("log_roughness", "log_roughness_st")))$p.value
# }
# 
# if (grepl("mean_treecover", as.character(terms)) ==TRUE) {
#     mod_comps$treecover_coeff[r]=(tidy%>% filter(term=="mean_treecover"))$estimate
#     mod_comps$treecover_se[r]=(tidy%>% filter(term=="mean_treecover"))$std.error
#     mod_comps$treecover_p[r]=(tidy%>% filter(term=="mean_treecover"))$p.value
# }
# 
# if (grepl("mean_hfi", as.character(terms)) ==TRUE) {
#     mod_comps$mean_hfi_coeff[r]=(tidy%>% filter(term=="mean_hfi"))$estimate
#     mod_comps$mean_hfi_se[r]=(tidy%>% filter(term=="mean_hfi"))$std.error
#     mod_comps$mean_hfi_p[r]=(tidy%>% filter(term=="mean_hfi"))$p.value
# }
# 
# if (grepl("mean_road_cover", as.character(terms)) ==TRUE) {
#     mod_comps$road_cover_coeff[r]=(tidy%>% filter(term=="mean_road_cover"))$estimate
#     mod_comps$road_cover_se[r]=(tidy%>% filter(term=="mean_road_cover"))$std.error
#     mod_comps$road_cover_p[r]=(tidy%>% filter(term=="mean_road_cover"))$p.value  
# }
# 
# if (grepl("disruptfast", as.character(terms)) ==TRUE) {
#     mod_comps$movement_disruptfast_coeff[r]=(tidy%>% filter(term=="disruptfast1"))$estimate
#     mod_comps$movement_disruptfast_se[r]=(tidy%>% filter(term=="disruptfast1"))$std.error
#     mod_comps$movement_disruptfast_p[r]=(tidy%>% filter(term=="disruptfast1"))$p.value
# }
# 
# if (grepl("pursuit", as.character(terms)) ==TRUE) {
#     mod_comps$movement_pursuit_coeff[r]=(tidy%>% filter(term=="pursuit1"))$estimate
#     mod_comps$movement_pursuit_se[r]=(tidy%>% filter(term=="pursuit1"))$std.error
#     mod_comps$movement_pursuit_p[r]=(tidy%>% filter(term=="pursuit1"))$p.value
# }
# 
# if (grepl("slowwalking", as.character(terms)) ==TRUE) {
#     mod_comps$movement_slow_coeff[r]=(tidy%>% filter(term=="slowwalking1"))$estimate
#     mod_comps$movement_slow_se[r]=(tidy%>% filter(term=="slowwalking1"))$std.error
#     mod_comps$movement_slow_p[r]=(tidy%>% filter(term=="slowwalking1"))$p.value
# }
# 
# if (grepl("hunting_cooperativity", as.character(terms)) ==TRUE) {
#     mod_comps$cooperativity_solitary_coeff[r]=(tidy%>% filter(term=="hunting_cooperativityTypically Solitary"))$estimate
#     mod_comps$cooperativity_solitary_se[r]=(tidy%>% filter(term=="hunting_cooperativityTypically Solitary"))$std.error
#     mod_comps$cooperativity_solitary_p[r]=(tidy%>% filter(term=="hunting_cooperativityTypically Solitary"))$p.value
# 
# }
# 
# if (grepl("speed_est", as.character(terms)) ==TRUE) {
#     mod_comps$speed_est_coeff[r]=(tidy%>% filter(term%in%c("speed_est", "speed_est_st")))$estimate
#     mod_comps$speed_est_se[r]=(tidy%>% filter(term%in%c("speed_est", "speed_est_st")))$std.error
#     mod_comps$speed_est_p[r]=(tidy%>% filter(term%in%c("speed_est", "speed_est_st")))$p.value
# }
# 
# if (grepl("seasonality_dhi_gpp", as.character(terms)) ==TRUE) {
#     mod_comps$seasonality_gpp_coeff[r]=(tidy%>% filter(term%in%c("seasonality_dhi_gpp",
#                                                                  "seasonality_dhi_gpp_st")))$estimate
#     mod_comps$seasonality_gpp_se[r]=(tidy%>% filter(term%in%c("seasonality_dhi_gpp",
#                                                               "seasonality_dhi_gpp_st")))$std.error
#     mod_comps$seasonality_gpp_p[r]=(tidy%>% filter(term%in%c("seasonality_dhi_gpp",
#                                                              "seasonality_dhi_gpp_st")))$p.value
# }
# 
# if (grepl("inv_ess", as.character(terms)) ==TRUE) {
#     mod_comps$inv_ess_coeff[r]=(tidy%>% filter(term=="inv_ess"))$estimate
#     mod_comps$inv_ess_se[r]=(tidy%>% filter(term=="inv_ess"))$std.error
#     mod_comps$inv_ess_p[r]=(tidy%>% filter(term=="inv_ess"))$p.value
# }
# 
# if (grepl("log_hr", as.character(terms)) ==TRUE) {
#     mod_comps$log_homerange_coeff[r]=(tidy%>% filter(term%in%c("log_hr","log_hr_st")))$estimate
#     mod_comps$log_homerange_se[r]=(tidy%>% filter(term%in%c("log_hr","log_hr_st")))$std.error
#     mod_comps$log_homerange_p[r]=(tidy%>% filter(term%in%c("log_hr","log_hr_st")))$p.value
# }
# 
# mod_comps$aic[r] = AIC(FIT)
# mod_comps$bic[r] = BIC(FIT)

# uncomment/run the following lines based on whether you're dealing with density or probability

### clade difference for ridge density
EFF <- ranef2(FIT)$phylo
REST <- EFF$est # random effects
RCOV <- EFF$COV # error covariance

CANINE <- dat$clade=="canidae"
CANINE <- unique(dat$phylo[CANINE])
FELINE <- dat$clade=="felidae"
FELINE <- unique(dat$phylo[FELINE])

iCOV <- ctmm::pd.solve(RCOV)
dW.dl <- iCOV %*% cbind( names(REST) %in% CANINE , names(REST) %in% FELINE )
M <- rbind( colSums( dW.dl[CANINE,] ) , colSums( dW.dl[FELINE,] ) )
lambda <- c( solve(M) %*% c(1,-1) )
W <- c(dW.dl %*% lambda)
names(W) <- rownames(dW.dl)

DIFF <- c(W %*% REST)
VAR.DIFF <- c(W %*% RCOV %*% W)

# relative impact on ridge density (back transform)
# 15%--27% increase in ridge density

clade_diff <- exp(ctmm:::norm.ci(DIFF, VAR.DIFF))

#change it for logit back-calc
# inv.logit <- function(x) {
#   exp(x)/(1+exp(x))
# }
# 
# clade_diff <- inv.logit(ctmm:::norm.ci(DIFF, VAR.DIFF))
# 
mod_comps$clade_reldiff_low[r] <- clade_diff[[1]]
mod_comps$clade_reldiff_est[r] <- clade_diff[[2]]
mod_comps$clade_reldiff_high[r] <- clade_diff[[3]]

#results!
#z-test for significance of random slopes (test against hypothesis that difference of random slopes between clades=0)
z_stat <- (DIFF) / sqrt(VAR.DIFF)
#two sided p-test
p = 2*pnorm(q=z_stat, lower.tail=FALSE)

mod_comps$clade_reldiff_zval[r] <- z_stat
mod_comps$clade_reldiff_pval[r] <- p

### clade difference for ridge probabilty
# EFF <- ranef2(FIT)$phylo
# REST <- EFF$est # random effects
# RCOV <- EFF$COV # error covariance
# 
# CANINE <- dat$clade=="canidae"
# CANINE <- unique(dat$phylo[CANINE])
# FELINE <- dat$clade=="felidae"
# FELINE <- unique(dat$phylo[FELINE])
# 
# iCOV <- ctmm::pd.solve(RCOV)
# dW.dl <- iCOV %*% cbind( names(REST) %in% CANINE , names(REST) %in% FELINE )
# M <- rbind( colSums( dW.dl[CANINE,] ) , colSums( dW.dl[FELINE,] ) )
# lambda <- c( solve(M) %*% c(1,-1) )
# W <- c(dW.dl %*% lambda)
# names(W) <- rownames(dW.dl)
# 
# DIFF <- c(W %*% REST)
# VAR.DIFF <- c(W %*% RCOV %*% W)
# 
# #change it for logit back-calc
# inv.logit <- function(x) {
#   exp(x)/(1+exp(x))
# }
# 
# GRAD <- c(1, log_mass_st, pursuit, disruptfast, slowwalking,
#           log_hr_st, inv_ess, log_roughness_st, mean_treecover, mean_hfi, seasonality_dhi_gpp_st)
# 
# EST <- array(0,2) # clade vector (canine,feline)
# 
# EST[1] <- GRAD %*% FIT$beta + W[CANINE] %*% REST[CANINE]
# EST[2] <- GRAD %*% FIT$beta - W[FELINE] %*% REST[FELINE]
# 
# COV <- matrix( GRAD %*% FIT$vb %*% GRAD ,2,2) # clade matrix
# COV[1,1] <- COV[1,1] + W[CANINE] %*% RCOV[CANINE,CANINE] %*% W[CANINE]
# COV[1,2] <- COV[1,2] - W[CANINE] %*% RCOV[CANINE,FELINE] %*% W[FELINE]
# COV[2,1] <- COV[1,2]
# COV[2,2] <- COV[2,2] + W[FELINE] %*% RCOV[FELINE,FELINE] %*% W[FELINE]
# 
# #Then if you can create whatever comparison function you want (e.g.):
# Ratio <- function(x) {  inv.logit(x[1]) / inv.logit(x[2]) }
# 
# #Calculate the gradient
# grad <- numDeriv::grad(Ratio,EST)
# 
# #and the variance will be
# VAR <- grad %*% COV %*% grad
# SE <- sqrt(VAR)
# 
# clade_diff <- ctmm:::lognorm.ci(Ratio(EST),VAR)
# 
# mod_comps$clade_reldiff_low[r] <- clade_diff[[1]]
# mod_comps$clade_reldiff_est[r] <- clade_diff[[2]]
# mod_comps$clade_reldiff_high[r] <- clade_diff[[3]]
# 
# #results!
# #z-test for significance of random slopes (test against hypothesis that difference of random slopes between clades=0)
# z_stat <- (DIFF) / sqrt(VAR.DIFF)
# #two sided p-test
# p = 2*pnorm(q=z_stat, lower.tail=FALSE)
# 
# mod_comps$clade_reldiff_zval[r] <- z_stat
# mod_comps$clade_reldiff_pval[r] <- p



