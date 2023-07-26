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

#ridge <- ridge %>% filter(sp %ni% c("Vulpes bengalensis","Leopardus geoffroyi"))

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

###################################################
#panel a
subs="FULL" 
effsizes <- data.frame(subset=NA, var=NA, low=NA, est=NA, high=NA, type=NA, color=NA) %>%
  drop_na()

effsizes <- effsizes %>%
  rbind(c(subs, "Phylogenetic", append((exp(qnorm(1-0.05/2)*c(-1,1)*sqrt(FIT$sigma2[[1]])))-1, NA, after=1), 
          "Biological Variance", "Biological Variance")) %>%
  rbind(c(subs, "Individual", append((exp(qnorm(1-0.05/2)*c(-1,1)*sqrt(FIT$sigma2[[2]])))-1, NA, after=1), 
          "Biological Variance", "Biological Variance")) %>%
  rbind(c(subs, "Phylogenetic\n(Clade)", exp(ctmm:::norm.ci(DIFF, VAR.DIFF))-1, 
          "Model Coefficients and Effects", "Biological Variance"))

tidy <- broom::tidy(FIT) %>%
  mutate(ub = estimate + (std.error*qnorm(0.975)),
         lb = estimate - (std.error*qnorm(0.975)))

for (i in seq_along(tidy$term)) {
  
  if (tidy$term[[i]] == "speed_est") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Speed\n(m/s)", 1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "log_mass") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Mass (g)", (1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])))/exp(1), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "pursuit1") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Hunting Movement \n(Pursuit)", 
                        1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "disruptfast1") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Hunting Movement \n(Disruptive Fast)", 
                        1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "slowwalking1") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Hunting Movement \n(Slow Walking)", 
                        1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  
  if (tidy$term[[i]] == "log_hr") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Home Range \nArea (m²)", 
                        (1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])))/exp(1), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "log_roughness") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Terrain \nRoughness (m)", 
                        (1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])))/exp(1), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "log_hfi") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Human Footprint Index (%)", 
                        (1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])))/exp(1), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "mean_treecover") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Tree Cover (%)", 1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "mean_road_cover") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Road Cover (%)", 1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
  if (tidy$term[[i]] == "seasonality_dhi_gpp") {
    effsizes <- rbind(effsizes, 
                      c(subs, "Dynamic Habitat Index\n(Gross Primary Productivity)\n(kg C/m²))", 
                        1 - exp(c(tidy$lb[[i]], tidy$estimate[[i]], tidy$ub[[i]])), 
                        "Model Coefficients and Effects", "Model Coefficients and Effects"))
  }
  
}

colnames(effsizes) = c("subset", "var", "low", "est", "high", "type", "color")

panel_a <- effsizes %>%
  mutate(subset = factor(subset, levels = c("FULL", "dna only", "speed", 
                                            "untouched dat", "full year", "shared landscape"),
                         labels = c("Full", "DNA Only Tree", "Speed Included", 
                                    "No Preprocessing", "Year or Longer", "Shared Landscapes")),
         var = factor(var, levels=rev(unique(effsizes$var)))) %>% #flip list
  mutate_at(c("low","est","high"), function(x){parse_number(x)*100}) %>%
  ggplot(mapping=aes(x=est, y = var, color = color)) +
  geom_point() +
  geom_errorbar(mapping=aes(xmin=low, xmax=high)) +
  geom_vline(mapping=aes(xintercept=0), linetype="dashed") +
  labs(x="Effect Size (Percent Increase in Ridge Density)", y = NULL, color = NULL,
       subtitle="\nFull (N=1219, C=16, F=18)") +
  scale_color_manual(values=c("Model Coefficients and Effects"="#5ac85a","Biological Variance"="#D2AAF0")) +
  facet_grid(type~., scales="free", space="free") +
  ggforce::facet_col(type~., scales = 'free', space = 'free') +
  theme_bw() +
  theme(legend.direction = "horizontal",
        legend.position = "none",
        text=element_text(size=15))

plot(panel_a)

###################################################
#clade reldiff across various subsets
ridgedensreldiffs <- read_csv("final/mod/fits/file/here") #add cleaned version to repo

panel_b <- ridgedensreldiffs %>%
  filter(version!="same species diff landscape") %>%
  mutate_at(c("clade_reldiff_low", "clade_reldiff_est", "clade_reldiff_high"), 
            function (x) {return(x-1)}) %>%
  mutate(version = factor(version, levels = c("FULL", "dna only tree", "speed included", "untouched dat",
                                              "duration >1yr", "shared landscape"))) %>%
  ggplot(mapping=aes(x=version, y=clade_reldiff_est)) +
  geom_point(size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(mapping=aes(ymin=clade_reldiff_low, ymax=clade_reldiff_high),
                width = 0.5, position=position_dodge(width=0.5)) +
  theme_bw() +
  labs(x=NULL, y="Ridge Density Percent Difference \n(Canid/Felid)") +
  scale_x_discrete(breaks = c("FULL", "dna only tree", "speed included", "untouched dat",
                              "duration >1yr", "shared landscape"), 
                   labels = c("Full", "DNA Only\nTree", "Speed\nIncluded", "No\nPreprocessing",
                              "Year or\nLonger", "Shared\nLandscapes")) +
  geom_hline(mapping=aes(yintercept = 0), linetype="dashed") +
  geom_text(aes(label = c("N=1219\nC=16\nF=18", "N=1011\nC=16\nF=18", "N=1183\nC=15\nF=17",
                          "N=1064\nC=16\nF=18", "N=338\nC=13\nF=16", "N=216\nC=7\nF=7"),
                y = rep(c(0.05), times = 6)), size=3) +
  # geom_text(aes(label = c("N=1219\nC=16\nc", "N=1011",  "N=1183", "N=1064", "N=338", "N=216"), 
  #               y = rep(c(0.1), times = 6)), size=3) +
  scale_y_continuous(limits=c(0, 0.5), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), 
                     labels = c("0%","10%","20%","30%", "40%", "50%")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot(panel_b)


#combine panels
cowplot::plot_grid(panel_a, panel_b,
                   ncol=1,
                   labels=c("A)", "B)"),
                   rel_heights=c(0.6,0.4))

