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
mods <- ~ log_mass_st + pursuit + disruptfast + slowwalking + log_hr_st + inv_ess + log_roughness_st + mean_treecover + mean_hfi + seasonality_dhi_gpp_st

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

###################################################
#panel a
subs="full" #subset
effsizes <- data.frame(subset=NA, var=NA, low=NA, est=NA, high=NA, type=NA, color=NA) %>%
  drop_na()

get_coefs <- function(tidy) {
  
  x <- c(tidy$lb, tidy$estimate, tidy$ub)
  
  c <- abs(1-exp(x))*sign(x) 
  print(((abs(1-exp(x)))/exp(1))*sign(x))
  print(str_detect(tidy$term, "log"))
  
  return(c)
  
}


effsizes <- effsizes %>%
  rbind(c(subs, "Phylogenetic", append((exp(c(-1,1)*sqrt(FIT$sigma2[[1]])))-1, NA, after=1), 
          "Biological Variance", "Biological Variance")) %>%
  rbind(c(subs, "Individual", append((exp(c(-1,1)*sqrt(FIT$sigma2[[2]])))-1, NA, after=1), 
          "Biological Variance", "Biological Variance")) %>%
  rbind(c(subs, "Phylogenetic<br>(Clade)", exp(ctmm:::norm.ci(DIFF, VAR.DIFF))-1, 
          "Model Coefficients and Effects", "Biological Variance"))

tidy <- broom::tidy(FIT) %>%
  mutate(ub = estimate + (std.error*qnorm(0.975)),
         lb = estimate - (std.error*qnorm(0.975)))

for (i in seq_along(tidy$term)) {
  
  if (tidy$term[[i]] == "speed_est_st") {
    effsizes <- rbind(effsizes, 
                      #c(subs, "Speed\n(m/s)", 
                      c(subs, tidy$term[[i]], 
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  if (tidy$term[[i]] == "log_mass_st") {
    effsizes <- rbind(effsizes, 
                      c(subs, tidy$term[[i]], 
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  if (tidy$term[[i]] == "pursuit1") {
    effsizes <- rbind(effsizes, 
                      #c(subs, "Hunting Movement \n(Pursuit)", 
                      c(subs, tidy$term[[i]],
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  if (tidy$term[[i]] == "disruptfast1") {
    effsizes <- rbind(effsizes, 
                      #c(subs, "Hunting Movement \n(Disruptive Fast)", 
                      c(subs, tidy$term[[i]], 
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  if (tidy$term[[i]] == "slowwalking1") {
    effsizes <- rbind(effsizes, 
                      #c(subs, "Hunting Movement \n(Slow Walking)", 
                      c(subs, tidy$term[[i]],
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  
  if (tidy$term[[i]] == "log_hr_st") {
    effsizes <- rbind(effsizes, 
                      #c(subs, "Home Range Area (m²)", 
                      c(subs, tidy$term[[i]], 
                        (get_coefs(filter(tidy, term==tidy$term[[i]]))), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  
  if (tidy$term[[i]] == "log_roughness_st") {
    effsizes <- rbind(effsizes, 
                      #c(subs, "Terrain \nRoughness (m)", 
                      c(subs, tidy$term[[i]], 
                        (get_coefs(filter(tidy, term==tidy$term[[i]]))), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  # if (tidy$term[[i]] == "log_hfi") {
  #   effsizes <- rbind(effsizes, 
  #         c(subs, tidy$term[[i]], 
  #           (get_coefs(filter(tidy, term==tidy$term[[i]]))), 
  #           "Model Coefficient", "Model Coefficient"))
  # }
  
  if (tidy$term[[i]] == "mean_hfi") {
    effsizes <- rbind(effsizes, 
                      c(subs, tidy$term[[i]], 
                        (get_coefs(filter(tidy, term==tidy$term[[i]]))), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  if (tidy$term[[i]] == "mean_treecover") {
    effsizes <- rbind(effsizes, 
                      c(subs, tidy$term[[i]], 
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  if (tidy$term[[i]] == "mean_road_cover") {
    effsizes <- rbind(effsizes, 
                      c(subs, tidy$term[[i]], 
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
  if (tidy$term[[i]] == "seasonality_dhi_gpp_st") {
    effsizes <- rbind(effsizes, 
                      #c(subs, "Dynamic Habitat Index\n(Gross Primary Productivity)\n(kg C/m²))", 
                      c(subs, tidy$term[[i]], 
                        get_coefs(filter(tidy, term==tidy$term[[i]])), 
                        "Model Coefficient", "Model Coefficient"))
  }
  
}


colnames(effsizes) = c("subset", "var", "low", "est", "high", "type", "color")

#rename/reorganize
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

library(ggtext)

panel_a <- effsizes %>%
  filter(subset%in%c("full")) %>%
  mutate(subset = factor(subset, levels = c("full", "dna only", "speed inc", 
                                            "clean", "year", "same landscape"),
                         labels = c("Full", "DNA Only Tree", "Speed Included", 
                                    "No Preprocessing", "Year or Longer", "Shared Landscapes"))) %>%
  mutate(textcol = #phylo/rand effects
                   ifelse(var %in% c("Phylogenetic", "Individual", "Phylogenetic<br>(Clade)"), "**", 
                   #log transformed vars
                   ifelse(var %in% c("log_mass_st", "log_hr_st", "log_roughness_st"), "***",
                   #indicator vars
                   ifelse(var %in% c("pursuit1", "disruptfast1", "slowwalking1"), "_", 
                   #untransformed vars
                   ifelse(var %in% c("mean_hfi", "mean_treecover", 
                                     "seasonality_dhi_gpp_st"), "", NA))))) %>%
  mutate(label = paste(textcol, name, textcol, sep = ""),
         label = fct_reorder(label, rev(sort(as.character(label))))) %>%
  #mutate_at(c("low","est","high"), function(x){x*100}) %>%  #parse_number(x)*100}) %>%
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

plot(panel_a)

###################################################
#clade reldiff across various subsets
ridgedensreldiffs <- read_csv("mod_comps.csv")  %>% #summary statistics for each of the ridge subsets
  dplyr::select(c("version", "response_var", "terms",
                  "aic", "bic", "clade_reldiff_low", "clade_reldiff_est", "clade_reldiff_high"))

panel_b <- ridgedensreldiffs %>%
  filter(version %in% c("full", "dna only", "speed inc", "clean",
                        "duration >1 yr", "shared landscape")) %>%
  mutate_at(c("clade_reldiff_low", "clade_reldiff_est", "clade_reldiff_high"), 
            function (x) {return(x-1)}) %>%
  mutate(version = factor(version, levels = c("full", "dna only", "speed inc", 
                                              "clean", "duration > 1yr", "same landscape"))) %>%
  ggplot(mapping=aes(x=version, y=clade_reldiff_est)) +
  geom_point(size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(mapping=aes(ymin=clade_reldiff_low, ymax=clade_reldiff_high),
                width = 0.5, position=position_dodge(width=0.5)) +
  theme_bw() +
  labs(x=NULL, y="Ridge Density Percent Difference \n(Canid/Felid)") +
  scale_x_discrete(breaks = c("full", "dna only", "speed inc", "clean",
                              "duration >1 yr", "shared landscape"), 
                   labels = c("Full", "DNA Only\nTree", "Speed\nIncluded", "No\nPreprocessing",
                              "Year or\nLonger", "Shared\nLandscapes"))  +
  geom_hline(mapping=aes(yintercept = 0), linetype="dashed") +
  geom_text(aes(label = c("N=1239\nC=16\nF=18", "N=1032\nC=16\nF=18", "N=1201\nC=15\nF=17", 
                          "N=1070\nC=16\nF=18", "N=342\nC=13\nF=16", "N=219\nC=7\nF=7"),
                y = rep(c(0.55), times = 6)), size=3) +
  scale_y_continuous(limits=c(0, 0.6), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 
                     labels = c("0%","10%","20%","30%", "40%", "50%", "60%")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot(panel_b)


#combine panels
cowplot::plot_grid(panel_a, panel_b,
                   ncol=1,
                   labels=c("A)", "B)"),
                   rel_heights=c(0.6,0.4))

