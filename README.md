# Spatial memory, hunting strategies, and preferred travel routes in mammalian carnivores

<b>Authors:</b> W.F. Fagan and many, many, many others

<insert abstract here>

This repo contains the code and additional data required to reproduce the results in this manuscript. More detail on each file is provided below:

## R Files
- calculate_fits.R: Part 1 of ridge density calculations. Used to fit and save CTMM models for each individual. 
- getridgefromfits.R. Part 2 of ridge density calculations. Uses fit from calculate_fits.R to get ridge density for each individual.
- ranef.rma.mv.R: Hacked metafor function to return covariance matrix (by Chris Fleming)

figure scripts: contains R code to run analyses and reproduce all figures in text and supplements.

main figures: .png files of all figures in main text