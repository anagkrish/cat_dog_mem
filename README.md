# Wild canids and felids differ in their reliance on travel routeways

<b>Authors:</b> W.F. Fagan and many, many, many others

[![DOI](https://zenodo.org/badge/659077878.svg)](https://zenodo.org/doi/10.5281/zenodo.10038542)

Diverse factors, including environmental features and cognitive processes, can drive animals' movements and space use, with far-reaching implications. For example, repeated use of individual-level travel routeways (directionally constrained but imperfectly aligned routes), which results in spatial concentration of movement and activity, can shape encounter-based processes including predation, mate-finding, and disease transmission. However, how much variation in routeway usage exists across species remains unknown. By analyzing GPS movement tracks for 1239 range-resident mammalian carnivores-representing 16 canid and 18 felid species from six continents-we found strong evidence of a clade-level difference in species' reliance on repeatedly used travel routeways. Across the global dataset, tracked canids had a 15% (±7 CI) greater density of routeways within their home ranges than did felids, rising to 33% (±16 CI) greater in landscapes shared with tracked felids. Moreover, comparisons within species across landscapes revealed broadly similar home range routeway densities despite habitat differences. On average, canids also reused their travel routeways more intensively than did felids, with hunting strategies and spatial contexts also contributing to the intensity of routeway usage. Collectively, our results suggest that key aspects of carnivore routeway-usage have an evolutionary component. Striking interspecific and clade-level differences in carnivores' reliance on reused travel routeways within home ranges identify important ways in which the movement patterns of real-world predators depart from classical assumptions of predator-prey theory. Because such departures can drive key aspects of human-wildlife interactions and other encounter-based processes, continued investigations of the relationships between movement mechanisms and space use are critical.

This repo contains the code and additional data required to reproduce the results in this manuscript, as well as high resolution versions of all the figures. More detail on each file is provided below:

## R Files
- pglmm_model_fit.R: Simple workflow to fit the PGLMM used in this manuscript. Includes information on how to access the phylogeny used in this model.
- calculate_fits.R: Part 1 of ridge density calculations. Used to fit and save CTMM models for each individual. 
- getridgefromfits.R. Part 2 of ridge density calculations. Uses fit from calculate_fits.R to get ridge density for each individual.
- ranef.rma.mv.R: Hacked metafor function to return covariance matrix (by Dr. Chris Fleming)
- spatial_controls.R: Code to extract all spatial information from bounding box of each individuals' tracks. Information on how to obtain each spatial variable is included in this code.
- simulations.R & simulation_functions.R: Exploring the behavior of ridge estimation with both the CTMM and adehabitatHR packages on simulated tracks.

figure scripts: contains R code to run analyses and reproduce all figures in text and supplements.

main figures: .png files of all figures in main text
supp figures: .png files of all figures in supplementary text

## Additional Notes
Some individuals, particularly those that were resampled to a minimum sampling interval of 1 minute, had unstable model fits due to the use of a non-default setting for the ctmm.select function, level=0.95, that we used to speed up computation time. Therefore, we fit an additional error model to all individuals that were resampled. Only resampled individuals with initially IID and OU fits changed when fit with the error model, so we only included error models for these 24 individuals. See calculate_fits.R for more details.

## R Environment and Packages

Implementation of all modified PGLMMs and generation of all figures were conducted in the R statistical package (v.4.2.2 R Core Team 2022) using the tidyverse (v. 1.3.2), ctmm (v.1.2.1), metafor (v. 3.8-1), nlme (v. 3.1-160), ape (v. 5.7-1), ggtree (v. 3.6.2) and phytools (v. 1.3.2). Initial model fits and ridge calculations of all individuals were conducted on the HPCC Zaratan R environment, using R (4.1.1 R Core Team 2022) and the tidyverse (v. 1.3.0), ctmm (v. 1.1.1), sf (v. 1.0-14), sp (v. 2.1-1), grDevices (v. 4.1.1), raster (v. 3.6-20) and geosphere (1.5-14). The authors acknowledge the University of Maryland supercomputing resources (http://hpcc.umd.edu) made available for conducting the research reported in this paper. 

R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical
  Computing, Vienna, Austria. URL https://www.R-project.org/.

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A, Henry L, Hester J,
  Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K,
  Vaughan D, Wilke C, Woo K, Yutani H (2019). “Welcome to the tidyverse.” Journal of Open Source Software,
  4(43), 1686. doi:10.21105/joss.01686.

Fleming CH, Calabrese JM (2023). ctmm: Continuous-Time Movement Modeling.
  https://github.com/ctmm-initiative/ctmm, https://groups.google.com/g/ctmm-user.

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor package. Journal of Statistical
  Software, 36(3), 1-48. https://doi.org/10.18637/jss.v036.i03.

Pinheiro J, Bates D, R Core Team (2022). nlme: Linear and Nonlinear Mixed Effects Models. R package
  version 3.1-160, https://CRAN.R-project.org/package=nlme.

Paradis E, Schliep K (2019). “ape 5.0: an environment for modern phylogenetics and evolutionary analyses in
  R.” Bioinformatics, 35, 526-528. doi:10.1093/bioinformatics/bty633

Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam. ggtree: an R package for
  visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods
  in Ecology and Evolution 2017, 8(1):28-36. doi:10.1111/2041-210X.12628

Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things).
  Methods Ecol. Evol. 3 217-223. doi:10.1111/j.2041-210X.2011.00169.x

Pebesma E, Bivand R (2005). “Classes and methods for spatial data in R.” R News, 5(2), 9-13.
  https://CRAN.R-project.org/doc/Rnews/.

Bivand R, Pebesma E, Gomez-Rubio V (2013). Applied spatial data analysis with R, Second edition.
  Springer, NY. https://asdar-book.org/.

Pebesma, E., & Bivand, R. (2023). Spatial Data Science: With Applications in R. Chapman and Hall/CRC.
  https://doi.org/10.1201/9780429459016

Pebesma, E., 2018. Simple Features for R: Standardized Support for Spatial Vector Data. The R Journal 10
  (1), 439-446, https://doi.org/10.32614/RJ-2018-009

Hijmans R (2022). geosphere: Spherical Trigonometry. R package version 1.5-14,
  https://CRAN.R-project.org/package=geosphere.

Hijmans R (2023). raster: Geographic Data Analysis and Modeling. R package version 3.6-20,
  https://CRAN.R-project.org/package=raster.
  