Geostatistical inference under preferential sampling
================
Daniel Dinsdale & Matias Salibian
2018-04-30

Geostatistical inference under preferential sampling
----------------------------------------------------

This repository contains `R` code illustrating the computations in [Dinsdale, D. and Salibian-Barrera, M. (2018) Methods for preferential sampling in geostatistics](https://doi.org/10.1111/rssc.12286) *Journal of the Royal Statistical Society Series C*. We use [TMB](https://github.com/kaskr/adcomp/wiki) and [INLA](http://www.r-inla.org/). For the latter we adapted the scripts publicly available [here](http://www.r-inla.org/examples/case-studies/diggle09/simulation1), with contributions from an anonymous referee.

What is available in this repository is:

-   A script to generate synthetic data following the Poisson preferential sampling model [here](RFiles/dataExample.R)
-   A script to obtain parameter estimates using INLA in R [here](RFiles/INLAExample.R)
-   A script to obtain parameter estimates using TMB in R [here](RFiles/TMBExample.R) and [here](RFiles/TMBExample.cpp)
