Inference under preferential sampling
================
Daniel Dinsdale & Matias Salibian
2017-11-29

Geostatistical inference under preferential sampling
----------------------------------------------------

This repository contains `R` code illustrating the computations in *Dinsdale, D. and Salibian-Barrera, M. (2017) A note on geostatistical inference under preferential sampling (submitted)*. We use [TMB](https://github.com/kaskr/adcomp/wiki) and [INLA](http://www.r-inla.org/). For the latter we adapted the scripts publicly available [here](http://www.r-inla.org/examples/case-studies/diggle09/simulation1).

What is available in this repository is:

-   A script to generate synthetic data following the Poisson preferential sampling model [here](RFiles/dataExample.R)
-   A script to obtain parameter estimates using INLA in R [here](RFiles/INLAExample.R)
-   A script to obtain parameter estimates using TMB in R [here](RFiles/TMBExample.R) and [here](RFiles/TMBExample.cpp)
