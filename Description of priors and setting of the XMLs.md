#Settings and priors for Bayesian Inference of SARS-CoV-2 datasets.
-A maximum-likelihood analysis was done in IQtree2 using GTR+F+G and 1000 bootstraps, then the resulting tree was analyzed in TempEst for molecular clock. 
Clock rate was assesed to be close to 8x10^-4 subs/site/year, which was used as a prior for further Bayesian Inference analyses. 

-For the priors and setting of the BEAST analyses, 8 datasets composed of **B.1, B.1.111, B.1.1.348, B.1.1.7, P.1, B.1.526 and B.1.621** which were found to be the most prevalent lineages
withing the general Colombian dataset, a total of 1000 sequences were retrieved from GISAID and included in the analysis, sequences with more than 5% Ns were excluded from the dataset.

-Birth-death Skyline model was used for Re determination for each dataset in Beast v2.6.5 using Hasegawa-Kishino-Yano model with gamma categories set to empirical.
*note: Best-fit model was found to be GTR but given the size of the dataset and the lack of computing power, I decided to use HKY.
-Clock Model: **Strick clock** with **uniform** distribution prior, lower = 8x10^-4 s/s/y; upper 1x10^-3 s/s/y.
-Become Uninfectious Rate: priors set to a **exponential** distribution with a mean of 36 and Standard deviation of 0.0 (as used by Naveca et al in https://rdcu.be/cpM9S)
-Reproductive Number: a **LogNormal** distritribution was used, with µ and s parameters set at 0.8 and 0.5 respectively; 
 Re was allowed to change at equally time-spaced points in a tri-monthly manner using a BirthRateTimeChanges TreeSlicer object within the skylinetools package (or at least tried to).
-Sampling Proportion: a **beta** distribution was used with α and β parameters set at 1 and 100 respectively, for the analysis we assumed a Sampling proportion consisting of 2 dimensions, with a proportion of 0 before date of the oldest sample(TMRCA) by using TreeSlicer.
MCMC Lenght: 100 million for each dataset, ESS for all parameters was >200
