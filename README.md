# Details
Although the description of the repository is about generation of correlated variables this is part of a simulation study that I did in the past. The simulation study was intended to check wether generalized estimating ecuations (GEE) is better that generalized linear models (GLM) when data is longitudinal. For this purpose I performed a simulation study in which I generated longitudinal response data, that is, correlated intra-subject, and then fitted the models under study. The flow of the simulation is detailed in the report (simulation_study.pdf) but here a brief review:  


1. Define a linear model where the estimators are known in advanced and the response is correlated (Binomial or Poisson).
2. Fit a GLM and a GEEs for different correlation matrix.
3. Compute MSE, SD and confident interval for the estimations in step 2.
4. Check if introducing a correlation matrix performs better than don't introducing it.

Given that the hard task of the simulation was to implement the algorithms of correlated data generation I titled the repository that way. The algorithms were developed in the next papers:  

- Park, C.G., Park, T. y Shin, D.W. (1996). A simple method for generating
correlated binary variates. The American Statistician, 50, (4), 306-310.  
- Park, C.G. y Shin, D.W. (1998). An algorithm for generating correlated randon
variables in a class of infinitely divisible distributions. Journal Statist. Comput.
Simul., 61, 127-139.


A complete explanation of the simulation along with the obtained results and the explored papers can be found in the sumulation_study.pdf file.    


**Note**: the language of the study is spanish.


requirementelibrary(dplyr)
library(geepack)
library(readr)
library(igraph)