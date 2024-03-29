This file contains the code used for my Bayesian hierarchical analysis.


As indicated in the Appendix 1 of my master thesis, the hierarchical models are written in the Stan files:

- Model_1.stan contains our first hierarchical model, used for Model 1 of the master thesis.

- Model_2.stan contains our second hierarchical model, used for Model 2 and Model 2bis of the master thesis.


These models are applied to the data of the multifaceted program and their results are displayed with the files:

- Model_1.R for our Model 1.

- Model_2_no_baseline.R for our Model 2.

- Model_2_baseline.R for our Model 2bis.


In order to run a model, you just have to follow the steps indicated in the Appendix 1 of my master thesis (installing
Rstan and running the R file only).

The data file of Banerjee et al. (2015), necessary to run these models, is available at 
https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NHIXNT
We used in particular the household level data file "pooled_hh_postanalysis.dta", that can 
be obtained (in the folder ScienceDataRelease/data_modified) by running the Stata do files 
of the authors (that are in the folder ScienceDataRelease/dofiles).



