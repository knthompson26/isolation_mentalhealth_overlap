# Genetic and environmental overlap between social isolation and mental health symptoms from childhood to early adulthood 
[Katherine N Thompson](https://twitter.com/KTNThompson)

September, 2022

***

This repository holds code for for manuscript titled "Social isolation and poor mental health in young people: Testing genetic and environmental influences in a longitudinal cohort study". This paper is published in European Child & Adolescent Psychiatry and can be found here: https://doi.org/10.1007/s00787-024-02573-w 

We conducted two stages of analyses. First, we computed multivariate ACE models with a longitudinal Cholesky decomposition across ages 5, 7, 10, and 12 years to test for genetic and environmental influences on social isolation across childhood. Second, we computed an independent pathway model (IPM) to assess genetic and environmental influences on the overlap between social isolation and mental health symptoms. Third, I have  provided preliminary functions that you can use to compute a longitudinal IPM. All data from the E-Risk Study. 

Analyses for this project were conducted in R (Version 4.0.3). I have listed all the names of the scripts used below, with a brief explanation of what each script entails. 

***

**Cholesky decomposition of social isolation in childhood**

1. isolation_cholesky.Rmd
This script contains code for a saturated multivariate model, saturated submodels, ACE model with Cholesky decomposition and correlated factors solutions, AE model with Cholesky decomposition and correlated factors solution. Data are social isolation sum scores at ages 5, 7, 10, and 12 years from the E-Risk study. Interpretation and commentary is provided throughout script. 


***

**Genetic and environmental overlap between social isolation and mental health problems (IPM)**

1. isolation_mhealth_overlap_ACE_univariate.Rmd
This script contains univariate ACE models for all variables (social isolation, depression symptoms, conduct problems, psychotic experiences) at age 12 and 18.

2. isolation_mhealth_overlap_ACE_bivariate.Rmd 
This script contains bivariate ACE models per phenotype so assess changes from age 12 to  18. 

3. isolation_mhealth_overlap_ACE_sexdiffs.Rmd 
This script contains univariate ACE models with heterogeneity paths and scalars applied to check for sex differences for all variables. 

4. isolation_mhealth_overlap_IPM.Rmd
This script contains longitudinal independent pathway models for all variables with appropriate heterogeneity paths and scalars applied. 

5. isolation_mhealth_overlap_sensitivity.Rmd 
This script contains two sensitivity analyses: 1. Psychometric model for social isolation and conduct problems. 2. Sex differences for conduct problems at age 18.  

***

**IPM functions**

1. isolation_mhealth_functions.R 
This script contains output functions used in the above scripts to produce formatted output tables using path tracing mathematics: this includes heritability estimates, path estimates (standardised and squared), proportion of common and specific ACE effects per variable, and the percentage contribution of ACE to all bivariate phenotypic correlations in the model. 

2. ipm_functions.R
This script contains the IPM_ace() function that computes longitudinal independent pathway models at two time points without needing to specify the full twin modelling code. This function uses OpenMx functions without the additional details that is often the same across all twin modelling scripts. My goal is to expand this function to be for all twin models (univariate and multivariate) and then to create a package based on these functions. 

***
