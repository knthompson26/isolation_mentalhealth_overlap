# ---
#   title: "Functions for all scripts looking at the overlap between social isolation and mental health" 
# ---
  
# clear the global environment
remove(list = ls())

# Load packages
library(knitr)
library(psych)
library(OpenMx)
library(tidyr)
library(tidyverse)
library(dplyr)     # conflicts with tidyverse for e.g. rename and row_number

# Univariate 

## ACE estimates
ACE_esitmates_table <- function(data, variable, age){
  table <- as.tibble(data$CI) %>%
    mutate(ACE = c("h2", "c2", "e2")) %>%
    mutate(Variable = variable) %>%
    mutate(Age = age) %>%
    select(Variable, Age, ACE, lbound, estimate, ubound) %>%
    mutate(across(is.numeric, round, digits = 3))
  
  return(table)
}

###########################################################################################

# Sex diffs 
ACEhet_esitmates_table <- function(data, variable, age){
    table <- as.tibble(data$CI) %>%
      mutate(ACE = c("h2m", "c2m", "e2m", "h2f", "c2f", "e2f")) %>%
      mutate(Variable = variable) %>%
      mutate(Age = age) %>%
      select(Variable, Age, ACE, lbound, estimate, ubound) %>%
      mutate(across(is.numeric, round, digits = 3))
    
    return(table)
}

###########################################################################################

# Bivariate script

## ACE estimates for bivariate model - same variable at age 12 and 18
ACE_esitmates_table_crosstime <- function(data, variable){
  table <- as.tibble(data$CI) %>%
    mutate(ACE = c(rep("h2",2), rep("c2",2), rep("e2",2))) %>%
    mutate(Variable = variable) %>%
    mutate(Age = rep(c("12", "18"), 3)) %>%
    select(Variable, Age, ACE, estimate, lbound, ubound) %>% 
    mutate(across(is.numeric, round, digits = 3))
  
  return(table)
}

###########################################################################################

### ADD CHOLESKY ESTIMATES FUNCTION HERE

###########################################################################################

# IPM
# This function can be used to create output from a two time-point 4 variable IPM model. This model includes 2 common *ACE* factors, 1 represents shared A/C/E at all time points and all variables in the model, the second common factor represents shared A/C/E between variables at the second point in time. This output function can produce the following: 
# * "uniACEestimates" represents univariate ACE estimates. These are the model estimates h2, c2, and e2 for each variable incuded in the model. To call this output uniACEestimates = TRUE.
# * "pathestimates" represents common and specific factor path estimates for A/C/E. This is all the standardised and squared path estimates in the model, for common and specific factors of A, C, and E. To call this output pathestimates = TRUE.
# * "A_percent", "C_percent", and "E_percent" represent the percentage of total A, C, or E that are due to common or specific effects, respectively. To call this output A_percent = TRUE, C_percent = TRUE, E_percent = TRUE. 
# * "common_factor_contribution" represents the % of A/C/E that contributes to the phenotypic correlation between the variables in the model. This is estimated using path tracing rules, and done in 2x2 correlations between all combinations of the variables in the model. To call this output common_factor_contribution = TRUE. 
# Before running this function, you will need an OpenMx model fit using mxTryHard(model, intervals = TRUE), or mxRun(model, intervals = TRUE). The specification of this model will have to be consistent with the scripts in this project. I.e. the mxModel should be called "ACE", and the corresponding mxAlgebra parameters should be called "Rph", "h2", "c2", "e2". 
# Objects needed are:
# data = mxTryHard(model, intervals = TRUE). e.g. model_fit <- mxTryHard(model, intervals = TRUE), data = model_fit. 
# variable1 = first variable at time point 1.  e.g. variable1 = "Isolation12". This is *your* choice of naming for these variables. 
# variable2 = second variable at time point 1. e.g. variable2 = "Anxiety12"
# variable3 = first variable at time point 2.  e.g. variable3 = "Isolation18"
# variable4 = second variable at time point 2. e.g. variable4 = "Anxiety18"

# To estimate the contribution of genetic factors to the phenotypic correlation, we need to path trace the following:
# * cor(isolation12, anxiety12)   = sqrt(Isolation12_A1 * Anxiety12_A1)
# * cor(isolation12, isolation18) = sqrt(Isolation12_A1 * Isolation18_A1)
# * cor(isolation12, anxiety18)   = sqrt(Isolation12_A1 * Anxiety18_A1)
# * cor(anxiety12, isolation18)   = sqrt(Anxiety12_A1   * Isolation18_A1)
# * cor(isolation18, anxiety18)   = sqrt(Isolation18A1 * Anxiety18_A1) + sqrt(Isolation18_A2 * Anxiety18_A2)

IPM_esitmates_ACE <- function(data, variable1, variable2, variable3, variable4, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rph <- mxEval(ACE.Rph, data)
  
  # create table containing the univariate ACE estimates
  h2 <- as.data.frame(diag(round(mxEval(ACE.h2, data),3))) # heritability
  colnames(h2) <- c("h2")
  
  c2 <- as.data.frame(diag(round(mxEval(ACE.c2, data),3))) # shared environment 
  colnames(c2) <- c("c2")
  
  e2 <- as.data.frame(diag(round(mxEval(ACE.e2, data),3))) # unique environment
  colnames(e2) <- c("e2")
  
  uniACEestimates_table <- cbind(h2, c2, e2) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4)) %>%
    mutate(model = model) %>%
    select(model, variable, h2, c2, e2) %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table)} 
  
  # create data frame with all the path estimates and confidence intervals
  ACE_CI <- as.data.frame(summary$CI)
  ACE_CI <- ACE_CI %>%
    mutate(parameter = rownames(ACE_CI)) %>%
    as.tibble()
  
  # from this, create a table containing estimates for common factors (squared and standardised)
  ACE_CI.common <- ACE_CI %>%
    filter(grepl('ACE.stac2|ACE.stcc2|ACE.stec2', parameter)) %>% # squared and standardised estimates 
    mutate(ACE = c(rep("A1",4), rep("A2",4), rep("C1",4), rep("C2",4), rep("E1",4), rep("E2",4))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4),6)) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_CI.specific <- ACE_CI %>%
    filter(grepl('ACE.stas2|ACE.stcs2|ACE.stes2', parameter) & estimate>0) %>% # delete the empty cells in the matrix
    mutate(ACE = c(rep("A",4), rep("C",4), rep("E",4))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4),3)) %>%
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # combine the common and specific path estimates
  pathestimates_table <- rbind(ACE_CI.common, ACE_CI.specific)
  
  pathestimates_table_rounded <- pathestimates_table %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages <- pathestimates_table %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           -Total)
  if(A_percent == TRUE){print(A_percentages)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages <- pathestimates_table %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           -Total)
  
  if(C_percent == TRUE){print(C_percentages)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           -Total)
  
  if(E_percent == TRUE){print(E_percentages)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent <- pathestimates_table %>%
    filter(grepl('A1|A2', ACE)) %>%
    filter(ACE != "A2" | (variable != variable1 & variable != variable2)) %>% # remove common factor A2 for variable 1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rph[2,1]) %>% #square root common paths, * together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rph[4,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rph[3,2]) %>% 
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rph[4,2]) %>% 
    mutate(variable3_variable4 = (sqrt(variable3_A1 * variable4_A1) + sqrt(variable3_A2 * variable4_A2))/Rph[4,3]) %>% 
    select(variable1_variable2, variable1_variable3, variable1_variable4, variable2_variable3, variable2_variable4, variable3_variable4) %>%
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ac % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) 
  
  # Second, contribution of C to the phenotypic correlation
  Rc_contribution_percent <- pathestimates_table %>%
    filter(grepl('C1|C2', ACE)) %>%
    filter(ACE != "C2" | (variable != variable1 & variable != variable2)) %>% # remove common factor C2 for variable 1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename the variables so we can apply all equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_C1 * variable2_C1)/Rph[2,1]) %>%
    mutate(variable1_variable3 = sqrt(variable1_C1 * variable3_C1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_C1 * variable4_C1)/Rph[4,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_C1 * variable3_C1)/Rph[3,2]) %>% 
    mutate(variable2_variable4 = sqrt(variable2_C1 * variable4_C1)/Rph[4,2]) %>% 
    mutate(variable3_variable4 = (sqrt(variable3_C1 * variable4_C1) + sqrt(variable3_C2 * variable4_C2))/Rph[4,3]) %>% 
    select(variable1_variable2, variable1_variable3, variable1_variable4, variable2_variable3, variable2_variable4, variable3_variable4) %>%
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Cc % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) 
  
  # Third, contribution of E to the phenotypic correlation
  Re_contribution_percent <- pathestimates_table %>%
    filter(grepl('E1|E2', ACE)) %>%
    filter(ACE != "E2" | (variable != variable1 & variable != variable2)) %>% # remove common factor E2 for variable 1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename the variables so we can apply all equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_E1 * variable2_E1)/Rph[2,1]) %>%
    mutate(variable1_variable3 = sqrt(variable1_E1 * variable3_E1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_E1 * variable4_E1)/Rph[4,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_E1 * variable3_E1)/Rph[3,2]) %>% 
    mutate(variable2_variable4 = sqrt(variable2_E1 * variable4_E1)/Rph[4,2]) %>% 
    mutate(variable3_variable4 = (sqrt(variable3_E1 * variable4_E1) + sqrt(variable3_E2 * variable4_E2))/Rph[4,3]) %>% 
    select(variable1_variable2, variable1_variable3, variable1_variable4, variable2_variable3, variable2_variable4, variable3_variable4) %>%
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ec % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) 
  
  # bind the contribution of A, C, and E together into one table
  contribution_percent_table <- plyr::join_all(list(Ra_contribution_percent, Rc_contribution_percent, Re_contribution_percent),
                                               by = "phenotypic correlation")
  
  contribution_percent_table <- contribution_percent_table %>% 
    mutate(total = `Ac % contribution` + `Cc % contribution` + `Ec % contribution`) %>% 
    mutate(across(is.numeric, round, digits = 3))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}

###########################################################################################

# Same function for *AE* model
IPM_esitmates_AE <- function(data, variable1, variable2, variable3, variable4, model, uniACEestimates, pathestimates, A_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rph <- mxEval(ACE.Rph, data)
  
  # create table containing the univariate ACE estimates
  h2 <- as.data.frame(diag(round(mxEval(ACE.h2, data),3))) # heritability
  colnames(h2) <- c("h2")
  
  c2 <- as.data.frame(diag(round(mxEval(ACE.c2, data),3))) # shared environment - kept in as a check that C is set to 0
  colnames(c2) <- c("c2")
  
  e2 <- as.data.frame(diag(round(mxEval(ACE.e2, data),3))) # unique environment
  colnames(e2) <- c("e2")
  
  uniACEestimates_table <- cbind(h2, c2, e2) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4)) %>%
    mutate(model = model) %>%
    select(model, variable, h2, c2, e2) %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table)} 
  
  # create data frame with all the path estimates and confidence intervals
  AE_CI <- as.data.frame(summary$CI)
  AE_CI <- AE_CI %>%
    mutate(parameter = rownames(AE_CI)) %>%
    as.tibble()
  
  # from this, create a table containing estimates for common factors (squared and standardised)
  AE_CI.common <- AE_CI %>%
    filter(grepl('ACE.stac2|ACE.stec2', parameter)) %>% # squared and standardised estimates 
    mutate(AE = c(rep("A1",4), rep("A2",4), rep("E1",4), rep("E2",4))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4),4)) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, AE, variable, lbound, estimate, ubound) 
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  AE_CI.specific <- AE_CI %>%
    filter(grepl('ACE.stas2|ACE.stes2', parameter) & estimate>0) %>% # delete the empty cells in the matrix
    mutate(AE = c(rep("A",4), rep("E",4))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4),2)) %>%
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    select(model, factor, AE, variable, lbound, estimate, ubound) 
  
  # combine the common and specific path estimates
  pathestimates_table <- rbind(AE_CI.common, AE_CI.specific)
  
  pathestimates_table_rounded <- pathestimates_table %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages <- pathestimates_table %>% 
    filter(grepl('A1|A2|A', AE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = AE, values_from = estimate) %>%
    mutate(Total = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           -Total)
  if(A_percent == TRUE){print(A_percentages)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', AE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = AE, values_from = estimate) %>%
    mutate(Total = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           -Total)
  
  if(E_percent == TRUE){print(E_percentages)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent <- pathestimates_table %>%
    filter(grepl('A1|A2', AE)) %>%
    filter(AE != "A2" | (variable != variable1 & variable != variable2)) %>% # remove common factor A2 for variable 1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    pivot_wider(names_from = c(variable, AE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rph[4,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rph[3,2]) %>% 
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rph[4,2]) %>% 
    mutate(variable3_variable4 = (sqrt(variable3_A1 * variable4_A1) + sqrt(variable3_A2 * variable4_A2))/Rph[4,3]) %>% 
    select(variable1_variable2, variable1_variable3, variable1_variable4, variable2_variable3, variable2_variable4, variable3_variable4) %>%
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ac % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4)))
  
  # Second, contribution of E to the phenotypic correlation
  Re_contribution_percent <- pathestimates_table %>%
    filter(grepl('E1|E2', AE)) %>%
    filter(AE != "E2" | (variable != variable1 & variable != variable2)) %>% # remove common factor E2 for variable 1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename the variables so we can apply all equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    pivot_wider(names_from = c(variable, AE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_E1 * variable2_E1)/Rph[2,1]) %>%
    mutate(variable1_variable3 = sqrt(variable1_E1 * variable3_E1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_E1 * variable4_E1)/Rph[4,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_E1 * variable3_E1)/Rph[3,2]) %>% 
    mutate(variable2_variable4 = sqrt(variable2_E1 * variable4_E1)/Rph[4,2]) %>% 
    mutate(variable3_variable4 = (sqrt(variable3_E1 * variable4_E1) + sqrt(variable3_E2 * variable4_E2))/Rph[4,3]) %>% 
    select(variable1_variable2, variable1_variable3, variable1_variable4, variable2_variable3, variable2_variable4, variable3_variable4) %>%
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ec % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) 
  
  # bind the contribution of A, C, and E together into one table
  contribution_percent_table <- plyr::join_all(list(Ra_contribution_percent, Re_contribution_percent),
                                               by = "phenotypic correlation")
  
  contribution_percent_table <- contribution_percent_table %>% 
    mutate(total = `Ac % contribution` + `Ec % contribution`) %>% 
    mutate(across(is.numeric, round, digits = 3))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}

###########################################################################################

# # For 8 variable model 
# data = IPM_allvar_noanx_fit
# variable1 = "Isolation12"
# variable2 = "Depression12"
# variable3 = "Conduct12"
# variable4 = "Psychosis"
# variable5 = "Isolation18"
# variable6 = "Depression18"
# variable7 = "Conduct18"
# variable8 = "Psychosis18"
# model = "ACE"
# uniACEestimates = TRUE
# pathestimates = TRUE
# A_percent = TRUE
# C_percent = TRUE
# E_percent = TRUE
# common_factor_contribution = TRUE

IPM_esitmates_ACE_8var <- function(data, variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rph <- mxEval(ACE.Rph, data)
  
  # create table containing the univariate ACE estimates
  h2 <- as.data.frame(diag(round(mxEval(ACE.h2, data),3))) # heritability
  colnames(h2) <- c("h2")
  
  c2 <- as.data.frame(diag(round(mxEval(ACE.c2, data),3))) # shared environment 
  colnames(c2) <- c("c2")
  
  e2 <- as.data.frame(diag(round(mxEval(ACE.e2, data),3))) # unique environment
  colnames(e2) <- c("e2")
  
  uniACEestimates_table <- cbind(h2, c2, e2) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(model = model) %>%
    select(model, variable, h2, c2, e2) %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table)} 
  
  # create data frame with all the path estimates and confidence intervals
  ACE_CI <- as.data.frame(summary$CI)
  ACE_CI <- ACE_CI %>%
    mutate(parameter = rownames(ACE_CI)) %>%
    as.tibble()
  
  # from this, create a table containing estimates for common factors (squared and standardised)
  ACE_CI.common <- ACE_CI %>%
    filter(grepl('ACE.stac2|ACE.stcc2|ACE.stec2', parameter)) %>% # squared and standardised estimates 
    mutate(ACE = c(rep("A1",8), rep("A2",8), rep("C1",8), rep("C2",8), rep("E1",8), rep("E2",8))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8),6)) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_CI.specific <- ACE_CI %>%
    filter(grepl('ACE.stas2|ACE.stcs2|ACE.stes2', parameter) & estimate>0) %>% # delete the empty cells in the matrix
    mutate(ACE = c(rep("A",8), rep("C",8), rep("E",8))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8),3)) %>%
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # combine the common and specific path estimates
  pathestimates_table <- rbind(ACE_CI.common, ACE_CI.specific)
  
  pathestimates_table_rounded <- pathestimates_table %>%
    mutate(across(is.numeric, round, digits = 8))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages <- pathestimates_table %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           -Total)
  if(A_percent == TRUE){print(A_percentages)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages <- pathestimates_table %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           -Total)
  
  if(C_percent == TRUE){print(C_percentages)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           -Total)
  
  if(E_percent == TRUE){print(E_percentages)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent <- pathestimates_table %>%
    filter(grepl('A1|A2', ACE)) %>%
   filter(ACE != "A2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4)) %>% # remove common factor A2 for variable1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variable7"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variable8"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_A1 * variable5_A1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_A1 * variable6_A1)/Rph[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_A1 * variable7_A1)/Rph[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_A1 * variable8_A1)/Rph[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_A1 * variable5_A1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_A1 * variable6_A1)/Rph[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_A1 * variable7_A1)/Rph[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_A1 * variable8_A1)/Rph[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_A1 * variable4_A1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_A1 * variable5_A1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_A1 * variable6_A1)/Rph[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_A1 * variable7_A1)/Rph[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_A1 * variable8_A1)/Rph[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_A1 * variable5_A1)/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_A1 * variable6_A1)/Rph[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_A1 * variable7_A1)/Rph[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_A1 * variable8_A1)/Rph[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_A1 * variable6_A1) + sqrt(variable5_A2 * variable6_A2))/Rph[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_A1 * variable7_A1) + sqrt(variable5_A2 * variable7_A2))/Rph[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_A1 * variable8_A1) + sqrt(variable5_A2 * variable8_A2))/Rph[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_A1 * variable7_A1) + sqrt(variable6_A2 * variable7_A2))/Rph[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_A1 * variable8_A1) + sqrt(variable6_A2 * variable8_A2))/Rph[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_A1 * variable8_A1) + sqrt(variable7_A2 * variable8_A2))/Rph[8,7]) %>% # var 7 & 8
    select(variable1_variable2:variable7_variable8) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ac % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable7", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable8", variable8)))
    
  # Second, contribution of C to the phenotypic correlation
  Rc_contribution_percent <- pathestimates_table %>%
    filter(grepl('C1|C2', ACE)) %>%
    filter(ACE != "C2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4)) %>% # remove common factor C2 for variable 1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variable7"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variable8"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_C1 * variable2_C1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_C1 * variable3_C1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_C1 * variable4_C1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_C1 * variable5_C1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_C1 * variable6_C1)/Rph[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_C1 * variable7_C1)/Rph[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_C1 * variable8_C1)/Rph[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_C1 * variable3_C1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_C1 * variable4_C1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_C1 * variable5_C1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_C1 * variable6_C1)/Rph[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_C1 * variable7_C1)/Rph[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_C1 * variable8_C1)/Rph[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_C1 * variable4_C1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_C1 * variable5_C1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_C1 * variable6_C1)/Rph[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_C1 * variable7_C1)/Rph[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_C1 * variable8_C1)/Rph[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_C1 * variable5_C1)/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_C1 * variable6_C1)/Rph[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_C1 * variable7_C1)/Rph[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_C1 * variable8_C1)/Rph[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_C1 * variable6_C1) + sqrt(variable5_C2 * variable6_C2))/Rph[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_C1 * variable7_C1) + sqrt(variable5_C2 * variable7_C2))/Rph[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_C1 * variable8_C1) + sqrt(variable5_C2 * variable8_C2))/Rph[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_C1 * variable7_C1) + sqrt(variable6_C2 * variable7_C2))/Rph[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_C1 * variable8_C1) + sqrt(variable6_C2 * variable8_C2))/Rph[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_C1 * variable8_C1) + sqrt(variable7_C2 * variable8_C2))/Rph[8,7]) %>% # var 7 & 8
   select(variable1_variable2:variable7_variable8) %>%
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Cc % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable7", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable8", variable8))) 
  
  # Third, contribution of E to the phenotypic correlation
  Re_contribution_percent <- pathestimates_table %>%
    filter(grepl('E1|E2', ACE)) %>%
    filter(ACE != "E2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4)) %>% # remove common factor E2 for variable 1 and 2
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variable7"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variable8"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_E1 * variable2_E1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_E1 * variable3_E1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_E1 * variable4_E1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_E1 * variable5_E1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_E1 * variable6_E1)/Rph[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_E1 * variable7_E1)/Rph[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_E1 * variable8_E1)/Rph[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_E1 * variable3_E1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_E1 * variable4_E1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_E1 * variable5_E1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_E1 * variable6_E1)/Rph[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_E1 * variable7_E1)/Rph[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_E1 * variable8_E1)/Rph[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_E1 * variable4_E1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_E1 * variable5_E1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_E1 * variable6_E1)/Rph[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_E1 * variable7_E1)/Rph[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_E1 * variable8_E1)/Rph[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_E1 * variable5_E1)/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_E1 * variable6_E1)/Rph[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_E1 * variable7_E1)/Rph[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_E1 * variable8_E1)/Rph[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_E1 * variable6_E1) + sqrt(variable5_E2 * variable6_E2))/Rph[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_E1 * variable7_E1) + sqrt(variable5_E2 * variable7_E2))/Rph[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_E1 * variable8_E1) + sqrt(variable5_E2 * variable8_E2))/Rph[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_E1 * variable7_E1) + sqrt(variable6_E2 * variable7_E2))/Rph[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_E1 * variable8_E1) + sqrt(variable6_E2 * variable8_E2))/Rph[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_E1 * variable8_E1) + sqrt(variable7_E2 * variable8_E2))/Rph[8,7]) %>% # var 7 & 8
    select(variable1_variable2:variable7_variable8) %>%
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ec % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable7", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable8", variable8))) 
  
  # bind the contribution of A, C, and E together into one table
  contribution_percent_table <- plyr::join_all(list(Ra_contribution_percent, Rc_contribution_percent, Re_contribution_percent),
                                               by = "phenotypic correlation")
  
  contribution_percent_table <- contribution_percent_table %>% 
    mutate(total = `Ac % contribution` + `Cc % contribution` + `Ec % contribution`) %>% 
    mutate(across(is.numeric, round, digits = 5))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}

############################################
###            6 var 
###########################################

IPM_esitmates_ACE_6var <- function(data, variable1, variable2, variable3, variable4, variable5, variable6, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rph <- mxEval(ACE.Rph, data)
  
  # create table containing the univariate ACE estimates
  h2 <- as.data.frame(diag(round(mxEval(ACE.h2, data),3))) # heritability
  colnames(h2) <- c("h2")
  
  c2 <- as.data.frame(diag(round(mxEval(ACE.c2, data),3))) # shared environment 
  colnames(c2) <- c("c2")
  
  e2 <- as.data.frame(diag(round(mxEval(ACE.e2, data),3))) # unique environment
  colnames(e2) <- c("e2")
  
  uniACEestimates_table <- cbind(h2, c2, e2) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6)) %>%
    mutate(model = model) %>%
    select(model, variable, h2, c2, e2) %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table)} 
  
  # create data frame with all the path estimates and confidence intervals
  ACE_CI <- as.data.frame(summary$CI)
  ACE_CI <- ACE_CI %>%
    mutate(parameter = rownames(ACE_CI)) %>%
    as.tibble()
  
  # from this, create a table containing estimates for common factors (squared and standardised)
  ACE_CI.common <- ACE_CI %>%
    filter(grepl('ACE.stac2|ACE.stcc2|ACE.stec2', parameter)) %>% # squared and standardised estimates 
    mutate(ACE = c(rep("A1",6), rep("A2",6), rep("C1",6), rep("C2",6), rep("E1",6), rep("E2",6))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6),6)) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_CI.specific <- ACE_CI %>%
    filter(grepl('ACE.stas2|ACE.stcs2|ACE.stes2', parameter) & estimate>0) %>% # delete the empty cells in the matrix
    mutate(ACE = c(rep("A",6), rep("C",6), rep("E",6))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6),3)) %>%
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # combine the common and specific path estimates
  pathestimates_table <- rbind(ACE_CI.common, ACE_CI.specific)
  
  pathestimates_table_rounded <- pathestimates_table %>%
    mutate(across(is.numeric, round, digits = 8))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages <- pathestimates_table %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           -Total)
  if(A_percent == TRUE){print(A_percentages)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages <- pathestimates_table %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           -Total)
  
  if(C_percent == TRUE){print(C_percentages)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           -Total)
  
  if(E_percent == TRUE){print(E_percentages)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent <- pathestimates_table %>%
    filter(grepl('A1|A2', ACE)) %>%
    filter(ACE != "A2" | (variable != variable1 & variable != variable2 & variable != variable3)) %>% # remove common factor A2 for variable1, 2 and 3
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_A1 * variable5_A1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_A1 * variable6_A1)/Rph[6,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_A1 * variable5_A1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_A1 * variable6_A1)/Rph[6,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_A1 * variable4_A1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_A1 * variable5_A1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_A1 * variable6_A1)/Rph[6,3]) %>% 
    mutate(variable4_variable5 = (sqrt(variable4_A1 * variable5_A1) + sqrt(variable4_A2 * variable5_A2))/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = (sqrt(variable4_A1 * variable6_A1) + sqrt(variable4_A2 * variable6_A2))/Rph[6,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_A1 * variable6_A1) + sqrt(variable5_A2 * variable6_A2))/Rph[6,5]) %>% # var 5 & 6
    select(variable1_variable2:variable5_variable6) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ac % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) 
  
  # Second, contribution of C to the phenotypic correlation
  Rc_contribution_percent <- pathestimates_table %>%
    filter(grepl('C1|C2', ACE)) %>%
    filter(ACE != "C2" | (variable != variable1 & variable != variable2 & variable != variable3)) %>% # remove common factor A2 for variable1, 2 and 3
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_C1 * variable2_C1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_C1 * variable3_C1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_C1 * variable4_C1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_C1 * variable5_C1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_C1 * variable6_C1)/Rph[6,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_C1 * variable3_C1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_C1 * variable4_C1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_C1 * variable5_C1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_C1 * variable6_C1)/Rph[6,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_C1 * variable4_C1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_C1 * variable5_C1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_C1 * variable6_C1)/Rph[6,3]) %>% 
    mutate(variable4_variable5 = (sqrt(variable4_C1 * variable5_C1) + sqrt(variable4_C2 * variable5_C2))/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = (sqrt(variable4_C1 * variable6_C1) + sqrt(variable4_C2 * variable6_C2))/Rph[6,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_C1 * variable6_C1) + sqrt(variable5_C2 * variable6_C2))/Rph[6,5]) %>% # var 5 & 6
    select(variable1_variable2:variable5_variable6) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Cc % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) 
  
  # Third, contribution of E to the phenotypic correlation
  Re_contribution_percent <- pathestimates_table %>%
    filter(grepl('E1|E2', ACE)) %>%
    filter(ACE != "E2" | (variable != variable1 & variable != variable2 & variable != variable3)) %>% # remove common factor A2 for variable1, 2 and 3
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_E1 * variable2_E1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_E1 * variable3_E1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_E1 * variable4_E1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_E1 * variable5_E1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_E1 * variable6_E1)/Rph[6,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_E1 * variable3_E1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_E1 * variable4_E1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_E1 * variable5_E1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_E1 * variable6_E1)/Rph[6,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_E1 * variable4_E1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_E1 * variable5_E1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_E1 * variable6_E1)/Rph[6,3]) %>% 
    mutate(variable4_variable5 = (sqrt(variable4_E1 * variable5_E1) + sqrt(variable4_E2 * variable5_E2))/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = (sqrt(variable4_E1 * variable6_E1) + sqrt(variable4_E2 * variable6_E2))/Rph[6,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_E1 * variable6_E1) + sqrt(variable5_E2 * variable6_E2))/Rph[6,5]) %>% # var 5 & 6
    select(variable1_variable2:variable5_variable6) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ec % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) 
  
  # bind the contribution of A, C, and E together into one table
  contribution_percent_table <- plyr::join_all(list(Ra_contribution_percent, Rc_contribution_percent, Re_contribution_percent),
                                               by = "phenotypic correlation")
  
  contribution_percent_table <- contribution_percent_table %>% 
    mutate(total = `Ac % contribution` + `Cc % contribution` + `Ec % contribution`) %>% 
    mutate(across(is.numeric, round, digits = 5))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}

###########################################################
############             10 var
##########################################################

# data = IPM_allvar_fit
# variable1 = "Isolation12"
# variable2 = "Anxiety12"
# variable3 = "Depression12"
# variable4 = "Conduct12"
# variable5 = "Psychosis12"
# variable6 = "Isolation18"
# variable7 = "Anxiety18"
# variable8 = "Depression18"
# variable9 = "Conduct18"
# variable10 = "Psychosis18"
# model = "ACE"
# uniACEestimates = TRUE
# pathestimates = TRUE
# A_percent = TRUE
# C_percent = TRUE
# E_percent = TRUE
# common_factor_contribution = TRUE


IPM_esitmates_ACE_10var <- function(data, variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rph <- mxEval(ACE.Rph, data)
  
  # create table containing the univariate ACE estimates
  h2 <- as.data.frame(diag(round(mxEval(ACE.h2, data),3))) # heritability
  colnames(h2) <- c("h2")
  
  c2 <- as.data.frame(diag(round(mxEval(ACE.c2, data),3))) # shared environment 
  colnames(c2) <- c("c2")
  
  e2 <- as.data.frame(diag(round(mxEval(ACE.e2, data),3))) # unique environment
  colnames(e2) <- c("e2")
  
  uniACEestimates_table <- cbind(h2, c2, e2) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10)) %>%
    mutate(model = model) %>%
    select(model, variable, h2, c2, e2) %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table)} 
  
  # create data frame with all the path estimates and confidence intervals
  ACE_CI <- as.data.frame(summary$CI)
  ACE_CI <- ACE_CI %>%
    mutate(parameter = rownames(ACE_CI)) %>%
    as.tibble()
  
  # from this, create a table containing estimates for common factors (squared and standardised)
  ACE_CI.common <- ACE_CI %>%
    filter(grepl('ACE.stac2|ACE.stcc2|ACE.stec2', parameter)) %>% # squared and standardised estimates 
    mutate(ACE = c(rep("A1",10), rep("A2",10), rep("C1",10), rep("C2",10), rep("E1",10), rep("E2",10))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10),6)) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_CI.specific <- ACE_CI %>%
    filter(grepl('ACE.stas2|ACE.stcs2|ACE.stes2', parameter) & estimate>0) %>% # delete the empty cells in the matrix
    mutate(ACE = c(rep("A",10), rep("C",10), rep("E",10))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10),3)) %>%
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # combine the common and specific path estimates
  pathestimates_table <- rbind(ACE_CI.common, ACE_CI.specific)
  
  pathestimates_table_rounded <- pathestimates_table %>%
    mutate(across(is.numeric, round, digits = 8))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages <- pathestimates_table %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           -Total)
  if(A_percent == TRUE){print(A_percentages)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages <- pathestimates_table %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           -Total)
  
  if(C_percent == TRUE){print(C_percentages)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           -Total)
  
  if(E_percent == TRUE){print(E_percentages)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent <- pathestimates_table %>%
    filter(grepl('A1|A2', ACE)) %>%
    filter(ACE != "A2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4 & variable != variable5)) %>% # remove common factor A2 for variable1 to 5
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variable7"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variable8"))) %>% 
    mutate(across(variable, ~ str_replace(., variable9, "variable9"))) %>% 
    mutate(across(variable, ~ str_replace(., variable10, "variable10"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_A1 * variable5_A1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_A1 * variable6_A1)/Rph[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_A1 * variable7_A1)/Rph[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_A1 * variable8_A1)/Rph[8,1]) %>%
    mutate(variable1_variable9 = sqrt(variable1_A1 * variable9_A1)/Rph[9,1]) %>%
    mutate(variable1_variable10 = sqrt(variable1_A1 * variable10_A1)/Rph[10,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_A1 * variable5_A1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_A1 * variable6_A1)/Rph[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_A1 * variable7_A1)/Rph[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_A1 * variable8_A1)/Rph[8,2]) %>% 
    mutate(variable2_variable9 = sqrt(variable2_A1 * variable9_A1)/Rph[9,2]) %>% 
    mutate(variable2_variable10 = sqrt(variable2_A1 * variable10_A1)/Rph[10,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_A1 * variable4_A1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_A1 * variable5_A1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_A1 * variable6_A1)/Rph[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_A1 * variable7_A1)/Rph[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_A1 * variable8_A1)/Rph[8,3]) %>% 
    mutate(variable3_variable9 = sqrt(variable3_A1 * variable9_A1)/Rph[9,3]) %>% 
    mutate(variable3_variable10 = sqrt(variable3_A1 * variable10_A1)/Rph[10,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_A1 * variable5_A1)/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_A1 * variable6_A1)/Rph[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_A1 * variable7_A1)/Rph[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_A1 * variable8_A1)/Rph[8,4]) %>% 
    mutate(variable4_variable9 = sqrt(variable4_A1 * variable9_A1)/Rph[9,4]) %>% 
    mutate(variable4_variable10 = sqrt(variable4_A1 * variable10_A1)/Rph[10,4]) %>%  
    mutate(variable5_variable6 = sqrt(variable5_A1 * variable6_A1)/Rph[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = sqrt(variable5_A1 * variable7_A1)/Rph[7,5]) %>% 
    mutate(variable5_variable8 = sqrt(variable5_A1 * variable8_A1)/Rph[8,5]) %>% 
    mutate(variable5_variable9 = sqrt(variable5_A1 * variable9_A1)/Rph[9,5]) %>% 
    mutate(variable5_variable10 = sqrt(variable5_A1 * variable10_A1)/Rph[10,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_A1 * variable7_A1) + sqrt(variable6_A2 * variable7_A2))/Rph[7,6]) %>% # var 6 
    mutate(variable6_variable8 = (sqrt(variable6_A1 * variable8_A1) + sqrt(variable6_A2 * variable8_A2))/Rph[8,6]) %>% 
    mutate(variable6_variable9 = (sqrt(variable6_A1 * variable9_A1) + sqrt(variable6_A2 * variable9_A2))/Rph[9,6]) %>% 
    mutate(variable6_variable10 = (sqrt(variable6_A1 * variable10_A1) + sqrt(variable6_A2 * variable10_A2))/Rph[10,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_A1 * variable8_A1) + sqrt(variable7_A2 * variable8_A2))/Rph[8,7]) %>% # var 7
    mutate(variable7_variable9 = (sqrt(variable7_A1 * variable9_A1) + sqrt(variable7_A2 * variable9_A2))/Rph[9,7]) %>% 
    mutate(variable7_variable10 = (sqrt(variable7_A1 * variable10_A1) + sqrt(variable7_A2 * variable10_A2))/Rph[10,7]) %>% 
    mutate(variable8_variable9 = (sqrt(variable8_A1 * variable9_A1) + sqrt(variable8_A2 * variable9_A2))/Rph[9,8]) %>% # var 8
    mutate(variable8_variable10 = (sqrt(variable8_A1 * variable10_A1) + sqrt(variable8_A2 * variable10_A2))/Rph[10,8]) %>% 
    mutate(variable9_variable10 = (sqrt(variable9_A1 * variable10_A1) + sqrt(variable9_A2 * variable10_A2))/Rph[10,9]) %>% # var 9&10
    select(variable1_variable2:variable9_variable10) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ac % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable7", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable8", variable8))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable9", variable9))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable10", variable10)))
  
  # Second, contribution of C to the phenotypic correlation
  Rc_contribution_percent <- pathestimates_table %>%
    filter(grepl('C1|C2', ACE)) %>%
    filter(ACE != "C2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4 & variable != variable5)) %>% # remove common factor C2 for variable1 to 5
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variable7"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variable8"))) %>% 
    mutate(across(variable, ~ str_replace(., variable9, "variable9"))) %>% 
    mutate(across(variable, ~ str_replace(., variable10, "variable10"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_C1 * variable2_C1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_C1 * variable3_C1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_C1 * variable4_C1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_C1 * variable5_C1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_C1 * variable6_C1)/Rph[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_C1 * variable7_C1)/Rph[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_C1 * variable8_C1)/Rph[8,1]) %>%
    mutate(variable1_variable9 = sqrt(variable1_C1 * variable9_C1)/Rph[9,1]) %>%
    mutate(variable1_variable10 = sqrt(variable1_C1 * variable10_C1)/Rph[10,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_C1 * variable3_C1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_C1 * variable4_C1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_C1 * variable5_C1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_C1 * variable6_C1)/Rph[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_C1 * variable7_C1)/Rph[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_C1 * variable8_C1)/Rph[8,2]) %>% 
    mutate(variable2_variable9 = sqrt(variable2_C1 * variable9_C1)/Rph[9,2]) %>% 
    mutate(variable2_variable10 = sqrt(variable2_C1 * variable10_C1)/Rph[10,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_C1 * variable4_C1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_C1 * variable5_C1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_C1 * variable6_C1)/Rph[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_C1 * variable7_C1)/Rph[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_C1 * variable8_C1)/Rph[8,3]) %>% 
    mutate(variable3_variable9 = sqrt(variable3_C1 * variable9_C1)/Rph[9,3]) %>% 
    mutate(variable3_variable10 = sqrt(variable3_C1 * variable10_C1)/Rph[10,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_C1 * variable5_C1)/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_C1 * variable6_C1)/Rph[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_C1 * variable7_C1)/Rph[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_C1 * variable8_C1)/Rph[8,4]) %>% 
    mutate(variable4_variable9 = sqrt(variable4_C1 * variable9_C1)/Rph[9,4]) %>% 
    mutate(variable4_variable10 = sqrt(variable4_C1 * variable10_C1)/Rph[10,4]) %>%  
    mutate(variable5_variable6 = sqrt(variable5_C1 * variable6_C1)/Rph[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = sqrt(variable5_C1 * variable7_C1)/Rph[7,5]) %>% 
    mutate(variable5_variable8 = sqrt(variable5_C1 * variable8_C1)/Rph[8,5]) %>% 
    mutate(variable5_variable9 = sqrt(variable5_C1 * variable9_C1)/Rph[9,5]) %>% 
    mutate(variable5_variable10 = sqrt(variable5_C1 * variable10_C1)/Rph[10,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_C1 * variable7_C1) + sqrt(variable6_C2 * variable7_C2))/Rph[7,6]) %>% # var 6 
    mutate(variable6_variable8 = (sqrt(variable6_C1 * variable8_C1) + sqrt(variable6_C2 * variable8_C2))/Rph[8,6]) %>% 
    mutate(variable6_variable9 = (sqrt(variable6_C1 * variable9_C1) + sqrt(variable6_C2 * variable9_C2))/Rph[9,6]) %>% 
    mutate(variable6_variable10 = (sqrt(variable6_C1 * variable10_C1) + sqrt(variable6_C2 * variable10_C2))/Rph[10,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_C1 * variable8_C1) + sqrt(variable7_C2 * variable8_C2))/Rph[8,7]) %>% # var 7
    mutate(variable7_variable9 = (sqrt(variable7_C1 * variable9_C1) + sqrt(variable7_C2 * variable9_C2))/Rph[9,7]) %>% 
    mutate(variable7_variable10 = (sqrt(variable7_C1 * variable10_C1) + sqrt(variable7_C2 * variable10_C2))/Rph[10,7]) %>% 
    mutate(variable8_variable9 = (sqrt(variable8_C1 * variable9_C1) + sqrt(variable8_C2 * variable9_C2))/Rph[9,8]) %>% # var 8
    mutate(variable8_variable10 = (sqrt(variable8_C1 * variable10_C1) + sqrt(variable8_C2 * variable10_C2))/Rph[10,8]) %>% 
    mutate(variable9_variable10 = (sqrt(variable9_C1 * variable10_C1) + sqrt(variable9_C2 * variable10_C2))/Rph[10,9]) %>% # var 9&10
    select(variable1_variable2:variable9_variable10) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Cc % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable7", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable8", variable8))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable9", variable9))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable10", variable10)))
  
  # Third, contribution of E to the phenotypic correlation
  Re_contribution_percent <- pathestimates_table %>%
    filter(grepl('E1|E2', ACE)) %>%
    filter(ACE != "E2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4 & variable != variable5)) %>% # remove common factor E2 for variable1 to 5
    select(-ubound, -lbound, -factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variable1"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variable2"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variable3"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variable4"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variable5"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variable6"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variable7"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variable8"))) %>% 
    mutate(across(variable, ~ str_replace(., variable9, "variable9"))) %>% 
    mutate(across(variable, ~ str_replace(., variable10, "variable10"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    mutate(variable1_variable2 = sqrt(variable1_E1 * variable2_E1)/Rph[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_E1 * variable3_E1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_E1 * variable4_E1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_E1 * variable5_E1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_E1 * variable6_E1)/Rph[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_E1 * variable7_E1)/Rph[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_E1 * variable8_E1)/Rph[8,1]) %>%
    mutate(variable1_variable9 = sqrt(variable1_E1 * variable9_E1)/Rph[9,1]) %>%
    mutate(variable1_variable10 = sqrt(variable1_E1 * variable10_E1)/Rph[10,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_E1 * variable3_E1)/Rph[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_E1 * variable4_E1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_E1 * variable5_E1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_E1 * variable6_E1)/Rph[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_E1 * variable7_E1)/Rph[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_E1 * variable8_E1)/Rph[8,2]) %>% 
    mutate(variable2_variable9 = sqrt(variable2_E1 * variable9_E1)/Rph[9,2]) %>% 
    mutate(variable2_variable10 = sqrt(variable2_E1 * variable10_E1)/Rph[10,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_E1 * variable4_E1)/Rph[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_E1 * variable5_E1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_E1 * variable6_E1)/Rph[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_E1 * variable7_E1)/Rph[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_E1 * variable8_E1)/Rph[8,3]) %>% 
    mutate(variable3_variable9 = sqrt(variable3_E1 * variable9_E1)/Rph[9,3]) %>% 
    mutate(variable3_variable10 = sqrt(variable3_E1 * variable10_E1)/Rph[10,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_E1 * variable5_E1)/Rph[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_E1 * variable6_E1)/Rph[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_E1 * variable7_E1)/Rph[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_E1 * variable8_E1)/Rph[8,4]) %>% 
    mutate(variable4_variable9 = sqrt(variable4_E1 * variable9_E1)/Rph[9,4]) %>% 
    mutate(variable4_variable10 = sqrt(variable4_E1 * variable10_E1)/Rph[10,4]) %>%  
    mutate(variable5_variable6 = sqrt(variable5_E1 * variable6_E1)/Rph[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = sqrt(variable5_E1 * variable7_E1)/Rph[7,5]) %>% 
    mutate(variable5_variable8 = sqrt(variable5_E1 * variable8_E1)/Rph[8,5]) %>% 
    mutate(variable5_variable9 = sqrt(variable5_E1 * variable9_E1)/Rph[9,5]) %>% 
    mutate(variable5_variable10 = sqrt(variable5_E1 * variable10_E1)/Rph[10,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_E1 * variable7_E1) + sqrt(variable6_E2 * variable7_E2))/Rph[7,6]) %>% # var 6 
    mutate(variable6_variable8 = (sqrt(variable6_E1 * variable8_E1) + sqrt(variable6_E2 * variable8_E2))/Rph[8,6]) %>% 
    mutate(variable6_variable9 = (sqrt(variable6_E1 * variable9_E1) + sqrt(variable6_E2 * variable9_E2))/Rph[9,6]) %>% 
    mutate(variable6_variable10 = (sqrt(variable6_E1 * variable10_E1) + sqrt(variable6_E2 * variable10_E2))/Rph[10,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_E1 * variable8_E1) + sqrt(variable7_E2 * variable8_E2))/Rph[8,7]) %>% # var 7
    mutate(variable7_variable9 = (sqrt(variable7_E1 * variable9_E1) + sqrt(variable7_E2 * variable9_E2))/Rph[9,7]) %>% 
    mutate(variable7_variable10 = (sqrt(variable7_E1 * variable10_E1) + sqrt(variable7_E2 * variable10_E2))/Rph[10,7]) %>% 
    mutate(variable8_variable9 = (sqrt(variable8_E1 * variable9_E1) + sqrt(variable8_E2 * variable9_E2))/Rph[9,8]) %>% # var 8
    mutate(variable8_variable10 = (sqrt(variable8_E1 * variable10_E1) + sqrt(variable8_E2 * variable10_E2))/Rph[10,8]) %>% 
    mutate(variable9_variable10 = (sqrt(variable9_E1 * variable10_E1) + sqrt(variable9_E2 * variable10_E2))/Rph[10,9]) %>% # var 9&10
    select(variable1_variable2:variable9_variable10) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ec % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable1", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable2", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable3", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable4", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable5", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable6", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable7", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable8", variable8))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable9", variable9))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variable10", variable10)))
  
  # bind the contribution of A, C, and E together into one table
  contribution_percent_table <- plyr::join_all(list(Ra_contribution_percent, Rc_contribution_percent, Re_contribution_percent),
                                               by = "phenotypic correlation")
  
  contribution_percent_table <- contribution_percent_table %>% 
    mutate(total = `Ac % contribution` + `Cc % contribution` + `Ec % contribution`) %>% 
    mutate(across(is.numeric, round, digits = 5))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}
