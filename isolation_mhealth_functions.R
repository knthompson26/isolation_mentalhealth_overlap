# ---
#   title: "Functions for all scripts looking at the overlap between social isolation and mental health" 
# ---

# Load packages
library(knitr)
library(psych)
library(OpenMx)
library(tidyr)
library(tidyverse)
library(dplyr)     # conflicts with tidyverse for e.g. rename and row_number

###########################################################################################
###### Data prep 
###########################################################################################

# heat map
reorder_cor_matrix <- function(cor_matrix){  # reorder matrix
  dd <- as.dist((1-cor_matrix)/2)  # use correlation between variables as distance
  hc <- hclust(dd)
  cor_matrix <- cor_matrix[hc$order, hc$order]
}

get_upper_tri <- function(cor_matrix){  # Get upper triangle of the correlation matrix
  cor_matrix[lower.tri(cor_matrix)]<- NA
  return(cor_matrix)
}

###########################################################################################
###### Univariate 
###########################################################################################

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

ACE_esitmates_table_8varCI <- function(data, variables){
  table <- as.tibble(data$CI) %>%
    mutate(ACE = c(rep("h2", 8), rep("c2", 8), rep("e2", 8))) %>%
    mutate(Variable = rep(variables, 3)) %>%
    select(Variable, ACE, lbound, estimate, ubound) %>%
    mutate(across(is.numeric, round, digits = 3))
  
  return(table)
}

###########################################################################################
###### Sex differences
###########################################################################################

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
###### Bivariate 
###########################################################################################

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
###### Independent pathway model (IPM)
###########################################################################################

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


###########################################################################################
###### 8 variable IPM output
###########################################################################################

######## RENAME TO BE IPM_ace_estimates ????????????????????????

IPM_esitmates_ACE_8var <- function(data, variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rph <- mxEval(ACE.Rph, data)
  
  # create table containing the univariate ACE estimates
  h2 <- as.data.frame(diag(mxEval(ACE.h2, data))) # heritability
  colnames(h2) <- c("h2")
  
  c2 <- as.data.frame(diag(mxEval(ACE.c2, data))) # shared environment 
  colnames(c2) <- c("c2")
  
  e2 <- as.data.frame(diag(mxEval(ACE.e2, data))) # unique environment
  colnames(e2) <- c("e2")
  
  uniACEestimates_table <- cbind(h2, c2, e2) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(model = model) %>%
    select(model, variable, h2, c2, e2) 
  
  uniACEestimates_table_rounded <- uniACEestimates_table %>%
    mutate(across(is.numeric, round, digits = 5))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table_rounded)} 
  
  # create data frame with all the path estimates and confidence intervals
  ACE_CI <- as.data.frame(summary$CI)
  ACE_CI <- ACE_CI %>%
    mutate(parameter = rownames(ACE_CI)) %>%
    as_tibble()
  
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
    filter(parameter == "ACE.stas2[1,1]" |   # specific A - this is a long way to do this but couldn't get grepl to work on [1,1] etc
             parameter == "ACE.stas2[2,2]" | 
             parameter == "ACE.stas2[3,3]" | 
             parameter == "ACE.stas2[4,4]" | 
             parameter == "ACE.stas2[5,5]" | 
             parameter == "ACE.stas2[6,6]" | 
             parameter == "ACE.stas2[7,7]" | 
             parameter == "ACE.stas2[8,8]" |
             parameter == "ACE.stcs2[1,1]" |   # specific C
             parameter == "ACE.stcs2[2,2]" | 
             parameter == "ACE.stcs2[3,3]" | 
             parameter == "ACE.stcs2[4,4]" | 
             parameter == "ACE.stcs2[5,5]" | 
             parameter == "ACE.stcs2[6,6]" | 
             parameter == "ACE.stcs2[7,7]" | 
             parameter == "ACE.stcs2[8,8]" |
             parameter == "ACE.stes2[1,1]" |   # specific E
             parameter == "ACE.stes2[2,2]" | 
             parameter == "ACE.stes2[3,3]" | 
             parameter == "ACE.stes2[4,4]" | 
             parameter == "ACE.stes2[5,5]" | 
             parameter == "ACE.stes2[6,6]" | 
             parameter == "ACE.stes2[7,7]" | 
             parameter == "ACE.stes2[8,8]") %>%
    mutate(ACE = c(rep("A",8), rep("C",8), rep("E",8))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8),3)) %>%
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
    mutate(total_var = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(A1 + A2 + A, 3)) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           Total)
  
  h2_joined <- uniACEestimates_table %>%
    select(variable, h2) %>%
    mutate(across(h2, ~ round(.x*100, 3)))
  
  A_percentages_joined <- plyr::join_all(list(A_percentages, h2_joined),
                                         by = "variable") %>%
    select(variable, h2, `A1c %`, `A2c %`, `As %`, Total)
  
  if(A_percent == TRUE){print(A_percentages_joined)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages <- pathestimates_table %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(C1 + C2 + C, 3)) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           Total)
  
  c2_joined <- uniACEestimates_table %>%
    select(variable, c2) %>%
    mutate(across(c2, ~ round(.x*100, 3)))
  
  C_percentages_joined <- plyr::join_all(list(C_percentages, c2_joined),
                                         by = "variable") %>%
    select(variable, c2, `C1c %`, `C2c %`, `Cs %`, Total)
  
  if(C_percent == TRUE){print(C_percentages_joined)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(E1 + E2 + E, 3)) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           Total)
  
  e2_joined <- uniACEestimates_table %>%
    select(variable, e2) %>%
    mutate(across(e2, ~ round(.x*100, 3)))
  
  E_percentages_joined <- plyr::join_all(list(E_percentages, e2_joined),
                                         by = "variable") %>%
    select(variable, e2, `E1c %`, `E2c %`, `Es %`, Total)
  
  if(E_percent == TRUE){print(E_percentages_joined)}
  
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
    mutate(across(is.numeric, ~ round(.x*100, 3)))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}


###########################################################################################
###### 8 variable IPM output - NO CONFIDENCE INTERVALS
###########################################################################################

# You may want to run your scripts without calculating confidence intervals (intervals = FALSE) to save time. This output function will display the output when your model mxTryHard(model, intervals = FALSE). 

IPM_esitmates_ACE_8var_noCI <- function(data, variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rph <- mxEval(ACE.Rph, data)
  
  # create table containing the univariate ACE estimates
  h2 <- as.data.frame(diag(mxEval(ACE.h2, data))) # heritability
  colnames(h2) <- c("h2")
  
  c2 <- as.data.frame(diag(mxEval(ACE.c2, data))) # shared environment 
  colnames(c2) <- c("c2")
  
  e2 <- as.data.frame(diag(mxEval(ACE.e2, data))) # unique environment
  colnames(e2) <- c("e2")
  
  uniACEestimates_table <- cbind(h2, c2, e2) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(model = model) %>%
    select(model, variable, h2, c2, e2) 
  
  uniACEestimates_table_rounded <- uniACEestimates_table %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table_rounded)} 
  
  # create a table containing estimates for common factors (squared and standardised)
  ACE_stac2 <- as.data.frame(mxEval(ACE.stac2, data)) %>% 
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    select(variable, A1 = V1, A2 = V2) %>%
    pivot_longer(cols = starts_with("A"), names_to = "ACE", values_to = "estimate")
  
  ACE_stcc2 <- as.data.frame(mxEval(ACE.stcc2, data)) %>% 
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    select(variable, C1 = V1, C2 = V2) %>%
    pivot_longer(cols = starts_with("C"), names_to = "ACE", values_to = "estimate")
  
  ACE_stec2 <- as.data.frame(mxEval(ACE.stec2, data)) %>% 
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    select(variable, E1 = V1, E2 = V2) %>%
    pivot_longer(cols = starts_with("E"), names_to = "ACE", values_to = "estimate")
  
  ACE.common <- rbind(ACE_stac2, ACE_stcc2, ACE_stec2)
  
  ACE.common <- ACE.common %>% 
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, estimate)
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_stas2 <- as.data.frame(diag(mxEval(ACE.stas2, data))) 
  colnames(ACE_stas2) <- c("estimate")
  ACE_stas2 <- ACE_stas2 %>%
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(ACE = "A") 
  
  ACE_stcs2 <- as.data.frame(diag(mxEval(ACE.stcs2, data))) 
  colnames(ACE_stcs2) <- c("estimate")
  ACE_stcs2 <- ACE_stcs2 %>%
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(ACE = "C") 
  
  ACE_stes2 <- as.data.frame(diag(mxEval(ACE.stes2, data))) 
  colnames(ACE_stes2) <- c("estimate")
  ACE_stes2 <- ACE_stes2 %>%
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(ACE = "E")
  
  ACE.specific <- rbind(ACE_stas2, ACE_stcs2, ACE_stes2)
  
  ACE.specific <- ACE.specific %>% 
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, estimate)
  
  # combine the common and specific path estimates
  pathestimates_table <- rbind(ACE.common, ACE.specific)
  
  pathestimates_table_rounded <- pathestimates_table %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages <- pathestimates_table %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(A1 + A2 + A, 3)) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           Total)
  
  h2_joined <- uniACEestimates_table %>%
    select(variable, h2) %>%
    mutate(across(h2, ~ round(.x*100, 3)))
    
  A_percentages_joined <- plyr::join_all(list(A_percentages, h2_joined),
                                         by = "variable") %>%
    select(variable, h2, `A1c %`, `A2c %`, `As %`, Total)
  
  if(A_percent == TRUE){print(A_percentages_joined)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages <- pathestimates_table %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(C1 + C2 + C, 3)) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           Total)
  
  c2_joined <- uniACEestimates_table %>%
    select(variable, c2) %>%
    mutate(across(c2, ~ round(.x*100, 3)))
  
  C_percentages_joined <- plyr::join_all(list(C_percentages, c2_joined),
                                         by = "variable") %>%
    select(variable, c2, `C1c %`, `C2c %`, `Cs %`, Total)
  
  if(C_percent == TRUE){print(C_percentages_joined)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(E1 + E2 + E, 3)) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           Total)
  
  e2_joined <- uniACEestimates_table %>%
    select(variable, e2) %>%
    mutate(across(e2, ~ round(.x*100, 3)))
  
  E_percentages_joined <- plyr::join_all(list(E_percentages, e2_joined),
                                         by = "variable") %>%
    select(variable, e2, `E1c %`, `E2c %`, `Es %`, Total)
  
  if(E_percent == TRUE){print(E_percentages_joined)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent <- pathestimates_table %>%
    filter(grepl('A1|A2', ACE)) %>%
    filter(ACE != "A2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4)) %>% # remove common factor A2 for variable1 and 2
    select(-factor) %>%
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
    select(-factor) %>%
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
    select(-factor) %>%
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
    mutate(across(is.numeric, ~ round(.x*100, 3)))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}

















###########################################################################################
###### 8 variable IPM output - het for conduct included
###########################################################################################

IPM_esitmates_ACE_8var_het <- function(data, variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
  # create summary object of the openmx model
  summary <- summary(data, verbose = TRUE)
  
  # create matrix containing the model estimated phenotypic correlations
  Rphm <- mxEval(ACE.Rphm, data)
  Rphf <- mxEval(ACE.Rphf, data)
  
  # create table containing the univariate ACE estimates
  h2m <- as.data.frame(diag(mxEval(ACE.h2m, data))) # heritability for males
  colnames(h2m) <- c("h2m")
  h2f <- as.data.frame(diag(mxEval(ACE.h2f, data))) # heritability for females
  colnames(h2f) <- c("h2f")
  
  c2m <- as.data.frame(diag(mxEval(ACE.c2m, data))) # shared environment 
  colnames(c2m) <- c("c2m")
  c2f <- as.data.frame(diag(mxEval(ACE.c2f, data))) # shared environment 
  colnames(c2f) <- c("c2f")
  
  e2m <- as.data.frame(diag(mxEval(ACE.e2m, data))) # unique environment
  colnames(e2m) <- c("e2m")
  e2f <- as.data.frame(diag(mxEval(ACE.e2f, data))) # unique environment
  colnames(e2f) <- c("e2f")
  
  uniACEestimates_table_m <- cbind(h2m, c2m, e2m) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(model = model) %>%
    mutate(sex = "Male") %>%
    select(model, variable, sex, h2 = h2m, c2 = c2m, e2 = e2m) 
  
  uniACEestimates_table_f <- cbind(h2f, c2f, e2f) %>%           # bind all together into a table
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8)) %>%
    mutate(model = model) %>%
    mutate(sex = "Female") %>%
    select(model, variable, sex, h2 = h2f, c2 = c2f, e2 = e2f) 
  
  uniACEestimates_table <- rbind(uniACEestimates_table_m, uniACEestimates_table_f)
  
  uniACEestimates_table_rounded <- uniACEestimates_table %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(uniACEestimates == TRUE){print(uniACEestimates_table_rounded)} 
  
  # create data frame with all the path estimates and confidence intervals
  ACE_CI <- as.data.frame(summary$CI)
  ACE_CI <- ACE_CI %>%
    mutate(parameter = rownames(ACE_CI)) %>%
    as_tibble()
  
  # from this, create a table containing estimates for common factors (squared and standardised)
  ACE_CI.common_m <- ACE_CI %>%
    filter(grepl('ACE.stac2m|ACE.stcc2m|ACE.stec2m', parameter)) %>% # squared and standardised estimates 
    mutate(ACE = c(rep("A1",8), rep("A2",8), rep("C1",8), rep("C2",8), rep("E1",8), rep("E2",8))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8),6)) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    mutate(sex = "Male") %>%
    select(model, sex, factor, ACE, variable, lbound, estimate, ubound) 
  
  ACE_CI.common_f <- ACE_CI %>%
    filter(grepl('ACE.stac2f|ACE.stcc2f|ACE.stec2f', parameter)) %>% # squared and standardised estimates 
    mutate(ACE = c(rep("A1",8), rep("A2",8), rep("C1",8), rep("C2",8), rep("E1",8), rep("E2",8))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8),6)) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    mutate(sex = "Female") %>%
    select(model, sex, factor, ACE, variable, lbound, estimate, ubound) 
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_CI.specific_m <- ACE_CI %>%
    filter(parameter == "ACE.stas2m[1,1]" |   # specific A - this is a long way to do this but couldn't get grepl to work on [1,1] etc
             parameter == "ACE.stas2m[2,2]" | 
             parameter == "ACE.stas2m[3,3]" | 
             parameter == "ACE.stas2m[4,4]" | 
             parameter == "ACE.stas2m[5,5]" | 
             parameter == "ACE.stas2m[6,6]" | 
             parameter == "ACE.stas2m[7,7]" | 
             parameter == "ACE.stas2m[8,8]" |
             parameter == "ACE.stcs2m[1,1]" |   # specific C
             parameter == "ACE.stcs2m[2,2]" | 
             parameter == "ACE.stcs2m[3,3]" | 
             parameter == "ACE.stcs2m[4,4]" | 
             parameter == "ACE.stcs2m[5,5]" | 
             parameter == "ACE.stcs2m[6,6]" | 
             parameter == "ACE.stcs2m[7,7]" | 
             parameter == "ACE.stcs2m[8,8]" |
             parameter == "ACE.stes2m[1,1]" |   # specific E
             parameter == "ACE.stes2m[2,2]" | 
             parameter == "ACE.stes2m[3,3]" | 
             parameter == "ACE.stes2m[4,4]" | 
             parameter == "ACE.stes2m[5,5]" | 
             parameter == "ACE.stes2m[6,6]" | 
             parameter == "ACE.stes2m[7,7]" | 
             parameter == "ACE.stes2m[8,8]") %>%
    mutate(ACE = c(rep("A",8), rep("C",8), rep("E",8))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8),3)) %>%
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    mutate(sex = "Male") %>%
    select(model, sex, factor, ACE, variable, lbound, estimate, ubound) 
  
  ACE_CI.specific_f <- ACE_CI %>%
    filter(parameter == "ACE.stas2f[1,1]" |   # specific A - this is a long way to do this but couldn't get grepl to work on [1,1] etc
             parameter == "ACE.stas2f[2,2]" | 
             parameter == "ACE.stas2f[3,3]" | 
             parameter == "ACE.stas2f[4,4]" | 
             parameter == "ACE.stas2f[5,5]" | 
             parameter == "ACE.stas2f[6,6]" | 
             parameter == "ACE.stas2f[7,7]" | 
             parameter == "ACE.stas2f[8,8]" |
             parameter == "ACE.stcs2f[1,1]" |   # specific C
             parameter == "ACE.stcs2f[2,2]" | 
             parameter == "ACE.stcs2f[3,3]" | 
             parameter == "ACE.stcs2f[4,4]" | 
             parameter == "ACE.stcs2f[5,5]" | 
             parameter == "ACE.stcs2f[6,6]" | 
             parameter == "ACE.stcs2f[7,7]" | 
             parameter == "ACE.stcs2f[8,8]" |
             parameter == "ACE.stes2f[1,1]" |   # specific E
             parameter == "ACE.stes2f[2,2]" | 
             parameter == "ACE.stes2f[3,3]" | 
             parameter == "ACE.stes2f[4,4]" | 
             parameter == "ACE.stes2f[5,5]" | 
             parameter == "ACE.stes2f[6,6]" | 
             parameter == "ACE.stes2f[7,7]" | 
             parameter == "ACE.stes2f[8,8]") %>%
    mutate(ACE = c(rep("A",8), rep("C",8), rep("E",8))) %>%
    mutate(variable = rep(c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8),3)) %>%
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    mutate(sex = "Female") %>%
    select(model, sex, factor, ACE, variable, lbound, estimate, ubound) 
  
  # combine the common and specific path estimates
  pathestimates_table_m <- rbind(ACE_CI.common_m, ACE_CI.specific_m)
  pathestimates_table_f <- rbind(ACE_CI.common_f, ACE_CI.specific_f)
  
  pathestimates_table_rounded_m <- pathestimates_table_m %>%
    mutate(across(is.numeric, round, digits = 3))
  pathestimates_table_rounded_f <- pathestimates_table_f %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded_m)} 
  if(pathestimates == TRUE){print(pathestimates_table_rounded_f)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages_m <- pathestimates_table_m %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(A1 + A2 + A, 3)) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           Total)
  
  h2_joined_m <- uniACEestimates_table %>%
    filter(sex == "Male") %>%
    select(variable, h2) %>%
    mutate(across(h2, ~ round(.x*100, 3)))
  
  A_percentages_joined_m <- plyr::join_all(list(A_percentages_m, h2_joined_m),
                                         by = "variable") %>%
    mutate(sex = "Male") %>%
    select(sex, variable, h2, `A1c %`, `A2c %`, `As %`, Total)
  
  A_percentages_f <- pathestimates_table_f %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(A1 + A2 + A, 3)) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           Total)
  
  h2_joined_f <- uniACEestimates_table %>%
    filter(sex == "Female") %>%
    select(variable, h2) %>%
    mutate(across(h2, ~ round(.x*100, 3)))
  
  A_percentages_joined_f <- plyr::join_all(list(A_percentages_f, h2_joined_f),
                                           by = "variable") %>%
    mutate(sex = "Female") %>%
    select(sex, variable, h2, `A1c %`, `A2c %`, `As %`, Total)
  
  A_percentages_joined <- rbind(A_percentages_joined_m, A_percentages_joined_f)
  
  if(A_percent == TRUE){print(A_percentages_joined)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages_m <- pathestimates_table_m %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(C1 + C2 + C, 3)) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           Total)
  
  c2_joined_m <- uniACEestimates_table %>%
    filter(sex == "Male") %>%
    select(variable, c2) %>%
    mutate(across(c2, ~ round(.x*100, 3)))
  
  C_percentages_joined_m <- plyr::join_all(list(C_percentages_m, c2_joined_m),
                                         by = "variable") %>%
    mutate(sex = "Male") %>%
    select(sex, variable, c2, `C1c %`, `C2c %`, `Cs %`, Total)
  
  C_percentages_f <- pathestimates_table_f %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(C1 + C2 + C, 3)) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           Total)
  
  c2_joined_f <- uniACEestimates_table %>%
    filter(sex == "Female") %>%
    select(variable, c2) %>%
    mutate(across(c2, ~ round(.x*100, 3)))
  
  C_percentages_joined_f <- plyr::join_all(list(C_percentages_f, c2_joined_f),
                                           by = "variable") %>%
    mutate(sex = "Female") %>%
    select(sex, variable, c2, `C1c %`, `C2c %`, `Cs %`, Total)
  
  C_percentages_joined <- rbind(C_percentages_joined_m, C_percentages_joined_f)
  
  if(C_percent == TRUE){print(C_percentages_joined)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages_m <- pathestimates_table_m %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(E1 + E2 + E, 3)) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           Total)
  
  e2_joined_m <- uniACEestimates_table %>%
    filter(sex == "Male") %>%
    select(variable, e2) %>%
    mutate(across(e2, ~ round(.x*100, 3)))
  
  E_percentages_joined_m <- plyr::join_all(list(E_percentages_m, e2_joined_m),
                                         by = "variable") %>%
    mutate(sex = "Male") %>%
    select(sex, variable, e2, `E1c %`, `E2c %`, `Es %`, Total)
  
  E_percentages_f <- pathestimates_table_f %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-ubound, -lbound, -factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(total_var = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/total_var*100, 3))) %>%
    mutate(Total = round(E1 + E2 + E, 3)) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           Total)
  
  e2_joined_f <- uniACEestimates_table %>%
    filter(sex == "Female") %>%
    select(variable, e2) %>%
    mutate(across(e2, ~ round(.x*100, 3)))
  
  E_percentages_joined_f <- plyr::join_all(list(E_percentages_f, e2_joined_f),
                                         by = "variable") %>%
    mutate(sex = "Female") %>%
    select(sex, variable, e2, `E1c %`, `E2c %`, `Es %`, Total)
  
  E_percentages_joined <- rbind(E_percentages_joined_m, E_percentages_joined_f)
  
  if(E_percent == TRUE){print(E_percentages_joined)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent_m <- pathestimates_table_m %>%
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
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rphm[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rphm[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rphm[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_A1 * variable5_A1)/Rphm[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_A1 * variable6_A1)/Rphm[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_A1 * variable7_A1)/Rphm[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_A1 * variable8_A1)/Rphm[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rphm[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rphm[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_A1 * variable5_A1)/Rphm[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_A1 * variable6_A1)/Rphm[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_A1 * variable7_A1)/Rphm[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_A1 * variable8_A1)/Rphm[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_A1 * variable4_A1)/Rphm[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_A1 * variable5_A1)/Rphm[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_A1 * variable6_A1)/Rphm[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_A1 * variable7_A1)/Rphm[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_A1 * variable8_A1)/Rphm[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_A1 * variable5_A1)/Rphm[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_A1 * variable6_A1)/Rphm[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_A1 * variable7_A1)/Rphm[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_A1 * variable8_A1)/Rphm[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_A1 * variable6_A1) + sqrt(variable5_A2 * variable6_A2))/Rphm[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_A1 * variable7_A1) + sqrt(variable5_A2 * variable7_A2))/Rphm[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_A1 * variable8_A1) + sqrt(variable5_A2 * variable8_A2))/Rphm[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_A1 * variable7_A1) + sqrt(variable6_A2 * variable7_A2))/Rphm[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_A1 * variable8_A1) + sqrt(variable6_A2 * variable8_A2))/Rphm[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_A1 * variable8_A1) + sqrt(variable7_A2 * variable8_A2))/Rphm[8,7]) %>% # var 7 & 8
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
  Rc_contribution_percent_m <- pathestimates_table_m %>%
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
    mutate(variable1_variable2 = sqrt(variable1_C1 * variable2_C1)/Rphm[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_C1 * variable3_C1)/Rphm[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_C1 * variable4_C1)/Rphm[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_C1 * variable5_C1)/Rphm[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_C1 * variable6_C1)/Rphm[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_C1 * variable7_C1)/Rphm[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_C1 * variable8_C1)/Rphm[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_C1 * variable3_C1)/Rphm[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_C1 * variable4_C1)/Rphm[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_C1 * variable5_C1)/Rphm[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_C1 * variable6_C1)/Rphm[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_C1 * variable7_C1)/Rphm[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_C1 * variable8_C1)/Rphm[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_C1 * variable4_C1)/Rphm[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_C1 * variable5_C1)/Rphm[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_C1 * variable6_C1)/Rphm[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_C1 * variable7_C1)/Rphm[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_C1 * variable8_C1)/Rphm[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_C1 * variable5_C1)/Rphm[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_C1 * variable6_C1)/Rphm[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_C1 * variable7_C1)/Rphm[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_C1 * variable8_C1)/Rphm[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_C1 * variable6_C1) + sqrt(variable5_C2 * variable6_C2))/Rphm[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_C1 * variable7_C1) + sqrt(variable5_C2 * variable7_C2))/Rphm[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_C1 * variable8_C1) + sqrt(variable5_C2 * variable8_C2))/Rphm[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_C1 * variable7_C1) + sqrt(variable6_C2 * variable7_C2))/Rphm[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_C1 * variable8_C1) + sqrt(variable6_C2 * variable8_C2))/Rphm[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_C1 * variable8_C1) + sqrt(variable7_C2 * variable8_C2))/Rphm[8,7]) %>% # var 7 & 8
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
  Re_contribution_percent_m <- pathestimates_table_m %>%
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
    mutate(variable1_variable2 = sqrt(variable1_E1 * variable2_E1)/Rphm[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_E1 * variable3_E1)/Rphm[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_E1 * variable4_E1)/Rphm[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_E1 * variable5_E1)/Rphm[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_E1 * variable6_E1)/Rphm[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_E1 * variable7_E1)/Rphm[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_E1 * variable8_E1)/Rphm[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_E1 * variable3_E1)/Rphm[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_E1 * variable4_E1)/Rphm[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_E1 * variable5_E1)/Rphm[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_E1 * variable6_E1)/Rphm[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_E1 * variable7_E1)/Rphm[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_E1 * variable8_E1)/Rphm[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_E1 * variable4_E1)/Rphm[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_E1 * variable5_E1)/Rphm[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_E1 * variable6_E1)/Rphm[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_E1 * variable7_E1)/Rphm[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_E1 * variable8_E1)/Rphm[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_E1 * variable5_E1)/Rphm[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_E1 * variable6_E1)/Rphm[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_E1 * variable7_E1)/Rphm[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_E1 * variable8_E1)/Rphm[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_E1 * variable6_E1) + sqrt(variable5_E2 * variable6_E2))/Rphm[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_E1 * variable7_E1) + sqrt(variable5_E2 * variable7_E2))/Rphm[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_E1 * variable8_E1) + sqrt(variable5_E2 * variable8_E2))/Rphm[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_E1 * variable7_E1) + sqrt(variable6_E2 * variable7_E2))/Rphm[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_E1 * variable8_E1) + sqrt(variable6_E2 * variable8_E2))/Rphm[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_E1 * variable8_E1) + sqrt(variable7_E2 * variable8_E2))/Rphm[8,7]) %>% # var 7 & 8
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
  contribution_percent_table_m <- plyr::join_all(list(Ra_contribution_percent_m, Rc_contribution_percent_m, Re_contribution_percent_m),
                                               by = "phenotypic correlation")
  
  contribution_percent_table_m <- contribution_percent_table_m %>% 
    mutate(total = `Ac % contribution` + `Cc % contribution` + `Ec % contribution`) %>% 
    mutate(sex = "Male") %>%
    mutate(across(is.numeric, ~ round(.x*100, 3)))
  
  # FEMALE
  Ra_contribution_percent_f <- pathestimates_table_f %>%
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
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rphf[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rphf[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rphf[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_A1 * variable5_A1)/Rphf[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_A1 * variable6_A1)/Rphf[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_A1 * variable7_A1)/Rphf[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_A1 * variable8_A1)/Rphf[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rphf[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rphf[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_A1 * variable5_A1)/Rphf[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_A1 * variable6_A1)/Rphf[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_A1 * variable7_A1)/Rphf[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_A1 * variable8_A1)/Rphf[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_A1 * variable4_A1)/Rphf[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_A1 * variable5_A1)/Rphf[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_A1 * variable6_A1)/Rphf[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_A1 * variable7_A1)/Rphf[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_A1 * variable8_A1)/Rphf[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_A1 * variable5_A1)/Rphf[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_A1 * variable6_A1)/Rphf[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_A1 * variable7_A1)/Rphf[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_A1 * variable8_A1)/Rphf[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_A1 * variable6_A1) + sqrt(variable5_A2 * variable6_A2))/Rphf[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_A1 * variable7_A1) + sqrt(variable5_A2 * variable7_A2))/Rphf[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_A1 * variable8_A1) + sqrt(variable5_A2 * variable8_A2))/Rphf[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_A1 * variable7_A1) + sqrt(variable6_A2 * variable7_A2))/Rphf[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_A1 * variable8_A1) + sqrt(variable6_A2 * variable8_A2))/Rphf[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_A1 * variable8_A1) + sqrt(variable7_A2 * variable8_A2))/Rphf[8,7]) %>% # var 7 & 8
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
  Rc_contribution_percent_f <- pathestimates_table_f %>%
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
    mutate(variable1_variable2 = sqrt(variable1_C1 * variable2_C1)/Rphf[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_C1 * variable3_C1)/Rphf[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_C1 * variable4_C1)/Rphf[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_C1 * variable5_C1)/Rphf[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_C1 * variable6_C1)/Rphf[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_C1 * variable7_C1)/Rphf[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_C1 * variable8_C1)/Rphf[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_C1 * variable3_C1)/Rphf[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_C1 * variable4_C1)/Rphf[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_C1 * variable5_C1)/Rphf[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_C1 * variable6_C1)/Rphf[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_C1 * variable7_C1)/Rphf[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_C1 * variable8_C1)/Rphf[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_C1 * variable4_C1)/Rphf[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_C1 * variable5_C1)/Rphf[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_C1 * variable6_C1)/Rphf[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_C1 * variable7_C1)/Rphf[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_C1 * variable8_C1)/Rphf[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_C1 * variable5_C1)/Rphf[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_C1 * variable6_C1)/Rphf[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_C1 * variable7_C1)/Rphf[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_C1 * variable8_C1)/Rphf[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_C1 * variable6_C1) + sqrt(variable5_C2 * variable6_C2))/Rphf[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_C1 * variable7_C1) + sqrt(variable5_C2 * variable7_C2))/Rphf[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_C1 * variable8_C1) + sqrt(variable5_C2 * variable8_C2))/Rphf[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_C1 * variable7_C1) + sqrt(variable6_C2 * variable7_C2))/Rphf[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_C1 * variable8_C1) + sqrt(variable6_C2 * variable8_C2))/Rphf[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_C1 * variable8_C1) + sqrt(variable7_C2 * variable8_C2))/Rphf[8,7]) %>% # var 7 & 8
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
  Re_contribution_percent_f <- pathestimates_table_f %>%
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
    mutate(variable1_variable2 = sqrt(variable1_E1 * variable2_E1)/Rphf[2,1]) %>% #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable3 = sqrt(variable1_E1 * variable3_E1)/Rphf[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_E1 * variable4_E1)/Rphf[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_E1 * variable5_E1)/Rphf[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_E1 * variable6_E1)/Rphf[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_E1 * variable7_E1)/Rphf[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_E1 * variable8_E1)/Rphf[8,1]) %>%
    mutate(variable2_variable3 = sqrt(variable2_E1 * variable3_E1)/Rphf[3,2]) %>% # var 2 associated with all others
    mutate(variable2_variable4 = sqrt(variable2_E1 * variable4_E1)/Rphf[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_E1 * variable5_E1)/Rphf[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_E1 * variable6_E1)/Rphf[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_E1 * variable7_E1)/Rphf[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_E1 * variable8_E1)/Rphf[8,2]) %>% 
    mutate(variable3_variable4 = sqrt(variable3_E1 * variable4_E1)/Rphf[4,3]) %>% # var 3 associated with all others
    mutate(variable3_variable5 = sqrt(variable3_E1 * variable5_E1)/Rphf[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_E1 * variable6_E1)/Rphf[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_E1 * variable7_E1)/Rphf[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_E1 * variable8_E1)/Rphf[8,3]) %>% 
    mutate(variable4_variable5 = sqrt(variable4_E1 * variable5_E1)/Rphf[5,4]) %>% # var 4 associated with all others
    mutate(variable4_variable6 = sqrt(variable4_E1 * variable6_E1)/Rphf[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_E1 * variable7_E1)/Rphf[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_E1 * variable8_E1)/Rphf[8,4]) %>% 
    mutate(variable5_variable6 = (sqrt(variable5_E1 * variable6_E1) + sqrt(variable5_E2 * variable6_E2))/Rphf[6,5]) %>% # var 5 associated with all others
    mutate(variable5_variable7 = (sqrt(variable5_E1 * variable7_E1) + sqrt(variable5_E2 * variable7_E2))/Rphf[7,5]) %>% 
    mutate(variable5_variable8 = (sqrt(variable5_E1 * variable8_E1) + sqrt(variable5_E2 * variable8_E2))/Rphf[8,5]) %>% 
    mutate(variable6_variable7 = (sqrt(variable6_E1 * variable7_E1) + sqrt(variable6_E2 * variable7_E2))/Rphf[7,6]) %>% # var 6 associated with all others 
    mutate(variable6_variable8 = (sqrt(variable6_E1 * variable8_E1) + sqrt(variable6_E2 * variable8_E2))/Rphf[8,6]) %>% 
    mutate(variable7_variable8 = (sqrt(variable7_E1 * variable8_E1) + sqrt(variable7_E2 * variable8_E2))/Rphf[8,7]) %>% # var 7 & 8
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
  contribution_percent_table_f <- plyr::join_all(list(Ra_contribution_percent_f, Rc_contribution_percent_f, Re_contribution_percent_f),
                                                 by = "phenotypic correlation")
  
  contribution_percent_table_f <- contribution_percent_table_f %>% 
    mutate(total = `Ac % contribution` + `Cc % contribution` + `Ec % contribution`) %>% 
    mutate(sex = "Female") %>%
    mutate(across(is.numeric, ~ round(.x*100, 3)))
  
  contribution_percent_table <- rbind(contribution_percent_table_m, contribution_percent_table_f)
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}

###########################################################################################
###### 10 variable IPM output
########################################################################################### 

#### RENAME VARIABLE1, VARIABLE2, TO VARIABLEA, VARIABLEB 

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
    mutate(ACE = 
             case_when(
               grepl("stac", parameter) & grepl(",1", parameter) ~ "A1", # common A1
               grepl("stac", parameter) & grepl(",2", parameter) ~ "A2", # common A2
               grepl("stcc", parameter) & grepl(",1", parameter) ~ "C1", # common C1
               grepl("stcc", parameter) & grepl(",2", parameter) ~ "C2", # common C2
               grepl("stec", parameter) & grepl(",1", parameter) ~ "E1", # common E1
               grepl("stec", parameter) & grepl(",2", parameter) ~ "E2", # common E2
             )) %>%
    mutate(variable =
             case_when(
               grepl("1,", parameter) ~ variable1, 
               grepl("2,", parameter) ~ variable2, 
               grepl("3,", parameter) ~ variable3, 
               grepl("4,", parameter) ~ variable4, 
               grepl("5,", parameter) ~ variable5, 
               grepl("6,", parameter) ~ variable6, 
               grepl("7,", parameter) ~ variable7, 
               grepl("8,", parameter) ~ variable8, 
               grepl("9,", parameter) ~ variable9, 
               grepl("10,", parameter) ~ variable10, 
             )) %>%
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, lbound, estimate, ubound) 
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_CI.specific <- ACE_CI %>%
    filter(grepl('stas2|stcs2|stes2', parameter)) %>% # select all the specific variables
    filter(grepl('1,1|2,2|3,3|4,4|5,5|6,6|7,7|8,8|9,9|10,10', parameter)) %>% # delete the empty cells in the matrix
    filter(!grepl('1,10', parameter)) %>% # annoyingly the above includes the [1,10] in the matrix but will remove here
    mutate(ACE = 
             case_when(
               grepl("stas", parameter) ~ "A", # specific A
               grepl("stcs", parameter) ~ "C", # specific C
               grepl("stes", parameter) ~ "E", # specific E
             )) %>%
    mutate(variable =
             case_when(
               grepl("1,1", parameter) ~ variable1, 
               grepl("2,2", parameter) ~ variable2, 
               grepl("3,3", parameter) ~ variable3, 
               grepl("4,4", parameter) ~ variable4, 
               grepl("5,5", parameter) ~ variable5, 
               grepl("6,6", parameter) ~ variable6, 
               grepl("7,7", parameter) ~ variable7, 
               grepl("8,8", parameter) ~ variable8, 
               grepl("9,9", parameter) ~ variable9, 
               grepl("10,10", parameter) ~ variable10, 
             )) %>%
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
           Total)
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
           Total)
  
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
           Total)
  
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
    #square root common paths, times together, and divide my appropriate Rph
    mutate(variable1_variable2 = sqrt(variable1_A1 * variable2_A1)/Rph[2,1]) %>% 
    mutate(variable1_variable3 = sqrt(variable1_A1 * variable3_A1)/Rph[3,1]) %>%
    mutate(variable1_variable4 = sqrt(variable1_A1 * variable4_A1)/Rph[4,1]) %>%
    mutate(variable1_variable5 = sqrt(variable1_A1 * variable5_A1)/Rph[5,1]) %>%
    mutate(variable1_variable6 = sqrt(variable1_A1 * variable6_A1)/Rph[6,1]) %>%
    mutate(variable1_variable7 = sqrt(variable1_A1 * variable7_A1)/Rph[7,1]) %>%
    mutate(variable1_variable8 = sqrt(variable1_A1 * variable8_A1)/Rph[8,1]) %>%
    mutate(variable1_variable9 = sqrt(variable1_A1 * variable9_A1)/Rph[9,1]) %>%
    mutate(variable1_variable10 = sqrt(variable1_A1 * variable10_A1)/Rph[10,1]) %>%
    # var 2 associated with all others
    mutate(variable2_variable3 = sqrt(variable2_A1 * variable3_A1)/Rph[3,2]) %>% 
    mutate(variable2_variable4 = sqrt(variable2_A1 * variable4_A1)/Rph[4,2]) %>% 
    mutate(variable2_variable5 = sqrt(variable2_A1 * variable5_A1)/Rph[5,2]) %>% 
    mutate(variable2_variable6 = sqrt(variable2_A1 * variable6_A1)/Rph[6,2]) %>% 
    mutate(variable2_variable7 = sqrt(variable2_A1 * variable7_A1)/Rph[7,2]) %>% 
    mutate(variable2_variable8 = sqrt(variable2_A1 * variable8_A1)/Rph[8,2]) %>% 
    mutate(variable2_variable9 = sqrt(variable2_A1 * variable9_A1)/Rph[9,2]) %>% 
    mutate(variable2_variable10 = sqrt(variable2_A1 * variable10_A1)/Rph[10,2]) %>% 
    # var 3 associated with all others
    mutate(variable3_variable4 = sqrt(variable3_A1 * variable4_A1)/Rph[4,3]) %>% 
    mutate(variable3_variable5 = sqrt(variable3_A1 * variable5_A1)/Rph[5,3]) %>% 
    mutate(variable3_variable6 = sqrt(variable3_A1 * variable6_A1)/Rph[6,3]) %>% 
    mutate(variable3_variable7 = sqrt(variable3_A1 * variable7_A1)/Rph[7,3]) %>% 
    mutate(variable3_variable8 = sqrt(variable3_A1 * variable8_A1)/Rph[8,3]) %>% 
    mutate(variable3_variable9 = sqrt(variable3_A1 * variable9_A1)/Rph[9,3]) %>% 
    mutate(variable3_variable10 = sqrt(variable3_A1 * variable10_A1)/Rph[10,3]) %>% 
    # var 4 associated with all others
    mutate(variable4_variable5 = sqrt(variable4_A1 * variable5_A1)/Rph[5,4]) %>% 
    mutate(variable4_variable6 = sqrt(variable4_A1 * variable6_A1)/Rph[6,4]) %>% 
    mutate(variable4_variable7 = sqrt(variable4_A1 * variable7_A1)/Rph[7,4]) %>% 
    mutate(variable4_variable8 = sqrt(variable4_A1 * variable8_A1)/Rph[8,4]) %>% 
    mutate(variable4_variable9 = sqrt(variable4_A1 * variable9_A1)/Rph[9,4]) %>% 
    mutate(variable4_variable10 = sqrt(variable4_A1 * variable10_A1)/Rph[10,4]) %>% 
    # var 5 associated with all others
    mutate(variable5_variable6 = sqrt(variable5_A1 * variable6_A1)/Rph[6,5]) %>% 
    mutate(variable5_variable7 = sqrt(variable5_A1 * variable7_A1)/Rph[7,5]) %>% 
    mutate(variable5_variable8 = sqrt(variable5_A1 * variable8_A1)/Rph[8,5]) %>% 
    mutate(variable5_variable9 = sqrt(variable5_A1 * variable9_A1)/Rph[9,5]) %>% 
    mutate(variable5_variable10 = sqrt(variable5_A1 * variable10_A1)/Rph[10,5]) %>% 
    # var 6 associated with all others
    mutate(variable6_variable7 = (sqrt(variable6_A1 * variable7_A1) + sqrt(variable6_A2 * variable7_A2))/Rph[7,6]) %>% 
    mutate(variable6_variable8 = (sqrt(variable6_A1 * variable8_A1) + sqrt(variable6_A2 * variable8_A2))/Rph[8,6]) %>% 
    mutate(variable6_variable9 = (sqrt(variable6_A1 * variable9_A1) + sqrt(variable6_A2 * variable9_A2))/Rph[9,6]) %>% 
    mutate(variable6_variable10 = (sqrt(variable6_A1 * variable10_A1) + sqrt(variable6_A2 * variable10_A2))/Rph[10,6]) %>% 
    # var 7 associated with all others
    mutate(variable7_variable8 = (sqrt(variable7_A1 * variable8_A1) + sqrt(variable7_A2 * variable8_A2))/Rph[8,7]) %>% 
    mutate(variable7_variable9 = (sqrt(variable7_A1 * variable9_A1) + sqrt(variable7_A2 * variable9_A2))/Rph[9,7]) %>% 
    mutate(variable7_variable10 = (sqrt(variable7_A1 * variable10_A1) + sqrt(variable7_A2 * variable10_A2))/Rph[10,7]) %>% 
    # var 8 associated with all others
    mutate(variable8_variable9 = (sqrt(variable8_A1 * variable9_A1) + sqrt(variable8_A2 * variable9_A2))/Rph[9,8]) %>% 
    mutate(variable8_variable10 = (sqrt(variable8_A1 * variable10_A1) + sqrt(variable8_A2 * variable10_A2))/Rph[10,8]) %>% 
    # var 9 associated with 10
    mutate(variable9_variable10 = (sqrt(variable9_A1 * variable10_A1) + sqrt(variable9_A2 * variable10_A2))/Rph[10,9]) %>% 
    # select the newly created variables
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
    mutate(across(is.numeric, ~ round(.x*100, 3)))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}


###########################################################################################
###### 10 variable IPM output - NO CONFIDENCE INTERVALS
###########################################################################################

# You may want to run your scripts without calculating confidence intervals (intervals = FALSE) to save time. This output function will display the output when your model mxTryHard(model, intervals = FALSE). 

IPM_esitmates_ACE_10var_noCI <- function(data, variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10, model, uniACEestimates, pathestimates, A_percent, C_percent, E_percent, common_factor_contribution){ 
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
  
  # create a table containing estimates for common factors (squared and standardised)
  ACE_stac2 <- as.data.frame(mxEval(ACE.stac2, data)) %>% 
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10)) %>%
    select(variable, A1 = V1, A2 = V2) %>%
    pivot_longer(cols = starts_with("A"), names_to = "ACE", values_to = "estimate")
  
  ACE_stcc2 <- as.data.frame(mxEval(ACE.stcc2, data)) %>% 
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10)) %>%
    select(variable, C1 = V1, C2 = V2) %>%
    pivot_longer(cols = starts_with("C"), names_to = "ACE", values_to = "estimate")
  
  ACE_stec2 <- as.data.frame(mxEval(ACE.stec2, data)) %>% 
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10)) %>%
    select(variable, E1 = V1, E2 = V2) %>%
    pivot_longer(cols = starts_with("E"), names_to = "ACE", values_to = "estimate")
  
  ACE.common <- rbind(ACE_stac2, ACE_stcc2, ACE_stec2)
  
  ACE.common <- ACE.common %>% 
    mutate(factor = "Common") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, estimate)
  
  # create a table containing (non-zero) estimates for specific factors (squared and standardised)
  ACE_stas2 <- as.data.frame(diag(mxEval(ACE.stas2, data))) 
  colnames(ACE_stas2) <- c("estimate")
  ACE_stas2 <- ACE_stas2 %>%
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10)) %>%
    mutate(ACE = "A") 
  
  ACE_stcs2 <- as.data.frame(diag(mxEval(ACE.stcs2, data))) 
  colnames(ACE_stcs2) <- c("estimate")
  ACE_stcs2 <- ACE_stcs2 %>%
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10)) %>%
    mutate(ACE = "C") 
  
  ACE_stes2 <- as.data.frame(diag(mxEval(ACE.stes2, data))) 
  colnames(ACE_stes2) <- c("estimate")
  ACE_stes2 <- ACE_stes2 %>%
    mutate(variable = c(variable1, variable2, variable3, variable4, variable5, variable6, variable7, variable8, variable9, variable10)) %>%
    mutate(ACE = "E")

  ACE.specific <- rbind(ACE_stas2, ACE_stcs2, ACE_stes2)
  
  ACE.specific <- ACE.specific %>% 
    mutate(factor = "Specific") %>%
    mutate(model = model) %>%
    select(model, factor, ACE, variable, estimate)
  
  # combine the common and specific path estimates
  pathestimates_table <- rbind(ACE.common, ACE.specific)
  
  pathestimates_table_rounded <- pathestimates_table %>%
    mutate(across(is.numeric, round, digits = 3))
  
  if(pathestimates == TRUE){print(pathestimates_table_rounded)} 
  
  # create table that shows out of the total contribution of all A factors, what % is due to common and specific 
  A_percentages <- pathestimates_table %>% 
    filter(grepl('A1|A2|A', ACE)) %>%
    select(-factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = A1 + A2 + A) %>%
    mutate(across(A1:A, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `A1c %` = A1,
           `A2c %` = A2,
           `As %` = A,
           Total)
  if(A_percent == TRUE){print(A_percentages)}
  
  # create table that shows out of the total contribution of all C factors, what % is due to common and specific 
  C_percentages <- pathestimates_table %>%
    filter(grepl('C1|C2|C', ACE)) %>%
    select(-factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = C1 + C2 + C) %>%
    mutate(across(C1:C, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `C1c %` = C1,
           `C2c %` = C2,
           `Cs %` = C,
           Total)
  
  if(C_percent == TRUE){print(C_percentages)}
  
  # create table that shows out of the total contribution of all E factors, what % is due to common and specific 
  E_percentages <- pathestimates_table %>%
    filter(grepl('E1|E2|E', ACE)) %>%
    select(-factor) %>%
    pivot_wider(names_from = ACE, values_from = estimate) %>%
    mutate(Total = E1 + E2 + E) %>%
    mutate(across(E1:E, ~ round(.x/Total*100, 3))) %>%
    select(variable,
           `E1c %` = E1,
           `E2c %` = E2,
           `Es %` = E,
           Total)
  
  if(E_percent == TRUE){print(E_percentages)}
  
  # create a table using paths tracing to calculate the % of the phenotypic correlation that is due to the A, C, and E common factors. This essentially shows "how much of the association between variable1 and variable2 is due to shared genetic and environmental factors". 
  # First, contribution of A to the phenotypic correlation
  Ra_contribution_percent <- pathestimates_table %>%
    filter(grepl('A1|A2', ACE)) %>%
    filter(ACE != "A2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4 & variable != variable5)) %>% # remove common factor A2 for variable1 to 5
    select(-factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variablea"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variableb"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variablec"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variabled"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variablee"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variablef"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variableg"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variableh"))) %>% 
    mutate(across(variable, ~ str_replace(., variable9, "variablei"))) %>% 
    mutate(across(variable, ~ str_replace(., variable10, "variablej"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    #square root common paths, times together, and divide my appropriate Rph
    mutate(variablea_variableb = sqrt(variablea_A1 * variableb_A1)/Rph[2,1]) %>% 
    mutate(variablea_variablec = sqrt(variablea_A1 * variablec_A1)/Rph[3,1]) %>%
    mutate(variablea_variabled = sqrt(variablea_A1 * variabled_A1)/Rph[4,1]) %>%
    mutate(variablea_variablee = sqrt(variablea_A1 * variablee_A1)/Rph[5,1]) %>%
    mutate(variablea_variablef = sqrt(variablea_A1 * variablef_A1)/Rph[6,1]) %>%
    mutate(variablea_variableg = sqrt(variablea_A1 * variableg_A1)/Rph[7,1]) %>%
    mutate(variablea_variableh = sqrt(variablea_A1 * variableh_A1)/Rph[8,1]) %>%
    mutate(variablea_variablei = sqrt(variablea_A1 * variablei_A1)/Rph[9,1]) %>%
    mutate(variablea_variablej = sqrt(variablea_A1 * variablej_A1)/Rph[10,1]) %>%
    # var 2 associated with all others
    mutate(variableb_variablec = sqrt(variableb_A1 * variablec_A1)/Rph[3,2]) %>% 
    mutate(variableb_variabled = sqrt(variableb_A1 * variabled_A1)/Rph[4,2]) %>% 
    mutate(variableb_variablee = sqrt(variableb_A1 * variablee_A1)/Rph[5,2]) %>% 
    mutate(variableb_variablef = sqrt(variableb_A1 * variablef_A1)/Rph[6,2]) %>% 
    mutate(variableb_variableg = sqrt(variableb_A1 * variableg_A1)/Rph[7,2]) %>% 
    mutate(variableb_variableh = sqrt(variableb_A1 * variableh_A1)/Rph[8,2]) %>% 
    mutate(variableb_variablei = sqrt(variableb_A1 * variablei_A1)/Rph[9,2]) %>% 
    mutate(variableb_variablej = sqrt(variableb_A1 * variablej_A1)/Rph[10,2]) %>% 
    # var 3 associated with all others
    mutate(variablec_variabled = sqrt(variablec_A1 * variabled_A1)/Rph[4,3]) %>% 
    mutate(variablec_variablee = sqrt(variablec_A1 * variablee_A1)/Rph[5,3]) %>% 
    mutate(variablec_variablef = sqrt(variablec_A1 * variablef_A1)/Rph[6,3]) %>% 
    mutate(variablec_variableg = sqrt(variablec_A1 * variableg_A1)/Rph[7,3]) %>% 
    mutate(variablec_variableh = sqrt(variablec_A1 * variableh_A1)/Rph[8,3]) %>% 
    mutate(variablec_variablei = sqrt(variablec_A1 * variablei_A1)/Rph[9,3]) %>% 
    mutate(variablec_variablej = sqrt(variablec_A1 * variablej_A1)/Rph[10,3]) %>% 
    # var 4 associated with all others
    mutate(variabled_variablee = sqrt(variabled_A1 * variablee_A1)/Rph[5,4]) %>% 
    mutate(variabled_variablef = sqrt(variabled_A1 * variablef_A1)/Rph[6,4]) %>% 
    mutate(variabled_variableg = sqrt(variabled_A1 * variableg_A1)/Rph[7,4]) %>% 
    mutate(variabled_variableh = sqrt(variabled_A1 * variableh_A1)/Rph[8,4]) %>% 
    mutate(variabled_variablei = sqrt(variabled_A1 * variablei_A1)/Rph[9,4]) %>% 
    mutate(variabled_variablej = sqrt(variabled_A1 * variablej_A1)/Rph[10,4]) %>% 
    # var 5 associated with all others
    mutate(variablee_variablef = sqrt(variablee_A1 * variablef_A1)/Rph[6,5]) %>% 
    mutate(variablee_variableg = sqrt(variablee_A1 * variableg_A1)/Rph[7,5]) %>% 
    mutate(variablee_variableh = sqrt(variablee_A1 * variableh_A1)/Rph[8,5]) %>% 
    mutate(variablee_variablei = sqrt(variablee_A1 * variablei_A1)/Rph[9,5]) %>% 
    mutate(variablee_variablej = sqrt(variablee_A1 * variablej_A1)/Rph[10,5]) %>% 
    # var 6 associated with all others
    mutate(variablef_variableg = (sqrt(variablef_A1 * variableg_A1) + sqrt(variablef_A2 * variableg_A2))/Rph[7,6]) %>% 
    mutate(variablef_variableh = (sqrt(variablef_A1 * variableh_A1) + sqrt(variablef_A2 * variableh_A2))/Rph[8,6]) %>% 
    mutate(variablef_variablei = (sqrt(variablef_A1 * variablei_A1) + sqrt(variablef_A2 * variablei_A2))/Rph[9,6]) %>% 
    mutate(variablef_variablej = (sqrt(variablef_A1 * variablej_A1) + sqrt(variablef_A2 * variablej_A2))/Rph[10,6]) %>% 
    # var 7 associated with all others
    mutate(variableg_variableh = (sqrt(variableg_A1 * variableh_A1) + sqrt(variableg_A2 * variableh_A2))/Rph[8,7]) %>% 
    mutate(variableg_variablei = (sqrt(variableg_A1 * variablei_A1) + sqrt(variableg_A2 * variablei_A2))/Rph[9,7]) %>% 
    mutate(variableg_variablej = (sqrt(variableg_A1 * variablej_A1) + sqrt(variableg_A2 * variablej_A2))/Rph[10,7]) %>% 
    # var 8 associated with all others
    mutate(variableh_variablei = (sqrt(variableh_A1 * variablei_A1) + sqrt(variableh_A2 * variablei_A2))/Rph[9,8]) %>% 
    mutate(variableh_variablej = (sqrt(variableh_A1 * variablej_A1) + sqrt(variableh_A2 * variablej_A2))/Rph[10,8]) %>% 
    # var 9 associated with 10
    mutate(variablei_variablej = (sqrt(variablei_A1 * variablej_A1) + sqrt(variablei_A2 * variablej_A2))/Rph[10,9]) %>% 
    # select the newly created variables
    select(variablea_variableb:variablei_variablej) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ac % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablea", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableb", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablec", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variabled", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablee", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablef", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableg", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableh", variable8))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablei", variable9))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablej", variable10)))
  
  # Second, contribution of C to the phenotypic correlation
  Rc_contribution_percent <- pathestimates_table %>%
    filter(grepl('C1|C2', ACE)) %>%
    filter(ACE != "C2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4 & variable != variable5)) %>% # remove common factor A2 for variable1 to 5
    select(-factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variablea"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variableb"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variablec"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variabled"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variablee"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variablef"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variableg"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variableh"))) %>% 
    mutate(across(variable, ~ str_replace(., variable9, "variablei"))) %>% 
    mutate(across(variable, ~ str_replace(., variable10, "variablej"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    #square root common paths, times together, and divide my appropriate Rph
    mutate(variablea_variableb = sqrt(variablea_C1 * variableb_C1)/Rph[2,1]) %>% 
    mutate(variablea_variablec = sqrt(variablea_C1 * variablec_C1)/Rph[3,1]) %>%
    mutate(variablea_variabled = sqrt(variablea_C1 * variabled_C1)/Rph[4,1]) %>%
    mutate(variablea_variablee = sqrt(variablea_C1 * variablee_C1)/Rph[5,1]) %>%
    mutate(variablea_variablef = sqrt(variablea_C1 * variablef_C1)/Rph[6,1]) %>%
    mutate(variablea_variableg = sqrt(variablea_C1 * variableg_C1)/Rph[7,1]) %>%
    mutate(variablea_variableh = sqrt(variablea_C1 * variableh_C1)/Rph[8,1]) %>%
    mutate(variablea_variablei = sqrt(variablea_C1 * variablei_C1)/Rph[9,1]) %>%
    mutate(variablea_variablej = sqrt(variablea_C1 * variablej_C1)/Rph[10,1]) %>%
    # var 2 associated with all others
    mutate(variableb_variablec = sqrt(variableb_C1 * variablec_C1)/Rph[3,2]) %>% 
    mutate(variableb_variabled = sqrt(variableb_C1 * variabled_C1)/Rph[4,2]) %>% 
    mutate(variableb_variablee = sqrt(variableb_C1 * variablee_C1)/Rph[5,2]) %>% 
    mutate(variableb_variablef = sqrt(variableb_C1 * variablef_C1)/Rph[6,2]) %>% 
    mutate(variableb_variableg = sqrt(variableb_C1 * variableg_C1)/Rph[7,2]) %>% 
    mutate(variableb_variableh = sqrt(variableb_C1 * variableh_C1)/Rph[8,2]) %>% 
    mutate(variableb_variablei = sqrt(variableb_C1 * variablei_C1)/Rph[9,2]) %>% 
    mutate(variableb_variablej = sqrt(variableb_C1 * variablej_C1)/Rph[10,2]) %>% 
    # var 3 associated with all others
    mutate(variablec_variabled = sqrt(variablec_C1 * variabled_C1)/Rph[4,3]) %>% 
    mutate(variablec_variablee = sqrt(variablec_C1 * variablee_C1)/Rph[5,3]) %>% 
    mutate(variablec_variablef = sqrt(variablec_C1 * variablef_C1)/Rph[6,3]) %>% 
    mutate(variablec_variableg = sqrt(variablec_C1 * variableg_C1)/Rph[7,3]) %>% 
    mutate(variablec_variableh = sqrt(variablec_C1 * variableh_C1)/Rph[8,3]) %>% 
    mutate(variablec_variablei = sqrt(variablec_C1 * variablei_C1)/Rph[9,3]) %>% 
    mutate(variablec_variablej = sqrt(variablec_C1 * variablej_C1)/Rph[10,3]) %>% 
    # var 4 associated with all others
    mutate(variabled_variablee = sqrt(variabled_C1 * variablee_C1)/Rph[5,4]) %>% 
    mutate(variabled_variablef = sqrt(variabled_C1 * variablef_C1)/Rph[6,4]) %>% 
    mutate(variabled_variableg = sqrt(variabled_C1 * variableg_C1)/Rph[7,4]) %>% 
    mutate(variabled_variableh = sqrt(variabled_C1 * variableh_C1)/Rph[8,4]) %>% 
    mutate(variabled_variablei = sqrt(variabled_C1 * variablei_C1)/Rph[9,4]) %>% 
    mutate(variabled_variablej = sqrt(variabled_C1 * variablej_C1)/Rph[10,4]) %>% 
    # var 5 associated with all others
    mutate(variablee_variablef = sqrt(variablee_C1 * variablef_C1)/Rph[6,5]) %>% 
    mutate(variablee_variableg = sqrt(variablee_C1 * variableg_C1)/Rph[7,5]) %>% 
    mutate(variablee_variableh = sqrt(variablee_C1 * variableh_C1)/Rph[8,5]) %>% 
    mutate(variablee_variablei = sqrt(variablee_C1 * variablei_C1)/Rph[9,5]) %>% 
    mutate(variablee_variablej = sqrt(variablee_C1 * variablej_C1)/Rph[10,5]) %>% 
    # var 6 associated with all others
    mutate(variablef_variableg = (sqrt(variablef_C1 * variableg_C1) + sqrt(variablef_C2 * variableg_C2))/Rph[7,6]) %>% 
    mutate(variablef_variableh = (sqrt(variablef_C1 * variableh_C1) + sqrt(variablef_C2 * variableh_C2))/Rph[8,6]) %>% 
    mutate(variablef_variablei = (sqrt(variablef_C1 * variablei_C1) + sqrt(variablef_C2 * variablei_C2))/Rph[9,6]) %>% 
    mutate(variablef_variablej = (sqrt(variablef_C1 * variablej_C1) + sqrt(variablef_C2 * variablej_C2))/Rph[10,6]) %>% 
    # var 7 associated with all others
    mutate(variableg_variableh = (sqrt(variableg_C1 * variableh_C1) + sqrt(variableg_C2 * variableh_C2))/Rph[8,7]) %>% 
    mutate(variableg_variablei = (sqrt(variableg_C1 * variablei_C1) + sqrt(variableg_C2 * variablei_C2))/Rph[9,7]) %>% 
    mutate(variableg_variablej = (sqrt(variableg_C1 * variablej_C1) + sqrt(variableg_C2 * variablej_C2))/Rph[10,7]) %>% 
    # var 8 associated with all others
    mutate(variableh_variablei = (sqrt(variableh_C1 * variablei_C1) + sqrt(variableh_C2 * variablei_C2))/Rph[9,8]) %>% 
    mutate(variableh_variablej = (sqrt(variableh_C1 * variablej_C1) + sqrt(variableh_C2 * variablej_C2))/Rph[10,8]) %>% 
    # var 9 associated with 10
    mutate(variablei_variablej = (sqrt(variablei_C1 * variablej_C1) + sqrt(variablei_C2 * variablej_C2))/Rph[10,9]) %>% 
    # select the newly created variables
    select(variablea_variableb:variablei_variablej) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Cc % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablea", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableb", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablec", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variabled", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablee", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablef", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableg", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableh", variable8))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablei", variable9))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablej", variable10)))
  
  # Third, contribution of E to the phenotypic correlation
  Re_contribution_percent <- pathestimates_table %>%
    filter(grepl('E1|E2', ACE)) %>%
    filter(ACE != "E2" | (variable != variable1 & variable != variable2 & variable != variable3 & variable != variable4 & variable != variable5)) %>% # remove common factor A2 for variable1 to 5
    select(-factor) %>%
    mutate(across(variable, ~ str_replace(., variable1, "variablea"))) %>% # rename variables to apply equations
    mutate(across(variable, ~ str_replace(., variable2, "variableb"))) %>% 
    mutate(across(variable, ~ str_replace(., variable3, "variablec"))) %>% 
    mutate(across(variable, ~ str_replace(., variable4, "variabled"))) %>% 
    mutate(across(variable, ~ str_replace(., variable5, "variablee"))) %>% 
    mutate(across(variable, ~ str_replace(., variable6, "variablef"))) %>% 
    mutate(across(variable, ~ str_replace(., variable7, "variableg"))) %>% 
    mutate(across(variable, ~ str_replace(., variable8, "variableh"))) %>% 
    mutate(across(variable, ~ str_replace(., variable9, "variablei"))) %>% 
    mutate(across(variable, ~ str_replace(., variable10, "variablej"))) %>% 
    pivot_wider(names_from = c(variable, ACE), values_from = estimate) %>%
    #square root common paths, times together, and divide my appropriate Rph
    mutate(variablea_variableb = sqrt(variablea_E1 * variableb_E1)/Rph[2,1]) %>% 
    mutate(variablea_variablec = sqrt(variablea_E1 * variablec_E1)/Rph[3,1]) %>%
    mutate(variablea_variabled = sqrt(variablea_E1 * variabled_E1)/Rph[4,1]) %>%
    mutate(variablea_variablee = sqrt(variablea_E1 * variablee_E1)/Rph[5,1]) %>%
    mutate(variablea_variablef = sqrt(variablea_E1 * variablef_E1)/Rph[6,1]) %>%
    mutate(variablea_variableg = sqrt(variablea_E1 * variableg_E1)/Rph[7,1]) %>%
    mutate(variablea_variableh = sqrt(variablea_E1 * variableh_E1)/Rph[8,1]) %>%
    mutate(variablea_variablei = sqrt(variablea_E1 * variablei_E1)/Rph[9,1]) %>%
    mutate(variablea_variablej = sqrt(variablea_E1 * variablej_E1)/Rph[10,1]) %>%
    # var 2 associated with all others
    mutate(variableb_variablec = sqrt(variableb_E1 * variablec_E1)/Rph[3,2]) %>% 
    mutate(variableb_variabled = sqrt(variableb_E1 * variabled_E1)/Rph[4,2]) %>% 
    mutate(variableb_variablee = sqrt(variableb_E1 * variablee_E1)/Rph[5,2]) %>% 
    mutate(variableb_variablef = sqrt(variableb_E1 * variablef_E1)/Rph[6,2]) %>% 
    mutate(variableb_variableg = sqrt(variableb_E1 * variableg_E1)/Rph[7,2]) %>% 
    mutate(variableb_variableh = sqrt(variableb_E1 * variableh_E1)/Rph[8,2]) %>% 
    mutate(variableb_variablei = sqrt(variableb_E1 * variablei_E1)/Rph[9,2]) %>% 
    mutate(variableb_variablej = sqrt(variableb_E1 * variablej_E1)/Rph[10,2]) %>% 
    # var 3 associated with all others
    mutate(variablec_variabled = sqrt(variablec_E1 * variabled_E1)/Rph[4,3]) %>% 
    mutate(variablec_variablee = sqrt(variablec_E1 * variablee_E1)/Rph[5,3]) %>% 
    mutate(variablec_variablef = sqrt(variablec_E1 * variablef_E1)/Rph[6,3]) %>% 
    mutate(variablec_variableg = sqrt(variablec_E1 * variableg_E1)/Rph[7,3]) %>% 
    mutate(variablec_variableh = sqrt(variablec_E1 * variableh_E1)/Rph[8,3]) %>% 
    mutate(variablec_variablei = sqrt(variablec_E1 * variablei_E1)/Rph[9,3]) %>% 
    mutate(variablec_variablej = sqrt(variablec_E1 * variablej_E1)/Rph[10,3]) %>% 
    # var 4 associated with all others
    mutate(variabled_variablee = sqrt(variabled_E1 * variablee_E1)/Rph[5,4]) %>% 
    mutate(variabled_variablef = sqrt(variabled_E1 * variablef_E1)/Rph[6,4]) %>% 
    mutate(variabled_variableg = sqrt(variabled_E1 * variableg_E1)/Rph[7,4]) %>% 
    mutate(variabled_variableh = sqrt(variabled_E1 * variableh_E1)/Rph[8,4]) %>% 
    mutate(variabled_variablei = sqrt(variabled_E1 * variablei_E1)/Rph[9,4]) %>% 
    mutate(variabled_variablej = sqrt(variabled_E1 * variablej_E1)/Rph[10,4]) %>% 
    # var 5 associated with all others
    mutate(variablee_variablef = sqrt(variablee_E1 * variablef_E1)/Rph[6,5]) %>% 
    mutate(variablee_variableg = sqrt(variablee_E1 * variableg_E1)/Rph[7,5]) %>% 
    mutate(variablee_variableh = sqrt(variablee_E1 * variableh_E1)/Rph[8,5]) %>% 
    mutate(variablee_variablei = sqrt(variablee_E1 * variablei_E1)/Rph[9,5]) %>% 
    mutate(variablee_variablej = sqrt(variablee_E1 * variablej_E1)/Rph[10,5]) %>% 
    # var 6 associated with all others
    mutate(variablef_variableg = (sqrt(variablef_E1 * variableg_E1) + sqrt(variablef_E2 * variableg_E2))/Rph[7,6]) %>% 
    mutate(variablef_variableh = (sqrt(variablef_E1 * variableh_E1) + sqrt(variablef_E2 * variableh_E2))/Rph[8,6]) %>% 
    mutate(variablef_variablei = (sqrt(variablef_E1 * variablei_E1) + sqrt(variablef_E2 * variablei_E2))/Rph[9,6]) %>% 
    mutate(variablef_variablej = (sqrt(variablef_E1 * variablej_E1) + sqrt(variablef_E2 * variablej_E2))/Rph[10,6]) %>% 
    # var 7 associated with all others
    mutate(variableg_variableh = (sqrt(variableg_E1 * variableh_E1) + sqrt(variableg_E2 * variableh_E2))/Rph[8,7]) %>% 
    mutate(variableg_variablei = (sqrt(variableg_E1 * variablei_E1) + sqrt(variableg_E2 * variablei_E2))/Rph[9,7]) %>% 
    mutate(variableg_variablej = (sqrt(variableg_E1 * variablej_E1) + sqrt(variableg_E2 * variablej_E2))/Rph[10,7]) %>% 
    # var 8 associated with all others
    mutate(variableh_variablei = (sqrt(variableh_E1 * variablei_E1) + sqrt(variableh_E2 * variablei_E2))/Rph[9,8]) %>% 
    mutate(variableh_variablej = (sqrt(variableh_E1 * variablej_E1) + sqrt(variableh_E2 * variablej_E2))/Rph[10,8]) %>% 
    # var 9 associated with 10
    mutate(variablei_variablej = (sqrt(variablei_E1 * variablej_E1) + sqrt(variablei_E2 * variablej_E2))/Rph[10,9]) %>% 
    # select the newly created variables
    select(variablea_variableb:variablei_variablej) %>%  
    pivot_longer(cols = starts_with("variable"), names_to = "phenotypic correlation", values_to = "Ec % contribution") %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablea", variable1))) %>% # name the variables back 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableb", variable2))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablec", variable3))) %>% 
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variabled", variable4))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablee", variable5))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablef", variable6))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableg", variable7))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variableh", variable8))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablei", variable9))) %>%
    mutate(across(`phenotypic correlation`, ~ str_replace(., "variablej", variable10)))
  
  # bind the contribution of A, C, and E together into one table
  contribution_percent_table <- plyr::join_all(list(Ra_contribution_percent, Rc_contribution_percent, Re_contribution_percent),
                                               by = "phenotypic correlation")
  
  contribution_percent_table <- contribution_percent_table %>% 
    mutate(total = `Ac % contribution` + `Cc % contribution` + `Ec % contribution`) %>% 
    mutate(across(is.numeric, ~ round(.x*100, 3)))
  
  if(common_factor_contribution == TRUE){print(contribution_percent_table)}
}
