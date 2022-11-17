# ---
#   title: "Functions for all scripts looking at the overlap between social isolation and mental health" 
# ---

# clear the global environment - remove when happy
remove(list = ls())

# Load packages
library(knitr)
library(psych)
library(OpenMx)
library(tidyr)
library(tidyverse)
library(dplyr)     # conflicts with tidyverse for e.g. rename and row_number


###############################################
######## Independent pathway model ############
###############################################

# TESTING
# selvars_chosen_anx <- c("sisoe12norm_reg", "masce12norm_reg", "socisoe18norm_reg", "gadsxe18norm_reg", 
#                         "sisoy12norm_reg", "mascy12norm_reg", "socisoy18norm_reg", "gadsxy18norm_reg")
# 
# 
# IPM_ACE_2time(dataMZ = dat.twin.MZ,
#               dataDZ = dat.twin.DZ,
#               nv = 4,
#               selvars_chosen = selvars_chosen_anx,
#               model_name = "IPM_anxiety")

# The IPM_ACE_2time function runs an OpenMx independent pathway model on your data. To accurately specify this model, you will need longitudinal data consisting of several variables at two distinct time points. For example: anxiety and depression at age 12, and then measured again at age 18. 
# This IPM will estimate 1 set of ACE specific factors for each variable, as well as 2 sets of ACE common factors. All avriables load onto the first set of common ACE factors (e.g. genetic effects that are common across all variables regardless of the time point). Only variables introduced at the *new* time point load onto the second common factor (e.g. new genetic effects that occur at age 18 but did not exist at age 12).
# This function will estimate the model using OpenMx mxTryHard() which will run an extra 10 iterations of your model and pick the best solution. 
# This function is reliant on "OpenMx", "tidyverse" and "dplyr".
# The nv should be a multiple of ncf.

nv = 6
nfAc = 3

########### IPM function to run ACE model at two time points ############ THINK ABOUT NAME OF FUNCTION

IPM_ace <- function(dataMZ, dataDZ, nv, ncf, selvars_chosen, model_name){ 
  
  # number of variables
  nv = nv    					      # number of variables - e.g. 2 at age 12 and 2 at age 18
  ntv <- nv*2 				      # number of twin variables
  nlower <- nv*(nv+1)/2 		# number of free elements in an nv*nv lower matrix 
  
  # number of factors A (common)
  nfAc <- ncf # change dimension of A factor matrix Ac to have *2* common factors
  
  # Set if the parameters are estimated (TRUE) or not (FALSE) 
  if(nfAc == 2){Ac2Free <- c(rep(TRUE, nv), rep(FALSE, nv/nfAc), rep(TRUE, nv/nfAc))}
  if(nfAc == 3){Ac2Free <- c(rep(TRUE, nv), rep(FALSE, nv/nfAc), rep(TRUE, nv*2/nfAc), rep(FALSE, nv*2/nfAc), rep(TRUE, nv/nfAc))} # FOR THIS WOULD NEED TO WORK OUT IF IT BEHAVES LIKE CHOLESKY OR IF THIRD IS COMPLETELY SEPARATE TO 2ND TIME POINT 
  
  # Create start values for 2 common ACE Factors - set the start values for common 
  Ac2Values <- c(rep(.1, nv),  rep(0, nv/2), rep(.1, nv/2))  
  Cc2Values <- c(rep(.1, nv),  rep(0, nv/2), rep(.1, nv/2))
  Ec2Values <- c(rep(.1, nv),  rep(0, nv/2), rep(.1, nv/2))
  
  IPM <- mxModel(model_name,
                 mxModel("ACE",
                        # Common ace paths
                        mxMatrix("Full", nrow = nv, ncol = nfAc, free = Ac2Free, values = Ac2Values, name = "ac"),
                        mxMatrix("Full", nrow = nv, ncol = nfAc, free = Ac2Free, values = Cc2Values, name = "cc"),
                        mxMatrix("Full", nrow = nv, ncol = nfAc, free = Ac2Free, values = Ec2Values, name = "ec"),
                        
                        # Specific ace paths
                        mxMatrix("Diag", nrow = nv, ncol = nv, free = TRUE, values = 0.3,  name = "as"),
                        mxMatrix("Diag", nrow = nv, ncol = nv, free = TRUE, values = 0.3,  name = "cs"),
                        mxMatrix("Diag", nrow = nv, ncol = nv, free = TRUE, values = 0.3,  name = "es"),
                        
                        # ACE variance components
                        mxAlgebra(ac %*% t(ac) + as %*% t(as), name = "A"), 
                        mxAlgebra(cc %*% t(cc) + cs %*% t(cs), name = "C"),
                        mxAlgebra(ec %*% t(ec) + es %*% t(es), name = "E"),
                        mxAlgebra(A + C + E, name = "V"),
                        
                        # Calculate standard deviation for standardization
                        mxMatrix("Iden", nrow = nv, ncol = nv, name = "I"),
                        mxAlgebra(solve(sqrt(I*V)), name = "iSD"), # extracting the variances from the diagonal of the covariance matrix V 
                        
                        # Expected means vector - one grand mean
                        mxMatrix("Full", nrow = 1, ncol = ntv, free = TRUE, values = 0.01, name = "expMean"), 
                        
                        # Genetic and environmental correlations (correlated factor solution)
                        mxAlgebra(solve(sqrt(I*V)) %&% V, name = "Rph"),
                        mxAlgebra(solve(sqrt(I*A)) %&% A, name = "Ra"),
                        mxAlgebra(solve(sqrt(I*C)) %&% C, name = "Rc"),
                        mxAlgebra(solve(sqrt(I*E)) %&% E, name = "Re"),
                        
                        # Standardised components of variance
                        mxAlgebra(A/V, name = "h2"),
                        mxAlgebra(C/V, name = "c2"),
                        mxAlgebra(E/V, name = "e2"),
                        
                        # Standardise common parameters (unsquared)
                        mxAlgebra(iSD %*% ac, name = "stac"),
                        mxAlgebra(iSD %*% cc, name = "stcc"),
                        mxAlgebra(iSD %*% ec, name = "stec"),
                        
                        # Standardise common parameters (Squared)
                        mxAlgebra(stac*stac, name = "stac2"),
                        mxAlgebra(stcc*stcc, name = "stcc2"),
                        mxAlgebra(stec*stec, name = "stec2"),
                        
                        # Standardise specific parameters (unsquared)
                        mxAlgebra(iSD %*% as, name = "stas"),
                        mxAlgebra(iSD %*% cs, name = "stcs"),
                        mxAlgebra(iSD %*% es, name = "stes"),
                        
                        # Standardise specific parameters (Squared)
                        mxAlgebra(stas*stas, name = "stas2"),
                        mxAlgebra(stcs*stcs, name = "stcs2"),
                        mxAlgebra(stes*stes, name = "stes2"),	
                        
                        # MZ expected covariance matrix
                        mxAlgebra(rbind(cbind(A+C+E, A+C), cbind(A+C, A+C+E)), name = "expCovMZ"),
                        
                        # DZ expected covariance matrix
                        mxAlgebra(rbind(cbind(A+C+E, (0.5%x%A)+C), cbind((0.5%x%A)+C, A+C+E)), name = "expCovDZ")
                                 
                         ),
                         
                         mxModel("MZ",
                                 mxData(observed = dataMZ, type = "raw"),
                                 mxExpectationNormal(covariance = "ACE.expCovMZ", means = "ACE.expMean", dimnames = selvars_chosen),
                                 mxFitFunctionML()
                         ),
                         mxModel("DZ",
                                 mxData(observed = dataDZ, type = "raw"),
                                 mxExpectationNormal(covariance = "ACE.expCovDZ", means = "ACE.expMean", dimnames = selvars_chosen),
                                 mxFitFunctionML()
                         ),
                         
                         mxAlgebra(MZ.objective + DZ.objective, name = "m2LL"),
                         mxFitFunctionAlgebra("m2LL"),
                         mxCI(c("ACE.stac2", "ACE.stcc2", "ACE.stec2", "ACE.stas2", "ACE.stcs2", "ACE.stes2", "ACE.Rph", "ACE.h2", "ACE.c2", "ACE.e2"))
  )
  
  # assign the model to the global environment - save as the model_name specified
  assign(model_name, IPM, envir = .GlobalEnv)
  
  # run the specified model
  IPM_fit <- mxTryHard(IPM_anxiety, intervals = TRUE)
  
  # return the summary to see how the model ran
  return(summary(IPM_fit, verbose = TRUE)) 
}


########### IPM function to run ACE model at two time points - including scalar means ############

# Here, two extra bits of input are needed. 
# * sex_diff_list is a c() list containing either TRUE or FALSE depending on which variables you want the scalar to be applied to. The order of this will be the same as the selvars order you have already provided. e
# .g. if I wanted scalar applied to anxiety12 and depression18 in my example for both twins:
# selvars_chosen = c(anxiety12_twin1, depression12_twin1, anxiety18_twin1, depression18_twin1, 
#                    anxiety12_twin2, depression12_twin2, anxiety18_twin2, depression18_twin2) and 
# sex_diff_list = c(TRUE, FALSE, FALSE, TRUE, 
#                   TRUE, FALSE, FALSE, TRUE)
#
# * bounds is a list of the scalar names for these. All should start with "sc.." and follow with a number. This indicates the position in your list. In the same example, I would now set: 
# bounds = c("sc1", "sc3") to indicate that I want bounds for the 1st and 3rd variable in the sex_diff_list (the ones that are TRUE). This will be the same for both twins. 
# we will also need to specify separate data frames for the male and female MZ and DZ data. 

IPM_ACE_2time_sexmean <- function(dataMZm, dataMZf, dataDZ, dataDZ, nv, selvars_chosen, model_name, sex_diff_list, bounds){ 
  
  # number of variables
  nv = nv    					      # number of variables - e.g. 2 at age 12 and 2 at age 18
  ntv <- nv*2 				      # number of twin variables
  nlower <- nv*(nv+1)/2 		# number of free elements in an nv*nv lower matrix 
  
  # number of factors A (common)
  nfAc <- 2 # change dimension of A factor matrix Ac to have *2* common factors
  
  # Set if the parameters are estimated (TRUE) or not (FALSE) 
  Ac2Free   <- c(rep(TRUE, nv), rep(FALSE, nv/2), rep(TRUE, nv/2))  
  
  # Create start values for 2 common ACE Factors - set the start values for common 
  Ac2Values <- c(rep(.1, nv),  rep(0, nv/2), rep(.1, nv/2))  
  Cc2Values <- c(rep(.1, nv),  rep(0, nv/2), rep(.1, nv/2))
  Ec2Values <- c(rep(.1, nv),  rep(0, nv/2), rep(.1, nv/2))
  
  IPM <- mxModel(model_name,
                 mxModel("ACE",
                         # Common ace paths
                         mxMatrix("Full", nrow = nv, ncol = nfAc, free = Ac2Free, values = Ac2Values, name = "ac"),
                         mxMatrix("Full", nrow = nv, ncol = nfAc, free = Ac2Free, values = Cc2Values, name = "cc"),
                         mxMatrix("Full", nrow = nv, ncol = nfAc, free = Ac2Free, values = Ec2Values, name = "ec"),
                         
                         # Specific ace paths
                         mxMatrix("Diag", nrow = nv, ncol = nv, free = TRUE, values = 0.3,  name = "as"),
                         mxMatrix("Diag", nrow = nv, ncol = nv, free = TRUE, values = 0.3, name = "cs"),
                         mxMatrix("Diag", nrow = nv, ncol = nv, free = TRUE, values = 0.3,  name = "es"),
                         
                         # ACE variance components
                         mxAlgebra(ac %*% t(ac) + as %*% t(as), name = "A"), 
                         mxAlgebra(cc %*% t(cc) + cs %*% t(cs), name = "C"),
                         mxAlgebra(ec %*% t(ec) + es %*% t(es), name = "E"),
                         mxAlgebra(A + C + E, name = "V"),
                         
                         # Calculate standard deviation for standardization
                         mxMatrix("Iden", nrow = nv, ncol = nv, name = "I"),
                         mxAlgebra(solve(sqrt(I*V)), name = "iSD"), # extracting the variances from the diagonal of the covariance matrix V 
                         
                         # Expected means vector - one grand mean
                         mxMatrix("Full", nrow = 1, ncol = ntv, free = TRUE, values = 0.01, name = "expMean"), 
                         
                         # Genetic and environmental correlations (correlated factor solution)
                         mxAlgebra(solve(sqrt(I*V)) %&% V, name = "Rph"),
                         mxAlgebra(solve(sqrt(I*A)) %&% A, name = "Ra"),
                         mxAlgebra(solve(sqrt(I*C)) %&% C, name = "Rc"),
                         mxAlgebra(solve(sqrt(I*E)) %&% E, name = "Re"),
                         
                         # Standardised components of variance
                         mxAlgebra(A/V, name = "h2"),
                         mxAlgebra(C/V, name = "c2"),
                         mxAlgebra(E/V, name = "e2"),
                         
                         # Standardise common parameters (unsquared)
                         mxAlgebra(iSD %*% ac, name = "stac"),
                         mxAlgebra(iSD %*% cc, name = "stcc"),
                         mxAlgebra(iSD %*% ec, name = "stec"),
                         
                         # Standardise common parameters (Squared)
                         mxAlgebra(stac*stac, name = "stac2"),
                         mxAlgebra(stcc*stcc, name = "stcc2"),
                         mxAlgebra(stec*stec, name = "stec2"),
                         
                         # Standardise specific parameters (unsquared)
                         mxAlgebra(iSD %*% as, name = "stas"),
                         mxAlgebra(iSD %*% cs, name = "stcs"),
                         mxAlgebra(iSD %*% es, name = "stes"),
                         
                         # Standardise specific parameters (Squared)
                         mxAlgebra(stas*stas, name = "stas2"),
                         mxAlgebra(stcs*stcs, name = "stcs2"),
                         mxAlgebra(stes*stes, name = "stes2"),	
                         
                         # Scalar multiplier 
                         mxMatrix("Diag", nrow = ntv, ncol = ntv, free = sex_diff_list,
                                  values = 1,
                                  label  = c("sc1", "sc2", "sc3", "sc4", "sc5", "sc6", 
                                             "sc1", "sc2", "sc3", "sc4", "sc5", "sc6"),
                                  name   = "Scalar"),
                         
                         mxBounds(bounds, 0, NA), # creates bounds for these objects when = TRUE
                         
                         # MZ expected covariance matrix
                         mxAlgebra(Scalar %&% rbind(cbind(A+C+E, A+C), cbind(A+C, A+C+E)), name = "expCovMZm"), 
                         mxAlgebra(rbind(cbind(A+C+E, A+C), cbind(A+C, A+C+E)), name = "expCovMZf"),
                         
                         # DZ expected covariance matrix
                         mxAlgebra(Scalar %&% rbind(cbind(A+C+E, (0.5%x%A)+C), cbind((0.5%x%A)+C, A+C+E)), name = "expCovDZm"), 
                         mxAlgebra(rbind(cbind(A+C+E, (0.5%x%A)+C), cbind((0.5%x%A)+C, A+C+E)), name = "expCovDZf")
                         
                 ),
                 
                 mxModel("MZM",
                         mxData(observed = dataMZm, type = "raw"),
                         mxExpectationNormal(covariance = "ACE.expCovMZm", means = "ACE.expMeanM", dimnames = selvars_chosen),
                         mxFitFunctionML()
                 ),
                 mxModel("MZF",
                         mxData(observed = dataMZf, type = "raw"),
                         mxExpectationNormal(covariance = "ACE.expCovMZf", means = "ACE.expMeanF", dimnames = selvars_chosen),
                         mxFitFunctionML()
                 ),
                 mxModel("DZM",
                         mxData(observed = dataDZm, type = "raw"),
                         mxExpectationNormal(covariance = "ACE.expCovDZm", means = "ACE.expMeanM", dimnames = selvars_chosen),
                         mxFitFunctionML()
                 ),
                 mxModel("DZF",
                         mxData(observed = dataDZf, type = "raw"),
                         mxExpectationNormal(covariance = "ACE.expCovDZf", means = "ACE.expMeanF", dimnames = selvars_chosen),
                         mxFitFunctionML()
                 ),
                 
                 mxAlgebra(MZM.objective + DZM.objective + MZF.objective + DZF.objective, name = "m2LL"),
                 mxFitFunctionAlgebra("m2LL"),
                 mxCI(c("ACE.stac2", "ACE.stcc2", "ACE.stec2", "ACE.stas2", "ACE.stcs2", "ACE.stes2", "ACE.Rph", "ACE.h2", "ACE.c2", "ACE.e2"))
  )
  
  # assign the model to the global environment - save as the model_name specified
  assign(model_name, IPM, envir = .GlobalEnv)
  
  # run the specified model
  IPM_fit <- mxTryHard(IPM_anxiety, intervals = TRUE)
  
  # return the summary to see how the model ran
  return(summary(IPM_fit, verbose = TRUE)) 
}



