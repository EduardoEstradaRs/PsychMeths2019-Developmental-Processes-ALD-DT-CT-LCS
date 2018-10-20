# R code for the statistical models used in

# Estrada & Ferrer (2018). Discrete- and Continuous-Time Latent Change Score Modeling of
# Developmental Processes in Accelerated Cohort-Sequential Designs.
# Submitted for publication in Psychological Methods (2018)

library(OpenMx)

# dEmp is a data frame containing an empirical data set in wide format
# It includes repeated observations of the process y
# In this example, y was observed from age 5 to age 19 (variable names: oy0 - oy14)
# Variables age0 to age14 include the exact time for each observation in each case
# Missing values are coded as NA


## Estimation of Latent Change Score in Discrete Time (LCS-DT) with OpenMx -----------
# Adapted from:
# Ghisletta & McArdle (2012). Latent Curve Models and Latent Change Score Models Estimated in R.
# Structural Equation Modeling: A Multidisciplinary Journal, 19(4), 651–682.
# https://doi.org/10.1080/10705511.2012.713275

SEMmanifest <- paste0("oy", seq(0,14))
SEMlatents <- paste0("y", seq(0,14))
SEMdiffs <- paste0("d", seq(1,14))

LCSMx <- mxModel(name = "mxLCS_SEM", mxData(observed = dEmp, type="raw"),
                 type="RAM",
                 manifestVars = SEMmanifest,
                 latentVars = c("yInt", "ySlp", SEMlatents, SEMdiffs),
                 mxPath(from=SEMlatents, to=SEMmanifest,
                        arrows=1, free=FALSE, values=1),
                 mxPath(from=SEMlatents[1:14], to=SEMlatents[2:15],
                        arrows=1, free=FALSE, values=1),
                 mxPath(from=SEMdiffs, to=SEMlatents[2:15],
                        arrows=1, free=FALSE, values=1),
                 mxPath(from=SEMlatents[1:14], to=SEMdiffs,
                        arrows=1, free=TRUE, values=-.2, labels="b_y"),
                 mxPath(from="yInt", to="y0", arrows=1, free=FALSE, values=1),
                 mxPath(from="ySlp", to=SEMdiffs, arrows=1, free=FALSE, values=1),
                 mxPath(from="one", to=c("yInt","ySlp"), free=TRUE, values=c(12,7),
                        labels=c("yInMn","ySlMn")),
                 mxPath(from=c("yInt","ySlp",SEMmanifest),
                        arrows=2, free=TRUE, values=c(25,.7,rep(2, length(SEMmanifest) )),
                        labels=c("yInV","ySlV",rep("MerY", length(SEMmanifest) ))),
                 mxPath(from="yInt", to="ySlp", arrows=2, free=TRUE, values=3, labels="yInSlCv") )

LCS_DT <- mxRun(LCSMx)
#summary(LCS_DT)



## Estimation of State-Space Model in Continuous Time (SSM-CT) with OpenMx -----------
# Adapted from:
# Hunter (2018). State Space Modeling in an Open Source, Modular, SEM Environment.
# Structural Equation Modeling: A Multidisciplinary Journal, 25(2), 307–324.
# https://doi.org/10.1080/10705511.2017.1369354

## Restructure data to lists
cases <- seq(1:nrow(dEmp))
data2List <- function(x, dwork) {
  obsI <- unlist(dwork[x,paste0("oy",seq(0,14))])
  ageI <- unlist(dwork[x,paste0("age",seq(0,14))])
  dLI <- data.frame(y = obsI, age = ageI, ID = x, time = 0:14)
  dLI <- dLI[complete.cases(dLI),]
}
dataL <- lapply(cases, data2List, dwork = dEmp)

modNames_ct <- paste0("i", cases, "ODE")

## Specify invariant part of the OpenMx objects ##
## For more details on the meaning of each matrix in the SSM model,
## see Equations 4 and 5 in the main paper, and Equations 8 and 9 in Hunter (2018).

opmxL <- list()

# Autoregressive dynamics
opmxL$amat_ct <- mxMatrix(name = "A", "Full", 2, 2, free = c(T,F,F,F),
                          values = c(-.2,0,1,0),
                          dimnames = list( c("yIn", "ySl"), c("yIn", "ySl") ),
                          labels = c("b_y", NA,NA,NA),
                          lbound = c(-1, NA,NA,NA),
                          ubound = c(0, NA,NA,NA))


opmxL$bmat <- mxMatrix(name = "B", "Zero", 2, 1) # Input effects on the latent variables

# Factor loadings in the measurement model
opmxL$cmat <- mxMatrix(name = "C", "Full", 1, 2, free = FALSE,
                       values = c(1,0), 
                       dimnames = list( c("y"), c("yIn", "ySl") ),
                       labels = c(NA, NA)  )

opmxL$dmat <- mxMatrix("Zero", 1, 1, name = "D") # Input effects on the observed variables

opmxL$qmat <- mxMatrix("Zero", 2, 2, name = "Q") # Dynamic error (i.e., innovations)

opmxL$rmat <- mxMatrix("Diag", 1, 1, TRUE, 2,    # Measurement error
                       name = "R", labels = "MerY")

# Mean vector of latent variables at t=0
opmxL$xmat <- mxMatrix(name = "x0", "Full", 2, 1, free = TRUE,
                       values = c(12, 7),
                       labels = c("yInMn", "ySlMn"))

# Covariance matrix of latent variables at t=0
opmxL$pmat <- mxMatrix(name = "P0", "Symm", 2, 2, TRUE,
                       values = c(25, 3, .7),
                       labels = c("yInV", "yInSlCv", "ySlV"),
                       lbound = c(0, NA, 0))

opmxL$umat <- mxMatrix("Zero", 1, 1, name = "u") # Covariates

# Specification of the time index
opmxL$tmat <- mxMatrix('Full', 1, 1, name='time', labels='data.age')

opmxL$modL_ct <- with(opmxL, list(amat_ct, bmat, cmat, dmat,
                                  qmat, rmat, xmat, pmat,
                                  umat, tmat))

opmxL$expODE <- mxExpectationStateSpaceContinuousTime(A = "A", B = "B",
                                                      C = "C", D = "D",
                                                      Q = "Q", R = "R",
                                                      x0 = "x0", P0 = "P0",
                                                      u = "u", t = "time")

## Create multisubject model ##
genMxIndModels_ct <- function(x, dwork, modNames_ct) {
  DataSetForSubjectK <- dwork[[x]]
  indivmodels <- mxModel(name = modNames_ct[x],
                         opmxL$modL_ct,
                         opmxL$expODE,
                         mxFitFunctionML(),
                         mxData(DataSetForSubjectK, type ='raw')  )  }

indivmodels_ct <- lapply(cases, genMxIndModels_ct, dataL, modNames_ct)


## SSM Estimation ##
multiSubjODE <- mxModel(name = "MultiODE_CT", indivmodels_ct,
                        mxFitFunctionMultigroup(modNames_ct))
SSM_CT <- mxRun(multiSubjODE)  

#summary(SSM_CT)





