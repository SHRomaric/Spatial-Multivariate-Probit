# Simulations a Probit Bivariate-SAR process
rm(list=ls())
# Before starting, install package
#install.packages(c("RColorBrewer", "rgdal", "spdep", "readxl")) 
rm(list=ls())
# load the packages
require("rgdal")
require("spdep")
require("RColorBrewer") 
require("readxl")
require("ks")
require("Matrix")
require("MASS")
require("spatialreg")

# Library for sampling from Multivariate Normal distribution
require(mvtnorm)
# Library for sampling from Truncated Normal distribution
require(truncnorm)
# Library for sparse matrix
require(Matrix)
# Library for implementing the Gibbs sampler with the precision matrix  Hin package (Wilhelm and Manjunath, 2013)
require(tmvtnorm)
# Library for vec, vech, invvech functions
require(ks)
# 
require(Rcpp)
require(RcppArmadillo)
#
require(MCMCpack)
# Metropolis Hasting
require(mcmc)

# For reproducibility
set.seed(212)  
n <- 500
q <- 3
nq <- n*q

# Generate x data 
# design matrix with two standard normal variates as explanatory variables
x  <- cbind(intercept=1, x1=rnorm(n), x2=rnorm(n))
xl <- lapply(1:q, function(i) x)
X  <- bdiag(xl)
#X  <- as.matrix(bdiag(x,x,x))
# True values of regression coeffiecients beta
beta1 <- c(1, -2, 1.25) #c(2, -3, 1)
beta2 <- c(0,  -1, 0.5)
beta3 <- c(1,  0.5, -3)
Beta  <- rbind(beta1,beta2,beta3)

Beta0 <- as.matrix(c(beta1,beta2,beta3))
# True values of coeffiecients : Omega and Delta
R    <- invvech(c(1,0.5,0.6,1,0.0,1))
rho  <- diag(c(0.4,0.7,0.1))
# Simulation of bivariate latente variable
nb      <- knn2nb(knearneigh(cbind(x = rnorm(n), y = rnorm(n)), k = 10))
listw   <- nb2listw(nb, style = "W")
W       <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
I_nq    <- sparseMatrix(i = 1:nq, j = 1:nq, x = 1)
A       <- I_nq - kronecker(rho,W)
Epsilon <- mvrnorm(n, mu=rep(0,q), Sigma=R) 
Epsilon <- as.matrix(c(Epsilon))
Z0      <- X%*%Beta0 + Epsilon
Z       <- solve(qr(A), Z0)  
#system.time(Z1 <- solve(A, Z0))
Y <- as.double(Z>=0)
Y <- matrix(c(Y), byrow = FALSE, ncol=q)
colMeans(Y)
# Bayesian univariate estimation with 
require("spatialprobit")
sarprobit.fit1 <- sarprobit(Y[,1] ~ x - 1, W, ndraw = 2000, burn.in = 1000,   thinning = 1, m = 10, showProgress = TRUE)
summary(sarprobit.fit1)
# plot(sarprobit.fit1)
sarprobit.fit2 <- sarprobit(Y[,2] ~ x - 1, W, ndraw = 2000, burn.in = 1000,   thinning = 1, m = 10, showProgress = TRUE)
summary(sarprobit.fit2)
# plot(sarprobit.fit2)
sarprobit.fit3 <- sarprobit(Y[,3] ~ x - 1, W, ndraw = 2000, burn.in = 1000,   thinning = 1, m = 10, showProgress = TRUE)
summary(sarprobit.fit3)
# plot(sarprobit.fit3)

sourceCpp("code/MC_log_detCpp.cpp")
sourceCpp("code/MC_log_detCpp2.cpp")
source("code/autre_func.R")
source("code/msarprobit.R")
source("code/msarprobit_mceem.R")
source("code/msarprobit_Lawrence.R")
source("code/Draw_rho_q_metropolis.R")

beta  <- c(sarprobit.fit1$beta,sarprobit.fit2$beta,sarprobit.fit3$beta)
R     <- diag(1, q)
rho <- diag(c(sarprobit.fit1$rho,sarprobit.fit2$rho,sarprobit.fit3$rho))
z  <- c(as.vector(sarprobit.fit1$fitted.values),as.vector(sarprobit.fit2$fitted.values),as.vector(sarprobit.fit3$fitted.values))
library("profvis")
#profvis({
fit <- msar_probit_mcmc(Y, X, W, ndraw = 1000, burn.in = 2000, thinning = 1, 
                        prior = list(a1 = 1, a2 = 1, c = rep(0, ncol(X)),    T = diag(ncol(X)) * 1e+12, lflag = 0),
                        start = list(rho = rho, beta = beta, R = R, z=z),
                        scale=rep(1,ncol(Y)), m = 1, burn.in_MH=20, ndraw_MH=2, o=100, computeMarginalEffects = FALSE, showProgress = TRUE, pflag=TRUE,
                        mflag=0,maxit=100)
#})
fit <- msar_probit_mcmc_Lawr(Y=Y, x, W, ndraw = 1000, burn.in = 20000, thinning = 30, 
                             prior = list(a1 = 1, a2 = 1, c = rep(0, ncol(X)),    T = diag(ncol(X)) * 1e+12, lflag = 0),
                             start = list(rho = rho, beta = beta, R = R, z= z),
                             scale=rep(1,ncol(Y)), m = 1,burn.in_MH=20, ndraw_MH=1, o=100, computeMarginalEffects = FALSE, showProgress = TRUE, pflag=TRUE,
                             mflag=0, maxit=50)



fit <- msar_probit_mcmc(Y, X, W, ndraw = 100, burn.in = 400, thinning = 10, 
                        prior = list(a1 = 1.02, a2 = 1.02, c = rep(0, ncol(X)),    T = diag(ncol(X)) * 1e+12, lflag = 0),
                        start = list(rho = diag(0.75, ncol(Y)),beta = rep(0, ncol(X)), R = diag(ncol(Y)), z= rep(0,length(c(Y)))),
                        scale = rep(1,ncol(Y)), m = 1, burn.in_MH=10, ndraw_MH=1, o=25, computeMarginalEffects = FALSE, showProgress = TRUE) 

profvis({
fit <- msar_probit_mcem (Y, X, W, niter=10, burn.in = 50, ndraw = 50, 
                              start = list(rho = rho,beta = beta, R = R),
                              scale = rep(1, ncol(Y)), m = 1, o=100, showProgress = TRUE, pflag=TRUE, maxit=10) 
})

plot(x=1:3000, y=fit$B[,10], type="l")