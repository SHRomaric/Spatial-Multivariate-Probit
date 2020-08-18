################################################################################
#
# Bayesian estimation of Multivariate spatial autoregressive probit model (MSAR probit)
# using MCMC sampling
# 
# using code from Miguel Godinho de Matos <miguelgodinhomatos@cmu.edu> and Stefan Wilhelm <Stefan.Wilhelm@financial.com>
#
################################################################################

# estimated tr(W^i) for i=1,...,100
# see LeSage (2009), chapter 4, p.97/98
#
# "The expectation of the quadratic form u'Au equals tr(A).
#  Since u_i^2 follows a chi^2 distribution with one degree of freedom."
# Pre-calculate traces tr(W^i) i=1,...,100 for the x-impacts calculations
#
# @param W spatial weight matrix (n x n)
# @param o highest order of W^i = W^o
# @param iiter number of MCMC iterations (we run this 50 times)
# @return (n x o) matrix with tr(W^i) in each column, for i=1..,o

tracesWi <- function(W, o=100, iiter=50) {
  n     <- nrow(W)
  trW_i <- matrix( data=0, nrow=n, ncol=o )   # n x o
  u     <- t(mvrnorm(n=iiter, mu = rep(0,n), Sigma = diag(1,n)))#matrix(rnorm(n * iiter), nrow = n, ncol = iiter)   # (n x iiter)
  xx    <- u
  trW_i[,1] <- apply(u * as.matrix(xx), 1, sum)    # tr(W^0) = 1
  for(i in 2:o ){
    xx <- W %*% xx  # (n x iter)
    trW_i[,i] <- apply(u * as.matrix(xx), 1, sum)  # u'(W^i)u; sum across all iterations
  }
  trW_i <- colSums(trW_i / iiter)
  return(trW_i)
  
}
tracesWi_P <- function(W, o=100, iiter=50) {
  n     <- nrow(W)
  trW_i <- matrix( data=0, nrow=o, ncol=iiter )   # n x o
  u     <- matrix(rnorm(n * iiter), nrow = n, ncol = iiter)   # (n x iiter)
  xx    <- u
  #trW_i[,1] <- apply(u * as.matrix(xx), 1, sum)    # tr(W^0) = 1
  for(i in 1:o ){
    xx <- W %*% xx  # (n x iter)
    trW_i[i,] <- n*apply(u * as.matrix(xx), 2, sum) /apply(u * u, 2, sum)  # u'(W^i)u; sum across all iterations
  }
  trW_i <- rowMeans(trW_i)
  return(trW_i)
  
}
# PURPOSE: draw rho from conditional distribution p(rho | beta, z, y)
# ---------------------------------------------------
#  USAGE: 
#
#  results.rmin 
#  results.rmax maximum eigen value
#  results.time execution time
# ---------------------------------------------------
#where epe0 <- t(e0) %*% e0
#      eped <- t(ed) %*% ed
#      epe0d<- t(ed) %*% e0
#
# detval1 umbenannt in => rho_grid = grid/vector of rhos
# detval2 umbenannt in => lndet = vector of log-determinants
# detval1sq = rvec^2
# lnbprior = vector of log prior densities for rho (default: Beta(1,1))
#
draw_rho <- function (rho_grid, lndet, rho_gridsq, yy, epe0, eped, epe0d, 
                      rho, nmk, nrho, lnbprior, u) 
{
  # This is the scalar concentrated log-likelihood function value
  # lnL(rho). See eqn(3.7), p.48
  z <- epe0 - 2 * rho_grid * epe0d + rho_gridsq * eped
  z <- -nmk * log(z)
  den <- lndet + z + lnbprior          # vector of log posterior densities for rho vector
  # posterior density post(rho | data) \propto likelihood(data|rho,beta,z) * prior(rho)
  # log posterior density post(rho | data) \propto loglik(data|rho,beta,z) + log prior(rho)
  n <- nrho
  adj <- max(den)
  den <- den - adj                     # adjustieren/normieren der log density auf maximum 0; 
  x <- exp(den)                        # density von p(rho) --> pdf
  isum <- sum(yy * (x[2:n] - x[1:(n - 1)])/2)
  z <- abs(x/isum)
  den <- cumsum(z)
  rnd <- u * sum(z)
  ind <- which(den <= rnd)
  idraw <- max(ind)
  if (idraw > 0 && idraw < nrho) {
    results <- rho_grid[idraw]
  }
  else {
    results <- rho
  }
  return(results)
}

# faster update of matrix S = (I - kron(rho,W)) for new values of rho
#
# @param S template matrix of (I - kron(rho,W))
# @param ind indizes to replaced
# @param W spatial weights matrix W
# @return (I - kron(rho,W))
update_I_rW <- function(S, ind, rho, W) {
  S@x[ind] <- (-kronecker(rho,W))@x
  return(S)
}


# Bayesian estimation of Multivariate SAR probit model
#
# @param formula 


msar_probit_mcmc <- function (Y, X, W, ndraw = 1000, burn.in = 5000, thinning = 1, 
                        prior = list(a1 = 2, a2 = 2, c = rep(0, ncol(X)),    T = diag(ncol(X)) * 1e+12, lflag = 0),
                        start = list(rho = diag(0.75, ncol(Y)),beta = rep(0, ncol(X)), R = diag(ncol(Y)), z),
                        scale = rep(1, ncol(Y)), m = 1, burn.in_MH=5, ndraw_MH=1, o=100, computeMarginalEffects = FALSE, showProgress = TRUE, pflag=TRUE, mflag=0, maxit=20) 
{
  timet <- Sys.time()
  y  <- as.matrix(c(Y))
  q  <- ncol(Y)
  n  <- nrow(Y)
  n1 <- nrow(X)
  n2 <- nrow(W)
  k  <- ncol(X)
  nq <- n*q
  
  I_n  <- sparseMatrix(i = 1:n, j = 1:n, x = 1)
  I_nq <- sparseMatrix(i = 1:nq, j = 1:nq, x = 1)
  if (is.null(colnames(X))) 
    colnames(X) <- paste("x", 1:k, sep = "")
  if (length(c(which(y == 0), which(y == 1))) != length(y)) {
    stop("msarprobit: not all y-values are 0 or 1")
  }
  if (n1 != q*n2 && n1 != n) {
    stop("msarprobit: wrong size of spatial weight matrix W")
  }
  if (!inherits(W, "sparseMatrix") || any(diag(W) != 
                                          0)) {
    stop("msarprobit: spatial weights matrix W must be a sparse matrix with zeros in the main diagonal")
  }
  ind <- match(n, apply(X, 2, sum))
  if (is.na(ind)) {
    cflag <- 0
    p <- k
  } else if (ind == 1) {
    cflag <- 1
    p <- k - 1
  } else {
    stop("msarprobit: intercept term must be in first column of the X-matrix")
  }
  rho   <- start$rho
  beta  <- start$beta
  R     <- start$R
  c     <- rep(0, k)
  if (is.numeric(prior$c) && length(prior$c) == k) {
    c <- prior$c
  }
  if (is.matrix(prior$T) && ncol(prior$T) == k && isSymmetric(prior$T) && 
      det(prior$T) > 0) {
    T <- prior$T
  }  else {
    T <- diag(k) * 1e+12
  }
  Tinv <- solve(T)
  if (class(W) == "dgCMatrix") {
    I <- sparseMatrix(i = 1:nq, j = 1:nq, x = Inf)
    S <- (I - kronecker(rho , W))
    ind  <- which(is.infinite(S@x))
    ind2 <- which(!is.infinite(S@x))
    S@x[ind] <- 1
  } else {
    S <- I_nq - kronecker(rho , W)
  }
  Rinv      <- solve(R)
  kr_Rinv_I <- kronecker(Rinv, I_n)
  H         <- t(S) %*% kr_Rinv_I %*% S
  H         <- forceSymmetric(H)
  QR        <- qr(S)
  mu        <- solve(QR, X %*% beta)
  lower <- lapply(1:q, function(t) as.double(ifelse(Y[,t] > 0, 0, -Inf) ) )
  upper <- lapply(1:q, function(t) as.double(ifelse(Y[,t] > 0, Inf, 0)) )
  rmin  <- -1
  rmax  <- 1
  lflag <- 0
  if (is.numeric(prior$lflag) && lflag %in% c(0, 1, 2))
    lflag <- prior$lflag
  # tmp <- sar_lndet(lflag, W, rmin, rmax)
  # detval <- tmp$detval
  a1 <- 1
  a2 <- 1
  if (is.numeric(prior$a1)) 
    a1 <- prior$a1
  if (is.numeric(prior$a2)) 
    a2 <- prior$a2
  # lnbprior <- log(beta_prior(detval[, 1], a1, a2))
  # u <- runif(thinning * ndraw + burn.in)
  # nrho <- nrow(detval)
  # nmk <- (n - k)/2
  # rho_grid <- detval[, 1]
  # lndet <- detval[, 2]
  # rho_gridsq <- rho_grid * rho_grid
  # yy <- (rho_grid[2:nrho] + rho_grid[1:(nrho - 1)])
  B <- matrix(NA, ndraw, k + q)
  R_chain <- matrix(NA, ndraw, length(vech(R)))
  if (showProgress) {
    pb <- txtProgressBar(min = 0, max = (thinning * ndraw +
                                           burn.in), initial = 0, style = 3)
  }
  # Compute posterior variance of \beta
  tX  <- t(X)
  xpx <- tX %*% kr_Rinv_I %*% X
  xpx <- forceSymmetric(xpx)
  T2  <- solve(xpx + Tinv)
  # Simulation of \beta0
  betadraws <- rmvnorm(n = (burn.in + ndraw * thinning), mean = rep(0,k), sigma = diag(k))
  # Degree of freedom for R
  df_0 <- n+q+1
  # 
  direct   <- matrix(NA, ndraw, p)
  indirect <- matrix(NA, ndraw, p)
  total    <- matrix(NA, ndraw, p)
  
  zmean    <- rep(0, n)
  if (cflag == 0) {
    namesNonConstantParams <- colnames(X)
  }  else {
    namesNonConstantParams <- colnames(X)[-1]
  }
  colnames(total)    <- namesNonConstantParams
  colnames(direct)   <- namesNonConstantParams
  colnames(indirect) <- namesNonConstantParams
  if (computeMarginalEffects) {
    trW.i <- tracesWi(W, o = 100, iiter = 100)
  }
  z    <- rep(0, nq)
  ones <- rep(1, nq)
  trW_i <- trW(W, m=o, p=500)
  indice.v <- NULL
  Fin<- 0
  for(t in 1:q) {
    Debut <- Fin + 1
    Fin  <- t*n
    indice.v[[t]] <- Debut:(Fin)
  }
  
  
  for (i in (1 - burn.in):(ndraw * thinning)) {
    # Draw latent variable Z from its full conditional: z | \theta, \lambda, W, y, X using Gibbs
    S     <- I_nq - kronecker(rho , W) #   #update_I_rW(S, ind=ind2, rho, W)
    KRI_N <- kronecker(solve(R), I_n)
    H     <- t(S) %*% KRI_N %*% S
    QR    <- qr(S)
    mu    <- as.double(solve(QR, X %*% beta))
    sim_z <- TRUE  # check if z  is well simulated
    nbsim <- 0
    while( sim_z== TRUE && nbsim <= 10){
      if (m == 1) {
        # ztemp <- rtmvnorm.sparseMatrix(n = 1, mean = mu, 
        #                                H = H, lower = lower, upper = upper, burn.in.samples = m, 
        #                                start.value = z)
        
        for (t  in 1:q) {
          zt    <- z[indice.v[[t]]]
          z_t   <- z[-indice.v[[t]]]
          mut   <- mu[indice.v[[t]]]
          mu_t  <- mu[-indice.v[[t]]]
          Htt   <- H[indice.v[[t]],indice.v[[t]]]
          Ht_t  <- H[,indice.v[[t]]][-indice.v[[t]],]
          mut_t <- mut-solve(Htt) %*% t(Ht_t) %*% (z_t-mu_t)  
          ztemp <- rtmvnorm.sparseMatrix(n = 1, mean = mut_t, 
                                         H = Htt, lower = lower[[t]], upper = upper[[t]], burn.in.samples = m, 
                                         start.value = zt)
          z[indice.v[[t]]] <- ztemp
        }
        ztemp <- z
      }else {
        ztemp <- rtmvnorm.sparseMatrix(n = 1, mean = mu, 
                                       H = H, lower = lower, upper = upper, burn.in.samples = m)
      }
      if(is.nan(sum(ztemp))){
        sim_z <- TRUE
        #m     <- m+1
        nbsim <- nbsim + 1
        z    <- rep(0, nq)
      }else{
        sim_z <- FALSE
      }
    }
    z <- as.double(ztemp)
    # Draw variable \omega from the full conditional: \
    Sz     <- S %*% z
    Xbeta  <- X %*% beta
    e0     <- as.double(Sz - Xbeta)
    e0M    <- matrix(e0, byrow=FALSE, ncol = q)
    RS     <- t(e0M) %*% e0M
    # RSinv  <- solve(RS)
    # Rminv   <- rWishart(n=1, df=df_0, Sigma=RSinv)
    # R       <- solve(apply(Rminv, 1:2, mean))
    R <- RS/df_0
    
    
    # Draw variable \beta from its full conditional: \beta | z, X
    KRI_N <- kronecker(solve(R), I_n)
    xpx   <- tX %*% KRI_N %*% X 
    xpx   <- forceSymmetric(xpx)
    T2    <- solve(xpx + Tinv)
    KrSz  <- KRI_N %*% Sz
    c2    <- T2 %*% (tX %*% KrSz + Tinv %*% c)
    
   # beta  <- as.double(c2 +  mvrnorm(1, mu=rep(0,k), Sigma=T2)  )#as.double(c2 + colMeans( mvrnorm(100, mu=rep(0,k), Sigma=T2) ) )
    beta  <- as.double(c2 +  t(chol(T2))%*%betadraws[i+ burn.in,]  )

    KRI_N  <- kronecker(solve(R), I_n)
    Z      <- matrix(z, byrow = FALSE, ncol= nrow(rho))
    WZ     <- W%*% Z
    Xbeta  <- X%*% beta
    I_nq_z <- I_nq%*%z
    I_nq_zXbeta <- I_nq_z - Xbeta
    if(mflag==0){
      for(r in 1:ncol(rho)){
        init_param  <- rho[r,r]
        rho_res <- optim(par=init_param, 
                         fn=ludensity_rho_q_easy3, a=a1, q=r, rho=rho, KRI_N=KRI_N, trW_i=trW_i, o=o , W=W,Xbeta=Xbeta,  I_nq=I_nq, z=z,
                         #fn=ludensity_rho_q_easy2, a=a1, q=r, rho=rho, KRI_N=KRI_N, trW_i=trW_i, o=o , WZ=WZ, I_nq_zXbeta=I_nq_zXbeta,
                         method = "Brent",
                         hessian = TRUE,
                         lower = -1,
                         upper = 1,
                         control = list(maxit=maxit, fnscale=-1))
        rho[r,r] <-  rtruncnorm(n=1, a=-0.99, b=0.99, mean = rho_res$par, sd = sqrt(-1/(rho_res$hessian)))
        #rho[r,r] <-  rho_res$par
      }
    } else if(mflag == 1){
      #Methode 2
      lu  <- function(param){
        ludensity_rho(param, a1, beta, z, I_n, I_nq, KRI_N, W, X, R, type= "SAR")
      }
      init_param <- diag(rho)
      rho_res <- optim(par=init_param,
                       fn=lu,
                       method = "BFGS",
                       hessian = TRUE,
                       control = list(maxit=maxit, fnscale=-1))
      #rho  <-  diag(as.double(rtmvnorm(n=1, mean =rho_res$par,  sigma= (-solve(rho_res$hessian)) , lower= rep(-0.99, q), upper= rep(0.99, q))))
      rho <-  diag(rho_res$par)
    } else{
      # Metropolis Hasting
      for(r in 1:ncol(rho)){
        init_param  <- rho[r,r]
        draw_rho    <- draw_rho_q_metropolis(type="SAR", a1, n=ndraw_MH , q=r, rho, beta, z, I_n, I_nq, KRI_N, W, X, R,burn.in=burn.in_MH, start.value=init_param, c=scale[r])
        #rho[r,r]    <- mean(draw_rho$rho_t[(burn.in_MH):(ndraw_MH+burn.in_MH)])
        rho[r,r]    <- draw_rho$rho_t[(ndraw_MH+burn.in_MH)]
        scale[r]    <- draw_rho$c
      }
    }
    
    print(diag(rho))
    print(vech(R))
    print(beta)
    # Scaling parameters and z
    R_12     <- sqrt(diag(diag(R)))
    R_12inv  <- solve( R_12)
    z        <- as.double(kronecker(R_12inv,I_n)%*% as.matrix(z))
    R        <- R_12inv %*% R %*% R_12inv
    beta     <- as.double(t(R_12inv %*% matrix(beta, byrow = TRUE, nrow = q)))
    rho      <- R_12inv %*% rho %*% R_12
    
    
    print(diag(rho))
    print(vech(R))
    print(beta)
    ##
    if (i > 0) {
      if (thinning == 1) {
        ind <- i
      } else if (i%%thinning == 0) {
        ind <- i%/%thinning
      } else {
        next
      }
      B[ind, ] <- c(beta, diag(rho))
      R_chain[ind,] <- as.double(vech(R))
      zmean <- zmean + z
      if (computeMarginalEffects) {
        o <- 100
        rhovec <- rho^(0:(o - 1))
        if (cflag == 1) {
          beff <- beta[-1]
        }
        else if (cflag == 0) {
          beff <- beta
        }
        pdfz <- dnorm(as.numeric(mu))
        dd <- sparseMatrix(i = 1:n, j = 1:n, x = pdfz)
        dir <- as.double(t(pdfz) %*% trW.i %*% rhovec/n)
        avg_direct <- dir * beff
        avg_total <- mean(dd %*% qr.coef(QR, ones)) * 
          beff
        avg_indirect <- avg_total - avg_direct
        total[ind, ] <- avg_total
        direct[ind, ] <- avg_direct
        indirect[ind, ] <- avg_indirect
      }
    }
    if (showProgress) 
      setTxtProgressBar(pb, i + burn.in)
  }
  if (showProgress) 
    close(pb)
  beta <- colMeans(B)[1:k]
  rho  <- diag(colMeans(B)[(k + 1): (k+q)])
  R    <- invvech(colMeans(R_chain))
  S    <- (I_nq - kronecker(rho, W) )
  fitted.values <- solve(qr(S), X %*% beta)
  fitted.response <- as.numeric(fitted.values >= 0)
  results <- NULL
  results$time <- Sys.time() - timet
  results$nobs <- n
  results$nvar <- k
  results$y <- y
  results$Y <- Y
  results$zip <- 2*n - sum(y)
  results$beta <- colMeans(B)[1:k]
  results$rho <- colMeans(B)[(k + 1): (k+q)]
  results$R <- R
  results$R_chain <- R_chain
  results$coefficients <- colMeans(B)
  results$fitted.values <- fitted.values
  results$fitted.response <- fitted.response
  results$ndraw <- ndraw
  results$nomit <- burn.in
  results$a1 <- a1
  results$a2 <- a2
  results$rmax <- rmax
  results$rmin <- rmin
  results$tflag <- "plevel"
  results$lflag <- lflag
  results$cflag <- cflag
  #results$lndet <- detval
  results$names <- c(colnames(X), "rho")
  results$B <- B
  results$bdraw <- B[, 1:k]
  results$pdraw <- B[, (k + 1): (k+q)]
  results$total <- total
  results$direct <- direct
  results$indirect <- indirect
  results$W <- W
  results$X <- X
  class(results) <- "msarprobit"
  return(results)
}