# @param beta parameter vektor (k x 1)
# @param z imputed observed variables before truncation (n x 1)
# @param W spatial weight matrix (n x n)
# @param X design matrix (n x k)
density_rho_q <- function(rho_q, q, rho, beta, z, I_n, I_nq, W, X, R, type=c("SAR","SEM")) {
  rho[q,q] <- rho_q
  A <- (I_nq - kronecker(rho , W) )         # nq x nq
  if (type == "SAR") {
    S <- A %*% z - X %*% beta     # nq x 1; Residuals of SAR
  } else {
    S <- as.numeric(A %*% (z - X %*% beta))   # nq x 1; Residuals of SEM
  }
  
  val <- as.numeric(det(A) * exp(-0.5 * t(S) %*% kronecker(solve(R), I_n) %*% S))
  #as.vector(exp(log(det(A)) - 0.5 * t(S) %*% S))
  if(val==0 || is.nan(val) || is.na(val)) {
    val <- 1.e-16 
  }
  return(val)
}

ludensity_rho_q <- function(rho_q, a, q, rho, beta, z, I_n, I_nq, KRI_N, W, X, R, type=c("SAR","SEM"), trW_i, o) {
  if( rho_q < -0.9   || rho_q > 0.9){
    val <- - 1.e12
  }else{
    rho[q,q] <- rho_q
    #A <- (I_nq - kronecker(rho , W) )         # nq x nq
    #Aq <- (I_n - rho_q * W)
    if (type == "SAR") {
     Z <- matrix(z, byrow = FALSE, ncol= nrow(rho))
     S <- as.double(I_nq%*%z - vec(W %*% Z%*% t(rho)) - X %*% beta )
     #S <- A %*% z - X %*% beta     # nq x 1; Residuals of SAR
    } else {
      Z <- matrix((z - X %*% beta), byrow = FALSE, ncol= nrow(rho))
      S <- as.numeric(I_nq%*%(z - X %*% beta) - vec(W %*% Z%*% t(rho)) )
      #S <- as.numeric(A %*% (z - X %*% beta))   # nq x 1; Residuals of SEM
    }
    lprior <- log(beta_prior(rho_q,a,a))
    
    #ldet_A <- as.double(determinant(A, logarithm = TRUE)$modulus)
    #ldet_A <- as.double(determinant(A, logarithm = TRUE)$modulus) 
    trrho_i <- sapply(1:o, function(t) sum(diag(rho^o)))
    ldet_A <- -sum(trrho_i*trW_i/(1:o))
    val <- as.numeric( ldet_A - 0.5 * t(S) %*% KRI_N %*% S + lprior)
    #as.vector(exp(log(det(A)) - 0.5 * t(S) %*% S))
    if(val==0 || is.nan(val) || is.na(val)) {
      val <- -1.e12 
    }
  }
 
  return(val)
}
ludensity_rho <- function(param, a, beta, z, I_n, I_nq, KRI_N, W, X, R, type=c("SAR","SEM")) {
  if( any(param < -0.99   || param > 0.99)){
    val <- - 1.e12
  }else{
    rho <- diag(param)
    A <- (I_nq - kronecker(rho , W) )         # nq x nq
    #Aq <- (I_n - rho_q * W)
    if (type == "SAR") {
      Z <- matrix(z, byrow = FALSE, ncol= nrow(rho))
      S <- as.double(I_nq%*%z - vec(W %*% Z%*% t(rho)) - X %*% beta )
      #S <- A %*% z - X %*% beta     # nq x 1; Residuals of SAR
    } else {
      Z <- matrix((z - X %*% beta), byrow = FALSE, ncol= nrow(rho))
      S <- as.numeric(I_nq%*%(z - X %*% beta) - vec(W %*% Z%*% t(rho)) )
      #S <- as.numeric(A %*% (z - X %*% beta))   # nq x 1; Residuals of SEM
    }
    lprior <- sum(log(beta_prior(param,a,a)))
    
    #ldet_A <- as.double(determinant(A, logarithm = TRUE)$modulus)
    ldet_A <- as.double(determinant(A, logarithm = TRUE)$modulus) 
    val <- as.numeric( ldet_A - 0.5 * t(S) %*% KRI_N %*% S + lprior)
    #as.vector(exp(log(det(A)) - 0.5 * t(S) %*% S))
    if(val==0 || is.nan(val) || is.na(val)) {
      val <- -1.e12 
    }
  }
  
  return(val)
}


ludensity_rho_all <- function(param, a, beta, z, I_n, I_nq, KRI_N, W, X, R, type=c("SAR","SEM")) {
  if( any(param < -0.99   || param > 0.99)){
    val <- - 1.e12
  }else{
    rho <- invvec(param)
    A <- (I_nq - kronecker(rho , W) )         # nq x nq
    #Aq <- (I_n - rho_q * W)
    if (type == "SAR") {
      Z <- matrix(z, byrow = FALSE, ncol= nrow(rho))
      S <- as.double(I_nq%*%z - vec(W %*% Z%*% t(rho)) - X %*% beta )
      #S <- A %*% z - X %*% beta     # nq x 1; Residuals of SAR
    } else {
      Z <- matrix((z - X %*% beta), byrow = FALSE, ncol= nrow(rho))
      S <- as.numeric(I_nq%*%(z - X %*% beta) - vec(W %*% Z%*% t(rho)) )
      #S <- as.numeric(A %*% (z - X %*% beta))   # nq x 1; Residuals of SEM
    }
    lprior <- log(beta_prior(param,a,a))
    
    #ldet_A <- as.double(determinant(A, logarithm = TRUE)$modulus)
    ldet_A <- as.double(determinant(A, logarithm = TRUE)$modulus) 
    val <- as.numeric( ldet_A - 0.5 * t(S) %*% KRI_N %*% S + lprior)
    #as.vector(exp(log(det(A)) - 0.5 * t(S) %*% S))
    if(val==0 || is.nan(val) || is.na(val)) {
      val <- -1.e12 
    }
  }
  
  return(val)
}

ludensity_rho_q_easy <- function(a, rho_q, q, o, rho, trW_i , z_z, z_tXbeta, Xbeta2, I_nq, KRI_N, R, W) {
  if( rho_q < -0.9   || rho_q > 0.9){
    val <- - 1.e12
  }else{
    rho[q,q] <- rho_q
    A      <- (I_nq - kronecker(rho , W) )         # nq x nq
    e0     <- - 0.5 * sum(diag(solv(R) %*% (A %*% z_z %*% t(A) - A %*% z_tXbeta - t(A %*% z_tXbeta) + Xbeta2)))
    lprior <- log(beta_prior(rho_q,a,a))
    trrho_i <- sapply(1:o, function(t) sum(diag(rho^o)))
    ldet_A <- -sum(trrho_i*trW_i/(1:o))
    val <- as.numeric( ldet_A + e0 + lprior)
    if(val==0 || is.nan(val) || is.na(val)) {
      val <- -1.e12 
    }
  }
  
  return(val)
}

ludensity_rho_q_easy2 <- function(rho_q, a, q, rho, KRI_N, trW_i, o, WZ, I_nq_zXbeta) {
  rho[q,q] <- rho_q
  WZrho    <- tcrossprod(WZ, rho)
  S        <- as.double(I_nq_zXbeta - vec(WZrho) )
  lprior   <- log(beta_prior(rho_q,a,a))
  #trrho_i  <- sapply(1:o, function(t) sum(diag(rho^t)))
  trrho_i  <- rep(0, o)
  v <- diag(1,nrow(rho))
  for(i in 1:o){
    v <- rho%*%v
    trrho_i[i] <- sum(diag(v))
  }
  ldet_A   <- -sum(trrho_i*trW_i/(1:o))
  val      <- as.numeric( ldet_A - 0.5 * t(S) %*% KRI_N %*% S + lprior)
  if(val==0 || is.nan(val) || is.na(val)) {
    val <- -1.e12 
  }
  return(val)
}

ludensity_rho_q_easy3 <- function(rho_q, a, q, o, rho, trW_i , z, Xbeta, I_nq, KRI_N, W) {
  rho[q,q] <- rho_q
  A        <- (I_nq - kronecker(rho , W) )         # nq x nq
  e0       <- A %*% z - Xbeta
  epe      <- -0.5* t(e0) %*% KRI_N %*% e0
  
  lprior   <- log(beta_prior(rho_q,a,a))
  #trrho_i  <- sapply(1:o, function(t) sum(diag(rho^t)))
  trrho_i  <- rep(0, o)
  v        <- diag(1,nrow(rho))
  for(i in 1:o){
    v <- rho%*%v
    trrho_i[i] <- sum(diag(v))
  }
  
  ldet_A   <- -sum(trrho_i*trW_i/(1:o))
  val      <- as.numeric( ldet_A + epe + lprior)
  if(val==0 || is.nan(val) || is.na(val)) {
    val <- -1.e12 
  }
  
  return(val)
}
ludensity_rho_q_mcem <- function(rho_q, q, o, rho, KRI_N, trW_i, W, zpz, zpXb, XbpXb) {
  if( rho_q < -1   || rho_q > 1){
    val <- - 1.e12
  }else{
    rho[q,q] <- rho_q
    #trrho_i  <- sapply(1:o, function(t) sum(diag(rho^t)))
    trrho_i  <- rep(0, o)
    v        <- diag(1,nrow(rho))
    for(i in 1:o){
      v <- rho%*%v
      trrho_i[i] <- sum(diag(v))
    }
    ldet_A   <- -sum(trrho_i*trW_i/(1:o))
    S <- I_nq - kronecker(rho, W)
    SzpXb   <- S  %*% zpXb
    temp    <- S %*% tcrossprod(zpz, S)  - SzpXb - t(SzpXb) + XbpXb
    ltemp <- - 0.5*sum(diag(KRI_N %*% temp))
    val      <- as.numeric( ldet_A + ltemp)
    if(val==0 || is.nan(val) || is.na(val)) {
      val <- -1.e12 
    }
  }
  return(val)
}

ludensity_rho_q_mcem2 <- function(rho_q, q, o, rho, KRI_N, trW_i, W, z, Xbeta) {
  if( rho_q < -1   || rho_q > 1){
    val <- - 1.e12
  }else{
    rho[q,q] <- rho_q
    #trrho_i  <- sapply(1:o, function(t) sum(diag(rho^t)))
    trrho_i  <- rep(0, o)
    v        <- diag(1,nrow(rho))
    for(i in 1:o){
      v <- rho%*%v
      trrho_i[i] <- sum(diag(v))
    }
    ldet_A   <- -sum(trrho_i*trW_i/(1:o))
    S <- I_nq - kronecker(rho, W)
    Sz <- as.double(S %*% z )
    ltemp <- dmvnorm(Sz, mean = as.double(Xbeta), sigma = as.matrix(KRI_N), log = TRUE)
    val      <- as.numeric( ldet_A + ltemp)
    if(val==0 || is.nan(val) || is.na(val)) {
      val <- -1.e12 
    }
  }
  return(val)
}

# Metropolis-Hastings-Chain with tuned acceptance rate, see LeSage (2009)
# 
# proposal density is normal distribution
# g(rho* | rho_t) = rho_t + c * N(0, 1)
#
# @param type "SEM" or "SAR",
# @param n number of samples to draw
# @param start.value start value for rho
# @param c tuning parameter
draw_rho_q_metropolis <- function(type="SEM", a, n , q, rho, beta, z, I_n, I_nq, KRI_N, W, X, R, burn.in=100, start.value, c=1) {
  n <- n + burn.in
  u <- runif(n=n, 0, 1)  # samples for M-H-ratio
  s <- rnorm(n=n, 0, 1)  # realisations from proposal density
  rho_t <- numeric(n)    # vector for rho within the chain
  rho_t[1] <- start.value              # start value
  p_t      <- ludensity_rho_q(a, rho_t[1],q , rho, beta=beta, z=z, I_n=I_n, I_nq=I_nq, KRI_N=KRI_N, W=W, X=X, R=R, type=type)  # f(rho_t | beta, z)
    
  i    <- 2 # number of accepted proposals / length of chain
  acceptance_rate <- numeric(n) # running acceptance rate
  num_accepted <- 0
  while (i <= n) {
    
    # create proposal rho_p from proposal density g(rho_p | rho_{t-1}) ~ N(rho_{t-1}, c^2)
    rho_p <- rho_t[i-1] + c * s[i]   # proposal rho = letztes akzeptiertes Rho + normalverteilte Zufallsstreuung s
   # rho_p <- rtruncnorm(n=1, a=-0.99, b=0.99, mean = rho_t[i-1], sd = c)
     if( rho_p < -1 || rho_p > 1)
       rho_p <- rho_t[i-1]
    
    # Berechnung der Dichte f(rho_p) für das Proposal rho_p
    p_p <- ludensity_rho_q(a, rho_p, q , rho, beta=beta, z=z, I_n=I_n, I_nq=I_nq, KRI_N=KRI_N, W=W, X=X, R, type=type)       
    # SW: Berechnung der Dichte ist teuer, daher besser für ein Grid vorher rechnen?
    # Jain, Dichte p(rho|t,beta( hängt von z, beta ab. Aber zumindestens die Log-Determinante kann man vorher berechnen
    # Berechnung Dichte f(rho_{t-1}) für letzten Wert der Chain
    # p_t <- p(rho_t[i-1], beta=beta, z=z, W=W, X=X, type=type)  # SW: Kann ich die alte Dichte nicht merken, dann muss ich die nicht neu berechnen
    
    # Wegen Symmetrie der Normalverteilung als Proposal-Dichte g(rho_p | rho_t) = g(rho_t | rho_p)
    # vereinfacht sich das M-H-Ratio von [ f(rho_p) * g(rho_t | rho_p) ] / [ f(rho_t) * g(rho_p | rho_t) ]
    # auf f(rho_p) / f(rho_t)!
    # determine M-H-Ratio R(rho_p, rho_t) = min[ 1 ; p(rho_p) / p(rho_t) ]
    # see LeSage (2009), eqn (5.27)
    Rm <- min(1, exp(p_p - p_t))
    if (u[i] <= Rm) {    # see Givens/Hoeting, p.184, eqn (7.2)
      # accept proposal
      rho_t[i] <- rho_p
      p_t <- p_p  # save density
      num_accepted <- num_accepted + 1
    } else {
      # stay with the current value
      rho_t[i] <- rho_t[i-1]
    }
    acceptance_rate[i] <- num_accepted / i
    if (acceptance_rate[i] < 0.2) c <- c / 1.1   # Wenn Akzeptanzrate zu klein, dann verringere (?) Streuungsparameter "c"
    if (acceptance_rate[i] > 0.4) c <- c * 1.1
    i <- i + 1
  }
  return(list(rho_t=rho_t,acceptance_rate=acceptance_rate, c=c)) 
}