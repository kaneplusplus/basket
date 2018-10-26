
#' @export
boa.hpd <- function(x, alpha) {
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]    
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

#' @export
sigmoid <- function(psi) { 
  1/(1+exp(-psi))
}

#' @export
logit <- function(a0) { 
  log(a0) - log(1-a0)
}

#' @export
MEM.mat <- function(Indices, mod.mat, H) {
  M <- matrix(NA, H, H )
  diag(M) <- rep(1,dim(M)[1])
  for(i in 1:(length(Indices)-1)) { 
    M[(i+1):dim(M)[2], i] <- M[i, (i+1):dim(M)[2]] <- mod.mat[[i]][Indices[i],]   }
  M[dim(M)[2], i+1] <- M[i+1, dim(M)[2]] <- (0:1)[ Indices[length(Indices)] ]
  M
}

#' @export
rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}

#' @export
logMarg.Dens <- function(I, mod.mat, xvec, nvec, avec, bvec) {

  M <- MEM.mat(I, mod.mat, length(xvec))
  marg.vec <- rep(NA, dim(M)[1])

  #calculate the product portion of integrated marginal likelihood
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) 

  for(i in 1:dim(M)[1]) { 
    p.vec <- prod( prod.vec^(1-M[i,]) )
    marg.vec[i] <- p.vec * 
      ( beta(avec[i] + M[i,] %*% xvec, bvec[i] + M[i,] %*% (nvec-xvec)) / 
        beta(avec[i],bvec[i]) ) 
  } 
  sum(log(marg.vec))
}

#' @export
logMarg.DensSA <- function(M, mod.mat, xvec, nvec, avec, bvec) {

  marg.vec <- rep(NA, dim(M)[1])

  #calculate the product portion of integrated marginal likelihood
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) 

  for(i in 1:dim(M)[1]) { 
    p.vec <- prod( prod.vec^(1-M[i,]) )
    marg.vec[i] <- p.vec * 
      ( beta(avec[i] + M[i,] %*% xvec, bvec[i] + M[i,] %*% (nvec-xvec)) / 
        beta(avec[i], bvec[i]) ) * p.vec
  } 
  sum(log(marg.vec))
}

#' @export
logMarg.DensPE <- function(m, mod.size, xvec, nvec, avec, bvec, Lambda) {

  M <- diag(1, mod.size)
  u <- upper.tri(M)
  M[u] <- m
  M. <- t(M)
  l <- lower.tri(M.)
  M[l] <- M.[l]
  marg.vec <- rep(NA, dim(M)[1])

  #calculate the product portion of integrated marginal likelihood
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) 

  for(i in 1:dim(M)[1]) { 
    p.vec <- prod(prod.vec^(1-M[i,]))
    marg.vec[i] <- p.vec * 
      ( beta(avec[i] + M[i,] %*% xvec, bvec[i] + M[i,] %*% (nvec-xvec)) / 
        beta(avec[i],bvec[i]) ) * p.vec
  } 
  u <- upper.tri(M)
  -sum(log(marg.vec)) - Lambda*sum(M[u])
}


#' Simulated annealing method
#' @export
MEM_marginalLoop <- function(xvec, nvec, avec, bvec, INITIAL) {
  ###Function to calculate the MEM marginal densities for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  
  mod.mat <- list()
  k <- length(xvec)-1
  j <- 1

  while(k > 0) { 
    mod.mat[[j]] <- as.matrix(expand.grid( rep(list(c(0,1)), k) ))[
        order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 

    k <- k - 1
    j <- j + 1 
  }

  if( is.na(INITIAL[1]) ) { 
    M.init <- diag(1,length(xvec)) 
  } else { 
    M.init <- INITIAL 
  }
  
  MOld <- M.init  
  d <- length(xvec)   

  for(i in 1:(d-1)) {
    for(k in (i+1):d) {
      if (MOld[i,k]!=MOld[k,i]) {
        MOld[i,k]<-MOld[k,i]<-rbinom(1,1,0.5)
      }
    }
  }

  logDenOld <- logMarg.DensSA(MOld, mod.mat, xvec, nvec, avec, bvec)

  for(iter in 1:10000) {
    changed <- FALSE
    for( i in 1:(d-1)) {
      for(k in (i+1):d) {
        loc <- c(k,i) 
        MNew <- MOld
        MNew[loc[1],loc[2]] <- MNew[loc[2],loc[1]] <- 1 - MOld[loc[1],loc[2]]
        logDenNew <- logMarg.DensSA(MNew, mod.mat, xvec, nvec, avec, bvec)
        if (logDenNew > logDenOld) {
          MOld<-MNew
          changed<-TRUE
          logDenOld<-logDenNew
        }
      }
    }
    if (!changed) {
      break
    }
  }
  
  colnames(M.init) <- rownames(M.init) <- xvec
  colnames(MOld) <- rownames(MOld) <- xvec
  list(mod.mat = mod.mat, init=M.init, maximizer=MOld)
}


## Penalized Estimation with Simulated annealing method ##
library(GenSA)

#' @export
#' @importFrom GenSA GenSA
MEM_marginalPE <- function(xvec, nvec, avec, bvec, INITIAL, Lambda){
  ###Function to calculate the MEM marginal densities for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  
  mod.mat <- list()
  k <- length(xvec)-1
  j <- 1

  while(k > 0) { 
    mod.mat[[j]] <- (as.matrix(expand.grid( rep(list(c(0,1)), k) )))[
      order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
    k <- k - 1
    j <- j + 1 
  }

  if( is.na(INITIAL[1]) ) { 
    M.init <- diag(1,length(xvec)) 
  } else{ 
    M.init <- INITIAL 
  }
  
  sa.opt <- GenSA(par=M.init[upper.tri(M.init)], fn=logMarg.DensPE, 
    rep(0,sum(1:(length(xvec)-1))), rep(1,sum(1:(length(xvec)-1))), 
    mod.size=nrow(M.init), xvec=xvec, nvec=nvec, avec=avec, bvec=bvec, 
    Lambda=Lambda) 

  M.out <- diag(1, nrow(M.init))
  u <- upper.tri(M.out)
  M.out[u] <- sa.opt$par
  M. <- t(M.out)
  l <- lower.tri(M.)
  M.out[l] <- M.[l]

  colnames(M.init) <- rownames(M.init) <- xvec
  colnames(M.out) <- rownames(M.out) <- xvec
  list(mod.mat = mod.mat, init=M.init, maximizer=M.out)
}

#' @export
MEM_marginal <- function(xvec, nvec, avec, bvec) {
  ###Function to calculate the MEM marginal densities for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
  
  mod.mat <- list()
  k <- length(xvec) - 1
  j <- 1
  while(k > 0) { 
    mod.mat[[j]] <- (as.matrix(expand.grid( rep(list(c(0,1)), k) )))[
      order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
    k <- k - 1
    j <- j + 1 
  }
  H <- length(mod.mat)
  temp <- list()
  for(h in 1:(H-1)) { 
    temp[[h]] <- 1:dim(mod.mat[[h]])[1] 
  }
  temp[[H]] <- 1:2
  Mod.I <- as.matrix(expand.grid( temp ))

  # identify maximizer of marginal density ##
  log.Marg <- apply(Mod.I, MARGIN=1, FUN=logMarg.Dens, mod.mat, xvec, nvec, 
                    avec, bvec)
  o <- order(log.Marg,decreasing=TRUE)[1]

  MAX <- matrix(NA, length(xvec), length(xvec) )
  diag(MAX) <- rep(1,dim(MAX)[1])

  for(i in 1:(length(Mod.I[o,])-1)) { 
    MAX[(i+1):dim(MAX)[2], i] <- MAX[i, (i+1):dim(MAX)[2]] <- 
      mod.mat[[i]][Mod.I[o,i],] 
  }
  MAX[dim(MAX)[2], i+1] <- MAX[i+1, dim(MAX)[2]] <- 
    (0:1)[ Mod.I[o,length(Mod.I[o,])] ]
  colnames(MAX) <- rownames(MAX) <- xvec
  
  # Create list to return matrix which identifies which sources are included 
  # in each model and the marginal densities
  list(mod.mat = mod.mat, maximizer=MAX)
}

#' @export
MEM_modweight <- function(mod.mat, source.vec) {
  # Function to calculate model weights for each MEM given source inclusion 
  # probabilities applied to each row of the symmetric source inclusion prior 
  # probabilities
  # mod.mat: matrix which specifies which supplemental sources are included 
  #  for analysis of each subtype
  #source.mat: rows of symmetric matrix of source inclusion probabilties for MEM
  
  s.in <- source.vec #source inclusion probability
  s.ex <- 1-source.vec #source exclusion probability
  
  # Calculate model weights
  # TODO: change this to vapply
  q.mat <- sapply(1:dim(mod.mat)[1], 
    function(x) s.in^(mod.mat[x,]) * s.ex^(1-mod.mat[x,]))
  
  if( length(s.in) == 1 ) { 
    # if only one supplemental sources, q.mat is equal to q.vec for model 
    # weights
    q.vec <- q.mat 
  } else { 
    # if more than one supplemental source, take product of source 
    # probabilities for each model
    q.vec <- apply(q.mat, 2, prod) 
  }
  
  q.vec
}

#' @export
eval.Post <- function(p0, X, N, Omega, w){
  a=b=0.5 
  alph <- a + Omega%*%X
  beta <- b + (Omega%*%N-Omega%*%X)
  return( sum( (1-pbeta(p0, alph, beta))*w ) )
}

#' @export
gen.Post <- function(X, N, Omega){
  a <- b <- 0.5 
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  rbeta(1, alph, beta)
}

#' @export
ESS <- function(X,N,Omega) {
  a <- b <- 0.5 
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  alph + beta 
}

#' @export
samp.Post <- function(X, N, Omega, w) { 
  gen.Post(X, N, Omega[ which(rmultinom(1,1,w) == 1), ]) 
}

#' @export
post.HPD <- function(Data, pars, marg.M, alp) {
  U <- MEM.w(Data,pars,marg.M)
  out <- boa.hpd( 
    replicate(10000, samp.Post(Data$X, Data$N, U$models, U$weights[[1]]) ), 
    alp)

  K <- length(Data$X)

  for(j in 2:(K-1)) { 
    Ii <- c(j,1:(j-1),(j+1):K)
    out <- rbind(out, boa.hpd( 
      replicate(10000, samp.Post(Data$X[Ii], Data$N[Ii], U$models, 
                U$weights[[j]]) ), 
      alp))
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))

  out <- rbind(out, boa.hpd( 
    replicate(10000, samp.Post(Data$X[Ii], Data$N[Ii], U$models, 
              U$weights[[j]]) ), 
    alp))
  rownames(out) <- Data$X
  out
}

#' @export
post.ESS <- function(Data, pars, marg.M) {
  U <- MEM.w(Data,pars,marg.M)
  out <- U$weights[[1]]%*%ESS( Data$X, Data$N, U$models )
  K <- length(Data$X)
  for(j in 2:(K-1)){ 
    Ii <- c(j,1:(j-1),(j+1):K)
    out <- c(out, U$weights[[j]]%*%ESS(Data$X[Ii], Data$N[Ii], U$models) )
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))
  out <- out <- c(out, U$weights[[j]]%*%ESS(Data$X[Ii], Data$N[Ii], U$models) )
  names(out) <- Data$X
  out
}

#' @export
mem.Prior <- function(I, mod.mat, pr.Inclus) {
  M <- MEM.mat(I, mod.mat, nrow(pr.Inclus))
  mem <- M[ upper.tri(M) ]
  source.vec <- pr.Inclus[upper.tri(pr.Inclus)]
  s.in <- source.vec #source inclusion probability
  s.ex <- 1-source.vec #source exclusion probability
  prod( s.in^(mem) * s.ex^(1-mem) )
}

#' @export
mem.ij <- function(I, mod.mat, pr.Inclus,i,j) {
 MEM.mat(I, mod.mat, nrow(pr.Inclus))[i, j]
}

#' @export
ij.pep <- function(j, Mod.I, mod.mat, pr.Inclus, i, log.Marg, PRIOR) {
  I.1 <- apply( Mod.I, MARGIN=1, FUN=mem.ij, mod.mat, pr.Inclus, i, j)
  ( (exp(log.Marg) * PRIOR) %*% I.1 ) / sum( exp(log.Marg) * PRIOR ) 
}

#' @export
compPEP <- function(i, J, Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR) {
 sapply( (i+1):J, FUN=ij.pep, Mod.I, mod.mat, pr.Inclus, i, log.Marg, PRIOR) 
}

#' @export
mem.j <- function(I, mod.mat, pr.Inclus, j, m) {
 M <- MEM.mat(I, mod.mat, nrow(pr.Inclus)) 
 as.numeric(sum( as.numeric(M[j,] == m) ) == length(m)) 
}

#' @export
j.weight.1 <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
 I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
              c(1,mod.mat[[1]][i,]) )
 ( (exp(log.Marg) * PRIOR) %*% I.m) / sum( exp(log.Marg)*PRIOR ) 
}

#' @export
j.weight.j <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
  u <- mod.mat[[1]][i,]
  I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
               c( u[1:(j-1)],1,u[j:length(u)] ) )
  ( (exp(log.Marg) * PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR ) 
}

#' @export
j.weight.J <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
 I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
              c(mod.mat[[1]][i,],1) )
 ( (exp(log.Marg)*PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR ) 
}

#' @export
post.Weights <- function(j, J, Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR) {
  if (j==1) { 
    out <- sapply(1:nrow(mod.mat[[1]]), FUN=j.weight.1, Mod.I, mod.mat, 
                  pr.Inclus, j, log.Marg, PRIOR)
  } else if (j==J) { 
    out <- sapply(1:nrow(mod.mat[[1]]), FUN=j.weight.J, Mod.I, mod.mat, 
                  pr.Inclus, j, log.Marg, PRIOR)
 } else { 
    out <- sapply(1:nrow(mod.mat[[1]]), FUN=j.weight.j, Mod.I, mod.mat, 
                  pr.Inclus, j, log.Marg, PRIOR) }
 out
}

#' @export
MEM_FullInf <- function(xvec, nvec, avec, bvec, pr.Inclus, alp) {
  mod.mat <- list()
  k <- length(xvec)-1
  j <- 1

  while (k > 0) { 
    mod.mat[[j]] <- as.matrix(expand.grid( rep(list(c(0,1)), k) ))[
      order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
    k <- k - 1
    j <- j + 1 
  }
  H <- length(mod.mat)
  temp <- list()
  
  for (h in 1:(H-1)) { 
    temp[[h]] <- 1:dim(mod.mat[[h]])[1] 
  }
  temp[[H]] <- 1:2
  Mod.I <- as.matrix(expand.grid( temp ))

  ## identify maximizer of marginal density ##
  log.Marg <- apply(Mod.I, MARGIN=1, FUN=logMarg.Dens, mod.mat, xvec, nvec, 
                    avec, bvec)
  max.i <- order(log.Marg,decreasing=TRUE)[1]  
  MAX <- MEM.mat(Mod.I[max.i,], mod.mat, length(xvec))
  colnames(MAX) <- rownames(MAX) <- xvec

  ## compute prior + posterior MEM probabilities ##
  PRIOR <- apply( Mod.I, MARGIN=1, FUN=mem.Prior, mod.mat, pr.Inclus)
  POST <- ( exp(log.Marg) * PRIOR ) / sum( exp(log.Marg) * PRIOR )
  map.i <- order(POST,decreasing=TRUE)[1]
  MAP <- MEM.mat(Mod.I[map.i,], mod.mat, length(xvec))
  colnames(MAP) <- rownames(MAP) <- xvec

  ## Posterior Exchangeability Prob ##
  PEP <- matrix(NA, length(xvec), length(xvec)); colnames(PEP) <- rownames(PEP) <- xvec
  diag(PEP) <- 1
  for(i in 1:(nrow(PEP)-1)){ for(j in (i+1):ncol(PEP)){
   I.1 <- apply( Mod.I, MARGIN=1, FUN=mem.ij, mod.mat, pr.Inclus, i, j) #; print(sum(I.1))
   PEP[i,j] <- ( (exp(log.Marg)*PRIOR)%*%I.1)/sum( exp(log.Marg)*PRIOR )
  }}

  ## Marginal CDF and ESS ##
  pweights <- list()
  w <- rep(NA, nrow(mod.mat[[1]]))
  j <- 1

  for (i in 1:nrow(mod.mat[[1]])) {
    I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
                 c(1,mod.mat[[1]][i,]) )
   w[i] <- ( (exp(log.Marg) * PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR )
  }
  pweights[[j]] <- w

  for (j in 2:(length(xvec)-1)) {
    w <- rep(NA, nrow(mod.mat[[1]]))
    for(i in 1:nrow(mod.mat[[1]])) {
      u <- mod.mat[[1]][i,]
      I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
                   c( u[1:(j-1)],1,u[j:length(u)] ) )
      w[i] <- ( (exp(log.Marg)*PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR )
    }
    pweights[[j]] <- w
  }
  j <- length(xvec)

  for (i in 1:nrow(mod.mat[[1]])) {
    I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
                 c(mod.mat[[1]][i,],1) )
   w[i] <- ( (exp(log.Marg) * PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR )
  }
  pweights[[j]] <- w


  pESS <- CDF <- rep(NA, length(xvec))
  names(pESS) <- names(CDF) <- xvec
  HPD <- matrix(NA,2,length(xvec))
  rownames(HPD)=c("lower","upper")
  colnames(HPD) <- xvec
  marg.M <- MEM_marginal(xvec, nvec, avec, bvec)
  models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])

  CDF[1] <- eval.Post(p0, xvec, nvec, models, pweights[[1]])
  pESS[1] <- pweights[[1]] %*% ESS(xvec, nvec, models)
  HPD[,1] <- boa.hpd( 
    replicate(10000, samp.Post(xvec, nvec, models, pweights[[1]]) ), alp)
  K <- length(xvec)
  for (j in 2:(K-1)) { 
    Ii <- c(j,1:(j-1),(j+1):K)
    CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
    pESS[j] <- pweights[[j]]%*%ESS(xvec[Ii], nvec[Ii], models)
    HPD[,j] <- boa.hpd( 
      replicate(10000, samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), 
      alp)
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))
  CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
  pESS[j] <- pweights[[j]]%*%ESS(xvec[Ii], nvec[Ii], models)
  HPD[,j] <- boa.hpd( 
    replicate(10000, samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), 
    alp)

  list(mod.mat = mod.mat, maximizer=MAX, MAP=MAP, PEP=PEP, CDF=CDF, ESS=pESS, 
       HPD=HPD)
}

# START HERE

#' @export
MEM.w <- function(Data,pars,marg.M){
  #marg.M <- MEM_marginal(xvec=Data$X, nvec=Data$N, avec=rep(0.5,length(Data$N)), bvec=rep(0.5,length(Data$N)))
  pr.Inclus <- pars$UB*marg.M$maximizer; diag(pr.Inclus) <- rep(1,nrow(pr.Inclus)); pr.Inclus[which(pr.Inclus==0)] <- pars$LB
  weights <- list()
  for(i in 1:nrow(pr.Inclus)){ 
    weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], source.vec=pr.Inclus[i,][-i]) 
  }
  models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
  out <- list( models=models, weights=weights) #, alph=alph, beta=beta)
 return(out)
}

#' @export
MEM.cdf <- function(Data,pars,p0,marg.M){
  #marg.M <- MEM_marginal(xvec=Data$X, nvec=Data$N, avec=rep(0.5,length(Data$N)), bvec=rep(0.5,length(Data$N)))
  pr.Inclus <- pars$UB*marg.M$maximizer; diag(pr.Inclus) <- rep(1,nrow(pr.Inclus)); pr.Inclus[which(pr.Inclus==0)] <- pars$LB
  weights <- list() 
  for(i in 1:nrow(pr.Inclus)){ 
    weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], source.vec=pr.Inclus[i,][-i]) 
  }
  models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
  out <- eval.Post(p0, Data$X, Data$N, models, weights[[1]])
  K <- length(Data$X)
  for(j in 2:(K-1)){ Ii <- c(j,1:(j-1),(j+1):K)
                     out <- c(out, eval.Post(p0, Data$X[Ii], Data$N[Ii], models, weights[[j]]))
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))
  out <- c(out, eval.Post(p0, Data$X[Ii], Data$N[Ii], models, weights[[j]])) #; print(pr.Inclus)
  return(out)
}

#' @export
PostProb.mlab <- function(Data,pars,p0){
  pp.mem <- MEM.cdf(Data,pars,p0)
  if(is.na(pars$ESS)){ out <- sum( pars$m.BasketWeights*pp.mem )
  } else{ out <- rdirichlet(1, pars$m.BasketWeights*pars$ESS)*pp.mem }
  return(out)
}

#' @export
genObs.mlab <- function(hyp,N,pars,p0){
  D <- list(X=rbinom(rep(1,length(N)),N,hyp), N=N)
  return( PostProb.mlab(D,pars,p0) )
}

#' @export
checkThres.mlab <- function(THRES,null,N,pars,p0){ N.rep <- 5000
                                                   return( sum( replicate(N.rep, genObs.mlab(null,N,pars,p0) ) > THRES )/N.rep )
}

#' @export
CalibrateTypeIError.mlab <- function(null,Grid,N,pars,p0){
  typeIerror <- sapply(X=Grid,FUN=checkThres.mlab,null,N,pars,p0)
  return(list(Threshold=Grid, TypeIerror=typeIerror)) 
}

#' @export
Power.mlab <- function(ALT,thres,N,pars,p0){ N.rep <- 5000
                                             return( list( PostProbThreshold=thres, power=sum( replicate(N.rep, genObs.mlab(ALT,N,pars,p0) ) > thres )/N.rep ) )
}

#' @export
genObs.bask <- function(THRES,hyp,N,pars,p0){
  D <- list(X=rbinom(rep(1,length(N)),N,hyp), N=N)
  return( as.numeric( MEM.cdf(D,pars,p0) > THRES ) )
}

#' @export
checkThres.bask <- function(THRES,null,N,pars,p0){ N.rep <- 5000 ## Compute FWER ##
                                                   return(  1 - ( sum( rowSums( t( replicate(N.rep, genObs.bask(THRES,null,N,pars,p0) ) ) ) == 0 ) / N.rep ) )
}

#' @export
CalibrateFWER.bask <- function(null,Grid,N,pars,p0){
  FWER <- sapply(X=Grid,FUN=checkThres.bask,null,N,pars,p0)
  return(list(Threshold=Grid, FWER=FWER)) 
}

#' @export
Power.bask <- function(ALT,thres,N,pars,p0){ N.rep <- 5000
                                             return( list( PostProbThreshold=thres, 
                                                           power=colSums( t( replicate(N.rep, genObs.bask(thres,ALT,N,pars,p0) ) ) )/N.rep ) )
}
