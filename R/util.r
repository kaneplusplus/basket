
#' @importFrom foreach getDoParWorkers
num_workers <- function(factor = 1) {
  getDoParWorkers() * factor
}

boa.hpd <- function(x, alpha) {
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

sigmoid <- function(psi) { 
  1/(1+exp(-psi))
}

logit <- function(a0) { 
  log(a0) - log(1-a0)
}

MEM.mat <- function(Indices, mod.mat, H) {
  M <- matrix(NA, H, H )
  diag(M) <- rep(1,dim(M)[1])
  for(i in 1:(length(Indices)-1)) { 
    M[(i+1):dim(M)[2], i] <- M[i, (i+1):dim(M)[2]] <- mod.mat[[i]][Indices[i],]   }
  M[dim(M)[2], i+1] <- M[i+1, dim(M)[2]] <- (0:1)[ Indices[length(Indices)] ]
  M
}

#' @importFrom stats rgamma
rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}


logMarg.Dens <- function(I, mod.mat, xvec, nvec, avec, bvec) {
  M <- MEM.mat(I, mod.mat, length(xvec))
  marg.vec <- rep(NA, dim(M)[1])
  prod.vec <-
    beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  for (i in 1:dim(M)[1]) {
    p.vec <- prod(prod.vec ^ (1 - M[i, ]))
    marg.vec[i] <-
      (beta(avec[i] + M[i, ] %*% xvec, bvec[i] + M[i, ] %*% (nvec - xvec)) / beta(avec[i], bvec[i])) * p.vec
  } #return( list(mod = M, marg = prod(marg.vec)) )
  return(sum(log(marg.vec)))
}



logMarg.DensSA <- function(M, mod.mat, xvec, nvec, avec, bvec) {

  marg.vec <- rep(NA, dim(M)[1])

  #calculate the product portion of integrated marginal likelihood
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) 

  for(i in 1:dim(M)[1]) { 
    p.vec <- prod(prod.vec^(1-M[i,]))
    marg.vec[i] <-  
      ( beta(avec[i] + M[i,] %*% xvec, bvec[i] + M[i,] %*% (nvec-xvec)) / 
        beta(avec[i], bvec[i]) ) * p.vec
  } 
  sum(log(marg.vec))
}



logMarg.DensPE <-function(m, mod.size, xvec, nvec, avec, bvec, Lambda) {
    #M <- MEM.mat(I, mod.mat, length(xvec))
    M <- diag(1, mod.size)
    u <-
      upper.tri(M)
    M[u] <- m
    M. <- t(M)
    l <- lower.tri(M.)
    M[l] <- M.[l]
    marg.vec <- rep(NA, dim(M)[1])
    prod.vec <-
      beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
    for (i in 1:dim(M)[1]) {
      p.vec <- prod(prod.vec ^ (1 - M[i, ]))
      marg.vec[i] <-
        (beta(avec[i] + M[i, ] %*% xvec, bvec[i] + M[i, ] %*% (nvec - xvec)) / beta(avec[i], bvec[i])) * p.vec
    } #return( list(mod = M, marg = prod(marg.vec)) )
    u <- upper.tri(M)
    out <- sum(log(marg.vec)) - Lambda * sum(M[u])
    return(-out)
}


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

#' @importFrom itertools isplitRows
MEM_marginal <- function(xvec, nvec, avec, bvec) {
  ###Function to calculate the MEM marginal densities for binomial data with beta(alpha,beta) prior
  #xvec: vector with counts of those with event (primary source first, supplemental afterwards)
  #nvec: vector of total number of subjects
  #avec: vector of alpha parameters for beta priors for each source
  #bvec: vector of beta parameters for beta priors for each source
 
  mod.mat <- foreach(k = rev(seq_len(length(xvec) - 1))) %do% {
    m <- as.matrix(expand.grid( rep(list(c(0,1)), k) ))
    m[order(rowSums(m)),]
  }
 
  H <- length(mod.mat)
  temp <- foreach(h = seq_len(H-1)) %do% {
    seq_len(dim(mod.mat[[h]])[1])
  }
  temp[[H]] <- 1:2

  temp[[H]] <- 1:2
  Mod.I <- as.matrix(expand.grid( temp ))

  # identify maximizer of marginal density ##

  log.Marg <- foreach(mod_i = isplitRows(Mod.I, chunks = num_workers()), 
                      .combine = c) %dopar% {
    apply(mod_i, MARGIN=1, FUN=logMarg.Dens, mod.mat, xvec, nvec, avec, bvec)
  }
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

MEM_modweight <- function(mod.mat, source.vec) {
  # Function to calculate model weights for each MEM given source inclusion 
  # probabilities applied to each row of the symmetric source inclusion prior 
  # probabilities
  # mod.mat: matrix which specifies which supplemental sources are included 
  # for analysis of each subtype
  # source.mat: rows of symmetric matrix of source inclusion probabilties for MEM
  
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

#eval.Post <- function(p0, X, N, Omega, w){
#  a=b=0.5 
#  alph <- a + Omega %*% X
#  beta <- b + (Omega %*% N - Omega %*% X)
#  return( sum( (1-pbeta(p0, alph, beta))*w ) )
#}



ESS <- function(X,N,Omega,a,b) {
  #a <- b <- 0.5 
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  alph + beta 
}



post.HPD <- function(Data, pars, marg.M, alp) {
  U <- MEM.w(Data,pars,marg.M)
  out <- boa.hpd( 
    replicate(10000, samp.Post(Data$X, Data$N, U$models, U$weights[[1]]) ), 
    alp)

  K <- length(Data$X)

  out <- rbind(out,
    foreach(js = isplitVector(2:(K-1), chunks = num_workers()),
            .combine = rbind) %dopar% {
      foreach(j = js, .combine = rbind) %do% {
        Ii <- c(j,1:(j-1),(j+1):K)
        boa.hpd(replicate(10000, samp.Post(Data$X[Ii], Data$N[Ii], U$models,
                          U$weights[[j]])),
                alp)
      }
    })

#  for(j in 2:(K-1)) { 
#    Ii <- c(j,1:(j-1),(j+1):K)
#    out <- rbind(out, boa.hpd( 
#      replicate(10000, samp.Post(Data$X[Ii], Data$N[Ii], U$models, 
#                U$weights[[j]]) ), 
#      alp))
#  }
#  j <- j + 1
  j <- K
  Ii <- c(j,1:(j-1))

  out <- rbind(out, boa.hpd( 
    replicate(10000, samp.Post(Data$X[Ii], Data$N[Ii], U$models, 
              U$weights[[j]]) ), 
    alp))
  rownames(out) <- Data$X
  out
}

#' @importFrom itertools isplitVector
#' @importFrom foreach foreach %dopar%
post.ESS <- function(Data, pars, marg.M, shape1, shape2) {
  U <- MEM.w(Data,pars,marg.M)
  out <- U$weights[[1]]%*%ESS( Data$X, Data$N, U$models, shape1, shape2 )
  K <- length(Data$X)
  out <- c(out,
    foreach(js = isplitVector(2:(K-1), chunks = num_workers()),
            .combine = c) %dopar% {
      foreach(j = js, .combine = c) %do% {
        Ii <- c(j, 1:(j-1), (j+1):K)
        U$weights[[j]] %*% ESS(Data$X[Ii], Data$N[Ii], U$models, shape1, shape2)
      }
    })
#  for(j in 2:(K-1)){ 
#    Ii <- c(j,1:(j-1),(j+1):K)
#    out <- c(out, U$weights[[j]]%*%ESS(Data$X[Ii], Data$N[Ii], U$models) )
#  }
  j <- K
  Ii <- c(j,1:(j-1))
  out <- out <- c(out, U$weights[[j]]%*%ESS(Data$X[Ii], Data$N[Ii], U$models), shape1, shape2)
  names(out) <- Data$X
  out
}

mem.Prior <- function(I, mod.mat, pr.Inclus) {
  M <- MEM.mat(I, mod.mat, nrow(pr.Inclus))
  mem <- M[ upper.tri(M) ]
  source.vec <- pr.Inclus[upper.tri(pr.Inclus)]
  s.in <- source.vec #source inclusion probability
  s.ex <- 1-source.vec #source exclusion probability
  prod( s.in^(mem) * s.ex^(1-mem) )
}

mem.ij <- function(I, mod.mat, pr.Inclus,i,j) {
 MEM.mat(I, mod.mat, nrow(pr.Inclus))[i, j]
}

ij.pep <- function(j, Mod.I, mod.mat, pr.Inclus, i, log.Marg, PRIOR) {
  I.1 <- apply( Mod.I, MARGIN=1, FUN=mem.ij, mod.mat, pr.Inclus, i, j)
  ( (exp(log.Marg) * PRIOR) %*% I.1 ) / sum( exp(log.Marg) * PRIOR ) 
}

compPEP <- function(i, J, Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR) {
 sapply( (i+1):J, FUN=ij.pep, Mod.I, mod.mat, pr.Inclus, i, log.Marg, PRIOR) 
}

mem.j <- function(I, mod.mat, pr.Inclus, j, m) {
 M <- MEM.mat(I, mod.mat, nrow(pr.Inclus)) 
 as.numeric(sum( as.numeric(M[j,] == m) ) == length(m)) 
}

j.weight.1 <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
 I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
              c(1,mod.mat[[1]][i,]) )
 ( (exp(log.Marg) * PRIOR) %*% I.m) / sum( exp(log.Marg)*PRIOR ) 
}

j.weight.j <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
  u <- mod.mat[[1]][i,]
  I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
               c( u[1:(j-1)],1,u[j:length(u)] ) )
  ( (exp(log.Marg) * PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR ) 
}

j.weight.J <- function(i, Mod.I, mod.mat, pr.Inclus, j, log.Marg, PRIOR) {
 I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
              c(mod.mat[[1]][i,],1) )
 ( (exp(log.Marg)*PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR ) 
}

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

# TODO: What is this? Do we need it?
#MEM_FullInf <- function(xvec, nvec, avec, bvec, pr.Inclus, alp) {
#  mod.mat <- list()
#  k <- length(xvec)-1
#  j <- 1
#
#  while (k > 0) { 
#    mod.mat[[j]] <- as.matrix(expand.grid( rep(list(c(0,1)), k) ))[
#      order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
#    k <- k - 1
#    j <- j + 1 
#  }
#  H <- length(mod.mat)
#  temp <- list()
#  
#  for (h in 1:(H-1)) { 
#    temp[[h]] <- 1:dim(mod.mat[[h]])[1] 
#  }
#  temp[[H]] <- 1:2
#  Mod.I <- as.matrix(expand.grid( temp ))
#
#  ## identify maximizer of marginal density ##
#  log.Marg <- apply(Mod.I, MARGIN=1, FUN=logMarg.Dens, mod.mat, xvec, nvec, 
#                    avec, bvec)
#  max.i <- order(log.Marg,decreasing=TRUE)[1]  
#  MAX <- MEM.mat(Mod.I[max.i,], mod.mat, length(xvec))
#  colnames(MAX) <- rownames(MAX) <- xvec
#
#  ## compute prior + posterior MEM probabilities ##
#  PRIOR <- apply( Mod.I, MARGIN=1, FUN=mem.Prior, mod.mat, pr.Inclus)
#  POST <- ( exp(log.Marg) * PRIOR ) / sum( exp(log.Marg) * PRIOR )
#  map.i <- order(POST,decreasing=TRUE)[1]
#  MAP <- MEM.mat(Mod.I[map.i,], mod.mat, length(xvec))
#  colnames(MAP) <- rownames(MAP) <- xvec
#
#  ## Posterior Exchangeability Prob ##
#  PEP <- matrix(NA, length(xvec), length(xvec)); colnames(PEP) <- rownames(PEP) <- xvec
#  diag(PEP) <- 1
#  for(i in 1:(nrow(PEP)-1)){ for(j in (i+1):ncol(PEP)){
#   I.1 <- apply( Mod.I, MARGIN=1, FUN=mem.ij, mod.mat, pr.Inclus, i, j) #; print(sum(I.1))
#   PEP[i,j] <- ( (exp(log.Marg)*PRIOR)%*%I.1)/sum( exp(log.Marg)*PRIOR )
#  }}
#
#  ## Marginal CDF and ESS ##
#  pweights <- list()
#  w <- rep(NA, nrow(mod.mat[[1]]))
#  j <- 1
#
#  for (i in 1:nrow(mod.mat[[1]])) {
#    I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
#                 c(1,mod.mat[[1]][i,]) )
#   w[i] <- ( (exp(log.Marg) * PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR )
#  }
#  pweights[[j]] <- w
#
#  for (j in 2:(length(xvec)-1)) {
#    w <- rep(NA, nrow(mod.mat[[1]]))
#    for(i in 1:nrow(mod.mat[[1]])) {
#      u <- mod.mat[[1]][i,]
#      I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
#                   c( u[1:(j-1)],1,u[j:length(u)] ) )
#      w[i] <- ( (exp(log.Marg)*PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR )
#    }
#    pweights[[j]] <- w
#  }
#  j <- length(xvec)
#
#  for (i in 1:nrow(mod.mat[[1]])) {
#    I.m <- apply(Mod.I, MARGIN=1, FUN=mem.j, mod.mat, pr.Inclus, j, 
#                 c(mod.mat[[1]][i,],1) )
#   w[i] <- ( (exp(log.Marg) * PRIOR) %*% I.m ) / sum( exp(log.Marg) * PRIOR )
#  }
#  pweights[[j]] <- w
#
#
#  pESS <- CDF <- rep(NA, length(xvec))
#  names(pESS) <- names(CDF) <- xvec
#  HPD <- matrix(NA,2,length(xvec))
#  rownames(HPD)=c("lower","upper")
#  colnames(HPD) <- xvec
#  marg.M <- MEM_marginal(xvec, nvec, avec, bvec)
#  models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
#
#  CDF[1] <- eval.Post(p0, xvec, nvec, models, pweights[[1]])
#  pESS[1] <- pweights[[1]] %*% ESS(xvec, nvec, models)
#  HPD[,1] <- boa.hpd( 
#    replicate(10000, samp.Post(xvec, nvec, models, pweights[[1]]) ), alp)
#  K <- length(xvec)
#  for (j in 2:(K-1)) { 
#    Ii <- c(j,1:(j-1),(j+1):K)
#    CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
#    pESS[j] <- pweights[[j]]%*%ESS(xvec[Ii], nvec[Ii], models)
#    HPD[,j] <- boa.hpd( 
#      replicate(10000, samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), 
#      alp)
#  }
#  j <- j + 1
#  Ii <- c(j,1:(j-1))
#  CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
#  pESS[j] <- pweights[[j]]%*%ESS(xvec[Ii], nvec[Ii], models)
#  HPD[,j] <- boa.hpd( 
#    replicate(10000, samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), 
#    alp)
#
#  list(mod.mat = mod.mat, maximizer=MAX, MAP=MAP, PEP=PEP, CDF=CDF, ESS=pESS, 
#       HPD=HPD)
#}

MEM.w <- function(Data, pars, marg.M) {
  pr.Inclus <- pars$UB*marg.M$maximizer
  diag(pr.Inclus) <- rep(1,nrow(pr.Inclus))
  pr.Inclus[which(pr.Inclus==0)] <- pars$LB
  weights <- list()
  for(i in 1:nrow(pr.Inclus)) {
    weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], 
      source.vec = pr.Inclus[i,][-i]) 
  }
  models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
  list(models = models, weights = weights) 
}

#' @importFrom foreach foreach %dopar%
#' @importFrom itertools isplitVector
MEM.cdf <- function(Data, pars, p0, marg.M) {
  pr.Inclus <- pars$UB*marg.M$maximizer
  diag(pr.Inclus) <- rep(1,nrow(pr.Inclus))
  pr.Inclus[which(pr.Inclus==0)] <- pars$LB
  weights <- list() 
  for(i in 1:nrow(pr.Inclus)) {
    weights[[i]] <- MEM_modweight(mod.mat = marg.M$mod.mat[[1]], 
      source.vec = pr.Inclus[i,][-i]) 
  }
  models <- cbind(rep(1,dim(marg.M$mod.mat[[1]])[1]),marg.M$mod.mat[[1]])
  out <- eval.Post(p0, Data$X, Data$N, models, weights[[1]])
  K <- length(Data$X)
  out <- c(out,
    foreach(js = isplitVector(2:(K-1), chunks = num_workers()),
            .combine = c) %dopar% {
      foreach(j = js, .combine = c) %do% {
        Ii <- c(j,1:(j-1),(j+1):K)
        eval.Post(p0, Data$X[Ii], Data$N[Ii], models, weights[[j]])      
      }
    })
#  for(j in 2:(K-1)) { 
#    Ii <- c(j,1:(j-1),(j+1):K)
#    out <- c(out, eval.Post(p0, Data$X[Ii], Data$N[Ii], models, weights[[j]]))
#  }
#  j <- j + 1
  j <- K
  Ii <- c(j,1:(j-1))
  c(out, eval.Post(p0, Data$X[Ii], Data$N[Ii], models, weights[[j]])) 
}

PostProb.mlab <- function( Data, pars, p0) {

  pp.mem <- MEM.cdf(Data,pars,p0)

  if(is.na(pars$ESS)) { 
    out <- sum( pars$m.BasketWeights*pp.mem )
  } else{ 
    out <- rdirichlet(1, pars$m.BasketWeights*pars$ESS)*pp.mem 
  }
  out
}

genObs.mlab <- function(hyp, N, pars, p0) {
  D <- list(X=rbinom(rep(1,length(N)),N,hyp), N=N)
  PostProb.mlab(D,pars,p0)
}

checkThres.mlab <- function(THRES, null, N, pars, p0) { 
  N.rep <- 5000
  sum( replicate(N.rep, genObs.mlab(null,N,pars,p0) ) > THRES ) / N.rep 
}

CalibrateTypeIError.mlab <- function(null, Grid, N, pars, p0) {
  typeIerror <- sapply(X=Grid, FUN=checkThres.mlab, null, N, pars, p0)
  list(Threshold=Grid, TypeIerror=typeIerror)
}

Power.mlab <- function(ALT, thres, N, pars, p0) { 
  N.rep <- 5000
  list(PostProbThreshold=thres, 
       power=sum(replicate(N.rep, genObs.mlab(ALT,N,pars,p0)) > thres)/N.rep ) 
}

genObs.bask <- function(THRES, hyp, N, pars, p0) {
  D <- list(X=rbinom(rep(1,length(N)),N,hyp), N=N)
  as.numeric( MEM.cdf(D, pars, p0) > THRES ) 
}

checkThres.bask <- function(THRES, null, N, pars, p0) { 
  N.rep <- 5000 ## Compute FWER ##
  1 - (
    sum(rowSums(t(replicate(N.rep, genObs.bask(THRES, null, N, pars, p0))))==0)/
    N.rep ) 
}

CalibrateFWER.bask <- function(null, Grid, N, pars, p0){
  FWER <- sapply(X = Grid, FUN = checkThres.bask, null, N, pars, p0)
  list(Threshold = Grid, FWER = FWER)
}

Power.bask <- function(ALT, thres, N, pars, p0) { 
  N.rep <- 5000
  list(PostProbThreshold=thres, 
       power=colSums(t(replicate(N.rep, genObs.bask(thres,ALT,N,pars,p0)))) /
             N.rep)
}



####################################################################
#### Alternative Computation of ESS using HPD interval matching ####
####################################################################
#fit <- mem_full_bayes(responses = Data$X, size = Data$N, name = Data$Basket, p0 = 0.25, shape1 = 0.5, shape2 = 0.5, alpha = 0.05, call = NULL)


euc.dist <- function(x1, x2, w = c(1, 1)) {
  if (sum(is.na(x1)) > 1) {
    out <- Inf
  } else{
    out <- sqrt(sum(w * ((x1 - x2) ^ 2)))
  }
  return(out)
}


dist.beta.HPD <- function(ess, fit, alpha, jj) {
  al <- fit$mean_est[jj] * ess
  al <- max(1e-2, al)
  be <- ess - al
  be <- max(1e-2, be)
  return(euc.dist(fit$HPD[, jj], qbeta(
    c(alpha / 2, 1 - alpha / 2),
    al, be
  
  )))
}


dist.beta.HPDwid <- function(ess, fit, alpha, jj) {
  al <- fit$median_est[jj] * ess
  al <- max(1e-2, al)
  be <- ess - al
  be <- max(1e-2, be)
  return(abs((fit$HPD[2, jj] - fit$HPD[1, jj]) - (diff(
    qbeta(
      c(alpha / 2, 1 - alpha / 2),
      al, be
    )
  ))))
}


ESS.from.HPD.i <- function(jj, fit, alpha) {
  library(GenSA)
  opt <-
    GenSA::GenSA(
      par = 1,
      fn = dist.beta.HPD,
      lower = 0,
      upper = 10000000,
      #control=list(maxit=pars$DTW.maxit),
      fit = fit,
      alpha = alpha,
      jj = jj
    )
  return(opt$par)
}

ESS.from.HPDwid.i <- function(jj, fit, alpha) {
  library(GenSA)
  opt <-
    GenSA::GenSA(
      par = 1,
      fn = dist.beta.HPDwid,
      lower = 0,
      upper = 10000000,
      #control=list(maxit=pars$DTW.maxit),
      fit = fit,
      alpha = alpha,
      jj = jj
    )
  return(opt$par)
}


calc.ESS.from.HPD <- function(fit, alpha) {
  ## fit is list with median vec and HPD vec ##
  return(sapply(
    1:length(fit$mean_est),
    FUN = ESS.from.HPD.i,
    fit = fit,
    alpha = alpha
  ))
}


calc.ESS.from.HPDwid <- function(fit, alpha) {
  ## fit is list with median vec and HPD vec ##
  return(sapply(
    1:length(fit$median_est),
    FUN = ESS.from.HPDwid.i,
    fit = fit,
    alpha = alpha
  ))
}

#calc.ESS.from.HPD(fit, alpha=0.05)

#calc.ESS.from.HPDwid(fit, alpha=0.05)
####################################################################
####################################################################
####################################################################


####################################################################
########## MCMC for Bayes with Metropolis-Hastings #################
####################################################################

flip.MEM <- function(v, MOld, M){
  out <- MOld
  if( MOld[ which(M==v)[1] ]==1 ){ out[ which(M==v) ] <- 0
  } else{ out[ which(M==v) ] <- 1 }
  return(out)
}


mem.Prior.Mat <- function(M,mod.mat, pr.Inclus) {
  #M <- MEM.mat(I, mod.mat, nrow(pr.Inclus));
  mem <- M[upper.tri(M)]
  source.vec <- pr.Inclus[upper.tri(pr.Inclus)]
  s.in <- source.vec #source inclusion probability
  s.ex <- 1 - source.vec #source exclusion probability
  return(prod(s.in ^ (mem) * s.ex ^ (1 - mem)))
}

update.MH <-
  function(MOld,
           M,
           xvec,
           nvec,
           avec,
           bvec,
           mod.mat,
           Prior.EP) {
    K <- max(M, na.rm = TRUE) + 1
    v <-
      sample(1:(K - 1), (1:(K - 1))[which(rmultinom(1, 1, (((
        K - 1
      ):1) ^ 3) / sum(((
        K - 1
      ):1) ^ 3)) == 1)])
    
    MProp <- flip.MEM(v[1], MOld, M)
    if (length(v) > 1) {
      for (ii in 2:length(v)) {
        MProp <- flip.MEM(v[ii], MProp, M)
      }
    }
    
    rho <-
      exp(
        logMarg.DensSA(MProp, mod.mat, xvec, nvec, avec, bvec) + log(mem.Prior.Mat(MProp, mod.mat, Prior.EP)) -
          (
            logMarg.DensSA(MOld, mod.mat, xvec, nvec, avec, bvec) +  log(mem.Prior.Mat(MOld, mod.mat, Prior.EP))
          )
      )
    #cat(rho, rho1, xvec,nvec, avec, bvec)
    if (rho >= 1) {
      out <-
        MProp
    } else{
      if (rbinom(1, 1, rho) == 1) {
        out <- MProp
      } else{
        out <- MOld
      }
    }
    
    return(out)
  }

I.models <- function(hh, models, Samp) {
  K <- length(models[1, ])
  if (hh == 1) {
    Ii <- 1:K
  } else if (hh == K) {
    Ii <- c(K, 1:(K - 1))
  } else{
    Ii <- c(hh, 1:(hh - 1), (hh + 1):K)
  }
  return(which(apply(
    models,
    MARGIN = 1,
    FUN = function(x, t) {
      sum(x == t)
    },
    t = Samp[hh, Ii]
  ) == K))
}

models.Count <- function(Samp, models) {
  out <- matrix(0, nrow(models), ncol(models))
  u <-
    sapply(1:length(models[1, ]),
           FUN = I.models,
           models = models,
           Samp = Samp)
  for (i in 1:length(u)) {
    out[u[i], i] <- 1
  }
  return(out)
}



mem.PostProb <- function(model, method = "samples", fit) {
  ## Incorporate direction ##
  if (method == "mixture") {
    out <-
      eval.Post(
        model$p0[1],
        model$responses,
        model$size,
        model$models,
        model$pweights[[1]],
        model$shape1[1],
        model$shape2[1],
        alternative = model$alternative
      )
    K <- length(model$responses)
    for (j in 2:(K - 1)) {
      Ii <- c(j, 1:(j - 1), (j + 1):K)
      out <-
        c(
          out,
          eval.Post(
            model$p0[j],
            model$responses[Ii],
            model$size[Ii],
            model$models,
            model$pweights[[j]],
            model$shape1[j],
            model$shape2[j],
            alternative = model$alternative
          )
        )
    }
    j <- j + 1
    Ii <- c(j, 1:(j - 1))
    out <-
      c(
        out,
        eval.Post(
          model$p0[j],
          model$responses[Ii],
          model$size[Ii],
          model$models,
          model$pweights[[j]],
          model$shape1[j],
          model$shape2[j],
          alternative = model$alternative
        )
      )
  } else{
    if (model$alternative == "greater") {
      out <-
        sapply(
          1:ncol(fit$samples),
          FUN = function(j, x, t) {
            return(sum(x[, j] > t[j]) / length(x[, j]))
          },
          x = fit$samples,
          t = model$p0
        )
    } else{
      out <-
        sapply(
          1:ncol(fit$samples),
          FUN = function(j, x, t) {
            return(sum(x[, j] < t[j]) / length(x[, j]))
          },
          x = fit$samples,
          t = model$p0
        )
    }
  }
  names(out) <- model$name
  return(out)
}





#############################################################
### Added Shape and Direction to Posterior Prob #############
#############################################################
eval.Post <- function(p0, X, N, Omega, w, a, b, alternative = "greater") 
{
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  if (alternative == "greater") {
    out <- sum((1 - pbeta(p0, alph, beta)) * w)
  } else{
    out <- sum((pbeta(p0, alph, beta)) * w)
  }
  return(out)
}

#############################################################
### Fixed shape paramter specification in gen.Post() ########
#############################################################
gen.Post <- function(X, N, Omega, a, b) {
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  return(rbeta(1, alph, beta))
}

samp.Post <- function(X, N, Omega, w, a, b) {
  return(gen.Post(X, N, Omega[which(rmultinom(1, 1, w) == 1), ], a, b))
}

samplePostOneBasket <- function(index, model, num_samples = 10000)
{
  if (index == 1)
  {
    result <- replicate(
      num_samples,
      samp.Post(
        model$responses,
        model$size,
        model$models,
        model$pweights[[1]],
        model$shape1[1],
        model$shape2[1]
      )
    )
    return(result)
  }
  numBasket <- length(model$responses)
  if (index == numBasket)
  {
    Ii <- c(numBasket, 1:(numBasket - 1))
    result <- 
               replicate(
                 num_samples,
                 samp.Post(
                   model$responses[Ii],
                   model$size[Ii],
                   model$models,
                   model$pweights[[numBasket]],
                   model$shape1[numBasket],
                   model$shape2[numBasket]
                 )
               )
    return(result)
  }
  Ii <- c(index, 1:(index - 1), (index + 1):numBasket )
  result <- 
    replicate(
      num_samples,
      samp.Post(
        model$responses[Ii],
        model$size[Ii],
        model$models,
        model$pweights[[index]],
        model$shape1[index],
        model$shape2[index]
      )
    )
  return(result)
}

sample_posterior.full_bayes <- function(model, num_samples = 10000) {
  ret <- replicate(
    num_samples,
    samp.Post(
      model$responses,
      model$size,
      model$models,
      model$pweights[[1]],
      model$shape1[1],
      model$shape2[1]
    )
  )
  K <- length(model$responses)
  ret <- rbind(ret,
               foreach(j = 2:(K - 1), .combine = rbind) %do% {
                 Ii <- c(j, 1:(j - 1), (j + 1):K)
                 replicate(
                   num_samples,
                   samp.Post(
                     model$responses[Ii],
                     model$size[Ii],
                     model$models,
                     model$pweights[[j]],
                     model$shape1[j],
                     model$shape2[j]
                   )
                 )
               })
  j <- K
  Ii <- c(j, 1:(j - 1))
  ret <- rbind(ret,
               replicate(
                 num_samples,
                 samp.Post(
                   model$responses[Ii],
                   model$size[Ii],
                   model$models,
                   model$pweights[[j]],
                   model$shape1[j],
                   model$shape2[j]
                 )
               ))
  ret <- t(ret)
  dimnames(ret) <- list(NULL, model$name)
  ret
}

clusterComp <- function(basketRet) {
  PEP <- basketRet$PEP
  name <- basketRet$name
  p0 <- basketRet$p0
  allSamp <- basketRet$samples
  graph <-
    igraph::graph_from_adjacency_matrix(PEP,
                                        mode = "undirected",
                                        weighted = TRUE,
                                        diag = FALSE)
  result <- factor(cluster_louvain(graph, weights = E(graph)$weight)$membership)
  numClusters <- length(levels(result))

  sampleC <- list()
  cName <- c()
  clusterElement <- list()
  for (k in 1:numClusters) {
    rank <- which(result == levels(result)[k])
    
    cBasket <- name[rank]
    clusterElement[[k]] <- cBasket
    print(paste0("#### Cluster Assignment: ", k))
    print(rank)
    print(name[rank])
    cName <- c(cName, paste0("Cluster ", k))

    sampV <- as.vector(allSamp[, rank])
    sampleC[[k]] <- sampV
  }
  names(sampleC) <- cName
  
  ret <- basketRet
  ret$name <- cName
  ret$cluster <- clusterElement
  #browser()
  #ret$models <- models
  #ret$pweights <- pweights
  ret$samples <- sampleC
  
  p0Test <-unique(ret$p0)
  allCDF<-matrix(0, 0, numClusters)
  for (kk in 1:length(p0Test))
  {
    if (ret$alternative == "greater") {
      res <- unlist(lapply(
        1:numClusters,
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = sampleC,
        t = p0Test[kk]
      ))
    } else{
      res <- unlist(lapply(
        1:numClusters,
        FUN = function(j, x, t) {
          return(sum(x[[j]] > t) / length(x[[j]]))
        },
        x = sampleC,
        t = p0Test[kk]
      ))
    }
    allCDF <- rbind(allCDF, res)
  }
  colnames(allCDF) <- cName
  rownames(allCDF) <- p0Test
  ret$post.prob <- allCDF
  
  #ret$accept.rate <- (n.chg) / niter.MCMC
  ret$mean_est <- unlist(lapply(ret$samples, mean))
  ret$median_est <- unlist(lapply(ret$samples, median))
  #ret$PEP <- PEP
  ret$HPD <-matrix(unlist(
    lapply(
      ret$samples,
      #MARGIN = 2,
      FUN = boa.hpd,
      alpha = ret$alpha
    )), nrow = 2, byrow = FALSE)
  colnames(ret$HPD) <- cName
  row.names(ret$HPD) <-c("Lower Bound", "Upper Bound")
  ret$ESS <-
    calc.ESS.from.HPD(fit = ret, alpha = ret$alpha)
  names(ret$ESS) <- cName
  #ret$ESS2 <-
  #  calc.ESS.from.HPDwid(fit = ret, alpha = ret$alpha)
  #names(ret$ESS2) <- cName
  #ret$PostProb <- mem.PostProb(MODEL, method = "samples", fit = ret)
  class(ret) <- c("full_bayes", "exchangeability_model")
  return(ret)
}

