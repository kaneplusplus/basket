# arguments for bayes, empirical bayes, constrained empirical bayes.

runSim <- function(ITER, N, prob, pr.Inclus, p0, alp){

  seed <- round(runif(1, 1, 1000000000))
  set.seed(seed)
  D <- list(X = rbinom(N, N, prob), N = N)

  xvec <- D$X
  nvec <- D$N
  avec <- bvec <- rep(0.5, length(D$X))

  mod.mat <- list()
  k <- length(xvec)-1
  j <- 1
  while (k > 0) { 
    mod.mat[[j]] <- as.matrix(expand.grid(rep(list(c(0,1)), k)))[
      order(rowSums(as.matrix(expand.grid(rep(list(c(0,1)), k))))),] 
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
  log.Marg <- apply(Mod.I, MARGIN = 1, FUN = logMarg.Dens, mod.mat, xvec, nvec, 
                    avec, bvec)
  max.i <- order(log.Marg,decreasing = TRUE)[1]  
  MAX <- MEM.mat(Mod.I[max.i,], mod.mat, length(xvec))
  colnames(MAX) <- rownames(MAX) <- xvec

  ## compute prior + posterior MEM probabilities ##
  PRIOR <- apply(Mod.I, MARGIN = 1, FUN = mem.Prior, mod.mat, pr.Inclus)
  POST <- (exp(log.Marg) * PRIOR) / sum(exp(log.Marg)*PRIOR)
  map.i <- order(POST, decreasing = TRUE)[1]
  MAP <- MEM.mat(Mod.I[map.i, ], mod.mat, length(xvec))
  colnames(MAP) <- rownames(MAP) <- xvec

  ## Posterior Exchangeability Prob ##
  PEP. <- matrix(NA, length(xvec), length(xvec))
  colnames(PEP.) <- rownames(PEP.) <- xvec
  diag(PEP.) <- 1

#  pep <- foreach(i = 1:(nrow(PEP.)-1), .combine = list, .multicombine = TRUE)  %dopar% compPEP(i, length(xvec), Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR)
  pep <- foreach(i = 1:(nrow(PEP.)-1), .multicombine = TRUE) %dopar% {
    compPEP(i, length(xvec), Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR)
  }
  l <- lower.tri(PEP.); PEP.[l] <- unlist(pep); PEP <- t(PEP.)

  ## Marginal CDF and ESS ##
  pweights <- foreach( j = 1:length(xvec), .multicombine = TRUE) %dopar% {
    post.Weights(j, length(xvec), Mod.I, mod.mat, pr.Inclus, log.Marg, PRIOR)
  }

  pESS <- CDF <- rep(NA, length(xvec))
  names(pESS) <- names(CDF) <- xvec
  HPD <- matrix(NA,2,length(xvec))
  rownames(HPD)=c("lower","upper")
  colnames(HPD) <- xvec
  models <- cbind(rep(1,dim(mod.mat[[1]])[1]),mod.mat[[1]])

  CDF[1] <- eval.Post(p0, xvec, nvec, models, pweights[[1]])
  pESS[1] <- pweights[[1]] %*% ESS(xvec, nvec, models)
  HPD[,1] <- boa.hpd(replicate(10000, 
                               samp.Post(xvec, nvec, models, pweights[[1]]) ), 
                     alp)

  K <- length(xvec)
  for (j in 2:(K-1)) { 
    Ii <- c(j,1:(j-1),(j+1):K)
    CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
    pESS[j] <- pweights[[j]] %*% ESS(xvec[Ii], nvec[Ii], models)
    HPD[,j] <- boa.hpd(replicate(10000, 
      samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), alp)
  }
  j <- j + 1
  Ii <- c(j,1:(j-1))
  CDF[j] <- eval.Post(p0, xvec[Ii], nvec[Ii], models, pweights[[j]])
  pESS[j] <- pweights[[j]] %*% ESS(xvec[Ii], nvec[Ii], models)
  HPD[,j] <- boa.hpd(replicate(10000, 
    samp.Post(xvec[Ii], nvec[Ii], models, pweights[[j]]) ), alp)

  # TODO: Make this a parameter.
  if (!dir.exists("compare-OCs")) {
    dir.create("compare-OCs")
  }

  dput(list(mod.mat = mod.mat, maximizer = MAX, MAP = MAP, PEP = PEP, 
            CDF = CDF, ESS = pESS, HPD = HPD), 
       file = paste("compare-OCs/MEM-",ITER,".txt",sep=""))

  ###########
  ### SEM ###
  library(rjags)

  data <- list(J = 6, N = nvec, Y = xvec)
  inits <- list(mu = rep(0, 6), M = 0, inv.sig2 = 1) 
  m <- jags.model("../Utils/Berry-ig_model.txt", data, inits, n.chains=3, 
                  quiet=TRUE)
  update(m, 10000, progress.bar="none")
  S <- coda.samples(m, "pr", n.iter=10000, progress.bar="none")

  compSEM.PP <- function(j, S, p0) {
    return(sum(c(S[[1]][,j], S[[2]][,j], S[[3]][,j]) > p0) /
           length(c(S[[1]][,j], S[[2]][,j], S[[3]][,j]) ))
  }

  compSEM.HPD <- function(j, S, alp) {
    return( boa.hpd( c(S[[1]][,j],S[[2]][,j],S[[3]][,j]), alp) )
  }

  sem.PPROB <- foreach(j = 1:data$J, .combine = c, 
                       .multicombine = TRUE) %dopar% {
    compSEM.PP(j, S, p0)
  }

  sem.HPD <- foreach(j = 1:data$J, .combine = cbind, 
                     .multicombine = TRUE)  %dopar% {
    compSEM.HPD(j, S, alp)
  }
  names(sem.PPROB) <- colnames(sem.HPD) <- xvec

  dput(list(CDF = sem.PPROB, HPD = sem.HPD), 
       file = paste("compare-OCs/SEM-", ITER, ".txt", sep="") )

  return(ITER)
}

