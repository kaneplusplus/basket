
#############################################################
### Added Shape and Direction to Posterior Prob #############
#############################################################
eval.Post <- function(p0, X, N, Omega, w, a, b, alternative="greater"){
  alph <- a + Omega %*% X
  beta <- b + (Omega %*% N - Omega %*% X)
  if(alternative=="greater"){ out <- sum( (1-pbeta(p0, alph, beta))*w )
  } else{ out <- sum( (pbeta(p0, alph, beta))*w ) }
  return( out )
}

#############################################################
### Fixed shape paramter specification in gen.Post() ########
#############################################################
gen.Post <- function(X,N,Omega,a,b){ 
  alph <- a + Omega%*%X
  beta <- b + (Omega%*%N-Omega%*%X)
  return( rbeta(1, alph, beta) )
}

samp.Post <- function(X,N,Omega,w,a,b){ return( gen.Post(X,N,Omega[ which( rmultinom(1,1,w)==1 ), ], a, b) ) }

sample_posterior.full_bayes <- function(model, num_samples = 10000) {
  ret <- replicate(num_samples, samp.Post(model$responses, model$size,
                       model$models, model$pweights[[1]], model$shape1[1], model$shape2[1]))
  K <- length(model$responses)
  ret <- rbind(ret, 
    foreach(j = 2:(K-1), .combine = rbind) %do% {
      Ii <- c(j,1:(j-1),(j+1):K)
      replicate(num_samples, 
                samp.Post(model$responses[Ii], model$size[Ii], model$models, 
                          model$pweights[[j]], model$shape1[j], model$shape2[j]))
    })
  j <- K
  Ii <- c(j,1:(j-1))
  ret <- rbind(ret, 
    replicate(num_samples, 
              samp.Post(model$responses[Ii], model$size[Ii], model$models, 
                        model$pweights[[j]], model$shape1[j], model$shape2[j])))
  ret <- t(ret)
  dimnames(ret) <- list(NULL, model$name)
  ret
}

  mem.PostProb <- function(model, method="samples", fit){  ## Incorporate direction ##
   if(method=="mixture"){
    out <- eval.Post(model$p0[1], model$responses, model$size, model$models, model$pweights[[1]], model$shape1[1], model$shape2[1], alternative=model$alternative)
    K <- length(model$responses)
    for (j in 2:(K - 1)) {
        Ii <- c(j, 1:(j - 1), (j + 1):K)
        out <- c(out, eval.Post(model$p0[j], model$responses[Ii], model$size[Ii], model$models, model$pweights[[j]], model$shape1[j], model$shape2[j], alternative=model$alternative))
    }
    j <- j + 1
    Ii <- c(j, 1:(j - 1))
    out <- c(out, eval.Post(model$p0[j], model$responses[Ii], model$size[Ii], model$models, model$pweights[[j]], model$shape1[j], model$shape2[j], alternative=model$alternative))
   } else{
     if( model$alternative=="greater" ){ out <- sapply(1:ncol(fit$samples), FUN=function(j,x,t){ return(sum(x[,j]>t[j])/length(x[,j])) }, x=fit$samples, t=model$p0) 
     } else{ out <- sapply(1:ncol(fit$samples), FUN=function(j,x,t){ return(sum(x[,j]<t[j])/length(x[,j])) }, x=fit$samples, t=model$p0) }
   }
   names(out) <- model$name
   return(out)
  }

#############################################################
#### Check Log Marginal Density Functions ###################
#############################################################
logMarg.Dens <- function(I, mod.mat, xvec, nvec, avec, bvec){
  M <- MEM.mat(I, mod.mat, length(xvec))
  marg.vec <- rep(NA, dim(M)[1])
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  for(i in 1:dim(M)[1]){ p.vec <- prod( prod.vec^(1-M[i,])  )
                         marg.vec[i] <- (beta(avec[i] + M[i,]%*%xvec, bvec[i] + M[i,]%*%(nvec-xvec)) / beta(avec[i],bvec[i]) ) * p.vec
  } #return( list(mod = M, marg = prod(marg.vec)) )
  return( sum(log(marg.vec)) )
}

logMarg.DensSA <- function(M, mod.mat, xvec, nvec, avec, bvec){
  #M <- MEM.mat(I, mod.mat, length(xvec))
  marg.vec <- rep(NA, dim(M)[1])
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  for(i in 1:dim(M)[1]){ p.vec <- prod( prod.vec^(1-M[i,])  )
                         marg.vec[i] <- (beta(avec[i] + M[i,]%*%xvec, bvec[i] + M[i,]%*%(nvec-xvec)) / beta(avec[i],bvec[i]) ) * p.vec
  } #return( list(mod = M, marg = prod(marg.vec)) )
  return( sum(log(marg.vec)) )
}

logMarg.DensPE <- function(m, mod.size, xvec, nvec, avec, bvec, Lambda){
  #M <- MEM.mat(I, mod.mat, length(xvec))
  M <- diag(1, mod.size)
  u <- upper.tri(M); M[u] <- m; M. <- t(M); l <- lower.tri(M.); M[l] <- M.[l]
  marg.vec <- rep(NA, dim(M)[1])
  prod.vec <- beta(xvec + avec, nvec + bvec - xvec) / beta(avec, bvec) #calculate the product portion of integrated marginal likelihood
  for(i in 1:dim(M)[1]){ p.vec <- prod( prod.vec^(1-M[i,])  )
                         marg.vec[i] <- (beta(avec[i] + M[i,]%*%xvec, bvec[i] + M[i,]%*%(nvec-xvec)) / beta(avec[i],bvec[i]) ) * p.vec
  } #return( list(mod = M, marg = prod(marg.vec)) )
  u <- upper.tri(M)
  out <- sum(log(marg.vec)) - Lambda*sum(M[u])
  return(-out)
}
####################################################################
####################################################################
####################################################################

####################################################################
### Incorporate "method" argument for marginal max search ##########
####################################################################
#method="optimal" # default is optimal for baskets size <6 or 7; then use Loop by default #

# if(method=="optimal"){ marg.M <- MEM_marginal(Data$X, Data$N, rep(0.5,length(Data$X)), rep(0.5,length(Data$X))) 
# } else if(method=="PE"){ marg.M <- MEM_marginalPE(Data$X, Data$N, rep(0.5,length(Data$X)), rep(0.5,length(Data$X)), INITIAL=NA, Lambda=L.Penalty)
# } else{ marg.M <- MEM_marginalLoop(Data$X, Data$N, rep(0.5,length(Data$X)), rep(0.5,length(Data$X)), INITIAL=diag(1,length(Data$X))) }
####################################################################
####################################################################
####################################################################


####################################################################
#### Alternative Computation of ESS using HPD interval matching ####
####################################################################
#fit <- mem_full_bayes(responses = Data$X, size = Data$N, name = Data$Basket, p0 = 0.25, shape1 = 0.5, shape2 = 0.5, alpha = 0.05, call = NULL)

euc.dist <- function(x1, x2, w=c(1,1)){ 
 if( sum(is.na(x1))>1 ){ out <- Inf 
 } else{ out <- sqrt(sum( w*((x1 - x2) ^ 2) ) ) }
 return(out)
}

dist.beta.HPD <- function(ess, fit, alpha, jj){ return( euc.dist(fit$HPD[,jj], qbeta(c(alpha/2,1-alpha/2), fit$mean_est[jj]*ess, ess*(1-fit$mean_est[jj]))) ) }
dist.beta.HPDwid <- function(ess, fit, alpha, jj){ return( abs( (fit$HPD[2,jj]-fit$HPD[1,jj])-(diff( qbeta(c(alpha/2,1-alpha/2), fit$median_est[jj]*ess, ess*(1-fit$median_est[jj])))) ) ) }

ESS.from.HPD.i <- function(jj, fit, alpha){ library(GenSA)
 opt <- GenSA::GenSA(par=1, fn=dist.beta.HPD, lower=0, upper=10000000, #control=list(maxit=pars$DTW.maxit), 
                    fit=fit, alpha=alpha, jj=jj); return(opt$par) }

ESS.from.HPDwid.i <- function(jj, fit, alpha){ library(GenSA)
 opt <- GenSA::GenSA(par=1, fn=dist.beta.HPDwid, lower=0, upper=10000000, #control=list(maxit=pars$DTW.maxit), 
                    fit=fit, alpha=alpha, jj=jj); return(opt$par) }

calc.ESS.from.HPD <- function(fit, alpha){ ## fit is list with median vec and HPD vec ##
 return( sapply(1:length(fit$mean_est), FUN=ESS.from.HPD.i, fit=fit, alpha=alpha) ) }

calc.ESS.from.HPDwid <- function(fit, alpha){ ## fit is list with median vec and HPD vec ##
 return( sapply(1:length(fit$median_est), FUN=ESS.from.HPDwid.i, fit=fit, alpha=alpha) ) }

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

 mem.Prior.Mat <- function(M, mod.mat, pr.Inclus){
  #M <- MEM.mat(I, mod.mat, nrow(pr.Inclus)); 
  mem <- M[ upper.tri(M) ]
  source.vec <- pr.Inclus[upper.tri(pr.Inclus)]
  s.in <- source.vec #source inclusion probability
  s.ex <- 1-source.vec #source exclusion probability
  return( prod( s.in^(mem) * s.ex^(1-mem) ) )
 }

  update.MH <- function(MOld, M, xvec, nvec, avec, bvec, mod.mat, Prior.EP){
   K <- max(M, na.rm=TRUE)+1
   v <- sample(1:(K-1), (1:(K-1))[which( rmultinom(1, 1, (((K-1):1)^3)/sum( ((K-1):1)^3 )) == 1 )])
  
   MProp <- flip.MEM(v[1], MOld, M)
   if( length(v)>1 ){ for(ii in 2:length(v)){ MProp <- flip.MEM(v[ii], MProp, M) }}

   rho <- exp( logMarg.DensSA(MProp, mod.mat, xvec, nvec, avec, bvec) + log( mem.Prior.Mat(MProp, mod.mat, Prior.EP) ) - 
   ( logMarg.DensSA(MOld, mod.mat, xvec, nvec, avec, bvec) +  log( mem.Prior.Mat(MOld, mod.mat, Prior.EP) ) ) )

   if(rho>=1){ out <- MProp } else{ if(rbinom(1,1,rho)==1){ out <- MProp } else{ out <- MOld } }

   return(out)
  }

  I.models <- function(hh, models, Samp){
   K <- length(models[1,])
   if(hh==1){ Ii <- 1:K
   } else if(hh==K){ Ii <- c(K, 1:(K - 1))
   } else{ Ii <- c(hh, 1:(hh - 1), (hh + 1):K) }
   return( which( apply(models, MARGIN=1, FUN=function(x,t){ sum(x==t) }, t=Samp[hh,Ii]) == K ) )
  }

  models.Count <- function(Samp, models){
   out <- matrix(0, nrow(models), ncol(models))
   u <- sapply(1:length(models[1,]), FUN=I.models, models=models, Samp=Samp)
   for(i in 1:length(u)){ out[u[i],i] <- 1 }
   return(out)
  }

 boa.hpd <- function(x, alpha){
  n <- length(x)
  m <- max(1, ceiling(alpha * n))
  y <- sort(x)
  a <- y[1:m]
  b <- y[(n - m + 1):n]
  i <- order(b - a)[1]
  structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
 }

 mem_full_bayes_mcmc <- function(responses, size, name, null = 0.15, shape1 = 0.5, shape2 = 0.5, 
        Prior = diag(length(responses))/2 + matrix(0.5, nrow = length(responses), ncol = length(responses)),
        HPD.alpha = 0.05, alternative = "greater", niter.MCMC = 10000, Initial = NA, call = NULL){

    if (is.null(getDoParName())) {
        registerDoSEQ()
    }
    if (length(responses) != length(size)) {
        stop(paste("The length of the responses and size parameters", 
            "must be equal."))
    }
    if (length(shape1) == 1) {
        shape1 <- rep(shape1, length(responses))
    }
    if (length(shape2) == 1) {
        shape2 <- rep(shape2, length(responses))
    }
    if (length(null) == 1) {
        null <- rep(null, length(responses))
    }

  ### Produce sample space of MEM ###
  mod.mat <- list(); k <- length(responses)-1; j <- 1
  while(k > 0){ mod.mat[[j]] <- ( as.matrix(expand.grid( rep(list(c(0,1)), k) )) )[order(rowSums(as.matrix(expand.grid( rep(list(c(0,1)), k) )))),] 
                k <- k - 1; j <- j + 1 }

  ### Set initial MEM for MCMC ###
  if( is.na(Initial)[1] ){ M.init <- diag(1,length(responses))
  } else{ M.init <- Initial }
  
  MOld <- M.init
 
  ### Check symmetry of Initial MEM ###
  d <- length(responses)
  update.Init <- FALSE
  for( i in 1:(d-1))
    for(k in (i+1):d)
    {
      if (MOld[i,k]!=MOld[k,i])
      {
        MOld[i,k] <- MOld[k,i] <- rbinom(1,1,0.5)
        update.Init <- TRUE
      }
    }; if(update.Init){ M.init <- MOld }

  ### Create Map for Proposal Distribution ###
  M <- diag(NA, nrow(MOld)); K <- 1
  for(ii in 1:(nrow(MOld)-1)){ for(jj in (ii+1):ncol(MOld)){
   M[ii,jj] <- M[jj,ii] <- K; K <- K + 1
  }}

  ### Implement Metropolis-Hastings Alg ###
  n.chg <- 0
  models <- cbind(rep(1, dim(mod.mat[[1]])[1]), mod.mat[[1]])
  mweights <- matrix(0, nrow(models), length(responses)); colnames(mweights) <- name
  mem.Samp <- list( MOld ); mweights <- mweights + models.Count(Samp=mem.Samp[[1]], models=models)
  MAP.list <- list(mem.Samp[[1]]); MAP.count <- c(1)
  mem.Samp[[2]] <- update.MH(MOld, M, responses, size, shape1, shape2, mod.mat, Prior)
  mweights <- mweights + models.Count(Samp=mem.Samp[[2]], models=models)
  Samp.Sum <- mem.Samp[[1]] + mem.Samp[[2]]
  if( sum( mem.Samp[[2]] == mem.Samp[[1]] ) < length(mem.Samp[[2]]) ){ n.chg <- n.chg + 1 }
  i.Map <- which( unlist( lapply(MAP.list, FUN=function(old, new){ return( sum( new == old ) == length(new) ) }, new=mem.Samp[[2]]) ) )
  if( length(i.Map) > 0 ){ MAP.count[i.Map] <- MAP.count[i.Map] + 1 
  } else{ MAP.list[[length(MAP.list)+1]] <- mem.Samp[[2]]; MAP.count <- c(MAP.count,1) }
  for(KK in 3:niter.MCMC){
   mem.Samp[[KK]] <- update.MH(mem.Samp[[KK-1]], M, responses, size, shape1, shape2, mod.mat, Prior)
   mweights <- mweights + models.Count(Samp=mem.Samp[[KK]], models=models)
   Samp.Sum <- Samp.Sum + mem.Samp[[KK]]
   if( sum( mem.Samp[[KK]] == mem.Samp[[KK-1]] ) < length(mem.Samp[[KK-1]]) ){ n.chg <- n.chg + 1 }
   i.Map <- which( unlist( lapply(MAP.list, FUN=function(old, new){ return( sum( new == old ) == length(new) ) }, new=mem.Samp[[KK]]) ) )
   if( length(i.Map) > 0 ){ MAP.count[i.Map] <- MAP.count[i.Map] + 1 
   } else{ MAP.list[[length(MAP.list)+1]] <- mem.Samp[[KK]]; MAP.count <- c(MAP.count,1) }
  }

  ### Compute Posterior Model Weights ###
  pweights <- list(); for(KK in 1:ncol(mweights)){ pweights[[KK]] <- mweights[,KK]/niter.MCMC }

  ### List for post-processing ###
  MODEL <- list(responses=responses, size=size, name=name, shape1=shape1, shape2=shape2, models=models, 
                pweights=pweights, p0=null, alpha=HPD.alpha, alternative=alternative)

  ### Compute and output results ###
  PEP <- Samp.Sum/niter.MCMC; rownames(PEP) <- colnames(PEP) <- MODEL$name
  MAP <- MAP.list[[order(MAP.count, decreasing=TRUE)[1]]]; rownames(MAP) <- colnames(MAP) <- MODEL$name
 
  if (is.null(call)) {
      call <- match.call()
  }

  ret <- list(responses = responses, size = size, name = name, null = null, alternative=alternative,
              shape1 = shape1, shape2 = shape2, Prior = Prior, call = call)
  ret$mod.mat <- mod.mat
  ret$Inital <- M.init; rownames(ret$Inital) <- colnames(ret$Inital) <- MODEL$name
  ret$models <- models
  ret$pweights <- pweights
  ret$samples <- sample_posterior.full_bayes(MODEL)
  ret$accept.rate <- (n.chg)/niter.MCMC
  ret$mean_est <- colMeans(ret$samples)
  ret$median_est <- apply(ret$samples, 2, median)
  ret$PEP <- PEP
  ret$MAP <- MAP
  ret$HPD <- apply(ret$samples, MARGIN=2, FUN=boa.hpd, alpha=MODEL$alpha)
  ret$ESS <- calc.ESS.from.HPD(fit=ret, alpha=MODEL$alpha); names(ret$ESS) <- MODEL$name
  ret$ESS2 <- calc.ESS.from.HPDwid(fit=ret, alpha=MODEL$alpha); names(ret$ESS2) <- MODEL$name
  ret$PostProb <- mem.PostProb(MODEL, method="samples", fit=ret)
  class(ret) <- c("full_bayes", "exchangeability_model")
  return(ret)
 }

