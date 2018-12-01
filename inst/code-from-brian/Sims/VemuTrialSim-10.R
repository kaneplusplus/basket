
rm(list=ls())

### Seq Vermu Sim ###
source("../Utils/R-functions-v11.R")

########################################
## 1) Define Data and hyperparameters ##
########################################
Data <- dget("Vemu-data.txt") #list( X=round(runif(4,0,100)), N=rep(100,4) ) ## Vectors of Observed number of Responses (X) and Patients (N)
pars <- list(UB=1, 						## MMLE upper boundary
             LB=0) 						## MMLE lower boundary
p0 <- 0.15								## Null response rate for Posterior Probability Calculation

## EB optimal ##
 marg.M <- MEM_marginal(Data$X, Data$N, rep(0.5,length(Data$X)), rep(0.5,length(Data$X))); marg.M$maximizer
 round(MEM.cdf(Data,pars,p0,marg.M),3); round(post.ESS(Data,pars,marg.M),3); t(round(post.HPD(Data,pars,marg.M,0.05),3))

## constained EB optimal ##
 pars$UB <- 0.5
 marg.M <- MEM_marginal(Data$X, Data$N, rep(0.5,length(Data$X)), rep(0.5,length(Data$X))); marg.M$maximizer
 round(MEM.cdf(Data,pars,p0,marg.M),3); round(post.ESS(Data,pars,marg.M),3); t(round(post.HPD(Data,pars,marg.M,0.05),3))

## Indep ##
 pars$UB <- 0
 marg.M <- MEM_marginal(Data$X, Data$N, rep(0.5,length(Data$X)), rep(0.5,length(Data$X))); marg.M$maximizer
 round(MEM.cdf(Data,pars,p0,marg.M),3); round(post.ESS(Data,pars,marg.M),3); t(round(post.HPD(Data,pars,marg.M,0.05),3))

## Full Bayesian ##
 pr.Inclus <- matrix(0.5, length(Data$X), length(Data$X)); diag(pr.Inclus) <- 1
 allInf <- MEM_FullInf(Data$X, Data$N, rep(0.5,length(Data$X)), rep(0.5,length(Data$X)), pr.Inclus, 0.05)

## Full Bayesian PAR ##
 source("../R-functions-v11-par.R")
# source("runFB-Par.R")

## Berry Model Analysis ##
 library(rjags)
  data <- list(J=6, N=Data$N, Y=Data$X); inits <- list(mu=rep(0,6), M=0, inv.sig2=1) #; bugsData(data,"data.txt")
  m <- jags.model("../UtilsBerry-ig_model.txt", data, inits, n.chains=3, quiet=TRUE); update(m, 30000, progress.bar="none"); S <- coda.samples(m, "pr", n.iter=10000, progress.bar="none"); dput(data,"data.txt")
  PPROB <- c( sum( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) > p0 )/length( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) ),
              sum( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) > p0 )/length( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) ),
              sum( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) > p0 )/length( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) ),
              sum( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) > p0 )/length( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) ),
              sum( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) > p0 )/length( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) ),
              sum( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) > p0 )/length( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) ) )
 names(PPROB) <-  c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
 write.csv(PPROB,"Berry-result.csv")

 ########
 pars0 <- pars; pars0$UB=0; MEM.cdf(Data,pars0,p0,marg.M)
 pars1 <- pars; pars1$LB=pars1$UB=1; MEM.cdf(Data,pars1,p0,marg.M)

###################
## OC functions ###
###################
checkOC.MEM <- function(THRES, prob, N, pars, p0, L.Penalty=NA, type=NA){ 
 seed <- round(runif(1,1,1000000000)); set.seed(seed)
 D <- list(X=rbinom(N,N,prob),N=N)
 if(type=="opt"){ marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X))) 
 } else if(type=="PE"){ marg.M <- MEM_marginalPE(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=NA, Lambda=L.Penalty)
 } else{ marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$X))) }
 return( as.numeric( MEM.cdf(D,pars,p0,marg.M) > THRES ) ) 
}

checkOC.MEM.FB <- function(THRES, prob, N, pr.Inclus, p0){ 
 seed <- round(runif(1,1,1000000000)); set.seed(seed)
 D <- list(X=rbinom(N,N,prob),N=N)
 fb <- MEM_FullInf(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), pr.Inclus, 0.05)
 return( as.numeric( fb$CDF > THRES ) ) 
}

checkOC.SEM <- function(THRES, prob, N, p0){
 library(rjags)
 seed <- round(runif(1,1,1000000000)); set.seed(seed)
 D <- list(X=rbinom(N,N,prob),N=N)
  data <- list(J=6, N=D$N, Y=D$X); inits <- list(mu=rep(0,6), M=0, inv.sig2=1) #; bugsData(data,"data.txt")
  m <- jags.model("Berry-ig_model.txt", data, inits, n.chains=3, quiet=TRUE); update(m, 10000, progress.bar="none"); S <- coda.samples(m, "pr", n.iter=10000, progress.bar="none"); dput(data,"data.txt")
  PPROB <- c( sum( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) > p0 )/length( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) ),
              sum( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) > p0 )/length( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) ),
              sum( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) > p0 )/length( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) ),
              sum( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) > p0 )/length( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) ),
              sum( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) > p0 )/length( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) ),
              sum( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) > p0 )/length( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) ) )
 return( as.numeric( PPROB > THRES ) ) 
}

checkOC.MEM.fwer <- function(THRES, prob, N, pars, p0, L.Penalty=NA, type=NA){ 
 seed <- round(runif(1,1,1000000000)); set.seed(seed)
 D <- list(X=rbinom(N,N,prob),N=N)
 if(type=="opt"){ marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X))) 
 } else if(type=="PE"){ marg.M <- MEM_marginalPE(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=NA, Lambda=L.Penalty)
 } else{ marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$X))) }
 return( as.numeric( sum( as.numeric( MEM.cdf(D,pars,p0,marg.M) > THRES ))==1 ) ) 
}

gen.MEM.FB <- function(prob, N, pr.Inclus, p0){ 
 seed <- round(runif(1,1,1000000000)); set.seed(seed)
 D <- list(X=rbinom(N,N,prob),N=N)
 fb <- MEM_FullInf(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), pr.Inclus, 0.05)
 #return( as.numeric( sum( as.numeric( fb$CDF > THRES ))==1 ) ) 
 return(fb$CDF) 
}

gen.SEM <- function(prob, N, p0){
 library(rjags)
 seed <- round(runif(1,1,1000000000)); set.seed(seed)
 D <- list(X=rbinom(N,N,prob),N=N)
  data <- list(J=6, N=D$N, Y=D$X); inits <- list(mu=rep(0,6), M=0, inv.sig2=1) #; bugsData(data,"data.txt")
  m <- jags.model("Berry-ig_model.txt", data, inits, n.chains=3, quiet=TRUE); update(m, 10000, progress.bar="none"); S <- coda.samples(m, "pr", n.iter=10000, progress.bar="none"); dput(data,"data.txt")
  PPROB <- c( sum( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) > p0 )/length( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) ),
              sum( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) > p0 )/length( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) ),
              sum( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) > p0 )/length( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) ),
              sum( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) > p0 )/length( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) ),
              sum( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) > p0 )/length( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) ),
              sum( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) > p0 )/length( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) ) )
 #return( as.numeric( sum( as.numeric( PPROB > THRES ))==1 ) ) 
 return( PPROB ) 
}

checkESS <- function(prob, N, pars, L.Penalty=NA, type=NA){
 seed <- round(runif(1,1,1000000000)); set.seed(seed)
 D <- list(X=rbinom(N,N,prob),N=N)
 if(type=="opt"){ marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X))) 
 } else if(type=="PE"){ marg.M <- MEM_marginalPE(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=NA, Lambda=L.Penalty)
 } else{ marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$X))) }
 out <- post.ESS(D,pars,marg.M); names(out) <- prob
 return( out ) 
}

############
## SEM OCs #
############
nITER <- 100
for(i in 1:50){ T <- replicate(nITER, gen.SEM(prob=c(0.15,0.15,0.15,0.15,0.15,0.15), N=Data$N, p0=p0))
 write.csv(T, paste( "SEM-fwer-",i,".csv",sep="") )
}

############
## FB OCs ## 
############
 pr.Inclus <- matrix(0.5, length(Data$X), length(Data$X)); diag(pr.Inclus) <- 1
 nITER <- 20
 for(i in 1:50){ T <- replicate(nITER, gen.MEM.FB(prob=c(0.15,0.15,0.15,0.15,0.15,0.15), N=Data$N, pr.Inclus=pr.Inclus, p0=p0))
  write.csv(T, paste( "MEM-fwer-",i,".csv",sep="") )
 }


############
## SEM OCs #
############
nITER <- 100
for(i in 1:50){ T <- replicate(nITER, gen.SEM(prob=c(0.45,0.15,0.15,0.15,0.15,0.15), N=Data$N, p0=p0))
 write.csv(T, paste( "SEM-fwer-alt",i,".csv",sep="") )
}

############
## FB OCs ## 
############
 pr.Inclus <- matrix(0.5, length(Data$X), length(Data$X)); diag(pr.Inclus) <- 1
 nITER <- 20
 for(i in 1:50){ T <- replicate(nITER, gen.MEM.FB(prob=c(0.45,0.15,0.15,0.15,0.15,0.15), N=Data$N, pr.Inclus=pr.Inclus, p0=p0))
  write.csv(T, paste( "MEM-fwer-alt",i,".csv",sep="") )
 }




############
## EB OCs ##
############
 nITER <- 20000

 #nITER <- 3000;   PROB <- c(0.45,0.15,0.15,0.15,0.45,0.45)
 #pt <- seq(0.91,0.99,0.005); ub <- seq(0.05,0.5,0.02); R <- matrix(NA, length(pt), length(ub)); rownames(R)=pt; colnames(R)=ub
 #for(j in 1:dim(R)[1]){ for(i in 1:dim(R)[2]){
 # R[j,i] <- mean( (rowSums( replicate(nITER, checkOC.MEM(pt[j], prob=PROB, N=Data$N, pars=list(LB=0,UB=ub[i]), p0=p0, type="SA") ))/nITER)[2:4] )
 #}}

pt <- 0.97;  pars$UB <- 0.5

# MEM.w(Data,pars,marg.M); post.ESS(Data,pars,marg.M)
# MEM.w(Data,list(LB=0, UB=1),marg.M); post.ESS(Data,list(LB=0, UB=1),marg.M)

 rowSums( replicate(nITER, checkOC.MEM(pt, prob=c(0.45,0.15,0.15,0.15,0.45,0.45), N=Data$N, pars=pars, p0=p0, type="SA") ))/nITER
 rowSums( replicate(nITER, checkOC.MEM(pt, prob=c(0.45,0.45,0.45,0.45,0.45,0.45), N=Data$N, pars=pars, p0=p0, type="SA") ))/nITER
 rowSums( replicate(nITER, checkOC.MEM(pt, prob=c(0.15,0.15,0.15,0.15,0.15,0.15), N=Data$N, pars=pars, p0=p0, type="SA") ))/nITER

 sum( replicate(nITER, checkOC.MEM.fwer(pt, prob=c(0.15,0.15,0.15,0.15,0.15,0.15), N=Data$N, pars=pars, p0=p0, type="SA") ))/nITER

 ESSa <- replicate(3000, checkESS(prob=c(0.45,0.15,0.15,0.15,0.45,0.45), N=Data$N, pars=pars, type="SA") )
 dput(ESSa,"ESS-alt.txt"); rowMeans(ESSa)
 #0.80940 0.04985 0.03790 0.05865 0.72995 0.53320#
 #36.19100 32.73983 40.28117 31.93867 33.34833 29.79567

 rowSums( replicate(nITER, checkOC.MEM(pt, prob=c(0.15,0.15,0.15,0.15,0.15,0.15), N=Data$N, pars=pars, p0=p0, type="SA") ))/nITER
 ESSb <- replicate(3000, checkESS(prob=c(0.15,0.15,0.15,0.15,0.15,0.15), N=Data$N, pars=pars, type="SA") )
 dput(ESSb,"ESS-null.txt"); rowMeans(ESSb)
 #0.01065 0.01090 0.00980 0.01745 0.01000 0.01300
 #49.42000 44.37867 53.20200 43.10933 46.13417 42.73033

 ##################################
 ### Permutation Sequential Sim ###
 ##################################
 sim.MEM.trial <- function(DATA, blk, pars, p0, L.Penalty=NA, type="SA"){ 
  o <- sample(1:nrow(DATA),nrow(DATA))
  U <- DATA[o,] 

  STOP <- 0; K <- 6
  while( STOP==0 ){
   D <- list( X=c( sum( U[which(U[1:K,2]=="NSCLC"), 1] ), sum( U[which(U[1:K,2]=="CRC-V"), 1] ), sum( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  sum( U[which(U[1:K,2]=="BD"), 1] ), sum( U[which(U[1:K,2]=="ECD-LCH"), 1] ), sum( U[which(U[1:K,2]=="ATC"), 1] ) ),
              N=c( length( U[which(U[1:K,2]=="NSCLC"), 1] ), length( U[which(U[1:K,2]=="CRC-V"), 1] ), length( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  length( U[which(U[1:K,2]=="BD"), 1] ), length( U[which(U[1:K,2]=="ECD-LCH"), 1] ), length( U[which(U[1:K,2]=="ATC"), 1] ) ) )
   if( sum(D$N==0) < 1 ){ STOP <- 1 
   } else{ K <- K + 1 }
  }
  if(type=="opt"){ marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X))) 
  } else if(type=="PE"){ marg.M <- MEM_marginalPE(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N)), Lambda=L.Penalty)
  } else{ marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N))) }
 
  PPROB <- list(MEM.cdf(D,pars,p0,marg.M)); MEM <- list(marg.M$maximizer); STAGE <- K; Ess <- list(post.ESS(D,pars,marg.M)); J <- 1
  while( K < nrow(DATA) ){
   K <- min( K + blk, nrow(DATA) ); J <- J + 1
   D <- list( X=c( sum( U[which(U[1:K,2]=="NSCLC"), 1] ), sum( U[which(U[1:K,2]=="CRC-V"), 1] ), sum( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  sum( U[which(U[1:K,2]=="BD"), 1] ), sum( U[which(U[1:K,2]=="ECD-LCH"), 1] ), sum( U[which(U[1:K,2]=="ATC"), 1] ) ),
              N=c( length( U[which(U[1:K,2]=="NSCLC"), 1] ), length( U[which(U[1:K,2]=="CRC-V"), 1] ), length( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  length( U[which(U[1:K,2]=="BD"), 1] ), length( U[which(U[1:K,2]=="ECD-LCH"), 1] ), length( U[which(U[1:K,2]=="ATC"), 1] ) ) )
   if(type=="opt"){ marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X))) 
   } else if(type=="PE"){ marg.M <- MEM_marginalPE(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N)), Lambda=L.Penalty)
   } else{ marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N))) }
   PPROB[[J]] <- MEM.cdf(D,pars,p0,marg.M); MEM[[J]] <- marg.M$maximizer; STAGE <- c(STAGE,K); Ess[[J]] <- post.ESS(D,pars,marg.M)
  }

  return( list(PPROB=PPROB, MEM=MEM, STAGE=STAGE, ESS=Ess) )
 }

 USE <- as.data.frame( list( Y=c( rep(1,8),rep(0,11),  rep(0,10),  1,rep(0,25),  1,rep(0,7),  rep(1,6),rep(0,8),  1,1,rep(0,5) ),
              B=c(   rep("NSCLC",19), rep("CRC-V",10), rep("CRC-VC",26), rep("BD",8), rep("ECD-LCH",14), rep("ATC", 7) ) ) ); table(USE)

 #for(i in 1:50000){
 # R <- sim.MEM.trial(DATA=USE, blk=1, pars=list(LB=0, UB=0.5), p0=0.15, L.Penalty=NA, type="sa")
 # dput(R, paste("V-sim/R-",i,".txt",sep=""))
 #}

s <- 6:84; NReps <- rep(0,length(s))
Ess <- array(NA,c(6,length(s),100000))
FUT <- matrix(0,6,length(s)); EXC <- matrix(0,15,length(s)); colnames(FUT) <- names(NReps) <- colnames(EXC) <- s
 i=2; R <- dget(paste("V-sim/R-",i,".txt",sep=""))
 MEM <- R$MEM[[70]]; colnames(MEM) <- rownames(MEM) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
 TEMP <- MEM; for(i in 1:6){ for(j in 1:6){ TEMP[i,j] <- paste(colnames(MEM)[i]," - ",colnames(MEM)[j],sep="") }}
 temp <- upper.tri(TEMP); TEMP[temp]
rownames(EXC) <- TEMP[temp]

EXCp <- EXC; mESS <- Ess[,,1]

use <- c(1:100000) #[-c(3,4,5)]

for(ii in 1:length(use)){
 i <- use[ii]
 R <- dget(paste("V-sim/R-",i,".txt",sep=""))
 for(j in 1:length(R$STAGE)){
  u <- which(s==R$STAGE[j]); NReps[u] <- NReps[u] + 1
  EXC[,u] <- EXC[,u] + R$MEM[[j]][upper.tri(R$MEM[[j]])]
  Ess[,u,i] <- R$ESS[[j]]
 }
}

PPf <- 0.05; s <- 6:84; NReps. <- rep(0,length(s)); FUT. <- matrix(0,6,length(s)); colnames(FUT.) <- names(NReps.) <- s
for(ii in 1:length(use)){
 i <- use[ii]
 R <- dget(paste("V-sim/R-",i,".txt",sep=""))
 u <- which(s==R$STAGE[1]); NReps.[u:79] <- NReps.[u:79] + 1
 for(k in 1:6){
  j <- 1; GO <- 1
  while(GO){
   u <- which(s==R$STAGE[j])
   if( as.numeric( R$PPROB[[j]][k] < PPf )){
    FUT.[k,u:79] <- FUT.[k,u:79] + 1
    GO <- FALSE
   } else if( u==79 ){ GO <- FALSE 
   } else{ j <- j + 1
   }
  }
 }
}
rownames(FUT.) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")

## Early Stopping: pars$UB = 0.5 ##
FUT.[,78]/NReps.[78] 
# NSCLC   CRC.v   CRC.vc    BD    ED.LH    ATC   
#0.00664 0.79869 1.00000 0.16331 0.00576 0.01245

FUT <- FUT.; FUTp <-  matrix(0,6,length(s)); rownames(FUTp) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
for(j in 1:length(s)){  EXCp[,j] <- EXC[,j]/NReps[j]
 FUTp[,j] <- FUT[,j]/NReps[j]
 mESS[,j] <- rowMeans(Ess[,j,], na.rm = TRUE)
}; aESS <- rbind( mESS, colSums(mESS) )
rownames(aESS) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ","Total ") 

StageSum <- matrix(NA,6,6); rownames(StageSum) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
for(k in 1:6){ 
 StopFreq <- diff(FUT[k,])/diff(s) 
 SS <- rep(84, NReps[79] - sum( StopFreq ) )
 for(j in 2:78){ SS <- c(SS, rep(s[j], StopFreq[j-1]) ) }
 StageSum[k,] <- summary(SS)
}
colnames(StageSum) <- c("Min.", "1st Qu.",  "Median",    "Mean", "3rd Qu.",    "Max.")
write.csv(StageSum,"StageSum.csv")


library(fields)
BRK <- seq(min(mESS),max(mESS),length.out=65)
 Palette <- c("green","yellow","orange","red"); nbcol <- length(BRK)-1; cols <- colorRampPalette( Palette ); color <- cols(nbcol)
par(oma=c(9,1.7,0,3))
 image(6:84,1:6,as.matrix(t(mESS)), breaks=BRK, ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)
 for(i in 1:6){  mtext(rownames(aESS)[i],side=2,line=1,at=i,las=1,font=2)  }
 mtext("Basket",side=2,line=0.48,at=6.8,las=1,cex=1.4, font=2) 
 mtext("Enrollment Stage (n)", side=1, line=3.5, cex=2)
 mtext("Effective",side=4,line=-1.6,at=7.15,las=1,cex=1.4, font=2); mtext("Sample Size",side=4,line=-2.4,at=6.8,las=1,cex=1.4, font=2) 
 image.plot(range(c(mESS)), legend.only=TRUE, breaks=BRK, legend.mar=0,
  ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)

library(fields)
BRK <- seq(min(aESS[7,]),max(aESS[7,]),length.out=65)
 Palette <- c("navy","cyan","pink","red"); nbcol <- length(BRK)-1; cols <- colorRampPalette( Palette ); color <- cols(nbcol)
par(oma=c(15,1.7,0,3))
 image(6:84,1:2,as.matrix(t(rbind(aESS[7,],aESS[1,]))), breaks=BRK, ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)
 mtext("Total",side=2,line=1,at=1,las=1,font=2, cex=1.2)
 #mtext("Basket",side=2,line=0.48,at=6.8,las=1,cex=1.4, font=2) 
 #mtext("Enrollment Stage (n)", side=1, line=3.5, cex=2)
 mtext("Total Effective",side=4,line=-2,at=3,las=1,cex=1.4, font=2); mtext("Sample Size",side=4,line=-2.4,at=2.8,las=1,cex=1.4, font=2) 
 image.plot(range(c(aESS[7,])), legend.only=TRUE, breaks=BRK, legend.mar=0,
  ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)

library(fields)
BRK <- seq(0,1,length.out=65)
 Palette <- c("green","yellow","orange","red"); nbcol <- length(BRK)-1; cols <- colorRampPalette( Palette ); color <- cols(nbcol)
par(oma=c(1,1.7,0,3))
 image(6:84,1:6,as.matrix(t(FUTp)), breaks=BRK, ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)
 for(i in 1:6){  mtext(rownames(FUTp)[i],side=2,line=1,at=i,las=1,font=2)  }
 mtext("Basket",side=2,line=0.48,at=6.8,las=1,cex=1.4, font=2) 
 mtext("Enrollment Stage (n)", side=1, line=3.5, cex=2)
 mtext("Futility",side=4,line=-1,at=7.15,las=1,cex=1.4, font=2); mtext("Probability",side=4,line=-1.8,at=6.8,las=1,cex=1.4, font=2) 
 image.plot(range(c(FUTp)), legend.only=TRUE, breaks=BRK, legend.mar=0,
  ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)

win.metafile("Futility-plot-05.wmf")
par(oma=c(1,1.7,0,3))
 image(6:84,1:6,as.matrix(t(FUTp)), breaks=BRK, ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)
 for(i in 1:6){  mtext(rownames(FUTp)[i],side=2,line=1,at=i,las=1,font=2)  }
 mtext("Basket",side=2,line=0.48,at=6.8,las=1,cex=1.4, font=2) 
 mtext("Enrollment Stage (n)", side=1, line=3.5, cex=2)
 mtext("Futility",side=4,line=-1,at=7.15,las=1,cex=1.4, font=2); mtext("Probability",side=4,line=-1.8,at=6.8,las=1,cex=1.4, font=2) 
 image.plot(range(c(FUTp)), legend.only=TRUE, breaks=BRK, legend.mar=0,
  ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)
dev.off()

pdf("Futility-plot-05.pdf")
par(oma=c(1,1.7,0,3))
 image(6:84,1:6,as.matrix(t(FUTp)), breaks=BRK, ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)
 for(i in 1:6){  mtext(rownames(FUTp)[i],side=2,line=1,at=i,las=1,font=2)  }
 mtext("Basket",side=2,line=0.48,at=6.8,las=1,cex=1.4, font=2) 
 mtext("Enrollment Stage (n)", side=1, line=3.5, cex=2)
 mtext("Futility",side=4,line=-1,at=7.15,las=1,cex=1.4, font=2); mtext("Probability",side=4,line=-1.8,at=6.8,las=1,cex=1.4, font=2) 
 image.plot(range(c(FUTp)), legend.only=TRUE, breaks=BRK, legend.mar=0,
  ylab="",xlab="", xlim=c(6,84), col=color, cex.lab=1.9, cex.axis=1.8, yaxt="n", legend.cex=20.0)
dev.off()

##############################
## Sequential OC Simulation ##
##############################
 sim.MEM.trial2 <- function(PROB, blk, pars, p0, L.Penalty=NA, type="SA"){ 
  DATA <- as.data.frame( list( Y=c( rbinom(19,1,PROB[1]),  
                                    rbinom(10,1,PROB[2]),
                                    rbinom(26,1,PROB[3]), 
                                    rbinom(8, 1,PROB[4]), 
                                    rbinom(14,1,PROB[5]), 
                                    rbinom(7, 1,PROB[6]) ),
              B=c(   rep("NSCLC",19), rep("CRC-V",10), rep("CRC-VC",26), rep("BD",8), rep("ECD-LCH",14), rep("ATC", 7) ) ) ); table(DATA)

  o <- sample(1:nrow(DATA),nrow(DATA))
  U <- DATA[o,] 

  STOP <- 0; K <- 6
  while( STOP==0 ){
   D <- list( X=c( sum( U[which(U[1:K,2]=="NSCLC"), 1] ), sum( U[which(U[1:K,2]=="CRC-V"), 1] ), sum( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  sum( U[which(U[1:K,2]=="BD"), 1] ), sum( U[which(U[1:K,2]=="ECD-LCH"), 1] ), sum( U[which(U[1:K,2]=="ATC"), 1] ) ),
              N=c( length( U[which(U[1:K,2]=="NSCLC"), 1] ), length( U[which(U[1:K,2]=="CRC-V"), 1] ), length( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  length( U[which(U[1:K,2]=="BD"), 1] ), length( U[which(U[1:K,2]=="ECD-LCH"), 1] ), length( U[which(U[1:K,2]=="ATC"), 1] ) ) )
   if( sum(D$N<1) < 1 ){ STOP <- 1 
   } else{ K <- K + 1 }
  }
  if(type=="opt"){ marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X))) 
  } else if(type=="PE"){ marg.M <- MEM_marginalPE(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N)), Lambda=L.Penalty)
  } else{ marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N))) }
 
  PPROB <- list(MEM.cdf(D,pars,p0,marg.M)); MEM <- list(marg.M$maximizer); STAGE <- K; Ess <- list(post.ESS(D,pars,marg.M)); J <- 1
  while( K < nrow(DATA) ){
   K <- min( K + blk, nrow(DATA) ); J <- J + 1
   D <- list( X=c( sum( U[which(U[1:K,2]=="NSCLC"), 1] ), sum( U[which(U[1:K,2]=="CRC-V"), 1] ), sum( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  sum( U[which(U[1:K,2]=="BD"), 1] ), sum( U[which(U[1:K,2]=="ECD-LCH"), 1] ), sum( U[which(U[1:K,2]=="ATC"), 1] ) ),
              N=c( length( U[which(U[1:K,2]=="NSCLC"), 1] ), length( U[which(U[1:K,2]=="CRC-V"), 1] ), length( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  length( U[which(U[1:K,2]=="BD"), 1] ), length( U[which(U[1:K,2]=="ECD-LCH"), 1] ), length( U[which(U[1:K,2]=="ATC"), 1] ) ) )
   if(type=="opt"){ marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X))) 
   } else if(type=="PE"){ marg.M <- MEM_marginalPE(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N)), Lambda=L.Penalty)
   } else{ marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$N))) }
   PPROB[[J]] <- MEM.cdf(D,pars,p0,marg.M); MEM[[J]] <- marg.M$maximizer; STAGE <- c(STAGE,K); Ess[[J]] <- post.ESS(D,pars,marg.M)
  }

  return( list(PPROB=PPROB, MEM=MEM, STAGE=STAGE, ESS=Ess) )
 }

 #library(BRugs)
 library(rjags)
 sim.SEM.trial2 <- function(PROB, blk, p0){ 
  DATA <- as.data.frame( list( Y=c( rbinom(19,1,PROB[1]),  
                                    rbinom(10,1,PROB[2]),
                                    rbinom(26,1,PROB[3]), 
                                    rbinom(8, 1,PROB[4]), 
                                    rbinom(14,1,PROB[5]), 
                                    rbinom(7, 1,PROB[6]) ),
              B=c(   rep("NSCLC",19), rep("CRC-V",10), rep("CRC-VC",26), rep("BD",8), rep("ECD-LCH",14), rep("ATC", 7) ) ) ); table(DATA)

  o <- sample(1:nrow(DATA),nrow(DATA))
  U <- DATA[o,] 

  STOP <- 0; K <- 6
  while( STOP==0 ){
   D <- list( X=c( sum( U[which(U[1:K,2]=="NSCLC"), 1] ), sum( U[which(U[1:K,2]=="CRC-V"), 1] ), sum( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  sum( U[which(U[1:K,2]=="BD"), 1] ), sum( U[which(U[1:K,2]=="ECD-LCH"), 1] ), sum( U[which(U[1:K,2]=="ATC"), 1] ) ),
              N=c( length( U[which(U[1:K,2]=="NSCLC"), 1] ), length( U[which(U[1:K,2]=="CRC-V"), 1] ), length( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  length( U[which(U[1:K,2]=="BD"), 1] ), length( U[which(U[1:K,2]=="ECD-LCH"), 1] ), length( U[which(U[1:K,2]=="ATC"), 1] ) ) )
   if( sum(D$N<1) < 1 ){ STOP <- 1 
   } else{ K <- K + 1 }
  }

  data <- list(J=6, N=D$N, Y=D$X); inits <- list(mu=rep(0,6), M=0, inv.sig2=1) #; bugsData(data,"data.txt")
  #ig.fit <- BRugsFit(modelFile="Berry-ig_model.txt",data=data,inits=list(inits,inits,inits),numChains=3,para=c("pr"),nBurnin=3000,nIter=3000,DIC = FALSE, BRugsVerbose=FALSE)$Stats
  m <- jags.model("Berry-ig_model.txt", data, inits, n.chains=3, quiet=TRUE); update(m, 3000, progress.bar="none"); S <- coda.samples(m, "pr", n.iter=3000, progress.bar="none"); dput(data,"data.txt")
  PPROB <- list(c( sum( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) > p0 )/length( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) ),
              sum( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) > p0 )/length( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) ),
              sum( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) > p0 )/length( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) ),
              sum( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) > p0 )/length( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) ),
              sum( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) > p0 )/length( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) ),
              sum( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) > p0 )/length( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) ) ))
  STAGE <- K; J <- 1

  while( K < nrow(DATA) ){
   K <- min( K + blk, nrow(DATA) ); J <- J + 1
   D <- list( X=c( sum( U[which(U[1:K,2]=="NSCLC"), 1] ), sum( U[which(U[1:K,2]=="CRC-V"), 1] ), sum( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  sum( U[which(U[1:K,2]=="BD"), 1] ), sum( U[which(U[1:K,2]=="ECD-LCH"), 1] ), sum( U[which(U[1:K,2]=="ATC"), 1] ) ),
              N=c( length( U[which(U[1:K,2]=="NSCLC"), 1] ), length( U[which(U[1:K,2]=="CRC-V"), 1] ), length( U[which(U[1:K,2]=="CRC-VC"), 1] ),
                  length( U[which(U[1:K,2]=="BD"), 1] ), length( U[which(U[1:K,2]=="ECD-LCH"), 1] ), length( U[which(U[1:K,2]=="ATC"), 1] ) ) )
   data <- list(J=6, N=D$N, Y=D$X); inits <- list(mu=rep(0,6), M=0, inv.sig2=1) #; bugsData(data,"data.txt")
   #ig.fit <- BRugsFit(modelFile="Berry-ig_model.txt",data=data,inits=list(inits,inits,inits),numChains=3,para=c("pr"),nBurnin=3000,nIter=3000,DIC = FALSE, BRugsVerbose=FALSE)$Stats
   m <- jags.model("Berry-ig_model.txt", data, inits, n.chains=3, quiet=TRUE); update(m, 3000, progress.bar="none"); S <- coda.samples(m, "pr", n.iter=3000, progress.bar="none"); dput(data,"data.txt")
   PPROB[[J]] <- list(c( sum( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) > p0 )/length( c(S[[1]][,1],S[[2]][,1],S[[3]][,1]) ),
              sum( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) > p0 )/length( c(S[[1]][,2],S[[2]][,2],S[[3]][,2]) ),
              sum( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) > p0 )/length( c(S[[1]][,3],S[[2]][,3],S[[3]][,3]) ),
              sum( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) > p0 )/length( c(S[[1]][,4],S[[2]][,4],S[[3]][,4]) ),
              sum( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) > p0 )/length( c(S[[1]][,5],S[[2]][,5],S[[3]][,5]) ),
              sum( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) > p0 )/length( c(S[[1]][,6],S[[2]][,6],S[[3]][,6]) ) ))
   STAGE <- c(STAGE,K)
  }

  return( list(PPROB=PPROB, STAGE=STAGE) )
 }

 for(i in 1:1000){
  R <- sim.MEM.trial2(PROB=c(0.45,0.15,0.15,0.15,0.45,0.45), blk=10, pars=list(LB=0, UB=0.5), p0=0.15, L.Penalty=NA, type="sa")
   dput(R, paste("OC-sim/R-mem3-",i,".txt",sep=""))

  R <- sim.SEM.trial2(PROB=c(0.45,0.15,0.15,0.15,0.45,0.45), blk=10, p0=0.15)
   dput(R, paste("OC-sim/R-sem3-",i,".txt",sep=""))

  R <- sim.MEM.trial2(PROB=c(0.45,0.15,0.15,0.15,0.45,0.15), blk=10, pars=list(LB=0, UB=0.5), p0=0.15, L.Penalty=NA, type="sa")
   dput(R, paste("OC-sim/R-mem2-",i,".txt",sep=""))

  R <- sim.SEM.trial2(PROB=c(0.45,0.15,0.15,0.15,0.45,0.15), blk=10, p0=0.15)
   dput(R, paste("OC-sim/R-sem2-",i,".txt",sep=""))

  R <- sim.MEM.trial2(PROB=c(0.45,0.15,0.15,0.15,0.15,0.15), blk=10, pars=list(LB=0, UB=0.5), p0=0.15, L.Penalty=NA, type="sa")
   dput(R, paste("OC-sim/R-mem1-",i,".txt",sep=""))

  R <- sim.SEM.trial2(PROB=c(0.45,0.15,0.15,0.15,0.15,0.15), blk=10, p0=0.15)
   dput(R, paste("OC-sim/R-sem1-",i,".txt",sep=""))
 }


########################
## Poolability Sim #####
########################
library(ggplot2)
library(reshape2)
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

evalPool.bask <- function(prob,N){
  D <- list(X=rbinom(N,N,prob),N=N)
  #marg.M <- MEM_marginal(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)))
  marg.M <- MEM_marginalLoop(D$X, D$N, rep(0.5,length(D$X)), rep(0.5,length(D$X)), INITIAL=diag(1,length(D$X)))
  return( marg.M$maximizer )
}

################################
################################
### PEP, PRIOR, MAP, and MAX plots ####

## PEP ##
temp <- dget("PEP-05.txt"); rownames(temp) <- colnames(temp) <- 1:6

#diag(temp) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
upper_tri <- get_upper_tri(round(temp,2))
melted_cormat <- melt(upper_tri, na.rm = TRUE)

win.metafile("PEP2.wmf")
#pdf("PEP.pdf") 
#png("PEP.png",width=960,height=720)
 ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Exchangeable Freq") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  coord_fixed()  + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 8) +
  labs(x="",y="") +
  ggtitle("") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=25, hjust=0.5)) +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("1","2","3","4","5","6")) +
  theme(
    axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22), #element_text("Basket 1","Basket 2","Basket 3","Basket 4","Basket 5","Basket 6"), #element_blank(),
    axis.title.y = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()

## PRIOR ##
temp <- matrix(0.5,6,6); diag(temp) <- 1; rownames(temp) <- colnames(temp) <- 1:6

#diag(temp) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
upper_tri <- get_upper_tri(round(temp,2))
melted_cormat <- melt(upper_tri, na.rm = TRUE)

#win.metafile("Prior.wmf")
win.metafile("Prior.wmf") 
 ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Exchangeable Freq") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  coord_fixed()  + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 8, fontface = "bold") +
  labs(x="",y="") +
  ggtitle("") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=25, hjust=0.5)) +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("1","2","3","4","5","6")) +
  theme(
    axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22), #element_text("Basket 1","Basket 2","Basket 3","Basket 4","Basket 5","Basket 6"), #element_blank(),
    axis.title.y = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()

## MAP + MAX ##
temp <- dget("MAP-MAX.txt"); rownames(temp) <- colnames(temp) <- 1:6

#diag(temp) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
upper_tri <- get_upper_tri(round(temp,2))
melted_cormat <- melt(upper_tri, na.rm = TRUE)

win.metafile("MAP-MAX.wmf") 
 ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Exchangeable Freq") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  coord_fixed()  + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 8, fontface = "bold") +
  labs(x="",y="") +
  ggtitle("") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=25, hjust=0.5)) +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("1","2","3","4","5","6")) +
  theme(
    axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22), #element_text("Basket 1","Basket 2","Basket 3","Basket 4","Basket 5","Basket 6"), #element_blank(),
    axis.title.y = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal") +
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()


#### Simulation Based ####
##########################
### with observed prop ###
ITER=10000
R <- replicate(n=ITER, evalPool.bask( prob=Data$X/Data$N, N=Data$N ) )
 temp <- matrix(NA, 6, 6); diag(temp) <- 1
 for(i in 1:5){ for(j in (i+1):6 ){
  temp[i,j] <- sum( R[i,j,] )/ITER
 } }
dput(temp,"R_ObsProp.txt")

temp <- dget("R_ObsProp.txt")
#diag(temp) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
upper_tri <- get_upper_tri(round(temp,2))
melted_cormat <- melt(upper_tri, na.rm = TRUE)

win.metafile("Pool_ObsProp.wmf") 
 ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Exchangeable Freq") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  coord_fixed()  + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 7) +
  labs(x="Basket",y="Basket") +
  ggtitle("") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=25, hjust=0.5)) +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("1","2","3","4","5","6")) +
  theme(
    axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22), #element_text("Basket 1","Basket 2","Basket 3","Basket 4","Basket 5","Basket 6"), #element_blank(),
    axis.title.y = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()

pdf("temp3.pdf")
 ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Exchangeable Freq") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  coord_fixed()  + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 7) +
  labs(x="Basket",y="Basket") +
  ggtitle("") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=25, hjust=0.5)) +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("1","2","3","4","5","6")) +
  theme(
    axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22), #element_text("Basket 1","Basket 2","Basket 3","Basket 4","Basket 5","Basket 6"), #element_blank(),
    axis.title.y = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()

#######################################
#######################################
PROB <- c(0.4,0.15,0.15,0.15,0.4,0.4) #; evalPool.bask( prob=PROB, N=Data$N )

### with 0.15 and 0.4 ###
ITER=10000
R <- replicate(n=ITER, evalPool.bask( prob=PROB, N=Data$N ) )
 temp <- matrix(NA, 6, 6); diag(temp) <- 1
 for(i in 1:5){ for(j in (i+1):6 ){
  temp[i,j] <- sum( R[i,j,] )/ITER
 } }
dput(temp,"R_015v04.txt")

#diag(temp) <- c("NSCLC ","CRC.v ","CRC.vc","  BD  ","ED.LH "," ATC  ")
upper_tri <- get_upper_tri(round(temp,2))
melted_cormat <- melt(upper_tri, na.rm = TRUE)

win.metafile("Pool_015v04.wmf") 
 ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "black", high = "red", mid = "orange", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Exchangeable Freq") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0, 
                                   size = 22, hjust = 0, color="#666666"))+
  coord_fixed()  + 
  geom_text(aes(Var2, Var1, label = value), color = "white", size = 7) +
  labs(x="Basket",y="Basket") +
  ggtitle("") +
  theme(plot.title = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=25, hjust=0.5)) +
  scale_x_discrete(labels=c("1","2","3","4","5","6")) +
  scale_y_discrete(labels=c("1","2","3","4","5","6")) +
  theme(
    axis.title.x = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22), #element_text("Basket 1","Basket 2","Basket 3","Basket 4","Basket 5","Basket 6"), #element_blank(),
    axis.title.y = element_text(family = "Trebuchet MS", color="#666666", face="bold", size=22),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()




