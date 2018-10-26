
SumR <- function(REPS,path){
 setwd(path)
 completed=list.files(pattern="MEM-"); comp <- rep(0,REPS)
 for(j in 1:REPS){ comp[j] <- as.numeric(sum( completed==paste("MEM-",j,".txt",sep="") )==1) }; sum(comp==0)
 u <- which(comp==1); i=u[1]
 R1 <- dget(paste("MEM-",i,".txt",sep="")); R2 <- dget(paste("SEM-",i,".txt",sep=""))
 PP1 <- R1$CDF;   PP2 <- R2$CDF
 L1 <- R1$HPD[1,]; L2 <- R2$HPD[1,]
 U1 <- R1$HPD[2,]; U2 <- R2$HPD[2,]
 for(i in 2:length(u)){ 
  R1 <- dget(paste("MEM-",u[i],".txt",sep="")); R2 <- dget(paste("SEM-",u[i],".txt",sep=""))
  PP1 <- rbind(PP1, R1$CDF);    PP2 <- rbind(PP2, R2$CDF)
  L1 <- rbind(L1,   R1$HPD[1,]); L2 <- rbind(L2,  R2$HPD[1,])
  U1 <- rbind(U1, R1$HPD[2,]);    U2 <- rbind(U2, R2$HPD[2,])
 }
 return( list(PP1=PP1, PP2=PP2, L1=L1, L2=L2, U1=U1, U2=U2) )
}

## Indicator that CI excludes Null value ##  
Reject <- function(qnts,null){
 try1 <- qnts[1]>null
 try2 <- qnts[2]<null
 temp <- sum(try1,try2)
 if(temp < 2){
  out <- temp
 }
 else{
  out <- NA
  print("NA in Rejection")
 }
 return(out)
}


REPS <- 5000

#######################################################################
#R <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-null")
# setwd("D:/R/Proj/Stat/a-MEM/OneSample-sim")
# dput( R, "global-null.txt")
#
#R <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-alt1")
# setwd("D:/R/Proj/Stat/a-MEM/OneSample-sim")
# dput( R, "alt-1.txt")
#
#Ra <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-alt2a")
#Rb <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-alt2b")
# RR <- list( PP1=rbind(Ra$PP1,Rb$PP1), PP2=rbind(Ra$PP2,Rb$PP2), 
#              L1=rbind(Ra$L1,Rb$L1), 
#              L2=rbind(Ra$L2,Rb$L2), 
#              U1=rbind(Ra$U1,Rb$U1), 
#              U2=rbind(Ra$U2,Rb$U2) )
#  setwd("D:/R/Proj/Stat/a-MEM/OneSample-sim")
#  dput( RR, "alt-2.txt")
#
#R <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-alt3")
# setwd("D:/R/Proj/Stat/a-MEM/OneSample-sim")
# dput( R, "alt-3.txt")
#
#R <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-alt5")
# setwd("D:/R/Proj/Stat/a-MEM/OneSample-sim")
# dput( R, "alt-5.txt")
#
#
#Ra <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-alt6a")
#Rb <- SumR(REPS, "D:/R/Proj/Stat/a-MEM/OneSample-sim/compare-OCs-alt6b")
# RR <- list( PP1=rbind(Ra$PP1,Rb$PP1), PP2=rbind(Ra$PP2,Rb$PP2), 
#              L1=rbind(Ra$L1,Rb$L1), 
#              L2=rbind(Ra$L2,Rb$L2), 
#              U1=rbind(Ra$U1,Rb$U1), 
#              U2=rbind(Ra$U2,Rb$U2) )
#  setwd("D:/R/Proj/Stat/a-MEM/OneSample-sim")
#  dput( RR, "alt-6.txt")
##############################################

applyPost <- function(thres,X){ return( sum(X>thres)/length(X) ) }
applyFWEPost <- function(thres,X){ return( sum(sapply(X, FUN=chk.FWE, thres))/nrow(X) ) }
chk.FWE <- function(X,thres){ return( as.numeric( sum( as.numeric(X>thres) )>0 ) ) }

chk.Cov <- function(L,U,tar,K){
out <- rep(NA,K)
if(is.matrix(L)){ 
 for(i in 1:K){
    X <- cbind(L[,i],U[,i])
   out[i] <- sum( 1-apply( X=X, MARGIN=1, FUN=Reject, null=tar) )/nrow(X)
  }
 } else{
  X <- cbind(L,U)
  out <- sum( 1-apply( X=X, MARGIN=1, FUN=Reject, null=tar) )/nrow(X)
 }
 return(out)
}

chk.Wid <- function(L,U){
if(is.matrix(L)){ 
 out <- list()
 for(i in 1:ncol(U)){
  out[[i]] <- summary(U[,i]-L[,i])
 }
} else{
  out <- summary(U-L)
 }
 return(out)
}

plot.Int <- function(L,U,tar){
 hist(L,prob=TRUE,breaks=200,main="",xlab="",cex.axis=1.5,xlim=c(0,1)); abline(v=tar,lwd=2,col=2)
 hist(U,prob=TRUE,breaks=200,main="",xlab="",cex.axis=1.5,add=TRUE); abline(v=tar,lwd=2,col=2)
}

pThres.vec <- c(0, seq(0.7, 0.9, length.out=6), seq(0.92, 0.98, length.out=7), 0.99, 0.999, 0.9999, 0.99999, 1)

tab <- matrix(NA,4,4); rownames(tab) <- c("MEM-TypeIError","SEM-TypeIError","MEM-FWER","SEM-FWER")
thrs <- matrix(NA,4,4); rownames(thrs) <- c("MEM-TypeIError","SEM-TypeIError","MEM-FWER","SEM-FWER")

library(zoo)

### Posterior Prob Boxplots ###
library(ggplot2)

YLIM=c(0,1);   H <- dget("global-null.txt"); source("boxplot.R"); p
YLIM=c(0.65,1); H <- dget("alt-1.txt"); source("boxplot.R"); p
YLIM=c(0.75,1); H <- dget("alt-2.txt"); source("boxplot.R"); p
YLIM=c(0.8,1); H <- dget("alt-3.txt"); source("boxplot.R"); p
YLIM=c(0.9,1); H <- dget("alt-6.txt"); source("boxplot.R"); p
#############################################################
#############################################################


par(mfrow=c(2,4),mar=c(4,1,2,2),oma=c(2,4,4,1)); df.=9

########
H0 <- dget("global-null.txt"); H1 <- dget("alt-1.txt") 
 TyIEr.1 <- sapply( pThres.vec, FUN=applyPost, c(H0$PP1)); TyIEr.2 <-  sapply( pThres.vec, FUN=applyPost, c(H0$PP2))
 Power.1 <- sapply( pThres.vec, FUN=applyPost, H1$PP1[,1]); Power.2 <- sapply( pThres.vec, FUN=applyPost, H1$PP2[,1])
 source("plotROC.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[1:2,1] <- c(prd1[101],prd2[101]) 
 Cov.1 <- chk.Cov( c(H0$L1), c(H0$U1), 0.15 ); Cov.2 <- chk.Cov( c(H0$L2), c(H0$U2), 0.15 )
 Wid.1 <- chk.Wid( c(H0$L1), c(H0$U1) ); Wid.2 <- chk.Wid( c(H0$L2), c(H0$U2) )

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[1:2,1] <- c(prd1[101],prd2[101]) 

#########
H0 <- dget("global-null.txt"); H1 <- dget("alt-2.txt"); apply(H1$PP1,MARGIN=2,FUN=summary); v <- c(1,5)
 TyIEr.1 <- sapply( pThres.vec, FUN=applyPost, c(H0$PP1)); TyIEr.2 <-  sapply( pThres.vec, FUN=applyPost, c(H0$PP2))
 Power.1 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP1[,v])); Power.2 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP2[,v]))
 source("plotROC2.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[1:2,2] <- c(prd1[101],prd2[101]) 

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[1:2,2] <- c(prd1[101],prd2[101]) 

########
H0 <- dget("global-null.txt"); H1 <- dget("alt-3.txt"); apply(H1$PP1,MARGIN=2,FUN=summary); v <- c(1,5,6)
 TyIEr.1 <- sapply( pThres.vec, FUN=applyPost, c(H0$PP1)); TyIEr.2 <-  sapply( pThres.vec, FUN=applyPost, c(H0$PP2))
 Power.1 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP1[,v])); Power.2 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP2[,v]))
 source("plotROC2.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[1:2,3] <- c(prd1[101],prd2[101]) 

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[1:2,3] <- c(prd1[101],prd2[101]) 

########
df.=11
H0 <- dget("global-null.txt"); H1 <- dget("alt-6.txt"); apply(H1$PP1,MARGIN=2,FUN=summary); v <- c(1:6)
 TyIEr.1 <- sapply( pThres.vec, FUN=applyPost, c(H0$PP1)); TyIEr.2 <-  sapply( pThres.vec, FUN=applyPost, c(H0$PP2))
 Power.1 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP1[,v])); Power.2 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP2[,v]))
 source("plotROC2.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[1:2,4] <- c(prd1[101],prd2[101]) 

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[1:2,4] <- c(prd1[101],prd2[101]) 

##########
### FWE ##
##########
H0 <- dget("global-null.txt"); H1 <- dget("alt-1.txt") 
 TyIEr.1 <- sapply( pThres.vec, FUN=applyFWEPost, H0$PP1); TyIEr.2 <-  sapply( pThres.vec, FUN=applyFWEPost, H0$PP2)
 Power.1 <- sapply( pThres.vec, FUN=applyPost, H1$PP1[,1]); Power.2 <- sapply( pThres.vec, FUN=applyPost, H1$PP2[,1])
 source("plotROCf.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[3:4,1] <- c(prd1[101],prd2[101]) 

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[3:4,1] <- c(prd1[101],prd2[101]) 

#########
H0 <- dget("global-null.txt"); H1 <- dget("alt-2.txt"); apply(H1$PP1,MARGIN=2,FUN=summary); v <- c(1,5)
 TyIEr.1 <- sapply( pThres.vec, FUN=applyFWEPost, H0$PP1); TyIEr.2 <-  sapply( pThres.vec, FUN=applyFWEPost, H0$PP2)
 Power.1 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP1[,v])); Power.2 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP2[,v]))
 source("plotROC2f.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[3:4,2] <- c(prd1[101],prd2[101]) 

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[3:4,2] <- c(prd1[101],prd2[101]) 

########
H0 <- dget("global-null.txt"); H1 <- dget("alt-3.txt"); apply(H1$PP1,MARGIN=2,FUN=summary); v <- c(1,5,6)
 TyIEr.1 <- sapply( pThres.vec, FUN=applyFWEPost, H0$PP1); TyIEr.2 <-  sapply( pThres.vec, FUN=applyFWEPost, H0$PP2)
 Power.1 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP1[,v])); Power.2 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP2[,v]))
 source("plotROC2f.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[3:4,3] <- c(prd1[101],prd2[101]) 

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[3:4,3] <- c(prd1[101],prd2[101]) 

########
H0 <- dget("global-null.txt"); H1 <- dget("alt-6.txt"); apply(H1$PP1,MARGIN=2,FUN=summary); v <- c(1:6)
 TyIEr.1 <- sapply( pThres.vec, FUN=applyFWEPost, H0$PP1); TyIEr.2 <-  sapply( pThres.vec, FUN=applyFWEPost, H0$PP2)
 Power.1 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP1[,v])); Power.2 <- sapply( pThres.vec, FUN=applyPost, c(H1$PP2[,v]))
 source("plotROC2f.R"); points(TyIEr.1,Power.1,cex=2,pch=16); points(TyIEr.2,Power.2,cex=2,pch=17)
 tab[3:4,4] <- c(prd1[101],prd2[101]) 

 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,pThres.vec,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,pThres.vec,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 thrs[3:4,4] <- c(prd1[101],prd2[101]) 

Cov.tab <- matrix(NA,2,5); rownames(Cov.tab) <- c("Null","Alt")
R <- list(MEM.cov=Cov.tab, MEM.wid=Cov.tab, SEM.cov=Cov.tab, SEM.wid=Cov.tab)

########
H0 <- dget("global-null.txt")
 Cov.1 <- chk.Cov( c(H0$L1), c(H0$U1), 0.15 ); Cov.2 <- chk.Cov( c(H0$L2), c(H0$U2), 0.15 )
 Wid.1 <- chk.Wid( c(H0$L1), c(H0$U1) ); Wid.2 <- chk.Wid( c(H0$L2), c(H0$U2) )
 R$MEM.cov[1,1] <- Cov.1; R$SEM.cov[1,1] <- Cov.2
 R$MEM.wid[1,1] <- Wid.1[3]; R$SEM.wid[1,1] <- Wid.2[3]

j <- 2; H1 <- dget("alt-1.txt"); v0 <- c(2:6); v1=1 
 Cov.1 <-   chk.Cov( c(H1$L1[,v0]), c(H1$U1[,v0]), 0.15 )
 Cov.2 <-   chk.Cov( c(H1$L2[,v0]), c(H1$U2[,v0]), 0.15 )
 Wid.1 <-   chk.Wid( c(H1$L1[,v0]), c(H1$U1[,v0]) )
 Wid.2 <-   chk.Wid( c(H1$L2[,v0]), c(H1$U2[,v0]) )
 R$MEM.cov[1,j] <- Cov.1; R$SEM.cov[1,j] <- Cov.2
 R$MEM.wid[1,j] <- Wid.1[3]; R$SEM.wid[1,j] <- Wid.2[3]

 Cov.1 <-   chk.Cov( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 Cov.2 <-   chk.Cov( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 Wid.1 <-   chk.Wid( c(H1$L1[,v1]), c(H1$U1[,v1]) )
 Wid.2 <-   chk.Wid( c(H1$L2[,v1]), c(H1$U2[,v1]) )
 R$MEM.cov[2,j] <- Cov.1; R$SEM.cov[2,j] <- Cov.2
 R$MEM.wid[2,j] <- Wid.1[3]; R$SEM.wid[2,j] <- Wid.2[3]
 
j <- 3; H1 <- dget("alt-2.txt"); v0 <- c(2:4,6); v1=c(1,5) 
 Cov.1 <-   chk.Cov( c(H1$L1[,v0]), c(H1$U1[,v0]), 0.15 )
 Cov.2 <-   chk.Cov( c(H1$L2[,v0]), c(H1$U2[,v0]), 0.15 )
 Wid.1 <-   chk.Wid( c(H1$L1[,v0]), c(H1$U1[,v0]) )
 Wid.2 <-   chk.Wid( c(H1$L2[,v0]), c(H1$U2[,v0]) )
 R$MEM.cov[1,j] <- Cov.1; R$SEM.cov[1,j] <- Cov.2
 R$MEM.wid[1,j] <- Wid.1[3]; R$SEM.wid[1,j] <- Wid.2[3]

 Cov.1 <-   chk.Cov( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 Cov.2 <-   chk.Cov( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 Wid.1 <-   chk.Wid( c(H1$L1[,v1]), c(H1$U1[,v1]) )
 Wid.2 <-   chk.Wid( c(H1$L2[,v1]), c(H1$U2[,v1]) )
 R$MEM.cov[2,j] <- Cov.1; R$SEM.cov[2,j] <- Cov.2
 R$MEM.wid[2,j] <- Wid.1[3]; R$SEM.wid[2,j] <- Wid.2[3]

j <- 4; H1 <- dget("alt-3.txt"); v0 <- c(2:4); v1=c(1,5,6) 
 Cov.1 <-   chk.Cov( c(H1$L1[,v0]), c(H1$U1[,v0]), 0.15 )
 Cov.2 <-   chk.Cov( c(H1$L2[,v0]), c(H1$U2[,v0]), 0.15 )
 Wid.1 <-   chk.Wid( c(H1$L1[,v0]), c(H1$U1[,v0]) )
 Wid.2 <-   chk.Wid( c(H1$L2[,v0]), c(H1$U2[,v0]) )
 R$MEM.cov[1,j] <- Cov.1; R$SEM.cov[1,j] <- Cov.2
 R$MEM.wid[1,j] <- Wid.1[3]; R$SEM.wid[1,j] <- Wid.2[3]

 Cov.1 <-   chk.Cov( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 Cov.2 <-   chk.Cov( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 Wid.1 <-   chk.Wid( c(H1$L1[,v1]), c(H1$U1[,v1]) )
 Wid.2 <-   chk.Wid( c(H1$L2[,v1]), c(H1$U2[,v1]) )
 R$MEM.cov[2,j] <- Cov.1; R$SEM.cov[2,j] <- Cov.2
 R$MEM.wid[2,j] <- Wid.1[3]; R$SEM.wid[2,j] <- Wid.2[3]

j <- 5; H1 <- dget("alt-6.txt"); v0 <- NA; v1=c(1:6) 
 Cov.1 <-   chk.Cov( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 Cov.2 <-   chk.Cov( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 Wid.1 <-   chk.Wid( c(H1$L1[,v1]), c(H1$U1[,v1]) )
 Wid.2 <-   chk.Wid( c(H1$L2[,v1]), c(H1$U2[,v1]) )
 R$MEM.cov[2,j] <- Cov.1; R$SEM.cov[2,j] <- Cov.2
 R$MEM.wid[2,j] <- Wid.1[3]; R$SEM.wid[2,j] <- Wid.2[3]

dput(R,"Coverage-Width.txt")

### MEM ###
par(mfrow=c(2,5),mar=c(4,1,2,2),oma=c(2,8,4,1)); df.=9

H0 <- dget("global-null.txt")
 plot.Int( c(H0$L1), c(H0$U1), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 
 mtext("density", side=2, line=2.6, cex=1.1); 

H1 <- dget("alt-1.txt"); v0 <- c(2:6); v1=1 
 plot.Int( c(H1$L1[,v0]), c(H1$U1[,v0]), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-2.txt"); v0 <- c(2:4,6); v1=c(1,5) 
 plot.Int( c(H1$L1[,v0]), c(H1$U1[,v0]), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-3.txt"); v0 <- c(2:4); v1=c(1,5,6) 
 plot.Int( c(H1$L1[,v0]), c(H1$U1[,v0]), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

plot.new(); plot.new()

H1 <- dget("alt-1.txt"); v0 <- c(2:6); v1=1 
 plot.Int( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 
 mtext("density", side=2, line=2.6, cex=1.1); 

H1 <- dget("alt-2.txt"); v0 <- c(2:4,6); v1=c(1,5) 
 plot.Int( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-3.txt"); v0 <- c(2:4); v1=c(1,5,6) 
 plot.Int( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-6.txt"); v0 <- NA; v1=c(1:6) 
 plot.Int( c(H1$L1[,v1]), c(H1$U1[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 


### SEM ###
par(mfrow=c(2,5),mar=c(4,1,2,2),oma=c(2,8,4,1)); df.=9

H0 <- dget("global-null.txt")
 plot.Int( c(H0$L2), c(H0$U2), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 
 mtext("density", side=2, line=2.6, cex=1.1); 

H1 <- dget("alt-1.txt"); v0 <- c(2:6); v1=1 
 plot.Int( c(H1$L2[,v0]), c(H1$U2[,v0]), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-2.txt"); v0 <- c(2:4,6); v1=c(1,5) 
 plot.Int( c(H1$L2[,v0]), c(H1$U2[,v0]), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-3.txt"); v0 <- c(2:4); v1=c(1,5,6) 
 plot.Int( c(H1$L2[,v0]), c(H1$U2[,v0]), 0.15 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

plot.new(); plot.new()

H1 <- dget("alt-1.txt"); v0 <- c(2:6); v1=1 
 plot.Int( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 
 mtext("density", side=2, line=2.6, cex=1.1); 

H1 <- dget("alt-2.txt"); v0 <- c(2:4,6); v1=c(1,5) 
 plot.Int( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-3.txt"); v0 <- c(2:4); v1=c(1,5,6) 
 plot.Int( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

H1 <- dget("alt-6.txt"); v0 <- NA; v1=c(1:6) 
 plot.Int( c(H1$L2[,v1]), c(H1$U2[,v1]), 0.45 )
 mtext("HPD Interval Estimators", side=1, line=3, cex=1.2) 

