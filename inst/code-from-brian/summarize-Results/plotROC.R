 plot(-10,-10,xlab="", ylab="", main="", xlim=c(0,0.5), ylim=c(0,1), cex.axis=1.4)
 mtext("Power", side=2, line=2.8, cex=1.6); mtext("Type I Error", side=1, line=3, cex=1.6) 
 s <- seq(0,1,length.out=1000)
 sm <- smooth.spline(TyIEr.1,Power.1,df=df.); prd1 <- predict(sm,x=s)$y; prd1[prd1<0] <- 0; prd1[prd1>1] <- 1
 sm <- smooth.spline(TyIEr.2,Power.2,df=df.); prd2 <- predict(sm,x=s)$y; prd2[prd2<0] <- 0; prd2[prd2>1] <- 1
 for(i in 1:(length(prd1)-1)){ if( prd1[i+1]<prd1[i] ){ prd1[i+1] <- prd1[i] } }
 for(i in 1:(length(prd2)-1)){ if( prd2[i+1]<prd2[i] ){ prd2[i+1] <- prd2[i] } }
 AUC1 <- sum(diff(s[order(s[1:200])][1:200])*rollmean(prd1[order(s[1:200])][1:200],2))
 AUC2 <- sum(diff(s[order(s[1:200])][1:200])*rollmean(prd2[order(s[1:200])][1:200],2)) 
 lines(s, prd1, lwd=2, col=1, lty=1)
 lines(s, prd2, lwd=2, col=1, lty=3)
 legend("right",c(
                    paste("MEM: AUC = ",round(AUC1,2),sep=""),
                    paste(" SEM: AUC = ",round(AUC2,2),sep="")),lwd=c(2,2),lty=c(1,3),cex=1.6,bty="n")
