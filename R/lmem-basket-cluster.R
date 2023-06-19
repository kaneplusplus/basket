#Construct a similarity matrix based on the partition of baskets and 
#determine to pool or not to pool based on bayes factor

## partitions: results from running basket_cluster()
## a0, b0: beta priors for response rate
## x: vector of responses 
## n: vector of basket sizes

# output: list of
# similarity matrix
# parameters for posterior (beta) distribution: alpha, beta [vec]
# 

basket_cluster <- function(partitions, a0=1, b0=1, x, n, bf_thres=3.2){  
  R <- ncol(partitions) - 1
  ##the similarity/exchangeability matrix
  simi <- matrix(0, R, R)
  rownames(simi) <- colnames(simi) <- LETTERS[1:R]
  K <- nrow(partitions)
  for(k in 1:K){
    grp <- unlist(partitions[k, 1:R])
    mat <- NULL
    for(r in 1:R){
      mat <- rbind(mat, ifelse(grp[r]==grp, 1, 0))
    }
    simi <- simi + mat * partitions$pp_grp[k]
  }
  is_pool <- FALSE
  ##compute bayes factor of pool vs. not pool
  h1 <- which(apply(partitions[, 1:R], 1, 
                    function(x){ length(unique(x)) }) == R)
  bf <- (1-partitions$pp_grp[h1]) / (partitions$pp_grp[h1])
  
  if(bf > bf_thres) {  # pooling if bf > 3.2
    is_pool <- TRUE
    partitions <- partitions[-h1, ]  # exclude h1 (not pooling)
    # optimal partition with max pp
    member <- as.numeric(partitions[1, 1:R])  
  
    # local similarity matrix
    loc_s <- simi
    for(i in 1:R){
      for(j in 1:R){
       if(member[j]==member[i]){
          loc_s[i, j] <- simi[i, j]
        }else{
         loc_s[i, j] <- 0
        }
      }
    }
    #borrowing locally helps to control alpha, such
    #that a lonely basket without signal won't become rich
    #by excessive borrowing 
    avec <- a0 + loc_s %*% x
    bvec <- b0 + loc_s %*% (n-x)
  }
  else {
    avec <- a0 + x
    bvec <- b0 + (n-x)
    member <- 1:R
  }
  return(list("Smat"=simi, "alpha"=avec, "beta"=bvec, 
              "member"=member, "x"=x, "n"=n, "BF"=bf, "is_pool"=is_pool))
}
