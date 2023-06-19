#enumerate all possible partitions of baskets in a basket trial 
#and calculate the uncertainties related to each partition

##x: vector of responses in each basket
##n: vector of sample sizes in each basket
##beta prior for response rate: a0, b0 [default: a0 = b0 = 1]
##max_cl: maximum number of clusters a basket trial should be partitioned to
##if the number of baskets is small, max_cl should be chosen as the number of baskets
##if the number of baskets is large (>10), max_cl should be kept low (2, 3) to speed up computation
# Output: table of all possible partitions (represented by group indexes)
# with posterior probability of each partition (in descending order)

#' @importFrom stats aggregate
#' @importFrom partitions listParts
basket_part <- function(x, n, a0 = 1, b0 = 1, 
                        max_cl = 6){
  ##number of baskets(arms)
  R <- length(x)
  ## generate all the possible partitions 
  parts <- listParts(R)
  ## number of clusters
  n_cl <- sapply(parts, length)
  parts <- parts[n_cl <= max_cl]
  n_cl <- sapply(parts, length)
  K <- length(parts)  # total number of partitions
  
  dat <- NULL
  p <- NULL
  for(k in 1:K){
    grp <- 1:R
    for(j in 1:n_cl[k]){
      idx <- parts[[k]][[j]]
      grp[idx] <- j
    }
    dat <- rbind(dat, grp)
    
    S <- aggregate(x, by=list(grp), sum)$x
    N <- aggregate(n, by=list(grp), sum)$x
    #in priori, not to pool and pool are equally likely
    if(length(unique(grp)) == R){  # not pool: each trial forms a basket
     prior_k <- 0.5
    }else{  # pool
     prior_k <- (1-0.5)/(K-1)
    }
    #calculate marginal probs m(s_j)
    p[k] <- prod((beta(a0+S, b0+N-S)/beta(a0,b0)))*prior_k
  }
  dat <- data.frame(dat)
  names(dat) <- LETTERS[1:R]
  rownames(dat) <- NULL
  #marginal probs m(s_j)
  mp <- sum(p)
  ## calculate posterior prob of each grouping structure
  pp_grp <- p/mp
  dat$pp_grp <- pp_grp
  dat <- dat[order(-pp_grp), ]
  
  return(dat)
}
