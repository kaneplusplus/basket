# function to summarize posterior distribution for each basket

##test: object returned from running basket_cluster()
##q: defalt is to produce 95% HPD
##p0: the null response rate, can be a vector in practice
##basket_name: name for each basket
##basket_member: cluster membership of each basket

#' @importFrom pscl betaHPD
summary_basket <- function(test, q=0.95, p0,  
                          basket_name, basket_member = NULL){
  tab <- NULL
  for(i in seq_along(basket_name)) {
    med <- qbeta(p = 0.5, test$alpha[i], test$beta[i])
    alpha_hpd = test$alpha[i]+0.0001
    beta_hpd = test$beta[i]
    if (test$alpha[i] <= 1) {
      warning("alpha is less than one, using 1 instead.")
      test$alpha[i] = 1.0001
    }
    if (test$beta[i] <= 1) {
      warning("beta is less than one, using 1 instead.")
      test$beta[i] = 1.0001
    }
    ci <- betaHPD(test$alpha[i], test$beta[i], p = q)
    #ess <- sum(test$n*test$Smat[i, ])
    #ess <- sum(test$n*test$loc[i, ])
    ess <- test$alpha[i] + test$beta[i]
    pp <- 1-pbeta(q = p0, test$alpha[i], test$beta[i])
    tab <- rbind(tab, c(med, ci, ess, pp))
  }
  post_summary <- data.frame("Basket"=basket_name,
                             "Median"=tab[, 1],
                             "L2.5%"=tab[, 2], "U97.5%"=tab[, 3],
                             "ESS"=tab[, 4], "pp"=tab[, 5],
                             "Member"=basket_member)
  rownames(post_summary) <- NULL
  post_summary
}
