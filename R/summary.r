
summary.mem_exact <- function(object, ...) {
  # object <- MHResult1$cluster
  if (class(object$post.prob) == "matrix") {
    cdfS <- c()
    for (i in 1:dim(object$post.prob)[1]) {
      ss <- paste("Post.Prob ", rownames(object$post.prob)[i])
      cdfS <- c(cdfS, ss)
    }
  } else {
    cdfS <- "Post.Prob"
  }
  oct <- rbind(object$post.prob, object$HPD, object$ESS, object$mean_est, object$median_est)
  rownames(oct) <- c(cdfS, "HPD LB", "HPD HB", "ESS", "Mean", "Median")
  return(oct)
}

make_em_summary <- function(x) {
  ret <- list(basket = summary(x$basket), cluster = summary(x$cluster))
  class(ret) <- "mem_summary"
  ret
}

summary.exchangeability_model <- function(object, ...) {
  ret <- list(basket = summary(object$basket), 
              cluster = summary(object$cluster))
  class(ret) <- "mem_summary"
  ret
}


#' @importFrom cli rule
#' @importFrom crayon bold
print.mem_summary <- function(x, ...) {
  cat_line(rule(left = "The Basket Summary", col = "bold"))
  print(x$basket)
  cat_line(rule(left = "The Cluster Summary", col = "bold"))
  print(x$basket)
  invisible(x)
}

make_basket_summary <- function(object) {
   mm <- rbind(object$mean_est, object$median_est)
  rownames(mm) <- c("Mean", "Median")

  ret <- list(null = object$p0, post_prob = object$post.prob, 
              mm_resp = mm, hpd_signif = object$alpha, hpd = object$HPD, 
              ess = object$ESS)
  class(ret) <- "mem_basket_summary"
  ret
}

summary.mem_basket <- function(object, ...) {
  make_basket_summary(object)
}

#' @importFrom cli cat_line
#' @importFrom crayon bold
print.mem_basket_summary <- function(x, ...) {

  # Posterior Probability Response is Greater than the Null
  cat_line("\nPosterior Probability Reponse is ", 
           "Greater than Null of ", paste(x$null, collapse = ", "), 
           ":", col = "bold")
  print(round(x$post_prob, 3))

  # The Mean and Median Response Rates
  cat_line("\nMean and Median Response Rates:", col = "bold")
  print(round(x$mm_resp, 3))

  # The Highest Posterior Density Interval
  cat_line("\nHigest Posterior Density Interval with Significance of ", 
           x$hpd_signif, ":", col = "bold")
  print(round(x$hpd, 3))
  
  # The Effective Sample Size
  cat_line("\nThe Effective Sample Size:", col = bold)
  print(round(x$ess, 3))
  cat_line("")

  invisible(x)
}

make_cluster_summary <- function(object) {
   mm <- rbind(object$mean_est, object$median_est)
  rownames(mm) <- c("Mean", "Median")

  ret <- list(null = object$p0, post_prob = object$post.prob, 
              mm_resp = mm, hpd_signif = object$alpha, hpd = object$HPD, 
              ess = object$ESS)
  class(ret) <- "mem_basket_summary"
  ret
}

summary.mem_cluster <- function(object, ...) {
  make_cluster_summary(object)
}


#' @importFrom crayon bold
print.mem_basket_summary <- function(x, ...) {
  # Posterior Probability Response is Greater than the Null
  cat_line("\nPosterior Probability Reponse is ", 
           "Greater than Null of ", paste(x$null, collapse = ", "), 
           ":", col = "bold")
  print(round(x$post_prob, 3))

  # The Mean and Median Response Rates
  cat_line("\nMean and Median Response Rates:", col = "bold")
  print(round(x$mm_resp, 3))

  # The Highest Posterior Density Interval
  cat_line("\nHigest Posterior Density Interval with Significance of ", 
           x$hpd_signif, ":", col = "bold")
  print(round(x$hpd, 3))
  
  # The Effective Sample Size
  cat_line("\nThe Effective Sample Size:", col = bold)
  print(round(x$ess, 3))
  cat_line("")

  invisible(x)
}

