
#' @export
summary.exchangeability_model <- function(object, ...) {
  if (!is.na(object$cluster)[1]) {
    ret <- list(
      call = object$call,
      basket = summary(object$basket),
      cluster = summary(object$cluster)
    )
  } else {
    ret <- list(
      call = object$call,
      basket = summary(object$basket),
      cluster = NA
    )
  }

  class(ret) <- "mem_summary"
  ret
}


#' @importFrom cli rule
#' @importFrom crayon bold
#' @export
print.mem_summary <- function(x, ...) {
  cat("\n")
  cat_line(rule(left = "The MEM Model Call", col = "bold"))
  cat("\n")
  print(x$call)
  cat("\n")
  cat_line(rule(left = "The Basket Summary", col = "bold"))
  print(x$basket)

  if (!is.na(x$cluster)[1]) {
    cat_line(rule(left = "The Cluster Summary", col = "bold"))
    print(x$cluster)
  }
  invisible(x)
}

make_basket_summary <- function(object) {
  mm <- rbind(object$mean_est, object$median_est)
  rownames(mm) <- c("Mean", "Median")
  ret <- list(
    null = object$p0, post_prob = object$post_prob,
    mm_resp = mm, hpd_signif = object$alpha, hpd = object$hpd,
    ess = object$ess,
    alternative = object$alternative
  )
  class(ret) <- "mem_basket_summary"
  ret
}


#' @export
summary.mem_basket <- function(object, ...) {
  make_basket_summary(object)
}

#' @importFrom cli cat_line
#' @importFrom crayon bold
#' @export
print.mem_basket_summary <- function(x, ...) {
  # Posterior Probability Response is Greater than the Null
  cat_line("\nThe Null Response Rates (alternative is ", x$alternative, "):")
  if (is.matrix(x$post_prob)) {
    a <- round(x$post_prob, 3)
    rownames(a) <- paste("Posterior for null of", rownames(a))
    print(a)
  } else {
    a <- round(rbind(x$null, x$post_prob), 3)
    rownames(a) <- c("Null", "Posterior Prob")
    print(a)
  }

  # The Mean and Median Response Rates
  cat_line("\nPosterior Mean and Median Response Rates:")
  print(round(x$mm_resp, 3))

  # The Highest Posterior Density Interval
  cat_line(
    "\nHighest Posterior Density Interval with Coverage Probability ",
    1 - x$hpd_signif,
    ":"
  )
  print(round(x$hpd, 3))

  # The Effective Sample Size
  cat_line("\nPosterior Effective Sample Size:")
  em <- matrix(x$ess, nrow = 1, dimnames = list("", names(x$ess)))
  print(round(em, 3))
  cat_line("")

  invisible(x)
}

make_cluster_summary <- function(object) {
  mm <- rbind(object$mean_est, object$median_est)
  rownames(mm) <- c("Mean", "Median")
  ret <- list(
    null = object$p0,
    post_prob = object$post_prob,
    mm_resp = mm,
    hpd_signif = object$alpha,
    hpd = object$hpd,
    ess = object$ess,
    alternative = object$alternative,
    name = object$name,
    cluster = object$cluster
  )
  class(ret) <- "mem_cluster_summary"
  ret
}

#' @export
summary.mem_cluster <- function(object, ...) {
  make_cluster_summary(object)
}


#' @importFrom crayon bold
#' @importFrom cli cat_line cat_print
#' @export
print.mem_cluster_summary <- function(x, ...) {
  cat("\n")
  for (i in seq_along(x$cluster)) {
    cat(x$name[i])
    pm <- matrix(x$cluster[[i]], nrow = 1)
    colnames(pm) <- rep("", ncol(pm))
    rownames(pm) <- rep("", nrow(pm))
    print(pm)
  }
  print.mem_basket_summary(x, ...)
}
