

#' @importFrom tidyr gather
#' @importFrom tibble as_tibble
#' @export
plot_density <- function(x, ...) {
  UseMethod("plot_density")
}

#' @export
plot_density.default <- function(x, ...) {
  stop(paste("Don't know how to make a density plot with an object of type",
             class(x)))
}

#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual facet_grid
#' xlab ylab theme_minimal xlim geom_vline labeller label_wrap_gen
#' @export
plot_density.exchangeability_model <- function(x, ...) {
  dots <- list(...)
  Basket <- Density <- NULL
  d <- gather(as_tibble(x$samples), key = Basket, value = Density)

  if ("basket_levels" %in% names(dots)) {
    basket_levels <- dots$basket_levels
    d$Basket <- factor(d$Basket, levels = basket_levels)
  } else {
    d$Basket <- as.factor(d$Basket)
  }

  if ("basket_colors" %in% names(dots)) {
    basket_colors <- dots$basket_colors
  } else {
    basket_colors <- rep("black", length(x$responses))
#    if (length(levels(d$Basket)) == 1) {
#      basket_colors <- rainbow(3)[1]
#    }
#    else if (length(levels(d$Basket)) == 2) {
#      basket_colors <- rainbow(3)[-2]
#    } else {
#      basket_colors <- rainbow(length(levels(d$Basket)))
#    }
  }
  ggplot(d, aes(x = Density, fill = Basket)) + 
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = basket_colors, guide = FALSE) +
    facet_grid( Basket ~ ., labeller=label_wrap_gen(width = 10)) +
    xlab("") +
    ylab("Density") +
    xlim(0, 1) +
    geom_vline(xintercept = x$p0) +
    theme_minimal()
    
}

#' @export
plot_posterior_exchangeability <- function(x, ...) {
  UseMethod("plot_posterior_exchangeability")
}

#' @export
plot_posterior_exchangeability.default <- function(x, ...) {
  stop(paste("Don't know how to make a posterior exchangeability plot",
             "with an object of type", class(x)))
}

#' @importFrom GGally ggcorr
#' @export
plot_posterior_exchangeability.full_bayes <- function(x, ...) {
  mat <- x$PEP
  mat[is.na(mat)] <- 0
  diag(mat) <- diag(mat) / 2
  mat <- mat + t(mat)
  dimnames(mat) <- list(x$name, x$name)
  ggcorr(data = NULL, cor_matrix = mat, midpoint = 0.5, limits = c(0, 1), ...)
}

#' @export
plot_exchangeability <- function(x, ...) {
  UseMethod("plot_exchangeability")
}

#' @export
plot_exchangeability.default <- function(x, ...) {
  stop(paste("Don't know how to make a posterior exchangeability plot",
             "with an object of type", class(x)))
}

#' @importFrom GGally ggcorr
#' @importFrom stats runif
#' @export
plot_exchangeability.exchangeability_model <- function(x, ...) {
  mat <- x$maximizer
  dimnames(mat) <- list(x$name, x$name)
  if (sum(mat == 1) > 0) {
    mat[mat == 1] <- runif(sum(mat == 1), 0.99, 1)
  }
  if (sum(mat == 0) > 0) {
    mat[mat == 0] <- runif(sum(mat == 0), 0, 0.01)
  }
  ggcorr(data = NULL, cor_matrix = mat, midpoint = 0.5, limits = c(0, 1), ...)
}


