
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual facet_grid
#' xlab ylab theme_minimal
#' @importFrom tibble as_tibble
#' @export
plot.exchangeability_model<- function(x, y, ...) {
  dots <- list(...)
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
    if (length(levels(d$Basket)) == 1) {
      basket_colors <- rainbow(3)[1]
    }
    else if (length(levels(d$Basket)) == 2) {
      basket_colors <- rainbow(3)[-2]
    } else {
      basket_colors <- rainbow(length(levels(d$Basket)))
    }
  }
  ggplot(d, aes(x = Density, fill = Basket)) + 
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = basket_colors, guide = FALSE) +
    facet_grid( Basket ~ .) +
    xlab("") +
    ylab("Density") +
    theme_minimal()

}

