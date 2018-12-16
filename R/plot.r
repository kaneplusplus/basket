

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

#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal
#' theme element_text coord_fixed scale_x_discrete scale_y_discrete 
#' geom_text labs guides element_blank guide_colorbar
exchangeogram <- function(mat, low = "black", high = "red", mid = "orange",
                          expand = c(0.3, 0.3), text_size = 4,
                          legend_position = c(0.3, 0.9)) {

  if (!is.null(mat) && any(rownames(mat) != colnames(mat))) {
    stop("The matrix supplied must be symmetric in the values and names.")
  }

  mg <- as_tibble(mat) %>% 
    mutate(V1 = rownames(mat)) %>% 
    gather(key = V2, value = value, -V1) %>%
    na.omit

  label_pos <- tibble(x = seq_len(nrow(mat)) - 1, y = seq_len(nrow(mat)), 
                      label = colnames(mat))
  ggplot(mg, aes(V2, V1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = low, high = high, mid = mid,
                         midpoint = 0.5, limit = c(0, 1), space = "Lab",
                         name = "Exchangeability") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, vjust = 0,
                                     size = 10, hjust = 0, color="#666666"))+
    theme(axis.text.y = element_text(angle = 0, vjust = 0,
                                    size = 10, hjust = 0, color="#666666"))+
    coord_fixed() +
    scale_x_discrete(labels = NULL, expand = expand) +
    scale_y_discrete(labels = NULL) +
    geom_text(aes(V2, V1, label = value), data = mg, color = "white", 
              size = text_size) +
    geom_text(aes(x = x, y = y, label = label), data = label_pos, 
              inherit.aes = FALSE, size = text_size) +
    labs(x = "", y = "") +
    theme(panel.grid.major = element_blank(), panel.border = element_blank(),
          panel.background = element_blank(), axis.ticks = element_blank(),
          legend.position = legend_position, legend.direction = "horizontal") +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
           title.position = "top", title.hjust = 0.5))
}

#' @export
plot_posterior_exchangeability.full_bayes <- function(x, ...) {
  mat <- round(x$PEP, 3)
  if (!is.null(basket_name(x))) {
    dimnames(mat) <- list(basket_name(x), basket_name(x))
  } else {
    dn <- paste("Basket", seq_len(nrow(mat)))
    dimnames(mat) <- list(dn, dn)
  }
  exchangeogram(mat, ...)
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

#' @export
plot_exchangeability.exchangeability_model <- function(x, ...) {
  mat <- round(x$maximizer, 3)
  mat[lower.tri(mat)] <- NA
  if (!is.null(basket_name(x))) {
    dimnames(mat) <- list(x$name, x$name)
  } else {
    dn <- paste("Basket", seq_len(nrow(mat)))
    dimnames(mat) <- list(dn, dn)
  }
  exchangeogram(mat, ...)
}


