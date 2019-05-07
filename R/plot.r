
#' @title Plot the Response Densities in Basket Trials
#'
#' @description The MEM analysis calculates the probability of exchangeability
#' of baskets and clusters in a basket trial. This function creates density 
#' plots of the response rates of each basket or each cluster under the MEM 
#' design taking into acount the extent to which power can be borrowed from
#' similar trials.
#' @param x the exchangeability model.
#' @param ... other options. See Details for more information.
#' @details The ... options can be used to specify the colors of the response
#' density plot or, when plotting an object of class `exchangeability_model`
#' the type can be specified. In this case, the default is
#' `type = c("both", "basket", "cluster")`.
#' @examples
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(vemu_wide$responders,
#'                           vemu_wide$evaluable,
#'                           vemu_wide$baskets)
#'
#' plot_density(mem_analysis)
#' @importFrom tidyr gather
#' @importFrom tibble as_tibble
#' @importFrom gridExtra grid.arrange
#' @export
plot_density <- function(x, ...) {
  UseMethod("plot_density")
}

#' @importFrom crayon red
#' @export
plot_density.default <- function(x, ...) {
  stop(red(
    "Don't know how to make a density plot with an object of type",
    paste(class(x), collapse = ", "), "."))
}

#' @export
plot_density.exchangeability_model <- function(x, ...) {
  dots <- list(...)
  if ("type" %in% names(dots)) {
    ps <- dots$type
  } else {
    ps <- c("basket", "cluster")
  }
  plots <- lapply(ps, function(pt) plot_density(x[[pt]]))
  if (length(plots[[1]])) {
    plots[[1]]
  } else {
    suppressWarnings(do.call(grid.arrange, c(plots, ncol = length(plots))))
  }
}

#' @importFrom ggplot2 ggplot aes geom_density scale_fill_manual facet_grid
#' xlab ylab theme_minimal xlim geom_vline labeller label_wrap_gen
#' @export
plot_density.mem <- function(x, ...) {
  dots <- list(...)
  Basket <- Density <- p0 <- NULL
  if (length(x$p0) == 1 && length(x$size) > 1) {
    x$p0 <- rep(x$p0, length(x$size))
  }
  d <- gather(as_tibble(x$samples), key = Basket, value = Density)
  d$p0 <- x$p0

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
  }
  d$p0 <- x$p0[match(d$Basket, x$name)]
  d$basket_name <- paste0(d$Basket, " (p0=", d$p0, ")")
  ggplot(d, aes(x = Density, fill = Basket)) +
    geom_density(alpha = 0.7) +
    facet_grid( basket_name~ ., labeller = label_wrap_gen(width = 10)) +
    geom_vline(aes(xintercept = p0)) +
    scale_fill_manual(values = basket_colors, guide = FALSE) +
    xlab("") +
    ylab("Density") +
    xlim(0, 1) +
    theme_minimal()
}

#' @title Plot the Posterior Exchangeability of a Basket Trial
#'
#' @description The posterior exchangeability of the baskets in an an 
#' MEM analysis can be visualized via an exchangeogram using this function.
#' @param x the exchangeability model.
#' @param ... other options passed to ggplot2 to alter the visual 
#' characteristics of the plot. See Details for more information.
#' @details The `plot_pep` attempts to place the basket names to the
#' left of the main diagonal in a way that makes it easy to read. However,
#' for especially long basket names options are provided to ``fine tune''
#' the visualizations. These auxiliary options include:
#' \itemize{
#'  \item{[low] }{The color corresponding to a low degree of exchangeability.
#'  (Default "black")}
#'  \item{[high] }{The color corresponding to a high degree of exchangeability.
#'  (Default "red")}
#'  \item{[mid] }{The color corresponding to 50\% exchangeability.
#'  (Default "orange")}
#'  \item{[expand] }{The proportion to expand the viewport
#'  (Default expand = c(0.3, 0.3))}
#'  \item{[text_size] }{The text size. (Default 4)}
#'  \item{[legend_position] }{The legend position. 
#'    (Default legend_position = c(0.25, 0.8)}
#'  \item{[draw_legend] }{Should the legend be drawn? (Default TRUE)}
#'  \item{[basket_name_hoffset] }{The horizontal offset of the basket names..
#'  (Default 0)}
#'  \item{[basket_name_hjust] }{The basket name justification..
#'  (Default 1 - right justified)}
#' }
#' @examples
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(vemu_wide$responders,
#'                           vemu_wide$evaluable,
#'                           vemu_wide$baskets)
#'
#' plot_pep(mem_analysis)
#' @export
plot_pep <- function(x, ...) {
  UseMethod("plot_pep")
}

#' @importFrom crayon red
plot_pep.default <- function(x, ...) {
  stop(red(
    "Don't know how to make a posterior exchangeability plot",
    "with an object of type", 
    paste(class(x), collpase = ", "), "."))
}

#' @importFrom stats na.omit
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr mutate %>%
#' @importFrom tidyr gather
#' @importFrom crayon red
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal
#' theme element_text coord_fixed scale_x_discrete scale_y_discrete
#' geom_text labs guides element_blank guide_colorbar
exchangeogram <- function(mat, low = "black", high = "red", mid = "orange",
                          expand = c(0.3, 0.3), text_size = 4,
                          legend_position = c(0.25, 0.8), draw_legend = TRUE,
                          basket_name_hoffset = 0, basket_name_hjust = 1) {
  if (!is.null(mat) && any(rownames(mat) != colnames(mat))) {
    stop(red("The matrix supplied must be symmetric in the",
             "values and names."))
  }

  for (i in 1:dim(mat)[1])
    for (j in 1:dim(mat)[2])
    {
      if (i <= j) {
        next
      }
      mat[i, j] <- NA
    }

  V2 <- value <- V1 <- x <- y <- label <- NULL

  baskets <- unique(colnames(mat))
  mg <- as_tibble(mat) %>%
    mutate(V1 = rownames(mat)) %>%
    gather(key = V2, value = value, -V1) %>%
    na.omit() %>%
    mutate(
      V1 = factor(V1, levels = baskets),
      V2 = factor(V2, levels = baskets)
    )


  label_pos <- tibble(
    x = seq_len(nrow(mat)) - 1 + basket_name_hoffset,
    y = seq_len(nrow(mat)), label = colnames(mat)
  )
  tG <- ggplot(mg, aes(V2, V1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = low, high = high, mid = mid,
      midpoint = 0.5, limit = c(0, 1), space = "Lab",
      name = "Probability"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
      angle = 0, vjust = 0,
      size = 10, hjust = 0, color = "#666666"
    )) +
    theme(axis.text.y = element_text(
      angle = 0, vjust = 0,
      size = 10, hjust = 0, color = "#666666"
    )) +
    coord_fixed() +
    scale_x_discrete(labels = NULL, expand = expand) +
    scale_y_discrete(labels = NULL) +
    geom_text(aes(V2, V1, label = value),
      data = mg, color = "white",
      size = text_size
    ) + geom_text(aes(x = x, y = y, label = label),
      data = label_pos,
      inherit.aes = FALSE, size = text_size, hjust = basket_name_hjust
    ) +
    labs(x = "", y = "")
  if (draw_legend) {
    tG <- tG +
      theme(
        panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank(),
        legend.position = legend_position, legend.direction = "horizontal"
      ) +
      guides(fill = guide_colorbar(
        barwidth = 7, barheight = 1,
        title.position = "top", title.hjust = 0.5
      ))
  } else {
    tG <- tG +
      theme(
        panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank(),
        legend.position = legend_position + c(-110, -100), 
        legend.direction = "horizontal") +
      guides(fill = guide_colorbar(
        barwidth = 7, barheight = 1,
        title.position = "top", title.hjust = 0.5
      ))
  }

  return(tG)
}

#' @title Plot the Prior, MAP, and PEP of a Basket Trial
#' 
#' @description: Plot the Prior, MAP, and PEP Matrices
#' @param x the exchangeability model.
#' @param type the plot type that will be plotted.
#' @param ... other options. See Details for more information.
#' @export
plot_mem <- function(x, type = c("prior", "map", "pep"), ...) {
  UseMethod("plot_mem", x)
}

#' @export
plot_mem.default <- function(x, type = c("prior", "map", "pep"), ...) {
  stop(red(
    "Don't know how to make an MEM plot with an object of type", 
    paste(class(x), collpase = ", ", ".")))
}

#' @export
plot_mem.exchangeability_model <- 
  function(x, type = c("prior", "map", "pep"), ...) {
  plot_mem(x$basket, type, ...)
}

#' @export
plot_mem.mem <- function(x, type = c("prior", "map", "pep"), ...) {
  # library("gridExtra")
  numC <- 0
  allPlot <- list()
  if (any(type == "prior")) {
    mat <- round(x$prior, 3)
    if (!is.null(x$name)) {
      dimnames(mat) <- list(x$name, x$name)
    } else {
      dn <- paste("Basket", seq_len(nrow(mat)))
      dimnames(mat) <- list(dn, dn)
    }
    plot1 <- (exchangeogram(mat, draw_legend = FALSE, ...) +
      ggtitle("Prior") +
      theme(
        plot.title = element_text(
          # family = "Trebuchet MS",
          color = "#666666",
          face = "bold",
          size = 20,
          hjust = 0.35
        )
      ))
    allPlot <- c(allPlot, list(plot1))
    numC <- numC + 1
  }


  if (any(type == "map")) {
    mat <- round(x$MAP, 3)
    if (!is.null(x$name)) {
      dimnames(mat) <- list(x$name, x$name)
    } else {
      dn <- paste("Basket", seq_len(nrow(mat)))
      dimnames(mat) <- list(dn, dn)
    }
    plot2 <- (exchangeogram(mat, legend_position = c(0.45, -0.22), ...) +
      ggtitle("MAP") +
      theme(
        plot.title = element_text(
          # family = "Trebuchet MS",
          color = "#666666",
          face = "bold",
          size = 20,
          hjust = 0.35
        )
      ))
    allPlot <- c(allPlot, list(plot2))
    numC <- numC + 1
  }

  if (any(type == "pep")) {
    mat <- round(x$PEP, 3)
    if (!is.null(x$name)) {
      dimnames(mat) <- list(x$name, x$name)
    } else {
      dn <- paste("Basket", seq_len(nrow(mat)))
      dimnames(mat) <- list(dn, dn)
    }
    plot3 <- exchangeogram(mat, draw_legend = FALSE, ...) +
      ggtitle("Posterior Prob.") +
      theme(
        plot.title = element_text(
          color = "#666666",
          face = "bold",
          size = 20,
          hjust = 0.35
        )
      )
    allPlot <- c(allPlot, list(plot3))
    numC <- numC + 1
  }

  PLOTS <- allPlot
  # grid.arrange(plot2, ncol = 3, nrow = 1)
  do.call(grid.arrange, c(PLOTS, ncol = numC))
}

#' @export
plot_pep.mem <- function(x, ...) {
  mat <- round(x$PEP, 3)
  if (!is.null(x$name)) {
    dimnames(mat) <- list(x$name, x$name)
  } else {
    dn <- paste("Basket", seq_len(nrow(mat)))
    dimnames(mat) <- list(dn, dn)
  }
  exchangeogram(mat, ...) +
    ggtitle("Posterior Exchangeability Probability") +
    theme(plot.title = element_text(
      color = "#666666",
      face = "bold", size = 20, hjust = 0.35
    ))
}

#' @title Plot the Map Exchangeability of a Basket Trial
#'
#' @description The Maximum Aposteriori Probability (MAP) of an MEM is the
#' estimate of the exchangeability structure of a basket trial. This function
#' visualizes this matrix as an exchangeogram.
#' @param x the exchangeability model.
#' @param ... other options passed to ggplot2 to alter the visual 
#' @details The `plot_pep` attempts to place the basket names to the
#' left of the main diagonal in a way that makes it easy to read. However,
#' for especially long basket names options are provided to ``fine tune''
#' the visualizations. These auxiliary options include:
#' \itemize{
#'  \item{[low] }{The color corresponding to a low degree of exchangeability.
#'  (Default "black")}
#'  \item{[high] }{The color corresponding to a high degree of exchangeability.
#'  (Default "red")}
#'  \item{[mid] }{The color corresponding to 50\% exchangeability.
#'  (Default "orange")}
#'  \item{[expand] }{The proportion to expand the viewport
#'  (Default expand = c(0.3, 0.3))}
#'  \item{[text_size] }{The text size. (Default 4)}
#'  \item{[legend_position] }{The legend position. 
#'    (Default legend_position = c(0.25, 0.8)}
#'  \item{[draw_legend] }{Should the legend be drawn? (Default TRUE)}
#'  \item{[basket_name_hoffset] }{The horizontal offset of the basket names..
#'  (Default 0)}
#'  \item{[basket_name_hjust] }{The basket name justification..
#'  (Default 1 - right justified)}
#' }
#' @examples
#' # Create an MEM analysis of the Vemurafenib trial data.
#' data(vemu_wide)
#'
#' mem_analysis <- mem_exact(vemu_wide$responders,
#'                           vemu_wide$evaluable,
#'                           vemu_wide$baskets)
#'
#' plot_map(mem_analysis)
#' @export
plot_map <- function(x, ...) {
  UseMethod("plot_map")
}

#' @importFrom crayon red
#' @export
plot_map.default <- function(x, ...) {
  stop(red(
    "Don't know how to make a posterior exchangeability plot",
    "with an object of type", 
    paste(class(x), collpase = ", ", ".")))
}

#' @importFrom ggplot2 ggtitle element_text theme
#' @export
plot_map.mem <- function(x, ...) {
  mat <- round(x$MAP, 3)
  mat[lower.tri(mat)] <- NA
  if (!is.null(x$name)) {
    dimnames(mat) <- list(x$name, x$name)
  } else {
    dn <- paste("Basket", seq_len(nrow(mat)))
    dimnames(mat) <- list(dn, dn)
  }
  exchangeogram(mat, ...) +
    ggtitle("Maximum A Posteriori MEM") +
    theme(plot.title = element_text(
      # family = "Trebuchet MS",
      color = "#666666",
      face = "bold", size = 25, hjust = 0.5
    ))
}
