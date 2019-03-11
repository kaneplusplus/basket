
#' @title Get the Basketwise Posterior Exchangeability Matrix
#' 
#' @description TODO: WRITE THIS
#' @param x either an exchangeability model or basket object.
#' @examples
#' # TODO: WRITE THIS
#' @export
basket_pep <- function(x) {
  UseMethod("basket_pep", x)
}

#' @importFrom crayon red
basket_pep.default <- function(x) {
  stop(cat_line("Don't know how to extract posterior probability matrix from ",
                "an object of type ", class(x), "\n"))
}

basket_pep.mem_basket <- function(x) {
  x$PEP
}

basket_pep.exchangeability_model <- function(x) {
  basket_pep(x$basket)
}

#' @title Get the Basketwise Maximum A Posteriori Probability Matrix
#' @description TODO: WRITE THIS
#' @param x either an exchangeability model or basket object.
#' @examples
#' # TODO: WRITE THIS
#' @export
basket_map <- function(x) {
  UseMethod("basket_map", x)
} 

#' @importFrom crayon red
basket_map.default <- function(x) {
  stop(cat_line("Don't know how to extract maximum a posteriori probability ",
                "matrix from an object of type ", class(x), "\n"))
}

basket_map.mem_basket <- function(x) {
  x$PEP
}

basket_map.exchangeability_model <- function(x) {
  basket_map(x$basket)
}


#' @title Get the Clusterwise Posterior Exchangeability Matrix
#' 
#' @description TODO: WRITE THIS
#' @param x either an exchangeability model or basket object.
#' @examples
#' # TODO: WRITE THIS
#' @export
cluster_pep <- function(x) {
  UseMethod("cluster_pep", x)
}

#' @importFrom crayon red
cluster_pep.default <- function(x) {
  stop(cat_line("Don't know how to extract posterior probability matrix from ",
                "an object of type ", class(x), "\n"))
}

cluster_pep.mem_basket <- function(x) {
  x$PEP
}

cluster_pep.exchangeability_model <- function(x) {
  cluster_pep(x$basket)
}

#' @title Get the Clusterwise Maximum A Posteriori Probability Matrix
#' @description TODO: WRITE THIS
#' @param x either an exchangeability model or basket object.
#' @examples
#' # TODO: WRITE THIS
#' @export
cluster_map <- function(x) {
  UseMethod("cluster_map", x)
} 

#' @importFrom crayon red
cluster_map.default <- function(x) {
  stop(cat_line("Don't know how to extract maximum a posteriori probability ",
                "matrix from an object of type ", class(x), "\n"))
}

cluster_map.mem_basket <- function(x) {
  x$PEP
}

cluster_map.exchangeability_model <- function(x) {
  basket_map(x$basket)
}


