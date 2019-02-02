
#' @title Make the operating characteristics table
#'
#' @description From the input full_bayes class object, summarize the CDF, HPD, ESS, Mean, and Median results.
#' @param res the full_bayes class object..
#' 
#' @examples


#' @export
OCTable <- function(res)
{
  UseMethod("OCTable")
}

#' @export
OCTable.default <- function(res) 
{
  stop(paste("Object class is not handled in OCTable function",
             class(res)))
}

#' @export
OCTable.full_bayes <- function(res)
{
  #res <- MHResult1$clusterwise
  if (class(res$CDF) == "matrix")
  {
    cdfS <- c()
    for (i in 1:dim(res$CDF)[1])
    {
      ss <- paste("CDF ", rownames(res$CDF)[i])
      cdfS <- c(cdfS, ss)
    }
  } else {
    cdfS <- "CDF"
  }
  oct <- rbind(res$CDF, res$HPD, res$ESS, res$mean_est, res$median_est)
  rownames(oct) <- c(cdfS, "HPD LB", "HPD HB", "ESS",  "Mean", "Median")
  return(oct)
}


