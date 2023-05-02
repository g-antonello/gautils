#' Glimpse into a table-like object
#' A simple function to make a "head" of the data on both rows and column, of custom length
#' @param x An object coercible to \code{matrix} format
#' @param nrows How many rows to print (default = 8)
#' @param ncols How many columns to print (default = 8)
#'
#' @return a printed subset of the data

#' @export
#'
#' @examples
#' data(iris)
#' qklook(iris)
#' qklook(iris, 5,5)

qklook <- function(x, nrows=8, ncols=8){
  x_formatted <- as(x, "matrix")

  if(nrow(x) < nrows){
    nrows <- nrow(x)
  }
  if(ncol(x) < ncols){
    ncols <- ncol(x)
  }

  return(x[1:nrows, 1:ncols])
}
