#' Glimpse into a table-like object
#' A simple function to make a "head" of the data on both rows and column, of custom length
#' @param x H \code{data.frame} or \code{matrix-array} object
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

  if(nrows(x) < nrows){
    nrows <- nrows(x)
  }
  if(ncols(x) < ncols){
    ncols <- ncols(x)
  }

  return(x[1:nrows, 1:ncols])
}
