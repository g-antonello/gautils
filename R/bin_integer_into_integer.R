
#' Bin an integer/numeric variable
#'
#' This function bins decimal numbers to their lower or upper values.
#' This is particularly useful for semi-quantitative modelling, in case there is a preference for a set of values.
#'
#' @param x A `numeric` vector
#' @param multiples An `integer` or `numeric` vector to bin `x` by. For example, \code{multiples = 2} will return all values approximated to \code{0,2,4,6, ...}
#' @param include.lowest Do you want to include upper or lower values in the bin? eg: Take `multiples = 5`, the ranges will be [0,5), [5,10), ..., if FALSE, (0,5], (5,10], ...

#' @return A \code{integer} vector
#'
#' @export
#'
#' @examples
#' data(mtcars)
#' bin_integer_into_integer(mtcars$mpg, multiples = 10, include_highest= TRUE)
#' bin_integer_into_integer(mtcars$mpg, multiples = 5, include_highest= FALSE)
#'
bin_integer_into_integer <- function(x, multiples, include.lowest = FALSE){
  brks <- seq(0, max(x)+multiples, multiples)
  x_fact <- as.character(cut(x, breaks = brks,include.lowest = include.lowest)) %>%
    strsplit(split = "\\(|\\]") %>% sapply("[", 2) %>%
    strsplit(split = ",", fixed = TRUE) %>% sapply("[", 2) %>%
    as.numeric()

  return(x_fact)
}
