#' Add sample size to the value
#'
#' @param x A vector, of any type
#'
#' @return A \code{factor} with the label and how many times that label appeared in \code{x}

#' @export
#'
#' @examples
#'
#' data(mtcars)
#'
#' mtcars$cyl_n <- add_N_to_catVar(mtcars$cyl)
#'
#' # you can also pipe it with tidyverse
#' mtcars %>%
#' mutate(cyl_n = add_N_to_catVar(cyl))
#'
add_N_to_catVar <- function(x){
  x2 <- as.factor(x)

  tmp <- table(x)

  new_names <- paste0(names(tmp), " (n = ", tmp, ")")

  levels(x2) <- new_names

  return(x2)
}
