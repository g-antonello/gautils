#' Format P-values for tables
#'
#' Easy wrapper around format.pval to make readability of p-values better. It
#' converts only p-values below 0.00
#'
#' @param p A \code{numeric} vector of p-values to convert
#' @param threshold A \code{numeric} threshold under which to convert p-values
#' with scientific notation
#' @param scientific A \code{logical}, whether to use scientific notation.
#' (default = `TRUE`). NOTE: `scientific = FALSE` is basically meaningless
#'
#' @return A \code{character} vector with transformed p-values
#' @export
#'
#' @examples
#'
#' pvalues <- c(runif(5, 0, 0.001), runif(5, 0, 0.3), runif(5, 0, 1))
#'
#' format.pval2(pvalues)
#'
#' format.pval2(pvalues, threshold = 0.05)
#'
#' format.pval2(pvalues, threshold = 0.05, scientific = FALSE)
#'

format_pval2 <- function(p, threshold = 0.001, scientific = TRUE){
  return(
    ifelse(
      p < threshold,
      format.pval(p, scientific = scientific,digits = digits = 1),
      format.pval(p, digits = 3)
    )
  )
}


