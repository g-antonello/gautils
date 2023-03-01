#' calculate p-values in table1 objects
#'
#' Inspired by the table1 vignette itself (https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html#example-a-column-of-p-values), only to have it quick every time I want a table1 significance
#' 
#'
#' @param x a table in the `table1` set
#' @param ... extra parameters. so far, none implemented, that I know 
#'
#' @return a `table1` object
#' @export
#'
#' @examples
#' 
#' data(mtcars)
#' table1::table1( ~ vs + gear + carb |as.factor(cyl), 
#'   data  = mtcars,
#'   extra.col = list("significance" = table1.pvalues),
#'   overall = "All cars"
#'   )
#'   
table1.pvalues <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- kruskal.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g), correct = T)$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", ifelse(p >= 0.01, format(p, digits = 2), format(p, digits=2, scientific = TRUE))))
}
