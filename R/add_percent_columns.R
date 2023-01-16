#' Add percentages to multicolumn tables
#'
#' Takes a complex table object as input, and it calculates the percentage for each column
#'
#' @param tb a `table` R object
#'
#' @return a `matrix` with percent columns added next to each table column
#' @export
#'
#' @examples
#'
#' data(mtcars)
#' tb0 <- table(mtcars$cyl, mtcars$gear)
#' tb1 <- add_percent_columns(tb = tb0)
#'
#'
add_percent_columns <- function(tb){
  new_tb <- matrix(0,nrow = nrow(tb), ncol = ncol(tb)*2)
  rownames(new_tb) <- rownames(tb)
  colnames(new_tb) <- rep(c("abs", "%"), ncol(tb))
  for(i in 1:ncol(tb)){

    new_tb[,c(2*i-1, 2*i)] <- cbind(tb[,i], round(tb[,i]/sum(tb[,i])*100, digits = 1))

  }
  return(new_tb)
}
