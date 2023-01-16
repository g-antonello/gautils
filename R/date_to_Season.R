#' From date to its season
#'
#' Convert an object of Date class and turn into an ordered season factor, with Spring<Summer<Autumn<Winter
#'
#' @param date a \\code{Date} object
#'
#' @return An ordered \\code{factor}
#' @export
#'
#' @examples
#' dat <- c("2022-09-28", "2021-07-15", "2020-02-10", "2019-04-19")
#' season <- date_to_Season(date)
#'
#' levels(season)
#'
#' # want to have an unordered factor?
#' season <- factor(season, ordered = FALSE)
#' # or, with magrittr
#' # season %<>% factor(ordered = FALSE)
#'
date_to_Season <- function(date) {
  stopifnot(class(date) == "Date")
  scalarCheck <- function(date) {
    m <- as.POSIXlt(date)$mon + 1        # correct for 0:11 range
    d <- as.POSIXlt(date)$mday           # correct for 0:11 range
    if ((m == 3 &
         d >= 21) | (m == 4) | (m == 5) | (m == 6 & d < 21)) {
      r <- 1
    } else if ((m == 6 &
                d >= 21) | (m == 7) | (m == 8) | (m == 9 & d < 21)) {
      r <- 2
    } else if ((m == 9 &
                d >= 21) | (m == 10) | (m == 11) | (m == 12 & d < 21)) {
      r <- 3
    } else {
      r <- 4
    }
    r
  }

  res <- sapply(date, scalarCheck)
  res <-
    ordered(res, labels = c("Spring", "Summer", "Autumn", "Winter"))
  return(res)
}
