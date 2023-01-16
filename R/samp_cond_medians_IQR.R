#' Create Sample conditional means and standard deviation
#'
#' @param x the numeric/integer/double vector of values
#' @param g the factor or character vector to group values by
#' @param na.rm should
#'
#' @return a matrix object, containing median and IQR on two columns, and group names in rows
#' @export
#'
#' @examples
#'
#' data("gss_cat")
#' age_by_marital_status <- samp_cond_medians_IQR(gss_cat$age, g = gss_cat$marital)
#' age_by_marital_status
#'

samp_cond_medians_iqr <- function(x, g, na.rm = TRUE, fix_iqr = TRUE){

  median <- sapply(unique(g), function(i) median(x[grepl(i,g, fixed = T)], na.rm = na.rm))
  iqr <- sapply(unique(g), function(i) IQR(x[grepl(i,g, fixed = T)], na.rm = na.rm))
  names(median) <- unique(g)

  if(fix_iqr){
    iqr[is.na(iqr)] <- 0
  }

  return(data.frame(varCateg = unique(g), median = median, iqr = iqr))
}

