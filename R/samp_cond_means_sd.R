#' Create Sample conditional means and standard deviation
#'
#' @param x the numeric/integer/double vector of values
#' @param g the factor or character vector to group values by
#' @param na.rm should
#'
#' @return a matrix object, containing mean and standard deviation on two columns, and group names in rows
#' @export
#'
#' @examples
#'
#' data("gss_cat")
#' age_by_marital_status <- samp_cond_means_sd(gss_cat$age, g = gss_cat$marital)
#' age_by_marital_status
samp_cond_means_sd <- function(x, g, na.rm = TRUE, fix_sd = TRUE){

  mean <- sapply(unique(g), function(i) mean(x[grepl(i,g, fixed = T)], na.rm = na.rm))
  sd <- sapply(unique(g), function(i) sd(x[grepl(i,g, fixed = T)], na.rm = na.rm))
  names(mean) <- unique(g)

  if(fix_sd){
    sd[is.na(sd)] <- 0
  }

  return(data.frame(varCateg = unique(g), mean = mean, sd = sd))
}
