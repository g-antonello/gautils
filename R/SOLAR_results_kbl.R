#' Make a kable table of SOLAR results
#'
#' @param x An object obtained with \code{get_SOLAR_results} 
#' @param highlight_significant \code{numeric} A heritability (h2) P-value threshold under which you would like to highlight the row
#'
#' @import kableExtra
#' 
#' @return a \code{kableExtra} object. ready to print as .html
#' @export
#'
#' @examples
#' 
#' tmp <- get_SOLAR_results(path_to_h2_runs)
#' print(SOLAR_results_kbl(tmp, 0.01))

SOLAR_results_kbl <- function(x, highlight_significant = 0.05){
  
  x <- relocate(x, 
                h2_SE, 
                .before = h2_p) %>% 
    relocate(
      c2_SE,
      .before = c2_p
    )
  
  x2 <- x %>% 
    mutate(
      Estimate = paste0(round(h2, 3), " (", round(h2_SE, 3), ")"),
      Estimate2 = paste0(round(c2, 3), " (", round(c2_SE, 3), ")")
    ) %>% 
    select(trait, Estimate, h2_p, Estimate2, c2_p) %>% 
    mutate_if(is.numeric, format.pval, digits = 2) %>% 
    set_colnames(c(" ", rep(c("Estimate (Std.Err.)", "P-value"),2)))
  
  kable_basic <- kbl(x2) %>% 
    kable_styling() %>% 
    add_header_above(c(" " = 1, "h²" = 2, "c²" = 2))
  
  if(is.numeric(highlight_significant)){
    kable_basic %>% 
      row_spec(which(x$h2_p < highlight_significant), bold = T) %>% 
      row_spec(which(x$h2_p > highlight_significant), color = "lightgray") %>% 
      return()
    
  }else{
    return(kable_basic)
  }
}
