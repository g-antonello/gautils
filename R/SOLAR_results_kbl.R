#' Show clean SOLAR results
#'
#' Makes a kbl table from raw results obtained by by `get_SOLAR_results`, which is a useful function to retrieve multiple-trait heritability. But in principle it works even with just 1 trait
#'
#' @param x A \code{data.frame} object as retrieved as-is by `get_SOLAR_results`
#' @param padj_method A \code{character} with the adjustment required. 'none' is also accepted and crude P-values will be kept
#' @param complete.cases A \code{logical} saying whether you want to remove rows with 1 or more NAs. Default is `TRUE`
#'
#' @return a \code{kbl} table, with default settings
#' @export
#'
#' @examples
#'
#' # No microbiome SOLAR runs were stored as data.
#'

SOLAR_results_kbl <- function(x, padj_method = "BH", complete.cases = TRUE){

  if(complete.cases){
    x_complete <- x[complete.cases(x),]
  }

  if(padj_method == "none"){

    x2 <- x_complete %>%
      mutate(
        Estimate = paste0(round(h2, 3), " (", round(h2_SE, 3), ")"),
        Estimate2 = paste0(round(c2, 3), " (", round(c2_SE, 3), ")")
      ) %>%
      select(trait, Estimate, h2_p, Estimate2, c2_p) %>%
      mutate_if(is.numeric, format.pval, digits = 2) %>%
      as_tibble() %>%
      as.matrix() %>%
      set_colnames(c(" ", rep(c("Estimate (Std.Err.)", "P-value"),2)))

  } else{
    x_padj <- x_complete %>%
      mutate(
        h2_padj = p.adjust(h2_p, padj_method),
        c2_padj = p.adjust(c2_p, padj_method)


      )

    x2 <- x_padj %>%
      mutate(
        Estimate = paste0(round(h2, 3), " (", round(h2_SE, 3), ")"),
        Estimate2 = paste0(round(c2, 3), " (", round(c2_SE, 3), ")")
      ) %>%
      select(trait, Estimate, h2_padj, Estimate2, c2_padj) %>%
      mutate_if(is.numeric, format.pval, digits = 1) %>%
      as_tibble() %>%
      as.matrix() %>%
      set_colnames(c("Trait", rep(c("Estimate (Std.Err.)", paste0("Q-value", " (", padj_method, ")")),2)))
  }



  kable_basic <- kbl(x2) %>%
    kable_styling() %>%
    add_header_above(c(" " = 1, "h²" = 2, "c²" = 2))


  return(kable_basic)

}
