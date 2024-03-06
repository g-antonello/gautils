#' Quick wrapper to format PERMANOVA results
#'
#' @param adonis.out A \code{adonis2} output object as is
#' @param new_rownames A \code{character} vector containing prettier variable
#' names, which will be passed with set_rownames
#'
#' @return A \code{data.frame} Ready to be passed to kable/DT/flextable
#' 
#' @importFrom magrittr %>% set_rownames set_colnames
#' @importFrom dplyr mutate rename
#' @export
#'
#' @examples
#' library(gautils)
#' data(enterotype)
#' distMat <-  phyloseq::distance(microbiome::transform(enterotype, "compositional"), "bray")
#' 
#' adonisObj <- adonis2(distMat ~ SeqTech + Enterotype, by = "margin", permutations = 2000, data = microbiome::meta(enterotype), na.action = na.omit)
#' 
#' adonis_as_table <- adonis_to_table(adonisObj)
#' 
#' kableExtra::kbl(adonis_as_table, digits = 2) %>% 
#' kableExtra::kable_styling()

adonis_to_table <- function(adonis.out, new_rownames = NA){
  if(is.na(new_rownames)){
    tmp0 <- adonis.out %>% 
      as.data.frame() 
  } else{
    tmp0 <- adonis.out %>% 
      as.data.frame() %>% 
      set_rownames(c(new_rownames, c("Residual Variance", "Total"))) 
  }
  
  formatted_permanova_results <- tmp0 %>% 
    set_colnames(c("Df", "SumOfSqs", "R2", "F-stat", "pval")) %>% 
    mutate(R2 = R2*100, 
           pval = format.pval2(pval, threshold = 0.001, scientific = TRUE)) %>% 
    rename(
      `R2 (%)` = R2,
      `P-value` = pval
           )
  return(formatted_permanova_results)
}
