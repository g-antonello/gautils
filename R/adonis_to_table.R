#' Quick wrapper to format PERMANOVA results
#'
#' @param adonis.out A \code{adonis2} output object as is
#' @param new_rownames A \code{character} vector containing prettier variable
#' names, which will be passed with set_rownames
#' @param adjR2 A \code{logical}, should another column be calculated with the adjusted R squared? the thing is calculated with \code{vegan::RsquareAdj}
#'
#' @return A \code{data.frame} Ready to be passed to kable/DT/flextable
#'
#' @importFrom magrittr %>% set_rownames set_colnames
#' @importFrom dplyr mutate rename relocate any_of across contains
#' @importFrom vegan RsquareAdj
#' @export
#'
#' @examples
#'
#' library(gautils)
#' data(enterotype)
#' distMat <-  phyloseq::distance(microbiome::transform(enterotype, "compositional"), "bray")
#'
#' adonisObj <- adonis2(distMat ~ SeqTech + Enterotype, by = "margin", permutations = 2000, data = microbiome::meta(enterotype), na.action = na.omit)
#'
#' adonis_as_table <- adonis_to_table(adonisObj)
#' adonis.kbl <- kableExtra::kbl(adonis_as_table, digits = 2) %>%
#' kableExtra::kable_styling()
#'
#' adonis_as_table2 <- adonis_to_table(adonisObj, new_rownames = c("Sequencing Techonology", "Enterotype Code"))
#' adonis.kbl2 <- kableExtra::kbl(adonis_as_table2, digits = 2) %>%
#' kableExtra::kable_styling()
#'

adonis_to_table <- function(adonis.out, adjR2 = TRUE, new_rownames = NA){

  tmp0 <- as.data.frame(adonis.out)

  if(adjR2){
    tmp0 <- mutate(tmp0,
             R2.adj = RsquareAdj(R2, m = Df, n = tail(adonis.out$Df, 1))
             ) %>%
      relocate(R2.adj, .after = R2)
    tmp0[nrow(tmp0)-1, "R2.adj"] <- NA
  }

  if(!all(is.na(new_rownames))){

    tmp0 <- tmp0 %>%
      set_rownames(c(new_rownames, c("Residual Variance", "Total")))
  }

  lookup_colnames <- c(`Rsq (%)` = "R2", `Rsq adj. (%)` = "R2.adj",  `P-value` = "Pr(>F)", `F-stat` = "F")

  formatted_permanova_results <- tmp0 %>%
    mutate(across(contains("R2"))*100,
           `Pr(>F)` = format_pval2(`Pr(>F)`, threshold = 0.001, scientific = TRUE)
           ) %>%
    rename(any_of(lookup_colnames))

  return(formatted_permanova_results)
}
