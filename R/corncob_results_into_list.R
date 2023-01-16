
#' Extract corncob results
#'
#' An easy wrangling function to extract all results from a corncob differential abundance experiment. NB: so far two things are left to test:
#' (1) the order of the formulas might be important, and the last term is what this function will extract. (2) Extraction of results for regression against a continuous variable is not available yet
#'
#' @param res A `corncob differentialTest` result object. \\textbf{NB: only analyses with `phyloseq`objects as input are supported for now}.
#' @param trait A `character` vector with the metadata column name to extract the results from.
#' @param class
#'
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom stringr str_split
#' @importFrom microbiome meta
#'
#' @return A `list` of `data.frame`s with all contrasts applicable against the reference. NB: other pairwise contrasts are not made
#'
#' @export
#'
#' @examples
#' library(corncob)
#' # phyloseq example
#' data(soil_phylum_small)
#' da_analysis <- differentialTest(formula = ~ DayAmdmt,
#'                                phi.formula = ~ DayAmdmt,
#'                                formula_null = ~ 1,
#'                               phi.formula_null = ~ DayAmdmt,
#'                                test = "Wald", boot = FALSE,
#'                                data = soil_phylum_small,
#'                                fdr_cutoff = 0.05)
#'
#' da_analysis_wrangled <- corncob_results_into_list(da_analysis, class = "factor")
#'

corncob_results_into_list <- function(res, class = "factor"){

  trait0 <- str_split(as.character(res$significant_models[[1]]$formula), "\\+ ")[[3]]
  trait <- trait0[length(trait0)]

  classTrait <- sum(agrepl(trait, rownames(res$significant_models[[1]]$coefficients)))

  if(classTrait > 1){

    trait_vec <- as.factor(meta(res$data)[[trait]])

    step1 <- lapply(res$all_models, "[[", "coefficients") %>%
      lapply(as.data.frame) %>%
      lapply(rownames_to_column, "coef")

    names(step1) <- taxa_names(res$data)

    step2 <- bind_rows(step1, .id = "taxon")

    step3 <- sapply(levels(trait_vec)[2:length(levels(trait_vec))], function(levl) filter(step2, grepl(levl, coef)), USE.NAMES = T, simplify = F)
    step4 <- lapply(step3, left_join,
                    as.data.frame(res$p) %>% rownames_to_column() %>%  set_names(c("taxon", "p_nominal"))) %>%
      lapply(left_join,
             as.data.frame(res$p_fdr) %>% rownames_to_column() %>%  set_names(c("taxon", "padj"))) %>%
      lapply(rename, "SE" = "Std. Error") %>%
      lapply(select, -coef, -`Pr(>|t|)`)

    return(step4)
  } else{ # this means it's a numeric variable, but this chunk is not tested yet

    step1 <- lapply(res$all_models, "[[", "coefficients") %>%
      lapply(as.data.frame) %>%
      lapply(rownames_to_column, "coef")

    names(step1) <- taxa_names(res$data)

    step2 <- bind_rows(step1, .id = "taxon") %>%
      filter(agrepl(trait, coef))

    step3 <- mutate(step2,
                    pnominal = res$p,
                    padj = res$p_fdr) %>%
      rename("SE" = "Std. Error") %>%
      select(-coef, -`Pr(>|t|)`)

    return(step3)
  }
}
