
#' Drop NAs from selected variables
#'
#' Given a phyloseq object and a character vector, the function removes all NAs that are in at least one variable
#' @param physeq a phyloseq object
#' @param variables a character vector with variable names to remove NAs from
#' @param verbose A \code{logical}. If TRUE, it will print how many samples per variable it removed
#'
#' @return a clean phyloseq object
#'
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' set.seed(1111)
#' GlobalPatternsNA <- GlobalPatterns
#' GlobalPatternsNA@sam_data$newVar1 <- replicate(n = nsamples(GlobalPatternsNA), expr = sample(x =c(NA, 1,0), size = 1, prob = c(0.1, 0.4, 0.5)))
#' GlobalPatternsNA@sam_data$newVar2 <- replicate(n = nsamples(GlobalPatternsNA), expr = sample(x =c(NA, 1,0), size = 1, prob = c(0.1, 0.4, 0.5)))
#'
#' drop_NA_from_phyloseq_vars(physeq = GlobalPatternsNA, vars = c("newVar1", "newVar2"))
#' GlobalPatterns

phy_drop_NA_vars <- function(physeq, variables, verbose = TRUE){

  metadata_reduced <- microbiome::meta(physeq) %>%
    select(all_of(variables))

  if(verbose){
  nas <-  sapply(metadata_reduced, function(x) sum(is.na(x))) %>%
    as.data.frame() %>%
    set_names("tot_NAs")

  print(nas)
  }
  # get a logical vector to subset variables in those
  subset_logical <- complete.cases(metadata_reduced)

  physeq_cleaned <- speedyseq::filter_sample_data(physeq,subset_logical)
  physeq_cleaned <- prune_taxa(taxa_sums(physeq_cleaned) > 0,physeq_cleaned)

  return(physeq_cleaned)
}
