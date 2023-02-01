#' From DESeq2 results to a list of data.frames
#'
#' This function only works with categorical contrasts
#'
#' @param deseq_results a \code{DESeq2} results object
#' @param trait A \code{character} with the variable name to compare contrasts within
#' @param physeq_obj the phyloseq object you used, if any
#' @param taxon_lvls_added A \code{character} with the taxonomic levels to add to each results \code{data.frame}
#' @param sort_p.val a \code{logical} saying whether you want the `padj` column sorted or not
#'
#' @importFrom DESeq2 results
#'
#' @return A \code{list} of pairwise contrasts based on the factor levels of the \code{trait} variable
#'
#' @export
#'
#' @examples
#'
#' library(phyloseq)
#' library(DESeq2)
#' filepath <- system.file("extdata", "study_1457_split_library_seqs_and_mapping.zip", package="phyloseq")
#' kostic <- microbio_me_qiime(filepath)
#'
#' kostic = subset_samples(kostic, DIAGNOSIS != "None")
#'
#' diagdds = phyloseq_to_deseq2(kostic, ~ OSH_DIAGNOSIS)
#'
#' # calculate geometric means prior to estimate size factors
#'
#'  gMeans = apply(counts(diagdds), 1, function(x, na.rm=TRUE) exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x)))
#'
#' # do the glm modelling
#'
#' diagdds2 = estimateSizeFactors(diagdds, geoMeans = gMeans)
#' diagdds = DESeq(diagdds2, fitType="local")
#'
#' res_list_with_Phyloseq <- deseq_results_into_list(deseq_results = diagdds,
#'                         trait = "OSH_DIAGNOSIS",
#'                         physeq_obj = kostic,
#'                         taxon_lvls_added = c("Phylum", "Class", "Genus")
#'                         )
#'
#' res_list_withOut_Phyloseq <- deseq_results_into_list(deseq_results = diagdds,
#'                         trait = "OSH_DIAGNOSIS"
#'                         )


deseq_results_into_list <- function(deseq_results,
                                    trait,
                                    physeq_obj = NULL,
                                    taxon_lvls_added = "all",
                                    sort_p.val = TRUE
                                    ){

  combin <- (t(combn(levels(deseq_results[[trait]]),2)))

  results_list <- lapply(1:nrow(combin), function(i) results(deseq_results, contrast =c(trait, combin[i,1], combin[i,2])) %>%
           as.data.frame() %>%
           tibble::rownames_to_column("taxon")
         )

  names(results_list) <- paste(combin[,1], combin[,2], sep = "_vs_")

  if(!is.null(physeq_obj)) {
    if (!is.null(taxon_lvls_added)) {
      taxtable <- as.data.frame(tax_table(physeq_obj)) %>%
        rownames_to_column("taxon")

      if (any("all" == taxon_lvls_added)) {
        results_list <-
          lapply(results_list, left_join, taxtable, by = "taxon")

      }

      if (all("all" != taxon_lvls_added)) {
        taxtable_subset <- taxtable %>%
          select(all_of(c(taxon_lvls_added, "taxon")))

        results_list <- lapply(results_list,
                               left_join,
                               taxtable_subset,
                               by = "taxon")

      }
    }
  }

  # sort the p-value if requested, by default it does it

  if(sort_p.val){
    results_list <- lapply(results_list, arrange, padj)
  }

  return(results_list)

}

