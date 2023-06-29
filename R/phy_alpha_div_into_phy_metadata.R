#' Add alpha diversities into phyloseq metadata
#'
#' This function takes a phyloseq object, passes it throuh the `estimate_richness` function in `phyloseq` and returns the phyloseq object containing alpha diversity estimates as "alpha_[diversity_name]"
#'
#' @param physeq a `phyloseq` object
#' @param raref_depth the number of minimum reads to rarefy by
#' @param measures the alpha diversity measures, as implemented by the `phyloseq::estimate_richness` function
#' @param verbose if `TRUE`, will tell you what it is doing, in particular in the rarefaction step
#'
#' @import phyloseq
#' @import microbiome
#' @import tidyverse
#'
#' @return a `phyloseq` object with the alpha diversity estimates added to the metadata
#'
#' @export
#'
#' @examples
#'
#' data(GlobalPatterns)
#' phy_alpha_div_into_phy_metadata(GlobalPatterns,
#'      raref_depth = 5000,
#'      measures = c("Shannon", "Simpson"),
#'      verbose = TRUE)
#'

phy_alpha_div_into_phy_metadata <- function(physeq,
                                            raref_depth = NULL,
                                            measures,
                                            verbose = FALSE
){

  sample_names(physeq) <- paste0("tmp_",sample_names(physeq))

  if(is.null(raref_depth)){
    raref_depth = min(colSums(abundances(physeq)))
    if(verbose){
      message(paste("Rarefying all samples to", raref_depth, "counts each"))
    }
  }


  # check if some metrics have already been calculated
  measures <- measures[!(paste0("alpha_",measures) %in% colnames(meta(physeq)))]

  if(is_empty(measures)){
    return(physeq)
    if(verbose){
    stop("Alpha metrix already in the metadata, check them out")
    }else{
      stop()
    }
  }

  # get the estimates
  tmp <- rarefy_even_depth(physeq, raref_depth, verbose) %>%
    estimate_richness(measures = measures)
  # put alpha in front of variables for ease of subsetting and variable calling in the metadata
  colnames(tmp) <- paste0("alpha_", colnames(tmp))

  # get variable names to merge
  alpha_df <- rownames_to_column(tmp, "IDs")

  # create metadata to merge
  new_meta <- merge(meta(physeq) %>% rownames_to_column("IDs"),
                    alpha_df,
                    by = "IDs", all.x = TRUE, sort = FALSE) %>%
    column_to_rownames("IDs")

  physeq_final <- phy_substitute_metadata(physeq, new_meta)
  physeq_final@sam_data$merge_by <- NULL


  sample_names(physeq_final) <- gsub("tmp_", "", sample_names(physeq_final))
  return(physeq_final)
}

