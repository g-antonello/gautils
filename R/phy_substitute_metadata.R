#' Substitute a phyloseq metadata with a new data.frame
#'
#' Substitute metadata in phyloseq objects. the two data frames must have row names matching sample names. The sample order doesn't matter, because the functions will reorder the new metatadata rows, to match the sample order in the phyloseq object.
#'
#' @param physeq A phyloseq object contaning or not containing `sample_data`
#' @param new_metadata a data.frame object, with rownames containing `sample_names` of the `physeq` object
#'
#' @return a phyloseq object with the new metadata stored under the "sam_data" slot
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' exMeta <- microbiome::meta(GlobalPatterns)
#' newMeta <- exMeta %>% mutate(newVar = runif(n = nrow(exMeta), min =  0, max = 20))
#'
#' GlobalPatterns <- phy_substitute_metadata(physeq = GlobalPatterns, new_metadata = newMeta)
#'
phy_substitute_metadata <- function(physeq, new_metadata){
  if(any(!rownames(new_metadata)%in%sample_names(physeq))){
    stop ("not all row names of the new metadata are in the sample names of the phyloseq object")
  }
  # reorganize metadata
  new_metadata <- new_metadata[sample_names(physeq),]
  # remove old metadata
  physeq@sam_data <- NULL
  sample_data(physeq) <- sample_data(new_metadata)
  return(physeq)
}
