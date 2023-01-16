#' merge otu and metadata table
#'
#' @param physeq a phyloseq object with otu and metadata tables available
#'
#' @return a data.frame object, with otus as first columns, and metadata merged after them
#' @export
#'
#' @examples
#'
#' data(GlobalPatterns)
#'
#' otuMeta <- OtuMetaTable(GlobalPatterns)

phy_OtuMetaTable <- function(physeq){
  otu <- t(microbiome::abundances(physeq)) %>%
    as.data.frame() %>%
    rownames_to_column("sampIDs")

  meta <- microbiome::meta(physeq) %>%
    rownames_to_column("sampIDs")

  otuMeta <- merge(otu, meta, by = "sampIDs", sort = FALSE)

  return(otuMeta)
}
