
#' Agglomerate at higher taxonomic rank/level
#'
#' Function that wraps speedyseq's tax_glom function and renames taxa based on more readable names as seen in the taxonomy table
#'
#' @param physeq A \code{phyloseq} object
#' @param level The taxonomic rank to aggregate to. Make sure it matches exactly the column name in the taxonomy table (see \code{rank_names(physeq)})
#' @param make.unique whether you should make unique the names. Default is \code{TRUE}. The standard separation pattern is \code{_.}


#' @importFrom speedyseq tax_glom
#' @importFrom magrittr %<>%

#' @return An aggregated \code{phyloseq} object
#'
#' @export
#'
#' @examples
#'
#' library(gautils2)
#' data("GlobalPatterns")

#' gp_Family <- phy_tax_glom(GlobalPatterns, "Family", make.unique = FALSE)

#' # this code shoots a "duplicated names" error. But you can make these names unique artificially:
#'
#' gp_Family <- phy_tax_glom(GlobalPatterns, "Family", make.unique = TRUE)
#'
#'

phy_tax_glom <- function(physeq, level, prefix = NULL, make.unique = TRUE){
  phy_agglomerated <- tax_glom(physeq, level)

  if(make.unique){
    taxa_names(phy_agglomerated) <- make.unique(tax_table(phy_agglomerated)[, level], sep = "_.")
  }else{
    taxa_names(phy_agglomerated) <- tax_table(phy_agglomerated)[, level]
  }

  if(!is.null(prefix)){
    taxa_names(phy_agglomerated) %<>% paste0(prefix, .)
  }

  return(phy_agglomerated)
}
