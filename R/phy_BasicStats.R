#' Basic Phyloseq stats
#'
#' Generate a data frame with prevalence, mean abundance and variance of the taxa in a phyloseq object
#'
#' @param physeq Either a \code{phyloseq} or a  \code{otu.table} or a \code{matrix} object. If you have a matrix, make sure taxa are on the rows
#' @param transform Any transformation function allowed by phy_transform
#'
#' @return A \code{data.frame} object
#' @export
#'
#' @examples
#'
#' data(GlobalPatterns)
#' gp_genus <- tax_glom(GlobalPatterns, "Genus")
#' basic_statistics <- phy_BasicStats(gp_genus, transform = "clr")
#'

phy_BasicStats <- function (physeq, transform) {

  preval <- prevalence(physeq)

  physeq_transf <- phy_transform(physeq, transform)

  otumat <- abundances(physeq_transf)

  mean_abunds <- apply(otumat, 1, mean)
  median_abunds <- apply(otumat, 1, median)

  variance_abunds <- apply(otumat, 1, var)


  basic_stats_tibble <- tibble(
    names(preval),
    preval,
    mean_abunds,
    median_abunds,
    variance_abunds)

  colnames(basic_stats_tibble) <- c("taxon",
                                    paste0("prevalence_", transform),
                                    paste0("meanAbund_", transform),
                                    paste0("medianAbund_", transform),
                                    paste0("varAbund_", transform))

  return(basic_stats_tibble)
}
