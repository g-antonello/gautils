#' Remove taxa that have counts summing to n
#' 
#' A simple wrapper function usable in a pipe that allows for quick pruning of a phyloseq object
#' 
#' @param physeq A \code{phyloseq} object
#' @param n A \code{numeric} number of the pruning threshold to filter out any taxon that sums to strictly less than that value
#' 
#' @return A \code{phyloseq} object with only taxa with \code{taxa_sums(physeq) > 0}
#' 
#' @export
#'
#' @examples
#' 
#' data(enterotype)
#' 
#' entero_noPyro <- subset_samples(enterotype, SeqTech != "Pyro454")
#' 
#' # before pruning
#' entero_noPyro
#' 
#' # after pruning
#' prune0taxa(entero_noPyro)
#' 


phy_pruneTaxa <- function(physeq, n = 0){
  return(prune_taxa(taxa_sums(physeq) > n, physeq))
}
