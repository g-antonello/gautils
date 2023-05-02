#' Transform microbiota counts
#'
#' Same as transform(), but also adding further transformations, like hyperbolic arcsine, inverse-rank normal, binary transformation.
#' Additionally it can calculate geometric mean of a sample, although this is not really the same thing
#'
#' @param physeq a \code{phyloseq} object
#' @param transform \code{character} all options available in "transform()", extending it with the \code{arcsinh} ("asinh"), \code{geometric mean}("gm_mean"), \code{Hellinger} ("hellinger"), \code{binary}, and \code{Rank Inverse Normal transf.}("irn")
#' @param binary_preval_thresh \code{numeric}, the prevalence below which you want to transform your taxa into binary traits. default is \code{0}, which will not split the otu table
#'
#' @importFrom RNOmni RankNorm
#' @importFrom microbiome abundances
#' @importFrom microbiome transform
#' @importFrom microbiome prevalence
#'
#' @return A \code{phyloseq} object with sample-wise transformed data. also, it will create a phyloseq with Taxa in rows!
#'
#' @export
#'
#' @examples
#'
#' data(GlobalPatterns)
#' GP_Genus <- tax_glom(GlobalPatterns, "Genus")
#' tax <- "295395"
#' transf <- "clr"
#' transformed_GP <- transform_microbiota_ga(GP_Genus,transform = transf)
#'
#'
#' plot(density(abundances(GP_Genus)[tax,]), main = paste("taxon n°", tax), col = "blue")
#' lines(density(abundances(transformed_GP)[tax,]), col = "red")
#' legend("topright", legend = c("Untransformed", transf), col = c("blue", "red"), pch = 18)
#'

phy_transform <- function (physeq, transform, binary_preval_thresh = 0)
{
  prevalences <- prevalence(physeq)
  ntaxa <- ntaxa(physeq)
  if (any(prevalences < binary_preval_thresh)) {
    otu_base <- abundances(physeq)
    otu_above <- otu_base[prevalences >= binary_preval_thresh,
    ]
    otu_below <- otu_base[prevalences < binary_preval_thresh,
    ]
    if (!(transform %in% c("asinh", "irn"))) {
      otu_above_transf <- t(transform(otu_table(otu_above,
                                                taxa_are_rows = T), transform))
    }
    if (transform == "asinh") {
      otu_above_transf <- apply(otu_above, 1, asinh)
    }
    if (transform == "irn") {
      otu_above_transf <- apply(otu_above, 2, function(x) x/(sum(x))) %>%
        apply(2, RNOmni::RankNorm)
    }
    otu_below_transf <- apply(otu_below, 1, function(x) ifelse(x >
                                                                 0, 1, 0))
    otu_transf_ready <- cbind(otu_above_transf, otu_below_transf)[sample_names(physeq),
                                                                  taxa_names(physeq)]
  }
  else {
    if (!(transform %in% c("asinh", "irn"))) {
      otu_transf_ready <- t(abundances(transform(physeq,
                                                 transform)))
    }
    if (transform == "asinh") {
      otu_transf_ready <- apply(abundances(physeq), 1,
                                asinh)
    }
    if (transform == "irn") {
      otu_transf_ready <- apply(abundances(physeq), 2, function(x) x/(sum(x))) %>%
        apply(2, RNOmni::RankNorm)
    }
  }
  physeq_transf <- physeq
  if (taxa_are_rows(physeq)) {
    otu_table(physeq_transf) <- otu_table(t(otu_transf_ready),
                                          taxa_are_rows = TRUE)
  }
  else {
    otu_table(physeq_transf) <- otu_table(otu_transf_ready,
                                          taxa_are_rows = FALSE)
  }
  return(physeq_transf)
}

