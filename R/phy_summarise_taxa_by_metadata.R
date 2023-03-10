#' Summarise taxa by metadata
#'
#' A quick dplyr function that allows to summarise (mean, median, sum, ...) all taxa for each category of a metadata column. So far tested only on discrete variables.
#'
#' @param physeq A \code{phyloseq} object
#' @param metadata_var A \code{character} with the name of the factor to summarise by. it must be accepted by the code{group_by} dplyr function
#' @param .fun a \code{function} object to call by the summarise function in dplyr. I am thinking \code{mean, median, sum, scale, ...}
#'
#' @return A \code{matrix} object with the metadata variable in columns, and taxa in rows
#' @export
#'
#' @examples
#'
#' data(enterotype)
#' matrix_for_plot <- phy_summarise_taxa_by_metadata(microbiome::transform(enterotype, "clr"), metadata_var = "Enterotype", .fun = mean)
#' pheatmap::pheatmap(matrix_for_plot)
#'
phy_summarise_taxa_by_metadata <- function(physeq, metadata_var, .fun = mean){
  phyloseq_summarised <- psmelt(physeq) %>%
    group_by(OTU, eval(parse(text = metadata_var))) %>%
    summarise(feature_summary = mean(Abundance)) %>%
    ungroup() %>%
    set_names(c("feature", "metadata_var", "feature_summary")) %>%
    reshape2::acast(value.var = "feature_summary", formula = metadata_var ~ feature) %>%
    t()
    return(phyloseq_summarised)
}
