#' Fancy MDS Plots
#'
#' The function wraps the phyloseq functions "distance", "ordinate", and "plot_ordination".
#' This function is useful to see the gradient of abundance of a taxon in the 2D ordinated dataset
#'
#' @param physeq A phyloseq object, **CAREFUL**: the function doesn't transform internally, so if you choose to calculate the distance matrix internally with "ordinate", make sure you transform the counts prior to function call (usually relative abundance is accepted)
#' @param dist distance object or character vector saying the distance metric to apply
#' @param method one of the methods allowed by "ordinate" in the phyloseq object
#' @param axes The axes to be plotted, default are the first 2, which should discriminate samples better
#' @param taxon the taxon name to color the samples by. The palette is "viridis_c"
#'
#' @importFrom grDevices chull
#' @import tidyverse
#' @importFrom microbiome abundances
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' data("enterotype")
#'
#' phy_betaPlotTaxon(physeq = enterotype,
#' dist = "bray",
#' method = "PCoA",
#' axes = 1:2,
#' taxon = "Bacteroides")

phy_betaPlotTaxon <- function(physeq,
                     dist = "bray",
                     method = "PCoA",
                     axes = 1:2,
                     taxon){
# generic preparatory phase, getting data ordinated

ord <- ordinate(physeq = physeq, method = method, distance = dist)

# retrieve data used to plot
basic_plot_data <- as.data.frame(
  cbind(
    plot_ordination(physeq = physeq,
                    ordination = ord,
                    axes = axes,
                    justDF = TRUE),
    t(abundances(physeq))
    )
  )
plot0 <- plot_ordination(physeq = physeq,
                         ordination = ord,
                         axes = axes,
                         color = NULL)
## new labels
prop_var <- c(plot0$labels$x,
              plot0$labels$y) %>%
  strsplit(., split = " ") %>%
  sapply("[[", 4)

colnames_new <- c(paste0("Axis ", axes[1], "  ",prop_var[1]),
                  paste0("Axis ", axes[2], "  ",prop_var[2])
                  )

colnames(basic_plot_data)[1:2] <- c("tmp1", "tmp2")

######## from now on I will plot case by case

final_plot <-  ggplot(basic_plot_data, aes_string(x = "tmp1",
                                           y = "tmp2",
                                           color = taxon))+
                          geom_point()+
                          scale_color_viridis_c()
final_plot$labels$x <- colnames_new[1]
final_plot$labels$y <- colnames_new[2]

return(final_plot)

}



