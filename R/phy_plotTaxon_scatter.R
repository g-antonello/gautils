#' Scatterplot of a taxon
#'
#' Make sure the "x" parameter calls a continuous variable
#'
#' @param physeq a phyloseq object
#' @param y the name of the taxon you want to plot (y axis)
#' @param x the name of the trait you want to plot (x axis)
#' @param colour another variable you may want to color your points by
#' @param transform a transformation, all those in the phy_transform are allowed
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' data("enterotype")
#' enterotype@sam_data$otherVar <- rnorm(n = nsamples(enterotype), mean = 5, sd = 1.5)
#'
#' phy_plotTaxon_scatter(physeq = enterotype,
#'                      y = "Rikenella",
#'                      x = "otherVar",
#'                      colour = "SeqTech",
#'                     transform = "clr")
#'


phy_plotTaxon_scatter <- function(physeq,
                              y,
                              x,
                              colour = NULL,
                              transform= NULL
                          ){

  if(!(class(physeq@sam_data[[x]]) %in%
       c("integer", "double", "numeric", "Date", "character", "factor"))){
    stop("x select is of a non-compatible format")
  }

  real_name <- y
  easier_name <- make.names(real_name)

  # transform otu counts
  if(!is.null(transform)){
    physeq_transf <- gautils2::phy_transform(physeq, transform = transform)
    taxa_names(physeq_transf) <- make.names(taxa_names(physeq_transf))
  }else{
    physeq_transf <- physeq
  }

  data_merged <- phy_OtuMetaTable(physeq_transf)

  scatpl <- ggplot(data = data_merged, aes_string(x = x, y = easier_name)) +
      geom_point(aes_string(colour=colour)) +
      labs(title = real_name,
           caption = paste0("transformation: ", transform),
           y = "Taxon Abundance"
        )

  return(scatpl)

}
