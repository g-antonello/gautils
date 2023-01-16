
#' Barplot taxa (beta)
#'
#' A not-fully-finished but "MVP" (minimum viable product) function that plots phyloseq data preferrably in relative abundance.
#'
#' @param physeq A phyloseq object
#' @param x The x axis variable, default is "Sample", which will plot one column for each sample
#' @param tax_fill The taxonomic level to fill bars by
#' @param title A character vector to entitle the plot with
#' @param facet A character vector of a metadata variable's name
#' @param transform A character vector of the transformation method to apply
#' @importFrom speedyseq psmelt
#'
#' @return A \code{ggplot2} plot, ready to be printed
#'
#' @export

#' @examples
#'
#' data("GlobalPatterns")
#'
#' testplotSample <- phy_barplotTaxa(physeq = GlobalPatterns,
#'                                  x = "Sample",
#'                                  tax_fill = "Genus",
#'                                  title = "By Sample",
#'                                  facet = "SampleType",
#'                                  transform = "compositional")
#'
#' testplotSample <- phy_barplotTaxa(physeq = GlobalPatterns,
#'                                  x = "SampleType",
#'                                tax_fill = "Genus",
#'                                title = "By Sample",
#'                                transform = "compositional")
#'

phy_barplotTaxa <-  function(physeq,
                        x = "Sample",
                        tax_fill = NULL,
                        title = NULL,
                        facet = NULL,
                        transform = "compositional"){

    if(x != "Sample"){
      if(!any(c(facet, x) %in% sample_variables(physeq))){
        absent_vars <- c(facet, x)[c(facet, x) %in% sample_variables(physeq)]
        stop(paste0("Metadata doesn't contain the following variables: ",  absent_vars))
      }
    }

    if(x == "Sample"){
      if(!(x %in% sample_variables(physeq))){
        stop(paste0("Cannot facet, missing variable: ",  x))
      }
    }


  y <- "Abundance"

  physeq@phy_tree <- NULL
  physeq_transf <- tax_glom(physeq, tax_fill) %>%
    phy_transform(transform = transform)

  ##### plotting phase
  mdf <- psmelt(physeq_transf)
  if(x == "Sample"){
  p <- ggplot(mdf, aes_string(x = x, y = y, fill = fill))+
    geom_bar(stat = "identity", position = "stack",  color = "transparent")+
    theme(axis.text.x = element_text(angle = -90, hjust = 0))+
    ggtitle(title)

  if (!is.null(facet)) {
    p <- p + facet_grid(cols = vars(eval(parse(text = facet))))
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)
    }

  }

  if(x != "Sample"){
    # make mean of all samples in each group dictated by x, which has to be the name of a column of the metadata
  summarized_mdf <- mdf %>%
      group_by(eval(parse(text = tax_fill)),eval(parse(text = x))
               ) %>%
      summarise(Abundance = mean(Abundance)) %>%
      set_names(c(tax_fill, x, "Abundance"))

  p <- ggplot(summarized_mdf, aes_string(x = x, y = "Abundance", fill = tax_fill))+
    geom_bar(stat = "identity", color = "transparent")+
    theme(axis.text.x = element_text(angle = -90, hjust = 0))+
    ggtitle(title)

  if (!is.null(facet)) {
    p <- p + facet_grid(cols = vars(eval(parse(text = facet))))
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  }

  return(p)
}
