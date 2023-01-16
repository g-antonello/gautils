#' Heatmap of microbiota features, NOT EXPORTED, because I haven't tested it yet!!
#'
#' Make a quick heatmap of selected features, also transforming internally
#'
#' @param physeq A phyloseq object
#' @param variable The variable to consider for clustering
#' @param physeq_transformation The transformation you want to apply to your data. Default is \code{identity}, which is no transformation
#' @param taxa A facultative list of taxa to subset the phyloseq by. Default is \code{all}
#' @param aggregation_samples If you want to perform variable-wise aggregation. the allowed aggregation methods here are \code{mean, median, none}. NB: \code{none} doesn't work yet
#' @param taxa_into_rows A \code{logical}. It asks whether you want taxa as rows or columns of the heatmap. make sure you do the scaling right, then, in the next parameter
#' @param scale The \code{scale} option in \code{pheatmap}, either scale by row or column. default is "row", as also \code{taxa_into_rows} default is TRUE
#' @param add_sample_size A \code{logical}. It asks whether you want to add the per-category sample size in your variable's labels
#' @param angle_col A \code{numeric/character} value of the degrees of tilting of labels. Justification is taken care of automatically
#' @param legend A \code{logical}. It asks whether you want to draw the legend or not
#' @param treeheight_row A \code{numeric} value for drawing the height of the row-wise tree clustering
#' @param treeheight_col A \code{numeric} value for drawing the height of the column-wise tree clustering
#' @param legend_labels A \code{character} vector of length 3 that describes the 3 labels to print on the legend. These labels will be printed as bottom, middle, top
#'
#' @importFrom microbiome meta
#' @importFrom pheatmap pheatmap
#' @importFrom speedyseq select_sample_data
#'
#' @return A \code{pheatmap} plot object. you can convert it into a ggplot2 object (to some extent) with \code{pheatmap_ggplot <- ggplotify::as.ggplot(returned_pheatmap_object)[[2]]}
#'
#'
#' @examples
#' data(enterotype)
#'
#' phy_heatmap(physeq = enterotype,
#'             variable = "SeqTech",
#'             aggregation_samples = "median",
#'             taxa_into_rows = TRUE,
#'             scale = "row",
#'             add_sample_size = TRUE,
#'             angle_col = 45,
#'             legend = FALSE)
#'

phy_heatmap <- function(physeq,
                        variable,
                        physeq_transformation = "identity",
                        taxa = "all",
                        aggregation_samples = "median",
                        taxa_into_rows = TRUE,
                        scale = "row",
                        add_sample_size = TRUE,
                        angle_col = "90",
                        legend = FALSE,
                        treeheight_row = 20,
                        treeheight_col = 20,
                        legend_labels = NULL) {

  if (any(taxa != "all")) {
    physeq <- prune_taxa(taxa_names(physeq) %in% taxa, physeq)
    taxa <- taxa_names(physeq)
  }

  # transform the microbiome data as requested
  physeq_transf <-
    phy_transform(physeq, transform = physeq_transformation)

  # by default,  the function will plot the number of samples per category
  if (add_sample_size){
    sample_n <- meta(physeq_transf)[[variable]] %>%
      table() %>%
      as.data.frame() %>%
      set_names(c(variable, "counts")) %>%
      mutate(variable_with_n  = paste0(eval(parse(text = variable)), " (n=", counts, ")"))

    physeq_transf <- phy_add_metadata_variables(physeq_transf, sample_n, by = variable)
    # reassign original "variable" object, so that we can write one single code afterwards
    variable <- "variable_with_n"
  }

  # split ASV table based on variable chosen
  if (aggregation_samples %in% c("median", "mean")) {
    if(length(variable)>1){
      stop("multiple variables are supported only in sample-wise plotting")
    }

    metadata_and_asvs <-  select_sample_data(physeq_transf, variable) %>%
      phy_OtuMetaTable() %>%
      select(-sampIDs) %>%
      split.data.frame(f = .[[variable]], drop = TRUE) %>% # split the data frame and keep names
      lapply(function(df)
        select(df,-matches(variable))) # remove leftover variable used to split

    if (aggregation_samples == "median") {
      aggregated_data <-
        lapply(metadata_and_asvs, function(df)
          apply(df, 2, median, na.rm = FALSE)) %>%
        bind_rows(.id = variable) %>%
        column_to_rownames(variable)

      if(any(colSums(aggregated_data) == 0)){
        warning("removing some taxa, since median values of all columns are 0")
        aggregated_data <- aggregated_data[, colSums(aggregated_data) != 0]
      }


    }
    if (aggregation_samples == "mean") {
      aggregated_data <-
        lapply(metadata_and_asvs, function(df)
          apply(df, 2, mean, na.rm = FALSE)) %>%
        bind_rows(.id = variable) %>%
        column_to_rownames(variable)
    }

    # transpose to have taxa in rows or columns, as preferred
    scale_options <- c("row", "column")
    if(taxa_into_rows){ # if the variable is TRUE
      data_ready_to_plot <- t(aggregated_data)
      scale <- scale_options[!(scale == scale_options)]

    }else{
      data_ready_to_plot <- aggregated_data
      if(scale != "none"){

      }
    }

    if(!is.null(legend_labels)){
    return(pheatmap(
      data_ready_to_plot,
      # don't show colnames, which with this settings are sample names
      treeheight_row = treeheight_row,
      treeheight_col = treeheight_col, # don't draw trees
      scale = "row",
      angle_col = angle_col,
      legend = legend,
      legend_labels = legend_labels))
    }else{
      return(pheatmap(
        data_ready_to_plot,
        # don't show colnames, which with this settings are sample names
        treeheight_row = treeheight_row,
        treeheight_col = treeheight_col, # don't draw trees
        scale = "row",
        angle_col = angle_col,
        legend = legend))
    }
  }


  if (aggregation_samples == "none") {
    data_ready_to_plot <- t(abundances(physeq_transf))
    color_variable <-
      select(microbiome::meta(physeq_transf), all_of(variable))

    if(!is.null(legend_labels)){
    return(pheatmap(
      data_ready_to_plot,
      annotation_col = color_variable,
      show_colnames = FALSE,
      # don't show colnames, which with this settings are sample names
      treeheight_row = treeheight_row,
      treeheight_col = treeheight_col, # don't draw trees
      scale = scale,
      angle_col = angle_col,
      legend = legend,
      legend_labels = legend_labels))
    } else{
      return(pheatmap(
        data_ready_to_plot,
        annotation_col = color_variable,
        show_colnames = FALSE,
        # don't show colnames, which with this settings are sample names
        treeheight_row = treeheight_row,
        treeheight_col = treeheight_col, # don't draw trees
        scale = scale,
        angle_col = angle_col,
        legend = legend))
    }
  }
}
