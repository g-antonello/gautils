#' Add metadata variables to a phyloseq
#'
#' Adds a set of metadata variables to the metadata of a phyloseq object. Make sure `rownames` contain the same names. Row ordering is taken care of internally
#'
#' @param physeq A `phyloseq` object
#' @param df A `data.frame` object, with rownames matching `sample_names` of the `phyloseq` object
#' @param by A `character` containing the variable to merge the two `sample_data` by.
#' @param verbose Should the function tell you the missing samples?
#' @export
#'
#' @importFrom phyloseq prune_samples
#' @importFrom microbiome meta
#' @importFrom microbiome abundances
#' @import tidyverse
#'
#' @examples
#' data(GlobalPatterns)
#' df_test <- data.frame(variable_new = rnorm(n = nsamples(GlobalPatterns), mean = 2, sd = 0.8))
#' rownames(df) <- samples_names(GlobalPatterns)
#' phy_add_metadata_variables(physeq=GlobalPatterns,
#'                            df = df_test)
#'
phy_add_metadata_variables <- function(physeq, df, by, verbose = FALSE){

  tmp_meta <- meta(physeq) %>%
    rownames_to_column("sampNamestmp")

  intersect_names <- intersect(tmp_meta[[by]],
                               df[[by]])

  if(length(intersect_names) == 0){stop("Sample names don't match. Did you pick the right column to merge by?")}

  # report which phyloseq samples are not in the intersection
  if(!all(tmp_meta[[by]] %in% intersect_names)){ # if not all sample names are in the intersection
    missing_samples <- tmp_meta[[by]][!(tmp_meta[[by]] %in% intersect_names)]

    cat("\n")
    cat(paste(length(missing_samples), "samples missing in physeq:\t"))
    cat(missing_samples, sep = " | ")
    cat("\n")
    # filter the original metadata
    tmp_meta <- filter(tmp_meta, eval(parse(text = by)) %in% intersect_names)
  }
  # report which df samples are not in the intersection
  if(!all(df[[by]] %in% intersect_names)){ # if not all sample names are in the intersection
    missing_samples <- df[[by]][!(df[[by]] %in% intersect_names)]

    if(verbose){
      cat("\n")
      cat(paste(length(missing_samples), "samples missing in physeq:\t"))
      cat(missing_samples, sep = " | ")
      cat("\n")
    }

    #filter the df
    df<- filter(df, eval(parse(text = by)) %in% intersect_names)
  }

  new_meta <- left_join(tmp_meta, df, by = by) %>%
    column_to_rownames("sampNamestmp")

  physeq_filt <- prune_samples(samples = sample_names(physeq) %in% rownames(new_meta), x = physeq)

  return(phy_substitute_metadata(physeq_filt, new_meta))
}


