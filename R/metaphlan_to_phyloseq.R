#' Create a phyloseq object from Metaphlan profiles
#'
#' Given the profiles' location, or a pre-loaded table, the function creates a taxonomic table and creates a phyloseq object. if metadata and phylogenetic tree are available, it will try to include those
#'
#' @param mpa Either a \code{character} vector containing the path to the text table to read, or the already-read object
#' @param metadata If in possession, add metadata to the phyloseq object. make sure the name and order of samples is correct for the mpa columns and metadata rows
#' @param version the metaphlan version, by default, version 4 is taken (to date it is the latest - Apr 2022). Versions below 3 are not supported
#' @param verbose \code{logical}, whether you want to see which samples did not match between the otu and metadata tables
#' @param tax_lvl \code{character}, with the taxonomic level (eg. Phylum, Class, ..., Genus, Species) that you want to have in the phyloseq
#'
#' @importFrom phyloseq phyloseq
#' @importFrom phyloseq sam_data
#' @importFrom phyloseq otu_table
#' @importFrom magrittr %>%
#' @importFrom tidyr separate
#'
#' @return A \code{phyloseq} object
#'
#' @export
#'
#' @examples
#' data("metaphlanData")
#'
#' meta_physeq <- metaphlan_to_phyloseq <- function(metaphlanData,
#' metadata = NULL,
#' phyloTree = NULL,
#' version = 4)
#'

metaphlan_to_phyloseq <- function(mpa,
                                  metadata = NULL,
                                  version = 4,
                                  verbose = TRUE,
                                  tax_lvl = "Species"){


  if(version == 4){
    if(is.character(mpa)){
      # load raw metaphlan data
      mpa <- data.table::fread(mpa) %>%
        as.data.frame()
    }
  }
    if(version == 3){
      if(is.character(mpa)){
        # load raw metaphlan data
        mpa <- data.table::fread(mpa,skip = 1) %>%
          as.data.frame()
      }
    }
  # find for each row, to which depth of taxonomy it arrives (as integers)
  tax_lengths <- mpa$clade_name %>%
    strsplit("|", fixed = T) %>%
    sapply(length)
  # remove first element,we want to keep it
  tax_lengths <- tax_lengths[-1]
  # get integer equivalent of the taxonomic level we want
  tax_lvl_int <- match(tax_lvl, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB"))
  # subset the otu table, so to keep only UNASSIGNED on top and the rows exactly to te taxonomic level
  otu_cleaned <- mpa[c(T, tax_lengths == tax_lvl_int), 3:ncol(mpa)]

  if(!is.null(metadata)){
    inters_names <- intersect(colnames(otu_cleaned), rownames(metadata))
    if(verbose){
      if(length(colnames(otu_cleaned)[!(colnames(otu_cleaned) %in% inters_names)])!= 0){
        cat("Metaphlan table samples lost: ")
        cat(colnames(otu_cleaned)[!(colnames(otu_cleaned) %in% inters_names)], sep = " ")
        cat("\n")
        cat("\n")
      }

      if(length(rownames(metadata)[!(rownames(metadata) %in% inters_names)]) != 0){
        cat("Metadata table samples lost: ")
        cat(rownames(metadata)[!(rownames(metadata) %in% inters_names)], sep = " ")
        cat("\n")
        cat("\n")
      }
    }

    otu_cleaned <- otu_cleaned[, inters_names]
    metadata_cleaned <- metadata[match(inters_names, rownames(metadata)),]

  }

  # create taxonomy table, filling empty columns with NA

  taxonomy_tab <- mpa[c(T, tax_lengths == tax_lvl_int), 1] %>%
    as.data.frame() %>%
    set_names("clade_name") %>%
    separate(col="clade_name", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "SGB")[1:tax_lvl_int],sep = "\\|", fill = "right") %>%
    as.matrix()
  # rename the UNKNOWN clade as all unknown, otherwise we can't assign rownames
  taxonomy_tab[1,] <- rep("UNKNOWN",tax_lvl_int)

  #assign rownames to the tax level chosen
  rownames(taxonomy_tab) <- taxonomy_tab[,tax_lvl_int]
  rownames(otu_cleaned) <- rownames(taxonomy_tab)

    profiles <- phyloseq(otu_table(as.matrix(otu_cleaned), taxa_are_rows = TRUE),
                         tax_table(taxonomy_tab),
                         sample_data(metadata_cleaned, errorIfNULL = FALSE)
    )

    return(profiles)
  }
