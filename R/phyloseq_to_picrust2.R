#' Phyloseq to Picrust input tables
#'
#' @param physeq \code{phyloseq} object
#' @param output.dir The directory path to save your data into. NB: the function can create it
#'
#' @importFrom microbiome abundances
#' @importFrom biomformat make_biom
#' @importFrom biomformat write_biom
#' @importFrom dplyr select
#'
#' @return a character vector saying where the files were stored
#' @export
#'
#' @examples
#'
#' #There are not example datasets with also unique sequences, and I cannot provide our data due to our strict data protection policy
#'

phyloseq_to_picrust2 <- function(physeq, output.dir){

  # create the folder if required
  dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

  # write biom-formated ASV table
  write_biom(x = make_biom(abundances(physeq)),
                         biom_file = paste0(output.dir, "/", deparse(quote(physeq)), "_asv_table.biom"))
  # write refseq into FASTA (actually, .fna) file
  phyloseq_to_fasta(physeq = physeq,
                    destination_path =file.path(output.dir),
                    filename = "refseqs.fna")

  cat(paste("\n", "PICRUSt2 input files saved into", output.dir, "\n"))
}
