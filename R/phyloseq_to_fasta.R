#' Save Refseqs of phyloseq into a .fasta
#'
#' Write reference sequences of a phyloseq object to a text file with fasta formatting
#'
#' @param physeq A \code{phyloseq} object
#' @param destination_path Write the file you want your fasta sequences to be written into. NB: the path must exist already
#' @param filename A \code{character}, the name of the file to write the files into. NB: the extension is not written automatically
#'
#' @return either a \code{character} with the sequences named after the taxon ID, or it writes them into the specified file in Fasta format
#'
#' @export
#'
#'
phyloseq_to_fasta <- function(physeq,
                              destination_path = NULL,
                              filename = NULL){

  # format option 1: sequences are somewhere in the taxa table
  if(is.null(physeq@refseq)){
    tx <- as.data.frame(physeq@tax_table@.Data)
    pos <- sapply(tx, is.nucleotide)
    if(!any(pos)){stop("there is no reference sequence in the taxonomic table")}
    seqs <- as.character(tx[, pos])
    names(seqs) <- rownames(tx)

    if (is.null(destination_path)){
      return(seqs)
      stop("no destination file for saving refseq data")
    }
  }
  if(!is.null(physeq@refseq)){ # format 2, there is a refseq
    seqs <- as.character(refseq(physeq))

    if(is.null(destination_path)){
      return(refseq(physeq))
      stop("no destination file for saving refseq data")
    }

  }

  # write the fasta file

  dir.create(destination_path, recursive=TRUE, showWarnings = FALSE)
  if(is.null(filename)){
    filename <- deparse(quote(physeq))
  }

  cat(file = file.path(destination_path, filename),
      paste(">", names(seqs), "\n", seqs, "\n", sep = ""),
      sep = "")

  cat(paste("\nfasta file written: '", tools::file_path_as_absolute(file.path(destination_path, filename)), "'", sep = ""))
}
