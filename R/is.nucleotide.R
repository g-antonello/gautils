#' Is a sequence made of nucleotides?
#'
#' @param x a character vector to test
#'
#' @return a logical, saying if the character is only mad of nucleotide-typic characters (capital and non-capital A, T, G, C, U, N). NB: it does not support other characters, like W or whatnot.
#' @export
#'
#' @examples
#'
#' is.nucleotide("ATTGCATA")
#' is.nucleotide("aTTgNGnt")
#' is.nucleotide("boyscout")
#'
is.nucleotide <- function(x){
  l <- sapply(strsplit(x, "*"), function(i) as.character(i) %in%
                c("a", "A", "c", "C", "g", "G", "t", "T", "u", "U", "n", "N"))
  s <- sapply(l, function(i) all(all(i == TRUE) ==TRUE))
  return(all(s==TRUE))

}
