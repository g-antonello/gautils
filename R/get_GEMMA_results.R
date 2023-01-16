#' Get Gemma heritability results
#'
#' From a directory populated with GEMMA associations, this will read the chip heritability estimates
#'
#' @param base_dir \code{character} with the GEMMA base directory
#' @param prefix \code{chracter} specifying the trait prefix used to make trait names easily parsable. The function strips those prefixces out
#'
#' @importFrom readr parse_number
#'
#' @return A `data.frame`
#' @export
#'

get_GEMMA_results <- function(base_dir, prefix = "x__"){

  files_list <- list.files(base_dir, pattern = "result.log.txt", recursive = TRUE, full.names = TRUE)
  traits_names <- list.files(base_dir, pattern = "result.log.txt", recursive = TRUE)
  traits_names <- sapply(strsplit(traits_names, "/"), "[[", 1)

  files_raw <- lapply(files_list, function(f) readLines(f))
  names(files_raw) <- traits_names

  h2_estimates <- parse_number(sapply(files_raw, function(x) grep("pve estimate", x, value = T)))
  h2_SE <- parse_number(sapply(files_raw, function(x) grep("se\\(pve\\)", x, value = T)))

  return(data.frame(trait = gsub(prefix, "", traits_names),
                    h2_gemma  = h2_estimates,
                    h2_gemma_SE = h2_SE))
}
