#' Get SOLAR results
#'
#' This is a quick extractor of a SOLAR run results
#'
#' @param out_dir A \code{character} specifying the base directory of the SOLAR analysis
#' @param prefix A \code{character} specifying the character prefix used in the generation of the traits. my default is `x__`, if none is used, then input `""`. This step simply remove that prefix from the trait names
#'
#' @importFrom readr parse_number
#'
#' @return A `data.frame` with estimates, P-values and so on
#'
#' @export
#'

get_SOLAR_results <- function(out_dir, prefix = "x__"){

  result_files <- list.files(out_dir, pattern = "polygenic.out", recursive = T, full.names = T)

  trait_names <- list.dirs(out_dir, recursive = F) %>%
    str_split("/",simplify = T) %>%
    .[,ncol(.)]

  files_read <- lapply(result_files, function(x) readLines(x))

  h2_and_p_line <- files_read %>%
    # extract line h2 estimates and signif
    lapply(function(f) grep("H2r is ", f, value = TRUE))

  h2_estimates <- lapply(h2_and_p_line, function(x) as.numeric(str_extract(x, "\\d+\\.*\\d* "))) %>%
    sapply(function(x) ifelse(is_empty(x),
                              NA,
                              x)) %>%
    Reduce(c, .)

  h2_signif <- lapply(h2_and_p_line, function(x) str_extract(x, "p = .*")) %>%
    sapply(function(x) ifelse(is_empty(x), NA, parse_number(x))) %>%
    Reduce(c, .)

  h2_SE <- sapply(files_read, function(x) grep("H2r Std. Error.*", x, value = T)) %>%
    sapply(function(x) ifelse(is_empty(str_extract(x,  " \\d+\\.*\\d*")),
                              NA,
                              str_extract(x,  " \\d+\\.*\\d*"))
    ) %>%
    Reduce(c, .) %>%
    as.numeric()

  ## household part
  c2_and_p_line <- files_read %>%
    # extract line c2 estimates and signif
    lapply(function(f) grep("C2 is ", f, value = TRUE))

  c2_estimates <- lapply(c2_and_p_line, function(x) ifelse(is_null(as.numeric(str_extract(x, " \\d+\\.*\\d* "))), NA, as.numeric(str_extract(x, " \\d+\\.*\\d* ")))) %>%
    Reduce(c, .)

  c2_signif <- lapply(c2_and_p_line, function(x) str_extract(x, "p = .*")) %>%
    sapply(function(x) ifelse(is_empty(x),
                              NA,
                              parse_number(x))
    )

  c2_SE <- sapply(files_read, function(x) grep("C2 Std. Error", x, value = T)) %>%
    sapply(function(x) ifelse(is_empty(x),
                              NA,
                              str_extract(x,  " \\d+\\.*\\d*"))
    ) %>%
    as.numeric()

  solar_results <-data.frame(trait = gsub(prefix, "", trait_names),
                             h2 = h2_estimates,
                             h2_p = h2_signif,
                             h2_SE = h2_SE,
                             c2 = c2_estimates,
                             c2_p = c2_signif,
                             c2_SE = c2_SE)

  return(solar_results)
}
