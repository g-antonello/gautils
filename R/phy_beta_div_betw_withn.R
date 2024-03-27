#' Beta diversity into group boxplots
#'
#' Function taking a phyloseq object or a diversity matrix and spits out group-wise interindividual diversity
#'
#' @param physeq a \code{phyloseq} object, important for the metadata it contains
#' @param dist either specify a metric (see \code{?phyloseq::distance} for more info), or pass a distance matrix either in the \code{matrix} or \code{dist} format
#' @param variab variable to segregate beta diversity values by. must be a \code{character} with a name in the \code{physeq}'s metadata
#' @param verbose should it print out the steps it is performing? (default is \code{FALSE})
#' @param btwn_pairwise should a third element be calculated and appended to the results? these will be pairwise comparisons of samples in each variable level
#' @param returnSampleNames \code{logical}, whether pairs of ids should be returned (default = `TRUE`)
#'
#' @importFrom magrittr %>% %<>%
#'
#' @return a \code{data.frame} object containing within and between diversities based on the variable specified. "between_pairwise" are calculated as "for each level of the variable categories, find diversities of all samples belonging to it against all other samples not belonging to it".
#' @export
#'
#'
#' @examples
#'
#' library(phyloseq)
#' data(GlobalPatterns)
#'
#' betadiversitiesFull <- phy_beta_div_betw_withn(physeq = GlobalPatterns, dist = "bray", variab = "Site", verbose = TRUE, btwn_pairwise = TRUE)
#'
#' betadiversitiesSimpler <- phy_beta_div_betw_withn(physeq = GlobalPatterns, dist = "bray", variab = "Site", verbose = TRUE, btwn_pairwise = FALSE)


phy_beta_div_betw_withn <- function (physeq,
                                     dist,
                                     variab,
                                     verbose = FALSE,
                                     btwn_pairwise = FALSE,
                                     returnSampleNames = TRUE) {
  if(!(variab %in% colnames(meta(physeq)))){
    stop("Variable not in the metadata")
  }

  if (is(dist,  "matrix") | is(dist,"dist")) {
    if (verbose) {
      message("Using pre-calculated distance matrix")
    }
    dist_mat <- as.matrix(dist)

    # filter samples to only those in common
    samples_in_common <- intersect(rownames(dist_mat), sample_names(physeq))

    dist_mat <- dist_mat[samples_in_common, samples_in_common]
    physeq <- filter_sample_data(physeq, sample_names(physeq) %in% samples_in_common)
  }

  if (is(dist,"character")) {
    if (verbose) {
      message(paste0("calculating the <", dist, "> beta diversity..."))
    }
    dist <- phyloseq::distance(microbiome::transform(physeq,
                                                     "compositional"), method = dist)
    dist_mat <- as.matrix(dist)
  }
  # subset the matrix only with sample names you have in both your objects


  # remove the symmetric part of the distance, plus the diagonal, which is uninformative
  dist_upper <- dist_mat
  dist_upper[lower.tri(dist_upper, diag = TRUE)] <- NA
  tmp <- meta(physeq) %>% rownames_to_column("samp_names") %>%
    select(all_of(c("samp_names", variab)))
  split_samp_names <- split.data.frame(tmp, tmp[[variab]],
                                       drop = TRUE)

  split_samp_names <- split_samp_names[sapply(split_samp_names, nrow)>= 2]
  rm(tmp)

  if (verbose) {
    message("Extracting within-group diversities...")
  }

  within_divs <- lapply(split_samp_names, function(s)
        # if the length of the vector is 1, it is impossible to get within diversities, so you might as well skip that
        dist_upper[s$samp_names,
                   s$samp_names] %>%
          reshape2::melt()
      ) %>%
    lapply(function(x) x[complete.cases(x),])

  if (verbose) {
    message("Extracting between-group diversities, one-vs-all")
  }

between_divs_all <-
  lapply(split_samp_names, function(s)
    dist_upper[s$samp_names,!(colnames(dist_upper) %in% s$samp_names)] %>%
      reshape2::melt()) %>%
  lapply(function(x) x[complete.cases(x),]) %>%
  .[sapply(., nrow) != 0]

between_divs_all <- between_divs_all[names(between_divs_all) %in% names(within_divs)]

if (btwn_pairwise) {
  if (verbose) {
    message("Extracting between-group diversities, one-vs-one")
  }
  samp_names_only <- purrr::transpose(split_samp_names)$samp_names
  combinations_levels <-
    combn(names(purrr::transpose(split_samp_names)$samp_names),
          2)
  between_divs_pairwise <- list()
  for (i in 1:ncol(combinations_levels)) {
    between_divs_pairwise[[paste(combinations_levels[1,
                                                     i], combinations_levels[2, i], sep = " vs ")]] <-
      dist_upper[samp_names_only[[combinations_levels[1,
                                                      i]]], samp_names_only[[combinations_levels[2,
                                                                                                 i]]]]
    between_divs_pairwise_clean <- lapply(between_divs_pairwise,
                                          function(dist)
                                            reshape2::melt(dist) %>%
                                            .[complete.cases(.),]) %>%
      .[sapply(., length) != 0]
  }
}

## final manipulation before returning the data.frame
if (btwn_pairwise) {
  betadiv_results <- rbind.data.frame(
    bind_rows(lapply(within_divs,
                     as.data.frame), .id = variab) %>% set_names(c(variab,
                                                                   "Sample1",
                                                                   "Sample2",
                                                                   "value")) %>% mutate(comparison = "within"),
    bind_rows(lapply(between_divs_all,
                     as.data.frame), .id = variab) %>% set_names(c(variab,
                                                                   "Sample1",
                                                                   "Sample2",
                                                                   "value")) %>% mutate(comparison = "between_all"),
    bind_rows(lapply(
      between_divs_pairwise_clean, as.data.frame
    ),
    .id = variab) %>% set_names(c(variab,
                                  "Sample1",
                                  "Sample2",
                                  "value")) %>%
      mutate(comparison = "between_pairwise")
  )
} else {

    betadiv_results <- rbind.data.frame(
      bind_rows(lapply(within_divs,
                       as.data.frame), .id = variab) %>% set_names(c(variab,
                                                                     "Sample1",
                                                                     "Sample2",
                                                                     "value")) %>%
        mutate(comparison = "within"),
      bind_rows(lapply(between_divs_all,
                       as.data.frame), .id = variab) %>% set_names(c(variab,
                                                                     "Sample1",
                                                                     "Sample2",
                                                                     "value")) %>%
      mutate(comparison = "between_all")
    )
  }

# final "if", asking if you want sample ID pairs back in the returned data.frame

if(!returnSampleNames){
  betadiv_results %<>% select(-Sample1, -Sample2)
}

return(betadiv_results[complete.cases(betadiv_results),])
}
