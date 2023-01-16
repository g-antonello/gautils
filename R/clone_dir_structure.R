#' Clone directory structure (Not working yet)
#'
#' Takes a directory and clones its structure elsewhere, but without copying files
#' @param start_dir a character vector with a path to an existing directory
#' @param end_dir a character vector with a path to an existing directory
#'
#' @importFrom stringr str_remove
#'
#' @return Nothing, it will simply copy the directory subtree
#'
#' @export
#'
#' @examples
#'
#' dir.create("~/tmp/tmp2/tmp2.1", recursive = TRUE)
#' dir.create("~/tmp/tmp3/tmp3.1", recursive = TRUE)
#'
#' dir.create("~/destination_dir")
#' clone_dir_structure(start_dir = "~/tmp/",
#'                     end_dir = "~/destination_dir/")
#'
#' dir()
#'
#'
clone_dir_structure <- function(start_dir, end_dir){

  if (!endsWith(start_dir, "/")){ # need to remove "/" if it is at the end of the input text
    start_dir <- paste(start_dir, "/", sep = "")
  } # strip that character from the end of the start_dir

  d <- list.dirs(start_dir)
  d_clean <- str_remove(d, start_dir)
  d_clean <- d_clean[nchar(d_clean) != 0]

  while(any(startsWith(d_clean, "/"))){
    d_clean[startsWith(d_clean, "/")] <- substr(x = d_clean[startsWith(d_clean, "/")],
                                                start = 2,
                                                stop = nchar(d_clean[startsWith(d_clean, "/")]))
  }

  end_directories <- paste(paste(end_dir, d_clean, sep = ""), "/", sep = "")
  sapply(end_directories, function(i) dir.create(i, recursive = TRUE))
}
