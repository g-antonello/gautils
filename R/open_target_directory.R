#' Open file/plot directory
#'
#' Given a PATH, either ending with "/" or a file name, this function opens a window to show you that directory.
#'
#' @param element_path The path, written unix-style. it has to be absolute, either from the root, or from the std. home ("~/")
#' @return A window opened on your file browsing element ("Finder" or "Explorer")
#' @export
#'
#' @examples
#' # in Unix and Unix-like
#' open_target_directory("/home/user/Documents/my_superImportantDoc.txt")
#' # in Windows
#' open_target_directory("C:/Users/user/Documents/my_superImportantDoc.txt")
#'
#'
open_target_directory <- function(element_path){

  if (.Platform$OS.type == "windows"){
    if(!(str_sub(element_path, start = nchar(element_path), end = nchar(element_path)) == "/")){ # if the directory name given ends with "/", then we can play with it without manipulation, otherwise, we need to do some managing
      s <- strsplit(element_path, "/")[[1]]
      s_sub <- s[1:length(s)-1]
      element_path_windows <- paste(s_sub, collapse = "\\") %>% paste("\\", sep = "")
    } else{ # if it ends with "/", just translate the "/" character
      element_path_windows <- gsub(pattern = "/", replacement = "\\\\", x = element_path)
    }

    shell.exec(element_path_windows)
  } else{ # if it is not windows
    if(!(str_sub(element_path, start = nchar(element_path), end = nchar(element_path)) == "/")){
      s <- strsplit(element_path, "/")[[1]]
      s_sub <- s[1:length(s)-1]
      element_path_unix <- paste(s_sub, collapse = "/") %>% paste("/", sep = "")
      system(paste(Sys.getenv("R_BROWSER"), element_path_unix))
    } else{
      system(paste(Sys.getenv("R_BROWSER"), element_path))
    }

  }
}
