#' Delete temporary or test objects in `ls()`
#'
#' This function deletes all items (object, function, ...) containing either 'tmp' or 'test' in the name. Also, if the name is only 1 letter long, it will be deleted
#'
#' @param extra_rules A \code{character} that defines extra patterns to delete (eg: extra_rules = c("temp","delete")). NB: RegEx are also allowed
#' @param verbose A \code{logical} whether you want the function to print the objects it removed
#'
#' @export
#'
#' @examples
#'
#' papaya <- "good"
#' first_attempt <- "eat"
#' second_attempt <- "digest"
#'
#' declutter_Envir(extra_rules = c("papaya", "attempt"))
#'
#' ls()
#' # you should only see digest

declutter.envir <- function(extra_rules = "test", verbose = FALSE){
  all_rules.t0 <- paste(extra_rules,"tmp", sep = "|")

  # get all names in the environment matching the patterns
  obj_to_del <- grep(all_rules.t0, ls(envir = parent.frame(1)), value = TRUE)
  #
  #obj_to_del_grep <- obj_to_del[2:length(obj_to_del)]

  obj_to_del_final <- c(obj_to_del,
                        ls(envir = parent.frame(1))[nchar(ls(envir = parent.frame(1))) ==1]
                        )

  rm(list = obj_to_del_final, envir = parent.frame(1))
  objs.t1 <- ls(envir = parent.frame(1))
  if(verbose){
    cat("elements removed:\n")
    cat(objs.t0[!(objs.t0 %in% objs.t1)], sep = "  ")
  }
}
