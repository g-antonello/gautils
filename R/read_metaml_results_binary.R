#' Read MetaML binary classification output
#'
#' So far it works only with binary classification using Random Forest, and it assumes the output to have top 25 features too.
#'
#' @param output_dir A \code{character} vector to locate the output directory (see the code you used to launch the MetaML algorithm)
#' @param output_prefix A \code{character} vector to know which name prefix you gave to your files
#'
#' @return A \code{list} containing all statistics divided into "full" and "top25" dataset
#'
#' @export
#'
#' @examples
#'
#' ##### say you run something like the following code:
#'
#' # python metaml/classification_thomas-manghi.py
#' #    input_data/data_ready.txt # input file with all variables needed in rows
#' #   output/Exp1 # output prefix
#' #   -l rf # type of learning (my function was tested only on rf)
#' #   -z s__ # prefix of microbiome- related variables
#' #   -mf 0.1 # number of ...?
#' #   -nt 1000 # number of trees
#' #    -nsl 10 # number of samples per tree leaf
#' #   --define 1:varible_of_interest:category_of_interest # the rest becomes 0
#' #   -c entropy # way to estimate the variable importance
#' #   -p 10 # number of folds for cross-validation
#' #   -r 10 # number of total runs. NB: each run is independent, while folds within a run are not independent
#' #    -cc 25 # re-do the analysis with this number of top features, append results on the same files
#'
#' ##### when this is done, go back to R and do:
#'
#' # my_ML_results <- read_metaml_results_binary("output/", "Exp1")
#'

read_metaml_results_binary <- function(output_dir, output_prefix){

  if(!endsWith(output_dir, "/")){

    output_dir <- paste0(output_dir, "/")

  }
  ##################################################################################
  # file 1, general summary

  general_summary_raw <- readLines(paste0(output_dir, output_prefix,".txt"))

  general_summary_ordered <- list("full_dataset" = list("n_samples" = str_split(general_summary_raw[min(grep("samples",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                                                        "n_features" = str_split(general_summary_raw[min(grep("features",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                                                        "n_runs" = str_split(general_summary_raw[min(grep("runs",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                                                        "n_crossValidations" = str_split(general_summary_raw[min(grep("cv_folds",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                                                        "Accuracy"  = str_split(general_summary_raw[min(grep("accuracy",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                                                        "f1" =  str_split(general_summary_raw[min(grep("f1",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                                                        "Precision" = str_split(general_summary_raw[min(grep("precision",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                                                        "Recall" = str_split(general_summary_raw[min(grep("recall",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                                                        "AUC" = str_split(general_summary_raw[min(grep("auc",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                                                        "confusion_matrix" = str_split(general_summary_raw[min(grep("confusion matrix",general_summary_raw)):(min(grep("confusion matrix",general_summary_raw))+1)],
                                                                                       "\t", simplify = T)[,-1] %>%
                                                          as.numeric() %>%
                                                          matrix(nrow = 2, ncol = 2, byrow = FALSE)
  ),
  "top 25 features" =  list("n_samples" = str_split(general_summary_raw[max(grep("samples",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                            "n_features" = str_split(general_summary_raw[max(grep("features",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                            "n_runs" = str_split(general_summary_raw[max(grep("runs",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                            "n_crossValidations" = str_split(general_summary_raw[max(grep("cv_folds",general_summary_raw))], "\t", simplify = T)[2] %>% as.double(),
                            "Accuracy"  = str_split(general_summary_raw[max(grep("accuracy",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                            "f1" =  str_split(general_summary_raw[max(grep("f1",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                            "Precision" = str_split(general_summary_raw[max(grep("precision",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                            "Recall" = str_split(general_summary_raw[max(grep("recall",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                            "AUC" = str_split(general_summary_raw[max(grep("auc",general_summary_raw))], "\t", simplify = T)[-1] %>% as.double(),
                            "confusion_matrix" = str_split(general_summary_raw[max(grep("confusion matrix",general_summary_raw)):(max(grep("confusion matrix",general_summary_raw))+1)],
                                                           "\t", simplify = T)[,-1] %>%
                              as.numeric() %>%
                              matrix(nrow = 2, ncol = 2, byrow = FALSE)
  )
  )



  # get the feature importance
  feature_importance_indx <- agrep("feature importance", general_summary_raw)
  feature_importance <- tibble(tmp = general_summary_raw[(feature_importance_indx+1):length(general_summary_raw)]) %>%
    separate(tmp, into = c("rank", "name", "avg_importance", "stdev_imporance"), sep = "\t")

  general_summary_ordered[["Feature Importance"]] <- feature_importance

  ##################################################################################
  ## file 2, ROC curve

  ROC <- readLines(paste0(output_dir, output_prefix,"_roccurve.txt"))

  split_ROCs_indx <- max(grep("#features",ROC, fixed = T))

  runs_and_folds_indx <- which(startsWith(ROC, "run/fold"))[which(which(startsWith(ROC, "run/fold")) < (split_ROCs_indx-1))]

  list_allFeatures <- list()

  for(i in runs_and_folds_indx){

    line <- ROC[i]

    run_fold <- str_split(ROC[i],"\t", simplify = TRUE)[2] %>%
      str_split("/", simplify = TRUE)

    data_frame_to_add <- cbind.data.frame(ROC[i+1] %>%
                                            str_split("\t", simplify = TRUE) %>%
                                            .[1:(length(.) - 1)],
                                          ROC[i+2] %>%
                                            str_split("\t", simplify = TRUE) %>%
                                            .[1:(length(.) - 1)],
                                          ROC[i+3] %>%
                                            str_split("\t", simplify = TRUE) %>%
                                            .[1:(length(.) - 1)]
    ) %>%
      set_names(c("var1", "var2", "var3")) %>%
      mutate_all(as.numeric)

    list_allFeatures[[run_fold[1]]][[run_fold[2]]] <- data_frame_to_add
  }

  df_allFeatures <- lapply(list_allFeatures, bind_rows, .id = "fold") %>%
    bind_rows(.id = "run") %>%
    as_tibble()

  general_summary_ordered$full_dataset$ROC <- df_allFeatures


  runs_and_folds_indx <- which(startsWith(ROC, "run/fold"))[which(which(startsWith(ROC, "run/fold")) > split_ROCs_indx)]

  list_top25Features <- list()

  for(i in runs_and_folds_indx){

    line <- ROC[i]

    run_fold <- str_split(ROC[i],"\t", simplify = TRUE)[2] %>%
      str_split("/", simplify = TRUE)

    data_frame_to_add <- cbind.data.frame(ROC[i+1] %>%
                                            str_split("\t", simplify = TRUE) %>%
                                            .[1:(length(.) - 1)],
                                          ROC[i+2] %>%
                                            str_split("\t", simplify = TRUE) %>%
                                            .[1:(length(.) - 1)],
                                          ROC[i+3] %>%
                                            str_split("\t", simplify = TRUE) %>%
                                            .[1:(length(.) - 1)]
    ) %>%
      set_names(c("var1", "var2", "var3")) %>%
      mutate_all(as.numeric)

    list_top25Features[[run_fold[1]]][[run_fold[2]]] <- data_frame_to_add
  }

  df_top25Features <- lapply(list_top25Features, bind_rows, .id = "fold") %>%
    bind_rows(.id = "run") %>%
    as_tibble()


  general_summary_ordered$`top 25 features`$ROC <- df_allFeatures

  ###################################################################################
  # file 3, estimations

  # part1, all features

  estimations_raw <- readLines(paste0(output_dir,
                                      output_prefix,
                                      "_estimations.txt"))

  features_split <- max(grep("#features", estimations_raw, fixed = TRUE))

  runs_and_folds_indx <- which(grepl("run/fold", estimations_raw))

  runs_and_folds_indx_all_feature <- runs_and_folds_indx[runs_and_folds_indx < features_split]

  list_allRuns <- list()

  suppressWarnings(for (i in runs_and_folds_indx_all_feature) {
    run_fold <-
      str_split(estimations_raw[i], "\t", simplify = TRUE)[2] %>%
      str_split("/", simplify = TRUE)

    run <- run_fold[1]
    fold <- run_fold[2]

    true_labs <-
      str_split(estimations_raw[i + 1], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]

    estim_labs <-
      str_split(estimations_raw[i + 2], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]

    estim_probs <-
      str_split(estimations_raw[i + 3], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]

    sample_indx <-
      str_split(estimations_raw[i + 4], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]


    list_allRuns[[run]][[fold]] <-
      cbind(true_labs, estim_labs, estim_probs, sample_indx)

  })

  general_summary_ordered$full_dataset$labels_estimations <- lapply(list_allRuns, lapply, as.data.frame) %>%
    lapply(bind_rows, .id = "fold") %>%
    bind_rows(.id = "run") %>%
    as_tibble()

  # part 2, top 25 features
  runs_and_folds_indx_25_feature <- runs_and_folds_indx[runs_and_folds_indx > features_split]

  list_allRuns <- list()

  suppressWarnings(for (i in runs_and_folds_indx_25_feature) {
    run_fold <-
      str_split(estimations_raw[i], "\t", simplify = TRUE)[2] %>%
      str_split("/", simplify = TRUE)

    run <- run_fold[1]
    fold <- run_fold[2]

    true_labs <-
      str_split(estimations_raw[i + 1], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]

    estim_labs <-
      str_split(estimations_raw[i + 2], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]

    estim_probs <-
      str_split(estimations_raw[i + 3], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]

    sample_indx <-
      str_split(estimations_raw[i + 4], "\t", simplify = TRUE) %>%
      parse_number() %>%
      .[!is.na(.)]


    list_allRuns[[run]][[fold]] <-
      cbind(true_labs, estim_labs, estim_probs, sample_indx)

  })

  general_summary_ordered$`top 25 features`$labels_estimation <- lapply(list_allRuns, lapply, as.data.frame) %>%
    lapply(bind_rows, .id = "fold") %>%
    bind_rows(.id = "run") %>%
    as_tibble()

  return(general_summary_ordered)

}
