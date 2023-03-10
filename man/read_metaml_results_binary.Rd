% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_metaml_results_binary.R
\name{read_metaml_results_binary}
\alias{read_metaml_results_binary}
\title{Read MetaML binary classification output}
\usage{
read_metaml_results_binary(output_dir, output_prefix)
}
\arguments{
\item{output_dir}{A \code{character} vector to locate the output directory (see the code you used to launch the MetaML algorithm)}

\item{output_prefix}{A \code{character} vector to know which name prefix you gave to your files}
}
\value{
A \code{list} containing all statistics divided into "full" and "top25" dataset
}
\description{
So far it works only with binary classification using Random Forest, and it assumes the output to have top 25 features too.
}
\examples{

##### say you run something like the following code:

# python metaml/classification_thomas-manghi.py
#    input_data/data_ready.txt # input file with all variables needed in rows
#   output/Exp1 # output prefix
#   -l rf # type of learning (my function was tested only on rf)
#   -z s__ # prefix of microbiome- related variables
#   -mf 0.1 # number of ...?
#   -nt 1000 # number of trees
#    -nsl 10 # number of samples per tree leaf
#   --define 1:varible_of_interest:category_of_interest # the rest becomes 0
#   -c entropy # way to estimate the variable importance
#   -p 10 # number of folds for cross-validation
#   -r 10 # number of total runs. NB: each run is independent, while folds within a run are not independent
#    -cc 25 # re-do the analysis with this number of top features, append results on the same files

##### when this is done, go back to R and do:

# my_ML_results <- read_metaml_results_binary("output/", "Exp1")

}
