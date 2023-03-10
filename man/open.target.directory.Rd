% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/open.target.directory.R
\name{open.target.directory}
\alias{open.target.directory}
\title{Open file/plot directory}
\usage{
\method{open}{target.directory}(element_path)
}
\arguments{
\item{element_path}{The path, written unix-style. it has to be absolute, either from the root, or from the std. home ("~/")}
}
\value{
A window opened on your file browsing element ("Finder" or "Explorer")
}
\description{
Given a PATH, either ending with "/" or a file name, this function opens a window to show you that directory.
}
\examples{
# in Unix (Linux or Mac)
open.target.directory("/home/user/Documents/my_superImportantDoc.txt")
# in Windows
open.target.directory("C:/Users/user/Documents/my_superImportantDoc.txt")


}
