% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adonis_to_table.R
\name{adonis_to_table}
\alias{adonis_to_table}
\title{Quick wrapper to format PERMANOVA results}
\usage{
adonis_to_table(adonis.out, adjR2 = TRUE, new_rownames = NA)
}
\arguments{
\item{adonis.out}{A \code{adonis2} output object as is}

\item{adjR2}{A \code{logical}, should another column be calculated with the adjusted R squared? the thing is calculated with \code{vegan::RsquareAdj}}

\item{new_rownames}{A \code{character} vector containing prettier variable
names, which will be passed with set_rownames}
}
\value{
A \code{data.frame} Ready to be passed to kable/DT/flextable
}
\description{
Quick wrapper to format PERMANOVA results
}
\examples{

library(gautils)
data(enterotype)
distMat <-  phyloseq::distance(microbiome::transform(enterotype, "compositional"), "bray")

adonisObj <- adonis2(distMat ~ SeqTech + Enterotype, by = "margin", permutations = 2000, data = microbiome::meta(enterotype), na.action = na.omit)

adonis_as_table <- adonis_to_table(adonisObj)
adonis.kbl <- kableExtra::kbl(adonis_as_table, digits = 2) \%>\%
kableExtra::kable_styling()

adonis_as_table2 <- adonis_to_table(adonisObj, new_rownames = c("Sequencing Techonology", "Enterotype Code"))
adonis.kbl2 <- kableExtra::kbl(adonis_as_table2, digits = 2) \%>\%
kableExtra::kable_styling()

}
