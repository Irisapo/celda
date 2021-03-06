% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/celda_CG.R
\docType{methods}
\name{clusterProbability,celda_CG-method}
\alias{clusterProbability,celda_CG-method}
\title{Conditional probabilities for cells and features from a Celda_CG
 model}
\usage{
\S4method{clusterProbability}{celda_CG}(counts, celdaMod, log = FALSE,
  ...)
}
\arguments{
\item{counts}{Integer matrix. Rows represent features and columns represent
cells. This matrix should be the same as the one used to generate
`celdaMod`.}

\item{celdaMod}{Celda object of class `celda_CG`.}

\item{log}{Logical. If FALSE, then the normalized conditional probabilities
will be returned. If TRUE, then the unnormalized log probabilities will be
returned. Default FALSE.}

\item{...}{Additional parameters.}
}
\value{
A list containging a matrix for the conditional cell and feature
 cluster probabilities.
}
\description{
Calculates the conditional probability of each cell belonging
 to each subpopulation given all other cell cluster assignments as well as
 each feature belonging to each module given all other feature cluster
 assignments in a `celda_CG()` result.
}
\examples{
data(celdaCGSim, celdaCGMod)
clusterProb <- clusterProbability(celdaCGSim$counts, celdaCGMod)
}
\seealso{
`celda_CG()` for clustering features and cells
}
