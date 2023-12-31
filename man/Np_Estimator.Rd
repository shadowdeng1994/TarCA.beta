% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Np_Estimator.R
\name{Np_Estimator}
\alias{Np_Estimator}
\title{Estimate Np using cell phylogeny for each sub-population.}
\usage{
Np_Estimator(Tree, Ann, ExTree = NULL, ReturnExTree = FALSE)
}
\arguments{
\item{Tree}{A tree file of class "phylo" with node labels (using node ID by default).}

\item{Ann}{A dataframe with TipLabel and TipAnn (tip labels on the tree file and corresponding cell annotations).}

\item{ExTree}{(optional) A ExTree file obtained from conv_ExTree.}

\item{ReturnExTree}{Output in ExTree format (default FALSE).}
}
\value{
A dataframe storing the inferred Np for each cell type.
}
\description{
Estimate Np using cell phylogeny for each sub-population.
}
\examples{
library("TarCA.beta")
load(system.file("Exemplar","Exemplar_TCA.RData",package = "TarCA"))
tmp.res <- Np_Estimator(
Tree = ExemplarData_1$Tree,
Ann = ExemplarData_1$Ann
)
}
