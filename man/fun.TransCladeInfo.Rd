% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TransCladeInfo.R
\name{fun.TransCladeInfo}
\alias{fun.TransCladeInfo}
\title{Transform clade infomation into a dataframe}
\usage{
fun.TransCladeInfo(CCChar)
}
\arguments{
\item{CCChar}{A character of clade infomation.}
}
\value{
A dataframe with the columns of clade size and frequency.
}
\description{
Transform clade infomation into a dataframe
}
\examples{
library(TarCA.beta)
tmp.cladesize <- "1 (12), 2 (1), 3 (2)"
fun.TransCladeInfo(tmp.cladesize)

}
