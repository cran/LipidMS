% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internalProcessing.R
\name{clustering}
\alias{clustering}
\title{EIC extraction based on previous partitions generated by \link{partitioning}}
\usage{
clustering(msobject, dmzagglom, drtclust, minpeak, mslevel, cE)
}
\arguments{
\item{msobject}{msobject generated by \link{partitioning}}

\item{dmzagglom}{mz tolerance for clusters}

\item{drtclust}{RT window for clusters}

\item{minpeak}{minimum number of measures to define a peak}

\item{mslevel}{info to access msobject}

\item{cE}{info to access msobject}
}
\value{
msobject
}
\description{
EIC extraction based on previous partitions generated by \link{partitioning}
}
\references{
Peak-picking algorithm has been imported from enviPick R-package:
https://cran.r-project.org/web/packages/enviPick/index.html
}
\author{
M Isabel Alcoriza-Balaguer <maribel_alcoriza@iislafe.es>
}
\keyword{internal}
