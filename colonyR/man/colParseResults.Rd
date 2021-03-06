% Generated by roxygen2 (4.0.1): do not edit by hand
\name{colParseResults}
\alias{colParseResults}
\title{Parse results from a Colony run}
\usage{
colParseResults(file, dir = "")
}
\arguments{
\item{file}{basename of the output file name given to Colony}

\item{dir}{directory where the run took place}
}
\value{
A list with two elements: \code{run} with the monitoring
  information and \code{pedigree} with the final pedigree
}
\description{
Parse the monitoring information and final best configuration
  for a Colony run
}

