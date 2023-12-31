% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sumformula_from_CAS.R
\name{sumformula_from_CAS}
\alias{sumformula_from_CAS}
\title{sumformula_from_CAS.}
\usage{
sumformula_from_CAS(x = NULL)
}
\arguments{
\item{x}{Vector of CAS IDs.}
}
\value{
A character vector of length input vector.
}
\description{
\code{sumformula_from_CAS} will try to get a chemical sum formula from a CAS ID.
}
\details{
tbd.
}
\examples{
\dontrun{
x <- readLines("C:/Users/jlisec/Documents/Francesco Russo/RECTOX/RECTOX_GC-EI-MS_CASRN")
sf <- sumformula_from_CAS(x = x)
}
}
