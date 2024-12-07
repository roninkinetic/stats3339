% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{poisson_greater}
\alias{poisson_greater}
\title{Poisson Probability: Greater Than}
\usage{
poisson_greater(lambda, k)
}
\arguments{
\item{lambda}{The mean number of occurrences in the given time period.}

\item{k}{The threshold number of events.}
}
\value{
The probability of more than (k) events occurring.
}
\description{
Calculates the probability of more than (k) events occurring in a Poisson distribution.
}
\examples{
poisson_greater(lambda = 3, k = 3)
}
