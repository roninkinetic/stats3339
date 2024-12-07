% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{hypergeo_exactly}
\alias{hypergeo_exactly}
\title{Hypergeometric Probability: Exactly}
\usage{
hypergeo_exactly(N, K, n, k)
}
\arguments{
\item{N}{The total population size.}

\item{K}{The total number of successes in the population.}

\item{n}{The number of draws.}

\item{k}{The number of successes in the sample.}
}
\value{
The probability of exactly \( k \) successes.
}
\description{
Calculates the probability of selecting exactly \( k \) successes in a hypergeometric distribution.
}
\examples{
hypergeo_exactly(N = 13, K = 6, n = 6, k = 2)
}
