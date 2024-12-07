% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{poisson_exact_3var}
\alias{poisson_exact_3var}
\title{Poisson Probability for Exact Events}
\usage{
poisson_exact_3var(rate, duration, k)
}
\arguments{
\item{rate}{The average rate of events per time unit (e.g., complaints per hour).}

\item{duration}{The time period during which events are counted (e.g., 8 hours).}

\item{k}{The exact number of events to calculate the probability for.}
}
\value{
The probability of observing exactly `k` events.
}
\description{
Calculates the probability of observing exactly `k` events in a Poisson distribution.
}
\examples{
poisson_exact_3var(rate = 0.6, duration = 8, k = 4)
}
