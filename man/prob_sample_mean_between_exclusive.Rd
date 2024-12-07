% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{prob_sample_mean_between_exclusive}
\alias{prob_sample_mean_between_exclusive}
\title{Probability of Sample Mean Between Two Thresholds (Exclusive)}
\usage{
prob_sample_mean_between_exclusive(
  lower_bound,
  upper_bound,
  population_mean,
  population_sd,
  sample_size
)
}
\arguments{
\item{lower_bound}{Lower threshold for the sample mean.}

\item{upper_bound}{Upper threshold for the sample mean.}

\item{population_mean}{Population mean (mu).}

\item{population_sd}{Population standard deviation (sigma).}

\item{sample_size}{Size of the sample (n).}
}
\value{
The probability that the sample mean is between the two thresholds, exclusive.
}
\description{
Calculates the probability that the sample mean is strictly between two thresholds, exclusive.
}
\examples{
prob_sample_mean_between_exclusive(lower_bound = 45,
upper_bound = 50, population_mean = 47, population_sd = 10.2, sample_size = 25)
}
