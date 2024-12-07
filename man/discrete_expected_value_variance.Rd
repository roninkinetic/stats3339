% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{discrete_expected_value_variance}
\alias{discrete_expected_value_variance}
\title{Expected Value and Variance of Discrete Probability Distribution}
\usage{
discrete_expected_value_variance(outcomes, probabilities)
}
\arguments{
\item{outcomes}{A numeric vector representing the outcomes of the random variable.}

\item{probabilities}{A numeric vector representing the probabilities of each outcome.}
}
\value{
A list containing the expected value, variance, and an explanation of the calculations.
}
\description{
Calculates the expected value and variance for a given discrete probability distribution.
}
\examples{
outcomes <- c(1, 2, 3, 4, 5, 6)
probabilities <- c(0.1, 0.2, 0.3, 0.3, 0, 0.1)
discrete_expected_value_variance(outcomes, probabilities)
}
