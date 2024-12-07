% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{cumulative_probability_greater}
\alias{cumulative_probability_greater}
\title{cumulative Probability Greater Than}
\usage{
cumulative_probability_greater(x, p, threshold)
}
\arguments{
\item{x}{A numeric vector of possible values of the random variable.}

\item{p}{A numeric vector of probabilities corresponding to the values in `x`.}

\item{threshold}{A numeric value specifying the threshold.}
}
\value{
A list containing the calculated probability and a detailed explanation.
}
\description{
Calculates the probability that a discrete random variable is greater than a given threshold.
}
\examples{
cumulative_probability_greater(x = c(0, 1, 2, 3), p = c(4, 5, 8, 3), threshold = 2)
}
