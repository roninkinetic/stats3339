% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{required_sample_size}
\alias{required_sample_size}
\title{Required Sample Size for Proportions}
\usage{
required_sample_size(confidence_level, width, probability)
}
\arguments{
\item{confidence_level}{The confidence level as a percentage (e.g., 95 for 95 percent confidence).}

\item{width}{The desired total width of the confidence interval (e.g., 0.05).}

\item{probability}{The estimated probability of success (e.g., 0.5 for a fair coin, 1/6 for a six-sided die).}
}
\value{
The required sample size and a detailed explanation.
}
\description{
Calculates the required sample size to achieve a specified confidence interval width for any probability-based scenario.
}
\examples{
required_sample_size(confidence_level = 95, width = 0.05, probability = 0.5)
required_sample_size(confidence_level = 95, width = 0.05, probability = 1/6)
}
