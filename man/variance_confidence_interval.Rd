% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{variance_confidence_interval}
\alias{variance_confidence_interval}
\title{Confidence Interval for Population Variance}
\usage{
variance_confidence_interval(sd, n, confidence_level)
}
\arguments{
\item{sd}{The population standard deviation.}

\item{n}{The sample size.}

\item{confidence_level}{The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).}
}
\value{
A list containing the lower and upper bounds of the confidence interval for the variance.
}
\description{
Calculates the confidence interval for the population variance given the sample variance, sample size, and confidence level.
}
\examples{
variance_confidence_interval(sd = 5, n = 30, confidence_level = 0.95)
}
