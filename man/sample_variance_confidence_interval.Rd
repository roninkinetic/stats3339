% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{sample_variance_confidence_interval}
\alias{sample_variance_confidence_interval}
\title{Confidence Interval for Sample Variance}
\usage{
sample_variance_confidence_interval(sample_variance, n, confidence_level)
}
\arguments{
\item{sample_variance}{The variance of the sample.}

\item{n}{The sample size.}

\item{confidence_level}{The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).}
}
\value{
A list containing the lower and upper bounds of the confidence interval for the sample variance.
}
\description{
Calculates the confidence interval for the sample variance given the sample variance, sample size, and confidence level.
}
\examples{
sample_variance_confidence_interval(sample_variance = 25, n = 30, confidence_level = 0.95)
}
