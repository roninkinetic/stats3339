% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{ci_mean}
\alias{ci_mean}
\title{Confidence Interval for a Mean}
\usage{
ci_mean(data, confidence_level = 0.95)
}
\arguments{
\item{data}{A numeric vector of sample data.}

\item{confidence_level}{The confidence level (e.g., 0.95 for 95 percent confidence).}
}
\value{
A list with intermediate calculations and the confidence interval:
  \item{mean}{The sample mean.}
  \item{standard_error}{The standard error of the mean.}
  \item{t_value}{The critical t-value.}
  \item{margin_of_error}{The margin of error.}
  \item{lower_bound}{The lower bound of the confidence interval.}
  \item{upper_bound}{The upper bound of the confidence interval.}
}
\description{
Calculates the confidence interval for a sample mean and shows the intermediate steps.
}
\examples{
ci_mean(c(5, 10, 15, 20, 25), 0.95)

}
