% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{generic_hypothesis_test}
\alias{generic_hypothesis_test}
\title{Generic Hypothesis Test with Rejection Region}
\usage{
generic_hypothesis_test(
  n,
  sample_stat,
  hypothesized_value = 0,
  standard_error,
  alpha = 0.05,
  test_type = "two-tailed"
)
}
\arguments{
\item{n}{Sample size, the number of observations.}

\item{sample_stat}{The sample statistic (e.g., mean or proportion).}

\item{hypothesized_value}{The hypothesized population parameter (default is 0).}

\item{standard_error}{The standard error of the sample statistic (calculated outside the function if needed).}

\item{alpha}{The significance level for the test (default is 0.05).}

\item{test_type}{The type of test ("two-tailed", "greater", "less").}
}
\value{
A list containing the test type, rejection region, test statistic, p-value, decision, and conclusion, with all work shown.
}
\description{
Performs a hypothesis test for a population mean or other parameter with a specified rejection region.
}
\details{
This function calculates the test statistic, rejection region, and p-value based on the test type.
It provides step-by-step output, including all relevant calculations and decision criteria.
}
\examples{
# Example: Two-tailed test for a mean
generic_hypothesis_test(n = 49, sample_stat = 415.7, hypothesized_value = 420,
                        standard_error = 10 / sqrt(49), alpha = 0.01, test_type = "two-tailed")
}
