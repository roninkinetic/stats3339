% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{hypothesis_test_proportion}
\alias{hypothesis_test_proportion}
\title{Hypothesis Test for Population Proportion (with Exact Binomial Fallback)}
\usage{
hypothesis_test_proportion(n, x, p0, alpha = 0.05, alternative = "two.sided")
}
\arguments{
\item{n}{Sample size, the number of observations in the sample.}

\item{x}{The number of successes or cases meeting the condition of interest.}

\item{p0}{The hypothesized population proportion (e.g., 0.10 for 10 percent).}

\item{alpha}{The significance level for the test (default is 0.05).}

\item{alternative}{A character string specifying the alternative hypothesis.
Use "less", "greater", or "two.sided".}
}
\value{
A list containing the test statistic (if applicable), p-value, decision, and conclusion, with all work shown.
}
\description{
Performs a hypothesis test for a population proportion using a z-test or exact binomial test when appropriate.
}
\examples{
hypothesis_test_proportion(n = 100, x = 91, p0 = 0.91, alpha = 0.01, alternative = "two.sided")
}
