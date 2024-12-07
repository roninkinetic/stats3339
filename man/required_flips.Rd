% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{required_flips}
\alias{required_flips}
\title{Required Sample Size for Coin-Flipping Experiment}
\usage{
required_flips(confidence_level, width)
}
\arguments{
\item{confidence_level}{The confidence level as a percentage (e.g., 95 for 95 percent confidence).}

\item{width}{The desired total width of the confidence interval (e.g., 0.05).}
}
\value{
The required number of flips and a detailed explanation.
}
\description{
Calculates the required number of flips to achieve a specified confidence interval width for a fair coin.
}
\examples{
required_flips(confidence_level = 95, width = 0.05)
}
