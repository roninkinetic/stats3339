% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{probability_union}
\alias{probability_union}
\title{Probability of Union of Two Events}
\usage{
probability_union(p_A, p_B, p_A_and_B)
}
\arguments{
\item{p_A}{Probability of event \( A \) (e.g., coronary care).}

\item{p_B}{Probability of event \( B \) (e.g., oncology care).}

\item{p_A_and_B}{Probability of the intersection of events \( A \) and \( B \) (e.g., both coronary and oncology care).}
}
\value{
A list containing the probability of the union \( P(A union B) \) and a detailed explanation.
}
\description{
Calculates the probability of the union of two events \( A \) and \( B \), accounting for overlap.
}
\examples{
probability_union(0.14, 0.26, 0.06)
}
