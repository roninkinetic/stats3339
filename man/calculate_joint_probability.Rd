% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{calculate_joint_probability}
\alias{calculate_joint_probability}
\title{Joint Probability P(X operator1 value1, Y operator2 value2)}
\usage{
calculate_joint_probability(
  f_xy,
  x_values,
  y_values,
  operator_x,
  value_x,
  operator_y,
  value_y
)
}
\arguments{
\item{f_xy}{A function representing the joint probability distribution \( f(x, y) \).}

\item{x_values}{A numeric vector of possible X values.}

\item{y_values}{A numeric vector of possible Y values.}

\item{operator_x}{A string specifying the condition for X ("<", "<=", "=", ">=", ">").}

\item{value_x}{The value to compare X against.}

\item{operator_y}{A string specifying the condition for Y ("<", "<=", "=", ">=", ">").}

\item{value_y}{The value to compare Y against.}
}
\value{
A list containing the probability and a detailed explanation of the steps.
}
\description{
Calculates the joint probability for specified conditions on \( X \) and \( Y \) using a given joint probability distribution.
}
\examples{
# Define the joint probability density function
f_xy <- function(x, y) (x + 2 * y) / 42
x_values <- c(0, 1, 2, 3)
y_values <- c(0, 1, 2)
# Calculate P(X <= 2, Y = 1)
calculate_joint_probability(f_xy, x_values, y_values, "<=", 2, "=", 1)
}
