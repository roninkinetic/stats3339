utils::globalVariables(c("x"))

#' @title Poisson Probability: Exactly
#' @description Calculates the probability of exactly (k) events occurring in a Poisson distribution.
#' @param lambda The mean number of occurrences in the given time period.
#' @param k The number of events for which the probability is calculated.
#' @importFrom stats dpois ppois
#' @return The probability of exactly (k) events occurring.
#' @examples
#' poisson_exact(lambda = 3, k = 3)
#' @export
poisson_exact <- function(lambda, k) {
  result <- dpois(k, lambda)
  explanation <- paste("The probability of exactly", k, "events when lambda =", lambda, "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Poisson Probability: Greater Than
#' @description Calculates the probability of more than (k) events occurring in a Poisson distribution.
#' @param lambda The mean number of occurrences in the given time period.
#' @param k The threshold number of events.
#' @importFrom stats dpois ppois
#' @return The probability of more than (k) events occurring.
#' @examples
#' poisson_greater(lambda = 3, k = 3)
#' @export
poisson_greater <- function(lambda, k) {
  result <- 1 - ppois(k, lambda)
  explanation <- paste("The probability of more than", k, "events when lambda =", lambda, "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Poisson Probability: Less Than
#' @description Calculates the probability of fewer than (k) events occurring in a Poisson distribution.
#' @param lambda The mean number of occurrences in the given time period.
#' @param k The threshold number of events (exclusive).
#' @importFrom stats dpois ppois
#' @return The probability of fewer than (k) events occurring.
#' @examples
#' poisson_less(lambda = 3, k = 3)
#' @export
poisson_less <- function(lambda, k) {
  result <- ppois(k - 1, lambda)
  explanation <- paste("The probability of fewer than", k, "events when lambda =", lambda, "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Poisson Probability: Between Two Values
#' @description Calculates the probability of a Poisson random variable being between two values (inclusive).
#' @param lambda The mean number of occurrences in the given time period.
#' @param k1 The lower bound of events (inclusive).
#' @param k2 The upper bound of events (inclusive).
#' @importFrom stats dpois ppois
#' @return The probability of the random variable being between (k1) and (k2).
#' @examples
#' poisson_between(lambda = 3, k1 = 2, k2 = 4)
#' @export
poisson_between <- function(lambda, k1, k2) {
  result <- ppois(k2, lambda) - ppois(k1 - 1, lambda)
  explanation <- paste("The probability of events being between", k1, "and", k2,
                       "inclusive, when lambda =", lambda, "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}


#' @title Discrete Probability: Exactly
#' @description Calculates the probability of exactly \( x \) in a discrete probability distribution.
#' @param prob_dist A named numeric vector where names represent the values of \( X \) and values represent \( P(X) \).
#' @param x The value for which the probability is calculated.
#' @return The probability of \( X = x \).
#' @examples
#' prob_dist <- c("0" = 0.1, "1" = 0.2, "2" = 0.3, "3" = 0.4)
#' discrete_exactly(prob_dist, x = 2)
#' @export
discrete_exactly <- function(prob_dist, x) {
  if (as.character(x) %in% names(prob_dist)) {
    result <- prob_dist[as.character(x)]
    explanation <- paste("The probability of exactly X =", x, "is", result)
    return(list(probability = result, explanation = explanation))
  } else {
    return(list(probability = 0, explanation = paste("X =", x, "is not in the distribution.")))
  }
}

#' @title Discrete Probability: Greater Than
#' @description Calculates the probability of \( X > x \) in a discrete probability distribution.
#' @param prob_dist A named numeric vector where names represent the values of \( X \) and values represent \( P(X) \).
#' @param x The threshold value for \( X > x \).
#' @return The probability of \( X > x \).
#' @examples
#' prob_dist <- c("0" = 0.1, "1" = 0.2, "2" = 0.3, "3" = 0.4)
#' discrete_greater_than(prob_dist, x = 1)
#' @export
discrete_greater_than <- function(prob_dist, x) {
  filtered_prob <- prob_dist[as.numeric(names(prob_dist)) > x]
  result <- sum(filtered_prob)
  explanation <- paste("The probability of X >", x, "is", result)
  return(list(probability = result, explanation = explanation))
}

#' @title Discrete Probability: Less Than
#' @description Calculates the probability of \( X < x \) in a discrete probability distribution.
#' @param prob_dist A named numeric vector where names represent the values of \( X \) and values represent \( P(X) \).
#' @param x The threshold value for \( X < x \).
#' @return The probability of \( X < x \).
#' @examples
#' prob_dist <- c("0" = 0.1, "1" = 0.2, "2" = 0.3, "3" = 0.4)
#' discrete_less_than(prob_dist, x = 2)
#' @export
discrete_less_than <- function(prob_dist, x) {
  filtered_prob <- prob_dist[as.numeric(names(prob_dist)) < x]
  result <- sum(filtered_prob)
  explanation <- paste("The probability of X <", x, "is", result)
  return(list(probability = result, explanation = explanation))
}

#' @title Discrete Probability: Between Two Values
#' @description Calculates the probability of x1 <= X <= x2 in a discrete probability distribution.
#' @param prob_dist A named numeric vector where names represent the values of \( X \) and values represent \( P(X) \).
#' @param x1 The lower bound for \( X \) (inclusive).
#' @param x2 The upper bound for \( X \) (inclusive).
#' @return The probability of x1 <= X <= x2.
#' @examples
#' prob_dist <- c("0" = 0.1, "1" = 0.2, "2" = 0.3, "3" = 0.4)
#' discrete_between(prob_dist, x1 = 1, x2 = 3)
#' @export
discrete_between <- function(prob_dist, x1, x2) {
  filtered_prob <- prob_dist[as.numeric(names(prob_dist)) >= x1 & as.numeric(names(prob_dist)) <= x2]
  result <- sum(filtered_prob)
  explanation <- paste("The probability of", x1, "<= X <=", x2, "is", result)
  return(list(probability = result, explanation = explanation))
}


#' @title Hypergeometric Probability: Exactly
#' @description Calculates the probability of selecting exactly \( k \) successes in a hypergeometric distribution.
#' @param N The total population size.
#' @param K The total number of successes in the population.
#' @param n The number of draws.
#' @param k The number of successes in the sample.
#' @return The probability of exactly \( k \) successes.
#' @examples
#' hypergeo_exactly(N = 13, K = 6, n = 6, k = 2)
#' @export
hypergeo_exactly <- function(N, K, n, k) {
  result <- dhyper(k, K, N - K, n)
  explanation <- paste("The probability of selecting exactly", k,
                       "successes from a population of", N,
                       "with", K, "successes and", n, "draws is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Hypergeometric Probability: Greater Than
#' @description Calculates the probability of selecting more than \( k \) successes in a hypergeometric distribution.
#' @param N The total population size.
#' @param K The total number of successes in the population.
#' @param n The number of draws.
#' @param k The threshold number of successes.
#' @return The probability of selecting more than \( k \) successes.
#' @examples
#' hypergeo_greater_than(N = 13, K = 6, n = 6, k = 2)
#' @export
hypergeo_greater_than <- function(N, K, n, k) {
  result <- sum(dhyper((k + 1):min(n, K), K, N - K, n))
  explanation <- paste("The probability of selecting more than", k,
                       "successes from a population of", N,
                       "with", K, "successes and", n, "draws is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Hypergeometric Probability: Less Than
#' @description Calculates the probability of selecting fewer than \( k \) successes in a hypergeometric distribution.
#' @param N The total population size.
#' @param K The total number of successes in the population.
#' @param n The number of draws.
#' @param k The threshold number of successes (exclusive).
#' @return The probability of selecting fewer than \( k \) successes.
#' @examples
#' hypergeo_less_than(N = 13, K = 6, n = 6, k = 3)
#' @export
hypergeo_less_than <- function(N, K, n, k) {
  result <- sum(dhyper(0:(k - 1), K, N - K, n))
  explanation <- paste("The probability of selecting fewer than", k,
                       "successes from a population of", N,
                       "with", K, "successes and", n, "draws is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Hypergeometric Probability: Between Two Values
#' @description Calculates the probability of selecting between \( k1 \) and \( k2 \) successes (inclusive) in a hypergeometric distribution.
#' @param N The total population size.
#' @param K The total number of successes in the population.
#' @param n The number of draws.
#' @param k1 The lower bound of successes (inclusive).
#' @param k2 The upper bound of successes (inclusive).
#' @return The probability of selecting between \( k1 \) and \( k2 \) successes.
#' @examples
#' hypergeo_between(N = 13, K = 6, n = 6, k1 = 1, k2 = 3)
#' @export
hypergeo_between <- function(N, K, n, k1, k2) {
  result <- sum(dhyper(k1:k2, K, N - K, n))
  explanation <- paste("The probability of selecting between", k1, "and", k2,
                       "successes inclusive, from a population of", N,
                       "with", K, "successes and", n, "draws is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}


#' @title Binomial Probability: Exactly
#' @description Calculates the probability of exactly \( k \) successes in a binomial distribution.
#' @param n The number of trials.
#' @param p The probability of success in a single trial.
#' @param k The number of successes for which the probability is calculated.
#' @importFrom stats dbinom dhyper pbinom
#' @return The probability of exactly \( k \) successes.
#' @examples
#' binomial_exactly(n = 15, p = 0.7, k = 10)
#' @export
binomial_exactly <- function(n, p, k) {
  result <- dbinom(k, n, p)
  explanation <- paste("The probability of exactly", k,
                       "successes in", n, "trials with success probability", p,
                       "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Binomial Probability: Greater Than
#' @description Calculates the probability of more than \( k \) successes in a binomial distribution.
#' @param n The number of trials.
#' @param p The probability of success in a single trial.
#' @param k The threshold number of successes.
#' @importFrom stats dbinom dhyper pbinom
#' @return The probability of more than \( k \) successes.
#' @examples
#' binomial_greater_than(n = 15, p = 0.7, k = 8)
#' @export
binomial_greater_than <- function(n, p, k) {
  result <- 1 - pbinom(k, n, p)
  explanation <- paste("The probability of more than", k,
                       "successes in", n, "trials with success probability", p,
                       "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Binomial Probability: Less Than
#' @description Calculates the probability of fewer than \( k \) successes in a binomial distribution.
#' @param n The number of trials.
#' @param p The probability of success in a single trial.
#' @param k The threshold number of successes (exclusive)
#' @importFrom stats dbinom dhyper pbinom
#' @return The probability of fewer than \( k \) successes.
#' @examples
#' binomial_less_than(n = 15, p = 0.7, k = 10)
#' @export
binomial_less_than <- function(n, p, k) {
  result <- pbinom(k - 1, n, p)
  explanation <- paste("The probability of fewer than", k,
                       "successes in", n, "trials with success probability", p,
                       "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Binomial Probability: Between Two Values
#' @description Calculates the probability of k1 <= X < k2 successes in a binomial distribution.
#' @param n The number of trials.
#' @param p The probability of success in a single trial.
#' @param k1 The lower bound for successes (inclusive).
#' @param k2 The upper bound for successes (exclusive).
#' @importFrom stats dbinom dhyper pbinom
#' @return The probability of successes between \( k1 \) and \( k2 \) (inclusive for \( k1 \), exclusive for \( k2 \)).
#' @examples
#' binomial_between(n = 15, p = 0.7, k1 = 8, k2 = 11)
#' @export
binomial_between <- function(n, p, k1, k2) {
  result <- sum(dbinom(k1:(k2 - 1), n, p))
  explanation <- paste("The probability of successes between", k1, "and", k2 - 1,
                       "inclusive in", n, "trials with success probability", p,
                       "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Adjusted Mean and Standard Deviation
#' @description Calculates the adjusted mean and standard deviation of a dataset after applying a linear transformation.
#' @param mean_original The original mean of the dataset.
#' @param sd_original The original standard deviation of the dataset.
#' @param a The adjustment to be added to each value.
#' @param b The scaling factor to be applied after the adjustment.
#' @return A list containing the adjusted mean and standard deviation.
#' @examples
#' adjusted_stats(mean_original = 5.5, sd_original = 1.8, a = 0.2, b = 1.4)
#' @export
adjusted_stats <- function(mean_original, sd_original, a, b) {
  # Adjusted mean
  mean_adjusted <- (mean_original + a) * b

  # Adjusted standard deviation
  sd_adjusted <- sd_original * b

  # Explanation
  explanation_mean <- paste("The adjusted mean is calculated as (mean_original + a) * b =",
                            "(", mean_original, "+", a, ")", "*", b, "=", round(mean_adjusted, 4))
  explanation_sd <- paste("The adjusted standard deviation is calculated as sd_original * b =",
                          sd_original, "*", b, "=", round(sd_adjusted, 4))

  # Return results
  return(list(
    mean_adjusted = mean_adjusted,
    sd_adjusted = sd_adjusted,
    explanation_mean = explanation_mean,
    explanation_sd = explanation_sd
  ))
}

#' @title Joint Probability: Exactly
#' @description Calculates the probability that the sum of two independent random variables equals \( k \).
#' @param dist1 A named numeric vector representing the probability distribution of the first variable.
#' @param dist2 A named numeric vector representing the probability distribution of the second variable.
#' @param k The target sum of the two random variables.
#' @return The probability that the sum of the two variables equals \( k \).
#' @examples
#' dist1 <- c("0" = 0.5, "1" = 0.1, "2" = 0.25, "3" = 0.15)
#' dist2 <- c("0" = 0.25, "1" = 0.45, "2" = 0.1, "3" = 0.2)
#' joint_exactly(dist1, dist2, k = 0)
#' @export
joint_exactly <- function(dist1, dist2, k) {
  sum_prob <- 0
  for (i in names(dist1)) {
    j <- k - as.numeric(i)
    if (as.character(j) %in% names(dist2)) {
      sum_prob <- sum_prob + dist1[i] * dist2[as.character(j)]
    }
  }
  explanation <- paste("The probability that the sum of the two variables equals", k, "is", round(sum_prob, 4))
  return(list(probability = sum_prob, explanation = explanation))
}

#' @title Joint Probability: Greater Than
#' @description Calculates the probability that the sum of two independent random variables is greater than \( k \).
#' @param dist1 A named numeric vector representing the probability distribution of the first variable.
#' @param dist2 A named numeric vector representing the probability distribution of the second variable.
#' @param k The threshold for the sum of the two random variables.
#' @return The probability that the sum of the two variables is greater than \( k \).
#' @examples
#' dist1 <- c("0" = 0.5, "1" = 0.1, "2" = 0.25, "3" = 0.15)
#' dist2 <- c("0" = 0.25, "1" = 0.45, "2" = 0.1, "3" = 0.2)
#' joint_greater_than(dist1, dist2, k = 2)
#' @export
joint_greater_than <- function(dist1, dist2, k) {
  max_value <- as.numeric(names(dist1)) + as.numeric(names(dist2))
  result <- sum(sapply(max_value, function(x) if (x > k) joint_exactly(dist1, dist2, x)$probability else 0))
  explanation <- paste("The probability that the sum of the two variables is greater than", k, "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Joint Probability: Less Than
#' @description Calculates the probability that the sum of two independent random variables is less than \( k \).
#' @param dist1 A named numeric vector representing the probability distribution of the first variable.
#' @param dist2 A named numeric vector representing the probability distribution of the second variable.
#' @param k The threshold for the sum of the two random variables.
#' @return The probability that the sum of the two variables is less than \( k \).
#' @examples
#' dist1 <- c("0" = 0.5, "1" = 0.1, "2" = 0.25, "3" = 0.15)
#' dist2 <- c("0" = 0.25, "1" = 0.45, "2" = 0.1, "3" = 0.2)
#' joint_less_than(dist1, dist2, k = 2)
#' @export
joint_less_than <- function(dist1, dist2, k) {
  result <- sum(sapply(0:(k - 1), function(x) joint_exactly(dist1, dist2, x)$probability))
  explanation <- paste("The probability that the sum of the two variables is less than", k, "is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Joint Probability: Between Two Values
#' @description Calculates the probability that the sum of two independent random variables is between \( k1 \) and \( k2 \).
#' @param dist1 A named numeric vector representing the probability distribution of the first variable.
#' @param dist2 A named numeric vector representing the probability distribution of the second variable.
#' @param k1 The lower bound for the sum (inclusive).
#' @param k2 The upper bound for the sum (inclusive).
#' @return The probability that the sum of the two variables is between \( k1 \) and \( k2 \).
#' @examples
#' dist1 <- c("0" = 0.5, "1" = 0.1, "2" = 0.25, "3" = 0.15)
#' dist2 <- c("0" = 0.25, "1" = 0.45, "2" = 0.1, "3" = 0.2)
#' joint_between(dist1, dist2, k1 = 2, k2 = 4)
#' @export
joint_between <- function(dist1, dist2, k1, k2) {
  result <- sum(sapply(k1:k2, function(x) joint_exactly(dist1, dist2, x)$probability))
  explanation <- paste("The probability that the sum of the two variables is between", k1, "and", k2, "inclusive is", round(result, 4))
  return(list(probability = result, explanation = explanation))
}

#' @title Median and Shape of Distribution
#' @description Calculates the median and describes the shape of a distribution from a dataset.
#' @param data A numeric vector representing the dataset.
#' @return A list containing the median and a description of the shape of the distribution.
#' @examples
#' data <- c(2.3, 2.4, 2.6, 2.7, 2.7, 2.8, 2.8, 2.8, 2.9, 3.2, 3.2,
#'           3.3, 3.4, 3.6, 3.7, 3.8, 3.9, 4.1, 4.2, 4.2, 4.3, 4.4,
#'           5.3, 5.5, 5.5, 5.6)
#' median_and_shape(data)
#' @export
median_and_shape <- function(data) {
  # Sort data
  sorted_data <- sort(data)

  # Median
  n <- length(sorted_data)
  if (n %% 2 == 0) {
    median <- (sorted_data[n / 2] + sorted_data[n / 2 + 1]) / 2
  } else {
    median <- sorted_data[ceiling(n / 2)]
  }

  # Shape description
  range_diff <- max(sorted_data) - min(sorted_data)
  skewness <- ifelse(mean(data) > median, "right-skewed", "left-skewed")
  shape <- paste("The data is", skewness, "with a range of", range_diff, ".")

  return(list(median = median, shape = shape))
}



#' @title Complete Probability Distribution Analysis
#' @description Solves for unknown probabilities and calculates mean, variance, and transformations for a discrete random variable.
#' @param x_values A numeric vector of the possible values of the random variable \( X \).
#' @param probabilities A numeric vector of probabilities corresponding to \( x_values \), with unknowns as `NA`.
#' @param given_mean The given mean of the random variable \( X \).
#' @return A list containing answers to all parts with detailed explanations and formulas.
#' @examples
#' x_values <- c(-5, 1, 2, 4, 7)
#' probabilities <- c(NA, 0.1, 0.3, NA, 0.3)
#' given_mean <- 2.2
#' analyze_distribution(x_values, probabilities, given_mean)
#' @export
analyze_distribution <- function(x_values, probabilities, given_mean) {
  # Part a: Solve for m and n
  total_prob <- sum(probabilities, na.rm = TRUE)
  remaining_prob <- 1 - total_prob
  indices_na <- which(is.na(probabilities))

  if (length(indices_na) != 2) {
    stop("Exactly two probabilities must be unknown (NA).")
  }

  # Calculate m and n based on given mean
  m_index <- indices_na[1]
  n_index <- indices_na[2]
  m_coeff <- x_values[m_index]
  n_coeff <- x_values[n_index]

  sum_known <- sum(x_values * probabilities, na.rm = TRUE)
  m <- (remaining_prob * n_coeff - given_mean + sum_known) / (n_coeff - m_coeff)
  n <- remaining_prob - m

  probabilities[m_index] <- m
  probabilities[n_index] <- n

  # Part b: Calculate P(X < 4)
  p_x_less_4 <- sum(probabilities[x_values < 4])

  # Part c: Variance of X
  mean_x <- given_mean
  variance_x <- sum(probabilities * (x_values - mean_x)^2)

  # Part d: Expected value of Y = 4X + 1
  mean_y <- 4 * mean_x + 1

  # Part e: Variance and standard deviation of Y = 4X + 1
  variance_y <- 4^2 * variance_x
  sd_y <- sqrt(variance_y)

  # Explanations with formulas
  explanations <- list(
    part_a = paste(
      "Part a:",
      "Using the total probability rule,",
      "P(X = x) = 1:",
      "m + n + sum of known probabilities =", remaining_prob,
      "\nUsing the mean formula, E(X) = sum(x * P(X)):",
      "m *", m_coeff, "+ n *", n_coeff, "+ sum of known contributions =", given_mean,
      "\nm =", round(m, 4), ", n =", round(n, 4)
    ),
    part_b = paste(
      "Part b:",
      "P(X < 4) = sum(P(X = x) for X < 4):",
      "P(X < 4) =", round(p_x_less_4, 4)
    ),
    part_c = paste(
      "Part c:",
      "Variance formula, Var(X) = sum(P(X = x) * (x - E(X))^2):",
      "Var(X) =", round(variance_x, 4)
    ),
    part_d = paste(
      "Part d:",
      "Expected value formula for Y = 4X + 1:",
      "E(Y) = 4 * E(X) + 1 = 4 *", mean_x, "+ 1 =", round(mean_y, 4)
    ),
    part_e = paste(
      "Part e:",
      "Variance formula for Y = 4X + 1:",
      "Var(Y) = 4^2 * Var(X) = 16 *", round(variance_x, 4), "=", round(variance_y, 4),
      "\nStandard deviation, SD(Y) = sqrt(Var(Y)) = sqrt(", round(variance_y, 4), ") =", round(sd_y, 4)
    )
  )

  # Return results
  return(list(
    probabilities = probabilities,
    explanations = explanations
  ))
}

#' Confidence Interval for a Proportion
#'
#' Calculates the confidence interval for a population proportion and shows the intermediate steps.
#'
#' @param successes The number of successes in the sample (integer).
#' @param sample_size The total size of the sample (integer).
#' @param confidence_level The confidence level (e.g., 0.95 for 95 percent confidence).
#' @importFrom stats qnorm
#' @return A list with intermediate calculations and the confidence interval:
#'   \item{proportion}{The sample proportion (p_hat).}
#'   \item{standard_error}{The standard error of the proportion.}
#'   \item{z_value}{The critical z-value.}
#'   \item{margin_of_error}{The margin of error.}
#'   \item{lower_bound}{The lower bound of the confidence interval.}
#'   \item{upper_bound}{The upper bound of the confidence interval.}
#'
#' @examples
#' ci_proportion(50, 100, 0.95)
#'
#' @export
ci_proportion <- function(successes, sample_size, confidence_level = 0.95) {
  # Calculate sample proportion
  p_hat <- successes / sample_size

  # Calculate standard error
  standard_error <- sqrt((p_hat * (1 - p_hat)) / sample_size)

  # Calculate z-value
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)

  # Calculate margin of error
  margin_of_error <- z_value * standard_error

  # Calculate confidence interval
  lower_bound <- p_hat - margin_of_error
  upper_bound <- p_hat + margin_of_error

  # Return detailed results
  list(
    proportion = p_hat,
    standard_error = standard_error,
    z_value = z_value,
    margin_of_error = margin_of_error,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  )
}

#' Confidence Interval for a Mean
#'
#' Calculates the confidence interval for a sample mean and shows the intermediate steps.
#'
#' @param data A numeric vector of sample data.
#' @param confidence_level The confidence level (e.g., 0.95 for 95 percent confidence).
#'
#' @return A list with intermediate calculations and the confidence interval:
#'   \item{mean}{The sample mean.}
#'   \item{standard_error}{The standard error of the mean.}
#'   \item{t_value}{The critical t-value.}
#'   \item{margin_of_error}{The margin of error.}
#'   \item{lower_bound}{The lower bound of the confidence interval.}
#'   \item{upper_bound}{The upper bound of the confidence interval.}
#' @importFrom stats qt sd
#'
#' @examples
#' ci_mean(c(5, 10, 15, 20, 25), 0.95)
#'
#' @export
ci_mean <- function(data, confidence_level = 0.95) {
  n <- length(data)
  mean_value <- mean(data)

  # Calculate standard error
  standard_error <- sd(data) / sqrt(n)

  # Calculate t-value
  alpha <- 1 - confidence_level
  t_value <- qt(1 - alpha / 2, df = n - 1)

  # Calculate margin of error
  margin_of_error <- t_value * standard_error

  # Calculate confidence interval
  lower_bound <- mean_value - margin_of_error
  upper_bound <- mean_value + margin_of_error

  # Return detailed results
  list(
    mean = mean_value,
    standard_error = standard_error,
    t_value = t_value,
    margin_of_error = margin_of_error,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  )
}

#' Confidence Interval for the Difference in Proportions
#'
#' Calculates the confidence interval for the difference between two population proportions and shows intermediate steps.
#'
#' @param successes1 The number of successes in the first sample (integer).
#' @param sample_size1 The size of the first sample (integer).
#' @param successes2 The number of successes in the second sample (integer).
#' @param sample_size2 The size of the second sample (integer).
#' @param confidence_level The confidence level (e.g., 0.95 for 95 percent confidence).
#' @return A list with intermediate calculations and the confidence interval:
#'   \item{proportion1}{The first sample proportion (p1).}
#'   \item{proportion2}{The second sample proportion (p2).}
#'   \item{standard_error}{The standard error of the difference in proportions.}
#'   \item{z_value}{The critical z-value.}
#'   \item{margin_of_error}{The margin of error.}
#'   \item{lower_bound}{The lower bound of the confidence interval.}
#'   \item{upper_bound}{The upper bound of the confidence interval.}
#'
#' @examples
#' ci_diff_proportions(30, 100, 20, 100, 0.95)
#'
#' @export
ci_diff_proportions <- function(successes1, sample_size1, successes2, sample_size2, confidence_level = 0.95) {
  p1 <- successes1 / sample_size1
  p2 <- successes2 / sample_size2

  # Calculate standard error
  standard_error <- sqrt((p1 * (1 - p1) / sample_size1) + (p2 * (1 - p2) / sample_size2))

  # Calculate z-value
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)

  # Calculate margin of error
  margin_of_error <- z_value * standard_error

  # Calculate difference and confidence interval
  diff <- p1 - p2
  lower_bound <- diff - margin_of_error
  upper_bound <- diff + margin_of_error

  # Return detailed results
  list(
    proportion1 = p1,
    proportion2 = p2,
    standard_error = standard_error,
    z_value = z_value,
    margin_of_error = margin_of_error,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  )
}

#' Confidence Interval for Variance
#'
#' Calculates the confidence interval for the variance of a normally distributed population.
#' Intermediate steps are included in the output.
#'
#' @param data A numeric vector of sample data.
#' @param confidence_level The confidence level (e.g., 0.95 for 95 percent confidence).
#' @importFrom stats var qchisq
#' @return A list with intermediate calculations and the confidence interval:
#'   \item{sample_variance}{The sample variance.}
#'   \item{chi_sq_lower}{The lower critical chi-squared value.}
#'   \item{chi_sq_upper}{The upper critical chi-squared value.}
#'   \item{lower_bound}{The lower bound of the confidence interval.}
#'   \item{upper_bound}{The upper bound of the confidence interval.}
#'
#' @examples
#' ci_variance(c(10, 12, 15, 18, 20), 0.95)
#'
#' @export
ci_variance <- function(data, confidence_level = 0.95) {
  n <- length(data)
  s2 <- var(data)
  alpha <- 1 - confidence_level
  chi_sq_lower <- qchisq(1 - alpha / 2, df = n - 1)
  chi_sq_upper <- qchisq(alpha / 2, df = n - 1)

  lower_bound <- (n - 1) * s2 / chi_sq_lower
  upper_bound <- (n - 1) * s2 / chi_sq_upper

  list(
    sample_variance = s2,
    chi_sq_lower = chi_sq_lower,
    chi_sq_upper = chi_sq_upper,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  )
}

#' Confidence Interval for Difference in Means
#'
#' Calculates the confidence interval for the difference in means of two independent samples.
#' Intermediate steps are included in the output.
#'
#' @param data1 A numeric vector of the first sample data.
#' @param data2 A numeric vector of the second sample data.
#' @param confidence_level The confidence level (e.g., 0.95 for 95 percent confidence).
#' @importFrom stats var
#' @return A list with intermediate calculations and the confidence interval:
#'   \item{mean_difference}{The difference in sample means.}
#'   \item{pooled_standard_error}{The pooled standard error of the difference.}
#'   \item{t_value}{The critical t-value.}
#'   \item{margin_of_error}{The margin of error.}
#'   \item{lower_bound}{The lower bound of the confidence interval.}
#'   \item{upper_bound}{The upper bound of the confidence interval.}
#'
#' @examples
#' ci_diff_means(c(15, 18, 21, 24), c(10, 12, 14, 16), 0.95)
#'
#' @export
ci_diff_means <- function(data1, data2, confidence_level = 0.95) {
  n1 <- length(data1)
  n2 <- length(data2)
  mean1 <- mean(data1)
  mean2 <- mean(data2)
  var1 <- var(data1)
  var2 <- var(data2)

  pooled_se <- sqrt(var1 / n1 + var2 / n2)
  df <- (pooled_se^4) / (((var1 / n1)^2 / (n1 - 1)) + ((var2 / n2)^2 / (n2 - 1)))

  alpha <- 1 - confidence_level
  t_value <- qt(1 - alpha / 2, df = df)

  margin_of_error <- t_value * pooled_se
  diff <- mean1 - mean2

  list(
    mean_difference = diff,
    pooled_standard_error = pooled_se,
    t_value = t_value,
    margin_of_error = margin_of_error,
    lower_bound = diff - margin_of_error,
    upper_bound = diff + margin_of_error
  )
}

#' Confidence Interval for Median
#'
#' Calculates the confidence interval for the median using bootstrapping.
#' Intermediate steps are included in the output.
#'
#' @param data A numeric vector of sample data.
#' @param confidence_level The confidence level (e.g., 0.95 for 95 percent confidence).
#' @param n_bootstrap The number of bootstrap samples (default is 1000).
#' @importFrom stats median quantile
#' @return A list with intermediate calculations and the confidence interval:
#'   \item{bootstrap_medians}{A vector of bootstrapped medians.}
#'   \item{lower_bound}{The lower bound of the confidence interval.}
#'   \item{upper_bound}{The upper bound of the confidence interval.}
#'
#' @examples
#' ci_median(c(5, 10, 15, 20, 25), 0.95, 1000)
#'
#' @export
ci_median <- function(data, confidence_level = 0.95, n_bootstrap = 1000) {
  alpha <- 1 - confidence_level
  bootstrap_medians <- replicate(n_bootstrap, median(sample(data, length(data), replace = TRUE)))

  lower_bound <- quantile(bootstrap_medians, probs = alpha / 2)
  upper_bound <- quantile(bootstrap_medians, probs = 1 - alpha / 2)

  list(
    bootstrap_medians = bootstrap_medians,
    lower_bound = lower_bound,
    upper_bound = upper_bound
  )
}

#' @title Confidence Interval for True Standard Deviation
#' @description Constructs a confidence interval for the population standard deviation based on the sample standard deviation.
#' @param sample_sd The sample standard deviation (\( s \)).
#' @param n The sample size.
#' @param confidence_level The desired confidence level (e.g., 0.96 for 96 percent confidence).
#' @importFrom stats qchisq
#' @return A list containing the lower and upper bounds of the confidence interval and an explanation of the calculations.
#' @examples
#' ci_true_sd(sample_sd = 0.81, n = 31, confidence_level = 0.96)
#' @export
ci_true_sd <- function(sample_sd, n, confidence_level) {
  # Degrees of freedom
  df <- n - 1

  # Alpha
  alpha <- 1 - confidence_level

  # Chi-square critical values
  chi_sq_lower <- qchisq(alpha / 2, df, lower.tail = FALSE)
  chi_sq_upper <- qchisq(1 - alpha / 2, df, lower.tail = FALSE)

  # Confidence interval bounds for variance
  variance_lower <- (df * sample_sd^2) / chi_sq_upper
  variance_upper <- (df * sample_sd^2) / chi_sq_lower

  # Confidence interval bounds for standard deviation
  sd_lower <- sqrt(variance_lower)
  sd_upper <- sqrt(variance_upper)

  # Explanation
  explanation <- paste(
    "To construct a confidence interval for the true standard deviation:",
    "\nDegrees of freedom (df) =", df,
    "\nChi-square critical values: Lower =", round(chi_sq_lower, 4), ", Upper =", round(chi_sq_upper, 4),
    "\nConfidence interval for variance: Lower =", round(variance_lower, 4), ", Upper =", round(variance_upper, 4),
    "\nConfidence interval for standard deviation: Lower =", round(sd_lower, 4), ", Upper =", round(sd_upper, 4)
  )

  return(list(
    lower_bound = sd_lower,
    upper_bound = sd_upper,
    explanation = explanation
  ))
}


#' Calculate the Confidence Interval for a New Confidence Level with Detailed Information
#'
#' This function calculates the confidence interval for a new confidence level based on a given original
#' confidence interval, sample size, and critical t-values for both the original and new confidence levels.
#' It uses the t-distribution to calculate the new margin of error and confidence interval.
#'
#' @param lower_original The lower bound of the original confidence interval.
#' @param upper_original The upper bound of the original confidence interval.
#' @param n The sample size.
#' @param confidence_original The original confidence level (e.g., 0.90 for 90 percent).
#' @param confidence_new The new confidence level (e.g., 0.99 for 99 percent).
#' @return A list containing the following:
#'   - `new_confidence_interval`: The new confidence interval as a vector.
#'   - `sample_mean`: The sample mean calculated from the original confidence interval.
#'   - `margin_error_original`: The margin of error for the original confidence interval.
#'   - `t_original`: The t-value for the original confidence level.
#'   - `t_new`: The t-value for the new confidence level.
#'   - `sample_standard_deviation`: The sample standard deviation.
#'   - `margin_error_new`: The margin of error for the new confidence interval.
#'
#' @examples
#' ci_to_new_ci(15.34, 36.66, 5, 0.90, 0.99)
#' @export
ci_to_new_ci <- function(lower_original, upper_original, n, confidence_original, confidence_new) {

  # Calculate the sample mean from the original confidence interval
  sample_mean <- (lower_original + upper_original) / 2

  # Calculate the margin of error from the original confidence interval
  margin_error_original <- (upper_original - lower_original) / 2

  # Calculate sample standard deviation (s) based on the margin of error and sample size
  t_original <- qt(1 - (1 - confidence_original) / 2, df = n - 1)  # t value for the original confidence level
  s <- margin_error_original * sqrt(n) / t_original

  # Calculate the margin of error for the new confidence interval
  t_new <- qt(1 - (1 - confidence_new) / 2, df = n - 1)  # t value for the new confidence level
  margin_error_new <- t_new * (s / sqrt(n))

  # Calculate the new confidence interval
  lower_new <- sample_mean - margin_error_new
  upper_new <- sample_mean + margin_error_new

  # Return the detailed results as a list
  return(list(
    new_confidence_interval = c(lower_new, upper_new),
    sample_mean = sample_mean,
    margin_error_original = margin_error_original,
    t_original = t_original,
    t_new = t_new,
    sample_standard_deviation = s,
    margin_error_new = margin_error_new
  ))
}

#' Calculate the required sample size for a X Percent confidence interval with a specified margin of error
#'
#' This function calculates the required sample size for a confidence interval of a population proportion
#' with a given margin of error and confidence level.
#'
#' @param margin_of_error The desired margin of error for the confidence interval (numeric).
#' @param confidence_level The desired confidence level (numeric, between 0 and 1).
#' @param p_estimate The estimated proportion (default is 0.5 for a fair coin).
#'
#' @return The required sample size (numeric).
#'
#' @examples
#' sample_size_ci(0.1, 0.9)
#' @export
sample_size_ci <- function(margin_of_error, confidence_level, p_estimate = 0.5) {
  z_score <- qnorm(1 - (1 - confidence_level) / 2)
  n <- (z_score^2 * p_estimate * (1 - p_estimate)) / (margin_of_error^2)
  return(ceiling(n))
}


#' @title Required Sample Size for Mean
#' @description Calculates the sample size required to estimate the mean of a normal distribution with a specified margin of error and confidence level.
#' @param sd The standard deviation of the population.
#' @param margin_of_error The desired margin of error (precision).
#' @param confidence_level The confidence level (e.g., 0.90 for 90 percent confidence).
#' @return A list containing the required sample size and an explanation of the calculation.
#' @examples
#' sample_size_mean(sd = 0.03, margin_of_error = 0.005, confidence_level = 0.90)
#' @export
sample_size_mean <- function(sd, margin_of_error, confidence_level) {
  # Validate inputs
  if (sd <= 0 || margin_of_error <= 0 || confidence_level <= 0 || confidence_level >= 1) {
    stop("Standard deviation, margin of error, and confidence level must be positive, and confidence level must be between 0 and 1.")
  }

  # Calculate the critical z-value
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)

  # Calculate the required sample size
  n <- (z_value * sd / margin_of_error)^2

  # Round up to the next whole number
  n <- ceiling(n)

  explanation <- paste(
    "To achieve a margin of error of +/-", margin_of_error,  # Replace ± with +/-
    "with a confidence level of", confidence_level, ":",
    "\n1. The critical z-value is", round(z_value, 4),
    "\n2. The required sample size is calculated as:",
    "\n   n = (z * sd / margin_of_error)^2",
    "\n   =", "(", round(z_value, 4), "*", sd, "/", margin_of_error, ")^2",
    "\n   =", n
  )


  return(list(sample_size = n, explanation = explanation))
}

#' @title Sample Size for Proportion Confidence Interval
#' @description Calculates the required sample size to achieve a specified confidence interval width for a population proportion.
#' @param p The estimated population proportion (e.g., 0.5 for a fair coin).
#' @param confidence_level The desired confidence level (e.g., 0.95 for 95 percent confidence).
#' @param width The desired total width of the confidence interval.
#' @return A list containing the required sample size and an explanation of the calculation.
#' @examples
#' sample_size_proportion(p = 0.5, confidence_level = 0.95, width = 0.05)
#' @export
sample_size_proportion <- function(p, confidence_level, width) {
  # Validate inputs
  if (p <= 0 || p >= 1) {
    stop("The proportion (p) must be between 0 and 1.")
  }
  if (confidence_level <= 0 || confidence_level >= 1) {
    stop("The confidence level must be between 0 and 1.")
  }
  if (width <= 0) {
    stop("The width must be a positive number.")
  }

  # Calculate critical z-value
  alpha <- 1 - confidence_level
  z <- qnorm(1 - alpha / 2)

  # Calculate required sample size
  margin_of_error <- width / 2
  n <- (z^2 * p * (1 - p)) / (margin_of_error^2)

  # Round up to the next whole number
  n <- ceiling(n)

  # Explanation
  explanation <- paste(
    "To achieve a confidence interval width of", width,
    "with a confidence level of", confidence_level, ":",
    "\n1. The critical z-value is", round(z, 4),
    "\n2. The margin of error is", margin_of_error,
    "\n3. The required sample size is calculated as:",
    "\n   n = (z^2 * p * (1 - p)) / (margin_of_error^2)",
    "\n   =", "(", round(z, 4), "^2 *", p, "*", (1 - p), ") /", margin_of_error^2,
    "\n   =", n
  )

  return(list(sample_size = n, explanation = explanation))
}

#' Mean and Standard Deviation for Z
#'
#' This function returns the mean and standard deviation for the standard normal random variable Z,
#' which is always 0 for the mean and 1 for the standard deviation.
#'
#' @return A list containing the mean and standard deviation of Z.
#' @importFrom stats pnorm
#' @examples
#' normal_mean_and_sd_Z()
#'
#' @export
normal_mean_and_sd_Z <- function() {
  list(mean = 0, sd = 1)
}

#' Plot the Distribution of Z
#'
#' This function plots the standard normal distribution with a mean of 0 and a standard deviation of 1.
#'
#' @param x_range The range of x values to plot the distribution.
#' @importFrom stats pnorm dnorm
#' @importFrom graphics curve
#' @return A plot of the standard normal distribution.
#'
#' @examples
#' normal_plot_standard_normal_distribution()
#'
#' @export
normal_plot_standard_normal_distribution <- function(x_range = c(-4, 4)) {
  curve(dnorm(x), from = x_range[1], to = x_range[2], col = "blue",
        main = "Standard Normal Distribution", xlab = "Z", ylab = "Density")
}

#' Find Probability P(Z < target)
#'
#' This function calculates the probability that the standard normal variable Z is less than a given target.
#'
#' @param target The value of Z for which to calculate the probability.
#' @importFrom stats pnorm
#' @return The probability P(Z < target).
#'
#' @examples
#' normal_prob_Z_less_than_target(1.2)
#'
#' @export
normal_prob_Z_less_than_target <- function(target) {
  pnorm(target)
}

#' Find Probability P(Z > target)
#'
#' This function calculates the probability that the standard normal variable Z is greater than a given target.
#'
#' @param target The value of Z for which to calculate the probability.
#' @importFrom stats pnorm
#' @return The probability P(Z > target).
#'
#' @examples
#' normal_prob_Z_greater_than_target(1.2)
#'
#' @export
normal_prob_Z_greater_than_target <- function(target) {
  1 - pnorm(target)
}

#' Find Probability P(lower < Z < upper)
#'
#' This function calculates the probability that the standard normal variable Z lies between two values, lower and upper.
#'
#' @param lower The lower bound for Z.
#' @param upper The upper bound for Z.
#'
#' @return The probability P(lower < Z < upper).
#' @importFrom stats pnorm
#' @examples
#' normal_prob_Z_between_target_values(-0.45, 1.96)
#'
#' @export
normal_prob_Z_between_target_values <- function(lower, upper) {
  pnorm(upper) - pnorm(lower)
}

#' Explanation of Type I and Type II Errors
#'
#' This function explains the concept of Type I and Type II errors in the context of hypothesis testing.
#' It provides an interpretation based on a given scenario, where the null hypothesis and alternative hypothesis are defined.
#'
#' @param null_hypothesis The null hypothesis.
#' @param alternative_hypothesis The alternative hypothesis.
#' @param error_type The type of error to explain: "Type I" or "Type II".
#' @return A string explaining the selected type of error.
#'
#' @examples
#' error_explanation(
#'   "The team will not get the first down",
#'   "The team will get the first down",
#'   "Type I"
#' )
#' @export
error_explanation <- function(null_hypothesis, alternative_hypothesis, error_type) {
  if (error_type == "Type I") {
    return(paste("Type I Error (False Positive): You reject the null hypothesis ('", null_hypothesis,
                 "') when in fact it is true. In this case, you would decide to go for the first down, but the team will not succeed.", sep = ""))
  } else if (error_type == "Type II") {
    return(paste("Type II Error (False Negative): You fail to reject the null hypothesis ('", null_hypothesis,
                 "') when in fact it is false. In this case, you would decide to punt, but the team would have succeeded in getting the first down.", sep = ""))
  } else {
    return("Invalid error type. Please choose either 'Type I' or 'Type II'.")
  }
}

#' Find c such that P(Z < c) = target_probability
#'
#' This function finds the value of c such that P(Z < c) equals a given probability.
#'
#' @param target_probability The target probability for which to find the value of c.
#' @importFrom stats qnorm
#' @return The value of c such that P(Z < c) = target_probability.
#'
#' @examples
#' find_c_for_probability_target(0.845)
#'
#' @export
find_c_for_probability_target <- function(target_probability) {
  qnorm(target_probability)
}

#' Find c such that P(Z > c) = target_probability
#'
#' This function finds the value of c such that P(Z > c) equals a given probability.
#'
#' @param target_probability The target probability for the upper tail.
#' @importFrom stats qnorm
#' @return The value of c such that P(Z > c) = target_probability.
#'
#' @examples
#' find_c_for_probability_greater_than_target(0.845)
#'
#' @export
find_c_for_probability_greater_than_target <- function(target_probability) {
  qnorm(1 - target_probability)
}

#' Find c such that P(-c < Z < c) = target_probability
#'
#' This function finds the value of c such that P(-c < Z < c) equals a given probability.
#' @importFrom stats qnorm
#' @param target_probability The desired cumulative probability for the two-tailed probability.
#'
#' @return The value of c such that P(-c < Z < c) = target_probability.
#'
#' @examples
#' find_c_for_two_tailed_probability_target(0.845)
#'
#' @export
find_c_for_two_tailed_probability_target <- function(target_probability) {
  qnorm((1 + target_probability) / 2)
}


#' Calculate Sample Mean
#'
#' This function calculates the mean of a given sample.
#'
#' @param data A numeric vector of sample data.
#'
#' @return The sample mean.
#'
#' @examples
#' sample_mean(c(5, 10, 15))
#'
#' @export
sample_mean <- function(data) {
  mean(data)
}

#' Calculate Standard Error of the Sample Mean
#'
#' This function calculates the standard error of the sample mean given the population standard deviation and sample size.
#'
#' @param sd The standard deviation of the population.
#' @param n The sample size.
#'
#' @return The standard error of the sample mean.
#'
#' @examples
#' sample_mean_sd(sd = 10, n = 25)
#'
#' @export
sample_mean_sd <- function(sd, n) {
  sd / sqrt(n)
}

#' Probability for Sample Mean
#'
#' Calculates the probability that the sample mean is less than or greater than a given target value.
#'
#' @param mu The population mean.
#' @param sd The population standard deviation.
#' @param n The sample size.
#' @param target The target value for the sample mean.
#' @param type The type of probability to calculate. Use "less" for P(X < target) or "greater" for P(X > target).
#'
#' @return A list containing the z-score and the calculated probability.
#'
#' @examples
#' sample_mean_probability(mu = 100, sd = 15, n = 30, target = 105, type = "less")
#'
#' @export
sample_mean_probability <- function(mu, sd, n, target, type = "less") {
  standard_error <- sd / sqrt(n)
  z_score <- (target - mu) / standard_error

  prob <- if (type == "less") {
    pnorm(z_score)
  } else if (type == "greater") {
    1 - pnorm(z_score)
  } else {
    stop("Invalid type. Use 'less' or 'greater'.")
  }

  list(z_score = z_score, probability = prob)
}

#' Confidence Interval for Sample Mean
#'
#' Calculates the confidence interval for the sample mean given the sample mean, population standard deviation, sample size, and confidence level.
#'
#' @param mean The sample mean.
#' @param sd The population standard deviation.
#' @param n The sample size.
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#'
#' @return A list containing the lower and upper bounds of the confidence interval.
#'
#' @examples
#' sample_mean_confidence_interval(mean = 20, sd = 5, n = 30, confidence_level = 0.95)
#'
#' @export
sample_mean_confidence_interval <- function(mean, sd, n, confidence_level) {
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)
  standard_error <- sd / sqrt(n)
  margin_of_error <- z_value * standard_error

  list(
    lower_bound = mean - margin_of_error,
    upper_bound = mean + margin_of_error
  )
}

#' Required Sample Size for Sample Mean
#'
#' Calculates the required sample size to estimate the sample mean with a specified margin of error and confidence level.
#'
#' @param sd The standard deviation of the population.
#' @param margin_of_error The desired margin of error.
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#'
#' @return The required sample size as an integer.
#'
#' @examples
#' sample_mean_required_size(sd = 10, margin_of_error = 2, confidence_level = 0.95)
#'
#' @export
sample_mean_required_size <- function(sd, margin_of_error, confidence_level) {
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)
  n <- ceiling((z_value * sd / margin_of_error)^2)

  n
}

#' Probability Between Two Values for Sample Mean
#'
#' Calculates the probability that the sample mean lies between two specified values.
#'
#' @param mu The population mean.
#' @param sd The population standard deviation.
#' @param n The sample size.
#' @param lower The lower bound for the sample mean.
#' @param upper The upper bound for the sample mean.
#'
#' @return The probability that the sample mean lies between the specified bounds.
#'
#' @examples
#' sample_mean_probability_between(mu = 100, sd = 15, n = 30, lower = 95, upper = 105)
#'
#' @export
sample_mean_probability_between <- function(mu, sd, n, lower, upper) {
  standard_error <- sd / sqrt(n)
  z_lower <- (lower - mu) / standard_error
  z_upper <- (upper - mu) / standard_error

  prob <- pnorm(z_upper) - pnorm(z_lower)

  prob
}

#' Z-Score for Sample Mean
#'
#' Calculates the z-score for a given sample mean.
#'
#' @param mean The sample mean.
#' @param mu The population mean.
#' @param sd The population standard deviation.
#' @param n The sample size.
#'
#' @return The z-score for the sample mean.
#'
#' @examples
#' sample_mean_z_score(mean = 105, mu = 100, sd = 15, n = 30)
#'
#' @export
sample_mean_z_score <- function(mean, mu, sd, n) {
  standard_error <- sd / sqrt(n)
  z_score <- (mean - mu) / standard_error

  z_score
}

#' Calculate Sample Proportion
#'
#' This function calculates the proportion of successes in a given sample.
#'
#' @param successes The number of successes in the sample (integer).
#' @param sample_size The total size of the sample (integer).
#'
#' @return The sample proportion.
#'
#' @examples
#' proportion(successes = 50, sample_size = 100)
#'
#' @export
proportion <- function(successes, sample_size) {
  successes / sample_size
}

#' Calculate Standard Error of the Proportion
#'
#' This function calculates the standard error of a proportion given the sample proportion and sample size.
#'
#' @param p_hat The sample proportion.
#' @param sample_size The total size of the sample (integer).
#'
#' @return The standard error of the proportion.
#'
#' @examples
#' proportion_sd(p_hat = 0.5, sample_size = 100)
#'
#' @export
proportion_sd <- function(p_hat, sample_size) {
  sqrt((p_hat * (1 - p_hat)) / sample_size)
}

#' Probability for Proportion
#'
#' Calculates the probability that the sample proportion is less than or greater than a given target value.
#'
#' @param p The population proportion.
#' @param sample_size The total size of the sample (integer).
#' @param target The target proportion for comparison.
#' @param type The type of probability to calculate. Use "less" for P(p_hat < target) or "greater" for P(p_hat > target).
#'
#' @return A list containing the z-score and the calculated probability.
#'
#' @examples
#' proportion_probability(p = 0.5, sample_size = 100, target = 0.6, type = "less")
#'
#' @export
proportion_probability <- function(p, sample_size, target, type = "less") {
  standard_error <- sqrt((p * (1 - p)) / sample_size)
  z_score <- (target - p) / standard_error

  prob <- if (type == "less") {
    pnorm(z_score)
  } else if (type == "greater") {
    1 - pnorm(z_score)
  } else {
    stop("Invalid type. Use 'less' or 'greater'.")
  }

  list(z_score = z_score, probability = prob)
}

#' Confidence Interval for Proportion
#'
#' Calculates the confidence interval for a sample proportion given the sample size and confidence level.
#'
#' @param successes The number of successes in the sample (integer).
#' @param sample_size The total size of the sample (integer).
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#'
#' @return A list containing the lower and upper bounds of the confidence interval.
#'
#' @examples
#' proportion_confidence_interval(successes = 50, sample_size = 100, confidence_level = 0.95)
#'
#' @export
proportion_confidence_interval <- function(successes, sample_size, confidence_level) {
  p_hat <- successes / sample_size
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)
  standard_error <- sqrt((p_hat * (1 - p_hat)) / sample_size)
  margin_of_error <- z_value * standard_error

  list(
    lower_bound = p_hat - margin_of_error,
    upper_bound = p_hat + margin_of_error
  )
}

#' Required Sample Size for Proportion
#'
#' Calculates the required sample size to estimate a proportion with a specified margin of error and confidence level.
#'
#' @param p The estimated population proportion.
#' @param margin_of_error The desired margin of error.
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#'
#' @return The required sample size as an integer.
#'
#' @examples
#' proportion_required_size(p = 0.5, margin_of_error = 0.05, confidence_level = 0.95)
#'
#' @export
proportion_required_size <- function(p, margin_of_error, confidence_level) {
  alpha <- 1 - confidence_level
  z_value <- qnorm(1 - alpha / 2)
  n <- ceiling((z_value^2 * p * (1 - p)) / margin_of_error^2)

  n
}

#' Probability Between Two Proportions
#'
#' Calculates the probability that the sample proportion lies between two specified values.
#'
#' @param p The population proportion.
#' @param sample_size The total size of the sample (integer).
#' @param lower The lower bound for the proportion.
#' @param upper The upper bound for the proportion.
#'
#' @return The probability that the sample proportion lies between the specified bounds.
#'
#' @examples
#' proportion_probability_between(p = 0.5, sample_size = 100, lower = 0.4, upper = 0.6)
#'
#' @export
proportion_probability_between <- function(p, sample_size, lower, upper) {
  standard_error <- sqrt((p * (1 - p)) / sample_size)
  z_lower <- (lower - p) / standard_error
  z_upper <- (upper - p) / standard_error

  prob <- pnorm(z_upper) - pnorm(z_lower)

  prob
}

#' Z-Score for Proportion
#'
#' Calculates the z-score for a given sample proportion.
#'
#' @param p_hat The sample proportion.
#' @param p The population proportion.
#' @param sample_size The total size of the sample (integer).
#'
#' @return The z-score for the sample proportion.
#'
#' @examples
#' proportion_z_score(p_hat = 0.6, p = 0.5, sample_size = 100)
#'
#' @export
proportion_z_score <- function(p_hat, p, sample_size) {
  standard_error <- sqrt((p * (1 - p)) / sample_size)
  z_score <- (p_hat - p) / standard_error

  z_score
}

#' One-Sample t-Test for Mean
#'
#' Conducts a one-sample t-test for the population mean.
#'
#' @param data A numeric vector of sample data.
#' @param mu The hypothesized population mean.
#' @param alternative The alternative hypothesis. Use "two.sided", "less", or "greater".
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#' @importFrom stats pt
#' @return A list containing the test statistic, p-value, and confidence interval.
#'
#' @examples
#' one_sample_t_test_mean(c(5, 10, 15), mu = 10, alternative = "two.sided", confidence_level = 0.95)
#'
#' @export
one_sample_t_test_mean <- function(data, mu, alternative = "two.sided", confidence_level = 0.95) {
  n <- length(data)
  sample_mean <- mean(data)
  sample_sd <- sd(data)
  standard_error <- sample_sd / sqrt(n)
  t_statistic <- (sample_mean - mu) / standard_error

  p_value <- switch(alternative,
                    "two.sided" = 2 * pt(-abs(t_statistic), df = n - 1),
                    "less" = pt(t_statistic, df = n - 1),
                    "greater" = 1 - pt(t_statistic, df = n - 1),
                    stop("Invalid alternative hypothesis")
  )

  alpha <- 1 - confidence_level
  t_critical <- qt(1 - alpha / 2, df = n - 1)
  margin_of_error <- t_critical * standard_error

  confidence_interval <- switch(alternative,
                                "two.sided" = c(sample_mean - margin_of_error, sample_mean + margin_of_error),
                                "less" = c(-Inf, sample_mean + margin_of_error),
                                "greater" = c(sample_mean - margin_of_error, Inf)
  )

  list(
    t_statistic = t_statistic,
    p_value = p_value,
    confidence_interval = confidence_interval
  )
}

#' Two-Sample t-Test for Means
#'
#' Conducts a two-sample t-test for the difference in population means.
#'
#' @param data1 A numeric vector of the first sample data.
#' @param data2 A numeric vector of the second sample data.
#' @param alternative The alternative hypothesis. Use "two.sided", "less", or "greater".
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#' @param equal_variance Logical, whether to assume equal variances (default is FALSE).
#' @importFrom stats pt
#' @return A list containing the test statistic, p-value, and confidence interval.
#'
#' @examples
#' two_sample_t_test_means(c(5, 10, 15),
#' c(8, 12, 18),
#' alternative = "two.sided",
#' confidence_level = 0.95,
#' equal_variance = FALSE)
#'
#' @export
two_sample_t_test_means <- function(data1, data2, alternative = "two.sided", confidence_level = 0.95, equal_variance = FALSE) {
  n1 <- length(data1)
  n2 <- length(data2)
  mean1 <- mean(data1)
  mean2 <- mean(data2)
  sd1 <- sd(data1)
  sd2 <- sd(data2)

  if (equal_variance) {
    pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
    standard_error <- pooled_sd * sqrt(1 / n1 + 1 / n2)
    df <- n1 + n2 - 2
  } else {
    standard_error <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
    df <- ((sd1^2 / n1) + (sd2^2 / n2))^2 / (((sd1^2 / n1)^2 / (n1 - 1)) + ((sd2^2 / n2)^2 / (n2 - 1)))
  }

  t_statistic <- (mean1 - mean2) / standard_error

  p_value <- switch(alternative,
                    "two.sided" = 2 * pt(-abs(t_statistic), df = df),
                    "less" = pt(t_statistic, df = df),
                    "greater" = 1 - pt(t_statistic, df = df),
                    stop("Invalid alternative hypothesis")
  )

  alpha <- 1 - confidence_level
  t_critical <- qt(1 - alpha / 2, df = df)
  margin_of_error <- t_critical * standard_error

  confidence_interval <- switch(alternative,
                                "two.sided" = c((mean1 - mean2) - margin_of_error, (mean1 - mean2) + margin_of_error),
                                "less" = c(-Inf, (mean1 - mean2) + margin_of_error),
                                "greater" = c((mean1 - mean2) - margin_of_error, Inf)
  )

  list(
    t_statistic = t_statistic,
    p_value = p_value,
    confidence_interval = confidence_interval
  )
}

#' One-Sample z-Test for Proportion
#'
#' Conducts a one-sample z-test for a population proportion.
#'
#' @param successes The number of successes in the sample.
#' @param sample_size The size of the sample.
#' @param p The hypothesized population proportion.
#' @param alternative The alternative hypothesis. Use "two.sided", "less", or "greater".
#'
#' @return A list containing the z-statistic and p-value.
#'
#' @examples
#' one_sample_z_test_proportion(successes = 50, sample_size = 100, p = 0.5, alternative = "two.sided")
#'
#' @export
one_sample_z_test_proportion <- function(successes, sample_size, p, alternative = "two.sided") {
  p_hat <- successes / sample_size
  standard_error <- sqrt(p * (1 - p) / sample_size)
  z_statistic <- (p_hat - p) / standard_error

  p_value <- switch(alternative,
                    "two.sided" = 2 * pnorm(-abs(z_statistic)),
                    "less" = pnorm(z_statistic),
                    "greater" = 1 - pnorm(z_statistic),
                    stop("Invalid alternative hypothesis")
  )

  list(
    z_statistic = z_statistic,
    p_value = p_value
  )
}

#' Two-Sample z-Test for Proportions
#'
#' Conducts a two-sample z-test for the difference in population proportions.
#'
#' @param successes1 The number of successes in the first sample.
#' @param sample_size1 The size of the first sample.
#' @param successes2 The number of successes in the second sample.
#' @param sample_size2 The size of the second sample.
#' @param alternative The alternative hypothesis. Use "two.sided", "less", or "greater".
#'
#' @return A list containing the z-statistic and p-value.
#'
#' @examples
#' two_sample_z_test_proportions(successes1 = 50,
#' sample_size1 = 100,
#' successes2 = 60,
#' sample_size2 = 120,
#' alternative = "two.sided")
#'
#' @export
two_sample_z_test_proportions <- function(successes1, sample_size1, successes2, sample_size2, alternative = "two.sided") {
  p1 <- successes1 / sample_size1
  p2 <- successes2 / sample_size2
  pooled_p <- (successes1 + successes2) / (sample_size1 + sample_size2)
  standard_error <- sqrt(pooled_p * (1 - pooled_p) * (1 / sample_size1 + 1 / sample_size2))
  z_statistic <- (p1 - p2) / standard_error

  p_value <- switch(alternative,
                    "two.sided" = 2 * pnorm(-abs(z_statistic)),
                    "less" = pnorm(z_statistic),
                    "greater" = 1 - pnorm(z_statistic),
                    stop("Invalid alternative hypothesis")
  )

  list(
    z_statistic = z_statistic,
    p_value = p_value
  )
}



#' @title Normal Distribution: Between Two Values
#' @description Calculates the probability that a normal random variable is between two values.
#' @param mean The mean of the normal distribution.
#' @param sd The standard deviation of the normal distribution.
#' @param lower The lower bound of the range.
#' @param upper The upper bound of the range.
#' @return The probability that the random variable is between the given bounds.
#' @examples
#' normal_between(mean = 0, sd = 1, lower = -1, upper = 1)
#' @export
normal_between <- function(mean, sd, lower, upper) {
  prob <- pnorm(upper, mean, sd) - pnorm(lower, mean, sd)
  explanation <- paste(
    "The probability that X is between", lower, "and", upper, "is:",
    "P(", lower, "< X <", upper, ") = P(X <", upper, ") - P(X <", lower, ") =",
    round(prob, 4)
  )
  return(list(probability = prob, explanation = explanation))
}

#' @title Normal Distribution: Greater Than
#' @description Calculates the probability that a normal random variable is greater than a given value.
#' @param mean The mean of the normal distribution.
#' @param sd The standard deviation of the normal distribution.
#' @param value The value above which the probability is calculated.
#' @return The probability that the random variable is greater than the given value.
#' @examples
#' normal_greater(mean = 0, sd = 1, value = 1.5)
#' @export
normal_greater <- function(mean, sd, value) {
  prob <- 1 - pnorm(value, mean, sd)
  explanation <- paste(
    "The probability that X is greater than", value, "is:",
    "P(X >", value, ") = 1 - P(X <", value, ") =",
    round(prob, 4)
  )
  return(list(probability = prob, explanation = explanation))
}

#' @title Normal Distribution: Less Than
#' @description Calculates the probability that a normal random variable is less than a given value.
#' @param mean The mean of the normal distribution.
#' @param sd The standard deviation of the normal distribution.
#' @param value The value below which the probability is calculated.
#' @return The probability that the random variable is less than the given value.
#' @examples
#' normal_less(mean = 0, sd = 1, value = -1.5)
#' @export
normal_less <- function(mean, sd, value) {
  prob <- pnorm(value, mean, sd)
  explanation <- paste(
    "The probability that X is less than", value, "is:",
    "P(X <", value, ") =", round(prob, 4)
  )
  return(list(probability = prob, explanation = explanation))
}

#' @title Normal Distribution: Exactly
#' @description Calculates the probability that a normal random variable equals a given value (theoretically zero for continuous distributions).
#' @param mean The mean of the normal distribution.
#' @param sd The standard deviation of the normal distribution.
#' @param value The value at which the probability is calculated.
#' @return The probability that the random variable equals the given value (always zero for continuous variables).
#' @examples
#' normal_exact(mean = 0, sd = 1, value = 0)
#' @export
normal_exact <- function(mean, sd, value) {
  explanation <- paste(
    "The probability that X equals", value, "is:",
    "P(X =", value, ") = 0 for continuous distributions."
  )
  return(list(probability = 0, explanation = explanation))
}




#' Find X Value for a Given Probability
#'
#' This function finds the value of `x` such that the probability of a normally
#' distributed random variable being greater than `x` is equal to the specified
#' probability. This is useful for determining threshold values in statistical analyses.
#'
#' The calculation involves the following steps:
#'
#' 1. Compute the z-score corresponding to the specified probability using the quantile function
#'    for the standard normal distribution:
#'    z = qnorm(1 - probability)
#'    This gives the z-value such that the area to the right of z in the standard normal
#'    distribution equals the specified probability.
#'
#' 2. Convert the z-score to the corresponding x value using the formula:
#'    x = mu + z * sigma
#'    where mu is the mean and sigma is the standard deviation of the distribution.
#'
#' @param mean The mean of the distribution (numeric).
#' @param sd The standard deviation of the distribution (numeric).
#' @param probability The target probability (numeric, must be between 0 and 1).
#'
#' @return The `x` value corresponding to the specified probability (numeric).
#'
#' @examples
#' find_x_value(100, 15, 0.05)
#'
#' @export
find_x_value <- function(mean, sd, probability) {
  # Calculate the z-score corresponding to the target probability
  z_value <- qnorm(1 - probability)

  # Calculate and return the corresponding x value
  mean + z_value * sd
}

#' @title Find X Value Where P(X < x)
#' @description Finds the value of `x` such that the probability of a normally distributed random variable being less than `x` equals the specified probability.
#' @param mean The mean of the distribution (numeric).
#' @param sd The standard deviation of the distribution (numeric).
#' @param probability The target cumulative probability (numeric, must be between 0 and 1).
#' @return The `x` value corresponding to the specified probability (numeric).
#' @examples
#' find_x_less(mean = 100, sd = 15, probability = 0.05)
#' @export
find_x_less <- function(mean, sd, probability) {
  # Calculate the z-score corresponding to the target cumulative probability
  z_value <- qnorm(probability)
  # Calculate and return the corresponding x value
  x_value <- mean + z_value * sd
  explanation <- paste(
    "The value of x such that P(X < x) =", probability, "is calculated as:",
    "x =", mean, "+", z_value, "*", sd, "=", round(x_value, 4)
  )
  return(list(x_value = x_value, explanation = explanation))
}

#' @title Find X Value Where P(X > x)
#' @description Finds the value of `x` such that the probability of a normally distributed random variable being greater than `x` equals the specified probability.
#' @param mean The mean of the distribution (numeric).
#' @param sd The standard deviation of the distribution (numeric).
#' @param probability The target probability (numeric, must be between 0 and 1).
#' @return The `x` value corresponding to the specified probability (numeric).
#' @examples
#' find_x_greater(mean = 100, sd = 15, probability = 0.05)
#' @export
find_x_greater <- function(mean, sd, probability) {
  # Use the find_x_value function for this case
  x_value <- find_x_value(mean, sd, probability)
  explanation <- paste(
    "The value of x such that P(X > x) =", probability, "is calculated using the formula:",
    "x =", mean, "+ z *", sd, "where z corresponds to P(Z > z) =", probability,
    "\nx =", round(x_value, 4)
  )
  return(list(x_value = x_value, explanation = explanation))
}

#' @title Find X Values Where P(lower < X < upper)
#' @description Finds two `x` values such that the probability of a normally distributed random variable being between them equals the specified probability.
#' @param mean The mean of the distribution (numeric).
#' @param sd The standard deviation of the distribution (numeric).
#' @param probability The target probability (numeric, must be between 0 and 1).
#' @return A list containing the lower and upper `x` values corresponding to the specified probability.
#' @examples
#' find_x_between(mean = 100, sd = 15, probability = 0.95)
#' @export
find_x_between <- function(mean, sd, probability) {
  # Calculate the z-scores corresponding to the middle probability
  z_value <- qnorm((1 + probability) / 2)
  x_upper <- mean + z_value * sd
  x_lower <- mean - z_value * sd
  explanation <- paste(
    "The values of x such that P(lower < X < upper) =", probability, "are calculated as:",
    "x_lower =", mean, "-", z_value, "*", sd, "=", round(x_lower, 4),
    "and x_upper =", mean, "+", z_value, "*", sd, "=", round(x_upper, 4)
  )
  return(list(x_lower = x_lower, x_upper = x_upper, explanation = explanation))
}

#' @title Find Symmetric X Values Around Mean
#' @description Finds two `x` values symmetric around the mean such that the probability of a normally distributed random variable being between them equals the specified probability.
#' @param mean The mean of the distribution (numeric).
#' @param sd The standard deviation of the distribution (numeric).
#' @param probability The target probability (numeric, must be between 0 and 1).
#' @return A list containing the lower and upper symmetric `x` values.
#' @examples
#' find_x_symmetric(mean = 100, sd = 15, probability = 0.95)
#' @export
find_x_symmetric <- function(mean, sd, probability) {
  # This is effectively the same as the between calculation
  result <- find_x_between(mean, sd, probability)
  explanation <- paste(
    "The symmetric values of x around the mean", mean, "are:",
    "x_lower =", round(result$x_lower, 4),
    "and x_upper =", round(result$x_upper, 4),
    "\nThese values ensure P(x_lower < X < x_upper) =", probability
  )
  return(list(x_lower = result$x_lower, x_upper = result$x_upper, explanation = explanation))
}

#' @title Verify Binomial Distribution
#' @description Verifies whether a scenario satisfies the conditions of a binomial distribution.
#' @param n The number of trials (e.g., group size).
#' @param p The probability of success on a single trial.
#' @return A logical value indicating whether it satisfies the binomial conditions and an explanation.
#' @examples
#' binomial_verify(n = 6, p = 0.9)
#' @export
binomial_verify <- function(n, p) {
  # Conditions for a binomial distribution
  independent_trials <- TRUE
  fixed_n <- n > 0
  success_probability <- p >= 0 && p <= 1

  is_binomial <- independent_trials && fixed_n && success_probability

  explanation <- paste(
    "Conditions for a binomial distribution:",
    "\n1. Fixed number of trials (n):", fixed_n,
    "\n2. Each trial is independent (assumed):", independent_trials,
    "\n3. Probability of success (p) is constant:", success_probability,
    "\nConclusion: This scenario", ifelse(is_binomial, "satisfies", "does not satisfy"), "the binomial conditions."
  )

  return(list(is_binomial = is_binomial, explanation = explanation))
}

#' @title Critical t* Value
#' @description Calculates the critical t* value for a given right-tail probability and degrees of freedom.
#' @param prob The right-tail probability (a value between 0 and 1).
#' @param df The degrees of freedom (a positive integer).
#' @return A list containing the critical t* value and an explanation of the calculation.
#' @examples
#' critical_t_value(prob = 0.10, df = 8)
#' @export
critical_t_value <- function(prob, df) {
  # Validate inputs
  if (prob <= 0 || prob >= 1) {
    stop("The probability (prob) must be between 0 and 1.")
  }
  if (df <= 0) {
    stop("The degrees of freedom (df) must be a positive integer.")
  }

  # Calculate t* value
  t_star <- qt(1 - prob, df)

  # Explanation
  explanation <- paste(
    "The critical t* value is calculated using the qt function:",
    "\n1. Input right-tail probability:", prob,
    "\n2. Degrees of freedom:", df,
    "\n3. The critical t* value corresponds to a cumulative probability of", 1 - prob,
    "\n   t* =", round(t_star, 4)
  )

  return(list(t_star = t_star, explanation = explanation))
}

#' @title Standard Error of Sample Proportion
#' @description Calculates the standard error of a sample proportion based on the proportion and sample size.
#' @param p_hat The sample proportion (a value between 0 and 1).
#' @param n The sample size (a positive integer).
#' @return A list containing the standard error and an explanation of the calculation.
#' @examples
#' proportion_standard_error(p_hat = 0.25, n = 200)
#' @export
proportion_standard_error <- function(p_hat, n) {
  # Validate inputs
  if (p_hat < 0 || p_hat > 1) {
    stop("The sample proportion (p_hat) must be between 0 and 1.")
  }
  if (n <= 0) {
    stop("The sample size (n) must be a positive integer.")
  }

  # Calculate standard error
  se <- sqrt((p_hat * (1 - p_hat)) / n)

  # Explanation
  explanation <- paste(
    "The standard error of the sample proportion is calculated as:",
    "\nSE = sqrt((p_hat * (1 - p_hat)) / n)",
    "\n   = sqrt((", p_hat, "*", 1 - p_hat, ") /", n, ")",
    "\n   =", round(se, 4)
  )

  return(list(standard_error = se, explanation = explanation))
}

#' @title Z-Score for a Percentile
#' @description Calculates the z-score corresponding to a given percentile under the standard normal curve.
#' @param percentile The percentile for which to find the corresponding z-score (value between 0 and 1).
#' @return A list containing the z-score and an explanation of the calculation.
#' @examples
#' z_score_for_percentile(percentile = 0.30) # 30th percentile
#' z_score_for_percentile(percentile = 0.95) # 95th percentile
#' @export
z_score_for_percentile <- function(percentile) {
  # Validate input
  if (percentile <= 0 || percentile >= 1) {
    stop("Percentile must be between 0 and 1 (exclusive).")
  }

  # Calculate z-score
  z_score <- qnorm(percentile)

  # Explanation
  explanation <- paste(
    "For the standard normal curve, the z-score corresponding to the",
    percentile * 100, "th percentile is:", round(z_score, 4)
  )

  return(list(z_score = z_score, explanation = explanation))
}


#' @title Find Z-Score for a Given Area
#' @description Finds the z-score corresponding to a given cumulative area under the standard normal curve.
#' @param area The cumulative area under the standard normal curve. For a left-tail test, this is the area to the left of the z-score.
#' @param tail Specify the tail of the test: "left", "right", or "two.sided".
#' @return A numeric value representing the z-score and an explanation.
#' @examples
#' z_score_for_area(area = 0.95, tail = "left") # 95% cumulative area to the left
#' z_score_for_area(area = 0.05, tail = "right") # 5% area to the right
#' z_score_for_area(area = 0.95, tail = "two.sided") # Two-sided 95% confidence level
#' @export
z_score_for_area <- function(area, tail = "left") {
  if (tail == "left") {
    # Z-Score for left-tail
    z <- qnorm(area)
    explanation <- paste("The z-score corresponding to a left-tail area of", area, "is:", round(z, 4))
  } else if (tail == "right") {
    # Z-Score for right-tail
    z <- qnorm(1 - area)
    explanation <- paste("The z-score corresponding to a right-tail area of", area, "is:", round(z, 4))
  } else if (tail == "two.sided") {
    # Z-Score for two-sided (divide area into two tails)
    alpha <- 1 - area
    z <- qnorm(1 - alpha / 2)
    explanation <- paste("The z-score corresponding to a two-sided area of", area,
                         "is +-", round(z, 4), "for a two-tailed test.")
  } else {
    stop("Invalid tail specification. Use 'left', 'right', or 'two.sided'.")
  }

  return(list(z_score = z, explanation = explanation))
}


#' @title Expected Value and Variance of Discrete Probability Distribution
#' @description Calculates the expected value and variance for a given discrete probability distribution.
#' @param outcomes A numeric vector representing the outcomes of the random variable.
#' @param probabilities A numeric vector representing the probabilities of each outcome.
#' @return A list containing the expected value, variance, and an explanation of the calculations.
#' @examples
#' outcomes <- c(1, 2, 3, 4, 5, 6)
#' probabilities <- c(0.1, 0.2, 0.3, 0.3, 0, 0.1)
#' discrete_expected_value_variance(outcomes, probabilities)
#' @export
discrete_expected_value_variance <- function(outcomes, probabilities) {
  # Validate input
  if (length(outcomes) != length(probabilities)) {
    stop("The 'outcomes' and 'probabilities' vectors must have the same length.")
  }
  if (abs(sum(probabilities) - 1) > 1e-6) {
    stop("The probabilities must sum to 1.")
  }

  # Calculate expected value
  expected_value <- sum(outcomes * probabilities)

  # Calculate variance
  expected_value_of_squares <- sum((outcomes^2) * probabilities)
  variance <- expected_value_of_squares - (expected_value^2)

  # Explanation
  explanation <- paste(
    "To calculate the expected value and variance of the discrete probability distribution:",
    "\n1. Expected value (E[X]) = Sum(outcome * probability)",
    "\n   =", paste0(outcomes, "*", probabilities, collapse = " + "), "=", round(expected_value, 4),
    "\n2. Variance (Var(X)) = Sum(outcome^2 * probability) - (E[X])^2",
    "\n   =", paste0(outcomes, "^2*", probabilities, collapse = " + "),
    "-", round(expected_value, 4), "^2 =", round(variance, 4)
  )

  return(list(expected_value = expected_value, variance = variance, explanation = explanation))
}

#' @title Test Selection for Comparing Means
#' @description Determines the most appropriate test for comparing the means of two independent samples
#' based on whether the population standard deviations are known.
#' @param population_sd_known Logical, whether population standard deviations are known (TRUE or FALSE).
#' @return A string indicating the most appropriate test and an explanation.
#' @examples
#' test_selection_for_means(population_sd_known = TRUE)
#' test_selection_for_means(population_sd_known = FALSE)
#' @export
test_selection_for_means <- function(population_sd_known) {
  if (population_sd_known) {
    test <- "Two-sample z-test for means"
    explanation <- paste(
      "Since the population standard deviations are known, a two-sample z-test for means is appropriate.",
      "This test compares the means of two independent groups assuming the standard deviations of the populations are known."
    )
  } else {
    test <- "Two-sample t-test for means"
    explanation <- paste(
      "Since the population standard deviations are unknown, a two-sample t-test for means is appropriate.",
      "This test compares the means of two independent groups using sample standard deviations."
    )
  }

  return(list(test = test, explanation = explanation))
}

#' @title Regression Coefficients Explanation
#' @description Explains the roles of \( b_0 \) (intercept) and \( b_1 \) (slope) in the regression equation \( y = b_0 + b_1x \).
#' @return A list containing descriptions for \( b_0 \) and \( b_1 \), and the key factor explaining the relationship.
#' @examples
#' regression_coefficients_explanation()
#' @export
regression_coefficients_explanation <- function() {
  b0_explanation <- paste(
    "b0 (intercept): Represents the value of y when x = 0.",
    "It indicates the starting point or baseline value of y in the absence of x."
  )

  b1_explanation <- paste(
    "b1 (slope): Represents the rate of change of y with respect to x.",
    "It quantifies how much y changes for a one-unit increase in x.",
    "A positive b1 indicates a positive relationship (y increases as x increases),",
    "while a negative b1 indicates a negative relationship (y decreases as x increases)."
  )

  key_relationship <- "b1 (slope) best explains the relationship between x and y because it directly measures the impact of x on y."

  return(list(
    b0_explanation = b0_explanation,
    b1_explanation = b1_explanation,
    key_relationship = key_relationship
  ))
}

#' @title Correlation Coefficient Explanation
#' @description Explains the meaning of the correlation coefficient ( r) when it is close to 1, -1, or 0.
#' @param r The correlation coefficient (numeric, between -1 and 1).
#' @return An explanation of the relationship between two variables based on ( r).
#' @examples
#' correlation_explanation(0.9)
#' correlation_explanation(0)
#' correlation_explanation(-0.8)
#' @export
correlation_explanation <- function(r) {
  if (r > 0.8 && r <= 1) {
    explanation <- paste("The correlation coefficient (r =", r, ") is close to 1, indicating a strong positive linear relationship.",
                         "As one variable increases, the other variable also increases.")
  } else if (r < -0.8 && r >= -1) {
    explanation <- paste("The correlation coefficient (r =", r, ") is close to -1, indicating a strong negative linear relationship.",
                         "As one variable increases, the other variable decreases.")
  } else if (r > -0.2 && r < 0.2) {
    explanation <- paste("The correlation coefficient (r =", r, ") is close to 0, indicating little or no linear relationship.",
                         "Changes in one variable are not linearly associated with changes in the other variable.")
  } else {
    explanation <- paste("The correlation coefficient (r =", r, ") indicates a moderate linear relationship.",
                         ifelse(r > 0, "This is a positive relationship.", "This is a negative relationship."))
  }

  return(explanation)
}



#' Solve Uniform Distribution Problems
#'
#' This function solves problems related to a uniform distribution, where the density curve is defined
#' between a lower bound (x = 1) and an upper bound (x = 8). It calculates the height of the uniform distribution curve,
#' the density function value at a given point, the percentage of observations in a specified range,
#' and the probability of specific values (including below, above, and equal).
#'
#' @param x A numeric value or a vector of numeric values specifying the bounds for the calculations.
#' @param lower_bound The lower bound of the uniform distribution (default is 1).
#' @param upper_bound The upper bound of the uniform distribution (default is 8).
#'
#' @return A list containing:
#'   - `height`: The height of the uniform distribution curve.
#'   - `density_function`: The value of the uniform density function at `x`.
#'   - `percent_in_range`: The percent of observations falling within the specified range.
#'   - `percent_below`: The percent of observations below the specified value.
#'   - `percent_above`: The percent of observations above the specified value.
#'   - `percent_equal`: The percent of observations equal to a specified value (always 0 for continuous distributions).
#'
#' @examples
#' uniform_distribution(2, 1, 8)
#' uniform_distribution(c(2, 5), 1, 8)  # Percent between 2 and 5
#' uniform_distribution(4, 1, 8)  # Percent below 4
#' uniform_distribution(6, 1, 8)  # Percent above 6
#' uniform_distribution(7, 1, 8)  # Percent equal to 7 (always 0)
#'
#' @export
uniform_distribution <- function(x, lower_bound, upper_bound) {
  # Calculate the height of the uniform distribution curve
  height <- 1 / (upper_bound - lower_bound)

  # Calculate density function value at x
  density_function <- ifelse(x >= lower_bound & x <= upper_bound, height, 0)

  # Calculate percent between two points if two values are provided
  if (length(x) == 2) {
    percent_in_range <- (x[2] - x[1]) * height * 100
  } else {
    percent_in_range <- NA
  }

  # Calculate percent below x
  if (!is.null(x)) {
    percent_below <- (x - lower_bound) * height * 100
  }

  # Calculate percent above x
  percent_above <- (upper_bound - x) * height * 100

  # Percent equal to x (always 0 for continuous distributions)
  percent_equal <- 0

  return(list(height = height,
              density_function = density_function,
              percent_in_range = percent_in_range,
              percent_below = percent_below,
              percent_above = percent_above,
              percent_equal = percent_equal))
}

#' @title Generic Hypothesis Test for Population Mean
#' @description Performs a hypothesis test to determine if the population mean differs from a specified value.
#' @param data A numeric vector containing the sample data.
#' @param mu The hypothesized population mean to test against.
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the test type, test statistic, p-value, confidence interval, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic, p-value, and confidence interval for a two-tailed hypothesis test on the population mean.
#' It assumes the data are normally distributed and uses a t-test when the population standard deviation is unknown.
#' @examples
#' data <- c(34,54,73,38,89,52,75,33,50,39,42,42,40,66,72,85,28,71,52,47,41,36,33,38,49,51,55,63,72,78)
#' hypothesis_test_mean(data = data, mu = 50, alpha = 0.05)
#' @export
hypothesis_test_mean <- function(data, mu, alpha = 0.05) {
  # Step (a): Determine the type of test
  test_type <- "Two-tailed test, since Ha: mu != hypothesized mean"

  # Step (b): Calculate sample statistics
  n <- length(data)
  sample_mean <- mean(data)
  sample_sd <- sd(data)
  standard_error <- sample_sd / sqrt(n)

  # Step (c): Calculate the test statistic
  t_stat <- (sample_mean - mu) / standard_error
  test_statistic_details <- paste(
    "Calculate t = (sample_mean - mu) / SE:",
    "sample_mean =", round(sample_mean, 4),
    ", mu =", mu,
    ", SE =", round(standard_error, 4),
    ", t =", round(t_stat, 4)
  )

  # Step (d): Determine the rejection region
  t_critical <- qt(1 - alpha / 2, df = n - 1)
  rejection_region <- paste("Reject H0 if t < -", round(t_critical, 4), " or t > ", round(t_critical, 4))

  # Step (e): Obtain the p-value
  p_value <- 2 * (1 - pt(abs(t_stat), df = n - 1))  # two-tailed p-value
  p_value_details <- paste("Two-tailed p-value for t =", round(t_stat, 4), "=> p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if (abs(t_stat) > t_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (abs(t_stat) > t_critical) {
    "There is significant evidence that the population mean differs from the hypothesized mean."
  } else {
    "There is not enough evidence to conclude that the population mean differs from the hypothesized mean."
  }

  # Step (g): Calculate the confidence interval
  margin_of_error <- t_critical * standard_error
  ci_lower <- sample_mean - margin_of_error
  ci_upper <- sample_mean + margin_of_error
  confidence_interval <- paste("Confidence interval: (", round(ci_lower, 4), ", ", round(ci_upper, 4), ")")

  # Return all details in the output
  return(list(
    test_type = test_type,
    test_statistic_details = test_statistic_details,
    rejection_region = rejection_region,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion,
    confidence_interval = confidence_interval
  ))
}

#' @title Hypothesis Test for Population Proportion
#' @description Performs a hypothesis test to determine if the population proportion differs from a specified value.
#' @param n Sample size, the number of patients or cases.
#' @param x The number of patients or cases with the outcome of interest (e.g., adverse symptoms).
#' @param p0 The hypothesized population proportion (e.g., 0.10 for 10 percent).
#' @param alpha The significance level for the test (default is 0.05).
#' @return A list containing the sampling distribution, test statistic, p-value, rejection region, decision, and conclusion, with all work shown.
#' @details
#' This function checks assumptions, calculates the test statistic, p-value, and provides a one-tailed or two-tailed hypothesis test on the population proportion.
#' @examples
#' hypothesis_test_proportion(n = 440, x = 23, p0 = 0.10, alpha = 0.05)
#' @export
hypothesis_test_proportion <- function(n, x, p0, alpha = 0.05) {
  # Step (a): Sampling distribution of the sample proportion
  sample_proportion <- x / n
  standard_error <- sqrt((p0 * (1 - p0)) / n)
  sampling_distribution <- paste("Sampling distribution: N(mean =", round(p0, 4),
                                 ", standard deviation =", round(standard_error, 4), ")")

  # Step (b1): Assumption check
  np0 <- n * p0
  n1_minus_p0 <- n * (1 - p0)
  assumption_check <- np0 >= 10 && n1_minus_p0 >= 10
  assumption_details <- paste(
    "Check np0 >= 10 and n(1 - p0) >= 10:",
    "np0 =", np0,
    ", n(1 - p0) =", n1_minus_p0,
    "=> Normal approximation is valid:", assumption_check
  )

  # Step (b2): State the hypotheses
  hypotheses <- list(
    H0 = paste("p =", p0, "(The proportion of patients with adverse symptoms is", p0, ")"),
    Ha = paste("p <", p0, "(The proportion of patients with adverse symptoms is less than", p0, ")")
  )

  # Step (c): Determine the rejection region for a left-tailed test
  z_critical <- qnorm(alpha)
  rejection_region <- paste("Reject H0 if z <", round(z_critical, 4))

  # Step (d): Calculate the test statistic
  z <- (sample_proportion - p0) / standard_error
  test_statistic_details <- paste(
    "Calculate z = (sample_proportion - p0) / SE:",
    "sample_proportion =", round(sample_proportion, 4),
    ", p0 =", p0,
    ", SE =", round(standard_error, 4),
    ", z =", round(z, 4)
  )

  # Step (e): Obtain the p-value for a one-tailed test
  p_value <- pnorm(z)
  p_value_details <- paste("One-tailed p-value for z =", round(z, 4), "p-value =", round(p_value, 4))

  # Step (f): Decision and Conclusion
  decision <- if (z < z_critical) "Reject H0" else "Do not reject H0"
  conclusion <- if (z < z_critical) {
    "There is strong evidence that fewer than 10 percent of patients taking this medication have adverse symptoms."
  } else {
    "There is not enough evidence to conclude that fewer than 10 percent of patients taking this medication have adverse symptoms."
  }

  # Return all details in the output
  return(list(
    sampling_distribution = sampling_distribution,
    assumption_details = assumption_details,
    hypotheses = hypotheses,
    rejection_region = rejection_region,
    test_statistic_details = test_statistic_details,
    p_value_details = p_value_details,
    decision = decision,
    conclusion = conclusion
  ))
}

#' @title Hypothesis Test for Variance
#' @description Performs a hypothesis test to determine if the population variance differs from a specified value.
#' @param data A numeric vector of sample data.
#' @param variance_0 The hypothesized population variance.
#' @param alpha The significance level for the test (default is 0.05).
#' @param alternative The alternative hypothesis. Use "two.sided", "greater", or "less".
#' @importFrom stats pchisq
#' @return A list containing the test statistic, degrees of freedom, p-value, rejection region, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic and p-value for a hypothesis test on the population variance.
#' It uses the chi-square distribution and assumes the data are normally distributed.
#' @examples
#' data <- c(34, 54, 73, 38, 89, 52, 75, 33, 50, 39)
#' hypothesis_test_variance(data, variance_0 = 25, alpha = 0.05, alternative = "two.sided")
#' @export
hypothesis_test_variance <- function(data, variance_0, alpha = 0.05, alternative = "two.sided") {
  n <- length(data)
  sample_variance <- var(data)
  test_statistic <- (n - 1) * sample_variance / variance_0
  df <- n - 1

  if (alternative == "two.sided") {
    chi_critical_low <- qchisq(alpha / 2, df)
    chi_critical_high <- qchisq(1 - alpha / 2, df)
    p_value <- 2 * min(pchisq(test_statistic, df), 1 - pchisq(test_statistic, df))
    rejection_region <- paste("Reject H0 if test statistic < ", round(chi_critical_low, 4),
                              " or > ", round(chi_critical_high, 4))
  } else if (alternative == "greater") {
    chi_critical <- qchisq(1 - alpha, df)
    p_value <- 1 - pchisq(test_statistic, df)
    rejection_region <- paste("Reject H0 if test statistic > ", round(chi_critical, 4))
  } else if (alternative == "less") {
    chi_critical <- qchisq(alpha, df)
    p_value <- pchisq(test_statistic, df)
    rejection_region <- paste("Reject H0 if test statistic < ", round(chi_critical, 4))
  } else {
    stop("Invalid alternative hypothesis. Use 'two.sided', 'greater', or 'less'.")
  }

  decision <- if ((alternative == "two.sided" && (test_statistic < chi_critical_low || test_statistic > chi_critical_high)) ||
                  (alternative == "greater" && test_statistic > chi_critical) ||
                  (alternative == "less" && test_statistic < chi_critical)) {
    "Reject H0"
  } else {
    "Do not reject H0"
  }

  conclusion <- if (decision == "Reject H0") {
    paste("There is significant evidence that the population variance",
          ifelse(alternative == "two.sided", "differs from", ifelse(alternative == "greater", "is greater than", "is less than")),
          variance_0, ".")
  } else {
    paste("There is not enough evidence to conclude that the population variance",
          ifelse(alternative == "two.sided", "differs from", ifelse(alternative == "greater", "is greater than", "is less than")),
          variance_0, ".")
  }

  return(list(
    test_statistic = round(test_statistic, 4),
    degrees_of_freedom = df,
    p_value = round(p_value, 4),
    rejection_region = rejection_region,
    decision = decision,
    conclusion = conclusion
  ))
}

#' @title Paired t-Test with Detailed Steps
#' @description Performs a paired t-test to compare the means of two related groups and displays detailed work.
#' @param data1 A numeric vector containing the first sample data.
#' @param data2 A numeric vector containing the second sample data.
#' @param alpha The significance level for the test (default is 0.05).
#' @return A detailed list showing the test statistic, degrees of freedom, p-value, confidence interval, decision, and work explanation.
#' @examples
#' paired_t_test(c(12, 15, 14), c(10, 13, 11), alpha = 0.05)
#' @export
paired_t_test <- function(data1, data2, alpha = 0.05) {
  # Calculate the differences
  differences <- data1 - data2
  n <- length(differences)
  mean_diff <- mean(differences)
  sd_diff <- sd(differences)
  se_diff <- sd_diff / sqrt(n)
  t_stat <- mean_diff / se_diff
  df <- n - 1
  p_value <- 2 * (1 - pt(abs(t_stat), df))
  t_critical <- qt(1 - alpha / 2, df)
  margin_of_error <- t_critical * se_diff
  ci <- c(mean_diff - margin_of_error, mean_diff + margin_of_error)
  decision <- ifelse(abs(t_stat) > t_critical, "Reject the null hypothesis", "Do not reject the null hypothesis")

  # Step-by-step explanation
  explanation <- paste(
    "Step 1: State the hypotheses.",
    "   Null hypothesis (H0): The mean difference is zero, mean_diff = 0.",
    "   Alternative hypothesis (Ha): The mean difference is not zero, mean_diff != 0.",
    "",
    "Step 2: Calculate the sample statistics.",
    "   Mean of the differences (mean_diff):", round(mean_diff, 4),
    "   Standard deviation of the differences (sd_diff):", round(sd_diff, 4),
    "   Sample size (n):", n,
    "",
    "Step 3: Calculate the standard error.",
    "   Standard error (SE): sd_diff / sqrt(n) =", round(sd_diff, 4), "/ sqrt(", n, ") =", round(se_diff, 4),
    "",
    "Step 4: Calculate the test statistic.",
    "   Test statistic (t): mean_diff / SE = (", round(mean_diff, 4), ") / (", round(se_diff, 4), ") =", round(t_stat, 4),
    "",
    "Step 5: Determine the rejection region.",
    "   Critical t-value for a two-tailed test:", round(t_critical, 4),
    "   Rejection region: Reject the null hypothesis if |t| >", round(t_critical, 4),
    "",
    "Step 6: Calculate the p-value.",
    "   P-value = 2 * P(T > |t|) =", round(p_value, 4),
    "",
    "Step 7: Calculate the confidence interval.",
    "   Confidence interval = mean_diff +/- t_critical * SE.",
    "   CI = (", round(ci[1], 4), ",", round(ci[2], 4), ")",
    "",
    "Step 8: Conclusion.",
    ifelse(abs(t_stat) > t_critical,
           "   Reject the null hypothesis. There is significant evidence that the mean difference is not zero.",
           "   Do not reject the null hypothesis. There is not enough evidence to conclude that the mean difference is not zero.")
  )

  # Return results
  return(list(
    test_statistic = round(t_stat, 4),
    degrees_of_freedom = df,
    p_value = round(p_value, 4),
    confidence_interval = round(ci, 4),
    conclusion = decision,
    explanation = explanation
  ))
}

#' @title One-Way ANOVA Test with Detailed Steps
#' @description Performs a one-way ANOVA test to compare the means of three or more groups and displays detailed work.
#' @param ... Numeric vectors representing the groups. Each group should be passed as a separate argument.
#' @param alpha The significance level for the test (default is 0.05).
#' @return A detailed list showing the F-statistic, degrees of freedom, p-value, and step-by-step explanation.
#' @examples
#' group1 <- c(12, 15, 14, 10, 13)
#' group2 <- c(18, 20, 17, 19, 16)
#' group3 <- c(22, 25, 24, 23, 26)
#' one_way_anova(group1, group2, group3, alpha = 0.05)
#' @export
one_way_anova <- function(..., alpha = 0.05) {
  # Combine groups into a list
  groups <- list(...)
  k <- length(groups)  # Number of groups
  n <- sum(sapply(groups, length))  # Total number of observations

  # Calculate group means and overall mean
  group_means <- sapply(groups, mean)
  overall_mean <- mean(unlist(groups))

  # Calculate between-group sum of squares (SSB)
  ssb <- sum(sapply(groups, function(g) length(g) * (mean(g) - overall_mean)^2))

  # Calculate within-group sum of squares (SSW)
  ssw <- sum(sapply(groups, function(g) sum((g - mean(g))^2)))

  # Degrees of freedom
  df_between <- k - 1
  df_within <- n - k

  # Mean squares
  ms_between <- ssb / df_between
  ms_within <- ssw / df_within

  # F-statistic
  f_stat <- ms_between / ms_within

  # P-value
  p_value <- pf(f_stat, df_between, df_within, lower.tail = FALSE)

  # Critical F-value
  f_critical <- qf(1 - alpha, df_between, df_within)

  # Decision
  decision <- ifelse(f_stat > f_critical, "Reject the null hypothesis", "Do not reject the null hypothesis")

  # Explanation
  explanation <- paste(
    "Step 1: State the hypotheses.",
    "   Null hypothesis (H0): All group means are equal.",
    "   Alternative hypothesis (Ha): At least one group mean is different.",
    "",
    "Step 2: Compute the sample statistics.",
    "   Overall mean:", round(overall_mean, 4),
    "   Group means:", paste(round(group_means, 4), collapse = ", "),
    "",
    "Step 3: Calculate the sum of squares.",
    "   Between-group sum of squares (SSB):", round(ssb, 4),
    "   Within-group sum of squares (SSW):", round(ssw, 4),
    "",
    "Step 4: Compute the degrees of freedom.",
    "   Degrees of freedom between groups:", df_between,
    "   Degrees of freedom within groups:", df_within,
    "",
    "Step 5: Compute the mean squares.",
    "   Mean square between groups (MSB): SSB / df_between =", round(ms_between, 4),
    "   Mean square within groups (MSW): SSW / df_within =", round(ms_within, 4),
    "",
    "Step 6: Calculate the F-statistic.",
    "   F = MSB / MSW =", round(f_stat, 4),
    "",
    "Step 7: Determine the rejection region.",
    "   Critical F-value for alpha =", alpha, ":", round(f_critical, 4),
    "   Rejection region: Reject H0 if F >", round(f_critical, 4),
    "",
    "Step 8: Calculate the p-value.",
    "   P-value =", round(p_value, 4),
    "",
    "Step 9: Conclusion.",
    ifelse(f_stat > f_critical,
           "   Reject the null hypothesis. There is significant evidence that at least one group mean is different.",
           "   Do not reject the null hypothesis. There is not enough evidence to conclude that the group means are different.")
  )

  # Return results
  return(list(
    f_statistic = round(f_stat, 4),
    degrees_of_freedom_between = df_between,
    degrees_of_freedom_within = df_within,
    p_value = round(p_value, 4),
    decision = decision,
    explanation = explanation
  ))
}



#' @title One-Way ANOVA Test with Detailed Steps
#' @description Performs a one-way ANOVA test to compare the means of three or more groups and displays detailed work.
#' @param ... Numeric vectors representing the groups. Each group should be passed as a separate argument.
#' @param alpha The significance level for the test (default is 0.05).
#' @importFrom stats pf qf
#' @return A detailed list showing the F-statistic, degrees of freedom, p-value, and step-by-step explanation.
#' @examples
#' group1 <- c(12, 15, 14, 10, 13)
#' group2 <- c(18, 20, 17, 19, 16)
#' group3 <- c(22, 25, 24, 23, 26)
#' one_way_anova(group1, group2, group3, alpha = 0.05)
#' @export
one_way_anova <- function(..., alpha = 0.05) {
  # Combine groups into a list
  groups <- list(...)
  k <- length(groups)  # Number of groups
  n <- sum(sapply(groups, length))  # Total number of observations

  # Calculate group means and overall mean
  group_means <- sapply(groups, mean)
  overall_mean <- mean(unlist(groups))

  # Calculate between-group sum of squares (SSB)
  ssb <- sum(sapply(groups, function(g) length(g) * (mean(g) - overall_mean)^2))

  # Calculate within-group sum of squares (SSW)
  ssw <- sum(sapply(groups, function(g) sum((g - mean(g))^2)))

  # Degrees of freedom
  df_between <- k - 1
  df_within <- n - k

  # Mean squares
  ms_between <- ssb / df_between
  ms_within <- ssw / df_within

  # F-statistic
  f_stat <- ms_between / ms_within

  # P-value
  p_value <- pf(f_stat, df_between, df_within, lower.tail = FALSE)

  # Critical F-value
  f_critical <- qf(1 - alpha, df_between, df_within)

  # Decision
  decision <- ifelse(f_stat > f_critical, "Reject the null hypothesis", "Do not reject the null hypothesis")

  # Explanation
  explanation <- paste(
    "Step 1: State the hypotheses.",
    "   Null hypothesis (H0): All group means are equal.",
    "   Alternative hypothesis (Ha): At least one group mean is different.",
    "",
    "Step 2: Compute the sample statistics.",
    "   Overall mean:", round(overall_mean, 4),
    "   Group means:", paste(round(group_means, 4), collapse = ", "),
    "",
    "Step 3: Calculate the sum of squares.",
    "   Between-group sum of squares (SSB):", round(ssb, 4),
    "   Within-group sum of squares (SSW):", round(ssw, 4),
    "",
    "Step 4: Compute the degrees of freedom.",
    "   Degrees of freedom between groups:", df_between,
    "   Degrees of freedom within groups:", df_within,
    "",
    "Step 5: Compute the mean squares.",
    "   Mean square between groups (MSB): SSB / df_between =", round(ms_between, 4),
    "   Mean square within groups (MSW): SSW / df_within =", round(ms_within, 4),
    "",
    "Step 6: Calculate the F-statistic.",
    "   F = MSB / MSW =", round(f_stat, 4),
    "",
    "Step 7: Determine the rejection region.",
    "   Critical F-value for alpha =", alpha, ":", round(f_critical, 4),
    "   Rejection region: Reject H0 if F >", round(f_critical, 4),
    "",
    "Step 8: Calculate the p-value.",
    "   P-value =", round(p_value, 4),
    "",
    "Step 9: Conclusion.",
    ifelse(f_stat > f_critical,
           "   Reject the null hypothesis. There is significant evidence that at least one group mean is different.",
           "   Do not reject the null hypothesis. There is not enough evidence to conclude that the group means are different.")
  )

  # Return results
  return(list(
    f_statistic = round(f_stat, 4),
    degrees_of_freedom_between = df_between,
    degrees_of_freedom_within = df_within,
    p_value = round(p_value, 4),
    decision = decision,
    explanation = explanation
  ))
}

#' @title Proportion Test with Detailed Steps
#' @description Performs a test for proportions and displays detailed work.
#' @param x The number of successes in the sample (or a vector of successes for two samples).
#' @param n The sample size (or a vector of sample sizes for two samples).
#' @param p Null hypothesis proportion (for one sample) or NULL for two-sample test (default is NULL).
#' @param alternative The alternative hypothesis ("two.sided", "greater", or "less").
#' @param alpha The significance level for the test (default is 0.05).
#' @importFrom stats prop.test
#' @return A list containing the test statistic, p-value, confidence interval, decision, and step-by-step explanation.
#' @examples
#' # One-sample test
#' prop_test(x = 45, n = 100, p = 0.5, alternative = "two.sided", alpha = 0.05)
#' # Two-sample test
#' prop_test(x = c(30, 40), n = c(100, 120), alternative = "two.sided", alpha = 0.05)
#' @export
prop_test <- function(x, n, p = NULL, alternative = "two.sided", alpha = 0.05) {
  if (is.null(p)) {
    # Two-sample test
    p_hat1 <- x[1] / n[1]
    p_hat2 <- x[2] / n[2]
    pooled_p <- sum(x) / sum(n)
    se <- sqrt(pooled_p * (1 - pooled_p) * (1 / n[1] + 1 / n[2]))
    z_stat <- (p_hat1 - p_hat2) / se
    df <- NULL
    confidence <- paste0(100 * (1 - alpha), "% confidence interval")
    ci <- prop.test(x, n, alternative = alternative, correct = FALSE)$conf.int
  } else {
    # One-sample test
    p_hat <- x / n
    se <- sqrt(p * (1 - p) / n)
    z_stat <- (p_hat - p) / se
    pooled_p <- NULL
    df <- NULL
    confidence <- paste0(100 * (1 - alpha), "% confidence interval")
    ci <- prop.test(x, n, p = p, alternative = alternative, correct = FALSE)$conf.int
  }

  p_value <- 2 * (1 - pnorm(abs(z_stat)))
  z_critical <- qnorm(1 - alpha / 2)

  decision <- ifelse(abs(z_stat) > z_critical, "Reject H0", "Do not reject H0")

  explanation <- paste(
    "Step 1: State the hypotheses.",
    if (!is.null(p)) {
      paste("   Null hypothesis (H0): p =", p, " (The population proportion equals the hypothesized value).")
    } else {
      "   Null hypothesis (H0): p1 = p2 (The proportions are equal in both groups)."
    },
    "   Alternative hypothesis (Ha):",
    if (alternative == "two.sided") {
      "The proportion(s) are different."
    } else if (alternative == "greater") {
      "The proportion(s) are greater than the hypothesized value."
    } else {
      "The proportion(s) are less than the hypothesized value."
    },
    "",
    "Step 2: Calculate sample statistics.",
    if (!is.null(p)) {
      paste("   Sample proportion (p_hat):", round(p_hat, 4))
    } else {
      paste("   Sample proportion 1 (p_hat1):", round(p_hat1, 4),
            "   Sample proportion 2 (p_hat2):", round(p_hat2, 4))
    },
    "   Standard error (SE):", round(se, 4),
    "",
    "Step 3: Calculate the test statistic.",
    "   Test statistic (z):", round(z_stat, 4),
    "",
    "Step 4: Determine the rejection region.",
    "   Critical z-value for alpha =", alpha, ":", round(z_critical, 4),
    "   Rejection region: Reject H0 if |z| >", round(z_critical, 4),
    "",
    "Step 5: Calculate the p-value.",
    "   P-value =", round(p_value, 4),
    "",
    "Step 6: Calculate the confidence interval.",
    paste("   ", confidence, ":", round(ci[1], 4), "to", round(ci[2], 4)),
    "",
    "Step 7: Conclusion.",
    ifelse(abs(z_stat) > z_critical,
           "   Reject the null hypothesis. There is significant evidence to support the alternative hypothesis.",
           "   Do not reject the null hypothesis. There is not enough evidence to support the alternative hypothesis.")
  )

  return(list(
    test_statistic = round(z_stat, 4),
    p_value = round(p_value, 4),
    confidence_interval = round(ci, 4),
    decision = decision,
    explanation = explanation
  ))
}

#' @title Fit Linear Regression Model
#' @description Fits a simple or multiple linear regression model and returns a summary along with detailed steps.
#' @param formula A formula object specifying the model (e.g., y ~ x1 + x2).
#' @param data A data frame containing the variables in the model.
#' @importFrom stats lm predict
#' @return A list containing the model, coefficients, residuals, R-squared, adjusted R-squared, and detailed explanation.
#' @examples
#' # Simple linear regression
#' fit_linear_model(mpg ~ wt, data = mtcars)
#' # Multiple linear regression
#' fit_linear_model(mpg ~ wt + hp, data = mtcars)
#' @export
fit_linear_model <- function(formula, data) {
  # Fit the model
  model <- lm(formula, data = data)
  summary_model <- summary(model)

  # Extract key components
  coefficients <- summary_model$coefficients
  residuals <- model$residuals
  r_squared <- summary_model$r.squared
  adj_r_squared <- summary_model$adj.r.squared

  # Explanation
  explanation <- paste(
    "Step 1: Fit the model using linear regression.",
    "   Model formula:", deparse(formula),
    "",
    "Step 2: Calculate the coefficients.",
    paste("   Coefficients:", paste(names(coefficients[, 1]), round(coefficients[, 1], 4), sep = " = ", collapse = ", ")),
    "",
    "Step 3: Evaluate the goodness of fit.",
    "   R-squared:", round(r_squared, 4),
    "   Adjusted R-squared:", round(adj_r_squared, 4),
    "",
    "Step 4: Analyze residuals.",
    "   Residuals indicate the differences between observed and predicted values.",
    "   Summary of residuals:",
    paste("   Min =", round(min(residuals), 4), ", Max =", round(max(residuals), 4), ", Mean =", round(mean(residuals), 4))
  )

  # Return results
  return(list(
    model = model,
    coefficients = coefficients,
    residuals = residuals,
    r_squared = round(r_squared, 4),
    adjusted_r_squared = round(adj_r_squared, 4),
    explanation = explanation
  ))
}

#' @title Evaluate Linear Model Predictions
#' @description Evaluates predictions from a linear regression model and provides metrics such as RMSE and MAE.
#' @param model A fitted linear model (from lm).
#' @param newdata A data frame containing new data for prediction.
#' @param actual A vector of actual values corresponding to the predictions.
#' @importFrom stats lm predict
#' @return A list containing predictions, RMSE, MAE, and detailed explanation.
#' @examples
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' linear_model_predictions(model, newdata = mtcars, actual = mtcars$mpg)
#' @export
linear_model_predictions <- function(model, newdata, actual) {
  # Generate predictions
  predictions <- predict(model, newdata)

  # Calculate evaluation metrics
  residuals <- actual - predictions
  rmse <- sqrt(mean(residuals^2))
  mae <- mean(abs(residuals))

  # Explanation
  explanation <- paste(
    "Step 1: Generate predictions using the fitted model.",
    "   Predictions are calculated using the model and new data.",
    "",
    "Step 2: Calculate residuals.",
    "   Residuals = Actual - Predicted.",
    "",
    "Step 3: Evaluate the prediction performance.",
    "   Root Mean Square Error (RMSE):", round(rmse, 4),
    "   Mean Absolute Error (MAE):", round(mae, 4),
    "",
    "Step 4: Interpretation.",
    "   Lower RMSE and MAE values indicate better predictive performance."
  )

  # Return results
  return(list(
    predictions = predictions,
    rmse = round(rmse, 4),
    mae = round(mae, 4),
    explanation = explanation
  ))
}

#' @title Linear Regression Explanation
#' @description Provides a step-by-step explanation of how linear regression works, using the fitted model.
#' @param model A fitted linear model (from lm).
#' @return A detailed explanation of the linear regression process and its results.
#' @examples
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' explain_linear_model(model)
#' @export
explain_linear_model <- function(model) {
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  r_squared <- summary_model$r.squared
  adj_r_squared <- summary_model$adj.r.squared

  explanation <- paste(
    "Step 1: Understand the model structure.",
    "   The model predicts the dependent variable based on the independent variables.",
    "",
    "Step 2: Examine the coefficients.",
    paste("   Intercept:", round(coefficients["(Intercept)", 1], 4)),
    paste("   Slope coefficients:", paste(names(coefficients[-1, 1]), round(coefficients[-1, 1], 4), collapse = ", ")),
    "",
    "Step 3: Evaluate the model fit.",
    "   R-squared:", round(r_squared, 4), "indicates the proportion of variance explained by the model.",
    "   Adjusted R-squared:", round(adj_r_squared, 4), "accounts for the number of predictors in the model.",
    "",
    "Step 4: Analyze residuals.",
    "   Residuals are the differences between observed and predicted values.",
    "   A good model will have residuals that are randomly distributed with no clear patterns."
  )

  return(explanation)
}

#' @title Generic Hypothesis Test with Rejection Region
#' @description Performs a hypothesis test for a population mean or other parameter with a specified rejection region.
#' @param n Sample size, the number of observations.
#' @param sample_stat The sample statistic (e.g., mean or proportion).
#' @param hypothesized_value The hypothesized population parameter (default is 0).
#' @param standard_error The standard error of the sample statistic (calculated outside the function if needed).
#' @param alpha The significance level for the test (default is 0.05).
#' @param test_type The type of test ("two-tailed", "greater", "less").
#' @return A list containing the test type, rejection region, test statistic, p-value, decision, and conclusion, with all work shown.
#' @details
#' This function calculates the test statistic, rejection region, and p-value based on the test type.
#' It provides step-by-step output, including all relevant calculations and decision criteria.
#' @examples
#' # Example: Two-tailed test for a mean
#' generic_hypothesis_test(n = 49, sample_stat = 415.7, hypothesized_value = 420,
#'                         standard_error = 10 / sqrt(49), alpha = 0.01, test_type = "two-tailed")
#' @export
generic_hypothesis_test <- function(n, sample_stat, hypothesized_value = 0, standard_error, alpha = 0.05, test_type = "two-tailed") {
  # Step (a): Determine the type of test
  if (test_type == "two-tailed") {
    critical_value <- qnorm(1 - alpha / 2)
    rejection_region <- paste("Reject H0 if z < -", round(critical_value, 4), " or z > ", round(critical_value, 4))
    p_value <- 2 * (1 - pnorm(abs((sample_stat - hypothesized_value) / standard_error)))
  } else if (test_type == "greater") {
    critical_value <- qnorm(1 - alpha)
    rejection_region <- paste("Reject H0 if z >", round(critical_value, 4))
    p_value <- 1 - pnorm((sample_stat - hypothesized_value) / standard_error)
  } else if (test_type == "less") {
    critical_value <- qnorm(alpha)
    rejection_region <- paste("Reject H0 if z <", round(critical_value, 4))
    p_value <- pnorm((sample_stat - hypothesized_value) / standard_error)
  } else {
    stop("Invalid test_type. Choose 'two-tailed', 'greater', or 'less'.")
  }

  # Step (b): Calculate the test statistic
  z <- (sample_stat - hypothesized_value) / standard_error
  test_statistic_details <- paste(
    "Calculate z = (sample_stat - hypothesized_value) / SE:",
    "sample_stat =", sample_stat,
    ", hypothesized_value =", hypothesized_value,
    ", SE =", round(standard_error, 4),
    ", z =", round(z, 4)
  )

  # Step (c): Decision based on rejection region
  decision <- if ((test_type == "two-tailed" && abs(z) > critical_value) ||
                  (test_type == "greater" && z > critical_value) ||
                  (test_type == "less" && z < critical_value)) {
    "Reject H0"
  } else {
    "Do not reject H0"
  }

  # Step (d): Conclusion
  conclusion <- if ((test_type == "two-tailed" && abs(z) > critical_value) ||
                    (test_type == "greater" && z > critical_value) ||
                    (test_type == "less" && z < critical_value)) {
    paste("There is significant evidence to reject the null hypothesis in favor of the alternative hypothesis.")
  } else {
    paste("There is not enough evidence to reject the null hypothesis.")
  }

  # Return all details
  return(list(
    test_type = test_type,
    rejection_region = rejection_region,
    test_statistic_details = test_statistic_details,
    p_value_details = paste("P-value =", round(p_value, 4)),
    decision = decision,
    conclusion = conclusion
  ))
}


#' @title Generic Standard Error Calculator
#' @description Calculates the standard error for different scenarios: single mean, single proportion, difference in means, and difference in proportions.
#' @param type The type of calculation ("mean", "proportion", "difference_means", "difference_proportions").
#' @param sd Standard deviation for a single mean (required for "mean" and "difference_means").
#' @param n Sample size for a single mean or proportion.
#' @param p Proportion (required for "proportion" and "difference_proportions").
#' @param sd1 Standard deviation of the first group (required for "difference_means").
#' @param n1 Sample size of the first group (required for "difference_means" and "difference_proportions").
#' @param sd2 Standard deviation of the second group (required for "difference_means").
#' @param n2 Sample size of the second group (required for "difference_means" and "difference_proportions").
#' @param p1 Proportion of the first group (required for "difference_proportions").
#' @param p2 Proportion of the second group (required for "difference_proportions").
#' @return The calculated standard error.
#' @examples
#' # Standard error for a mean
#' calculate_standard_error(type = "mean", sd = 10, n = 50)
#' # Standard error for a proportion
#' calculate_standard_error(type = "proportion", p = 0.5, n = 100)
#' # Standard error for the difference in means
#' calculate_standard_error(type = "difference_means", sd1 = 15, n1 = 30, sd2 = 20, n2 = 40)
#' # Standard error for the difference in proportions
#' calculate_standard_error(type = "difference_proportions", p1 = 0.6, n1 = 50, p2 = 0.5, n2 = 60)
#' @export
calculate_standard_error <- function(type, sd = NULL, n = NULL, p = NULL,
                                     sd1 = NULL, n1 = NULL, sd2 = NULL, n2 = NULL,
                                     p1 = NULL, p2 = NULL) {
  if (type == "mean") {
    # Standard error for a single mean
    if (is.null(sd) || is.null(n)) stop("Provide 'sd' and 'n' for type 'mean'.")
    se <- sd / sqrt(n)
    explanation <- paste(
      "Standard error for mean = sd / sqrt(n):",
      "sd =", sd, ", n =", n, "=> SE =", round(se, 4)
    )
  } else if (type == "proportion") {
    # Standard error for a single proportion
    if (is.null(p) || is.null(n)) stop("Provide 'p' and 'n' for type 'proportion'.")
    se <- sqrt((p * (1 - p)) / n)
    explanation <- paste(
      "Standard error for proportion = sqrt(p * (1 - p) / n):",
      "p =", p, ", n =", n, "=> SE =", round(se, 4)
    )
  } else if (type == "difference_means") {
    # Standard error for the difference in means
    if (is.null(sd1) || is.null(n1) || is.null(sd2) || is.null(n2)) stop("Provide 'sd1', 'n1', 'sd2', and 'n2' for type 'difference_means'.")
    se <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
    explanation <- paste(
      "Standard error for difference in means = sqrt((sd1^2 / n1) + (sd2^2 / n2)):",
      "sd1 =", sd1, ", n1 =", n1, ", sd2 =", sd2, ", n2 =", n2, "=> SE =", round(se, 4)
    )
  } else if (type == "difference_proportions") {
    # Standard error for the difference in proportions
    if (is.null(p1) || is.null(n1) || is.null(p2) || is.null(n2)) stop("Provide 'p1', 'n1', 'p2', and 'n2' for type 'difference_proportions'.")
    se <- sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2))
    explanation <- paste(
      "Standard error for difference in proportions = sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2)):",
      "p1 =", p1, ", n1 =", n1, ", p2 =", p2, ", n2 =", n2, "=> SE =", round(se, 4)
    )
  } else {
    stop("Invalid 'type'. Choose from 'mean', 'proportion', 'difference_means', or 'difference_proportions'.")
  }

  return(list(
    standard_error = round(se, 4),
    explanation = explanation
  ))
}

#' @title Generic Standard Error Calculator
#' @description Calculates the standard error for different scenarios: single mean, single proportion, difference in means, and difference in proportions.
#' @param type The type of calculation ("mean", "proportion", "difference_means", "difference_proportions").
#' @param sd Standard deviation for a single mean (required for "mean" and "difference_means").
#' @param n Sample size for a single mean or proportion.
#' @param p Proportion (required for "proportion" and "difference_proportions").
#' @param sd1 Standard deviation of the first group (required for "difference_means").
#' @param n1 Sample size of the first group (required for "difference_means" and "difference_proportions").
#' @param sd2 Standard deviation of the second group (required for "difference_means").
#' @param n2 Sample size of the second group (required for "difference_means" and "difference_proportions").
#' @param p1 Proportion of the first group (required for "difference_proportions").
#' @param p2 Proportion of the second group (required for "difference_proportions").
#' @return The calculated standard error.
#' @examples
#' # Standard error for a mean
#' calculate_standard_error(type = "mean", sd = 10, n = 50)
#' # Standard error for a proportion
#' calculate_standard_error(type = "proportion", p = 0.5, n = 100)
#' # Standard error for the difference in means
#' calculate_standard_error(type = "difference_means", sd1 = 15, n1 = 30, sd2 = 20, n2 = 40)
#' # Standard error for the difference in proportions
#' calculate_standard_error(type = "difference_proportions", p1 = 0.6, n1 = 50, p2 = 0.5, n2 = 60)
#' @export
calculate_standard_error <- function(type, sd = NULL, n = NULL, p = NULL,
                                     sd1 = NULL, n1 = NULL, sd2 = NULL, n2 = NULL,
                                     p1 = NULL, p2 = NULL) {
  if (type == "mean") {
    # Standard error for a single mean
    if (is.null(sd) || is.null(n)) stop("Provide 'sd' and 'n' for type 'mean'.")
    se <- sd / sqrt(n)
    explanation <- paste(
      "Standard error for mean = sd / sqrt(n):",
      "sd =", sd, ", n =", n, "=> SE =", round(se, 4)
    )
  } else if (type == "proportion") {
    # Standard error for a single proportion
    if (is.null(p) || is.null(n)) stop("Provide 'p' and 'n' for type 'proportion'.")
    se <- sqrt((p * (1 - p)) / n)
    explanation <- paste(
      "Standard error for proportion = sqrt(p * (1 - p) / n):",
      "p =", p, ", n =", n, "=> SE =", round(se, 4)
    )
  } else if (type == "difference_means") {
    # Standard error for the difference in means
    if (is.null(sd1) || is.null(n1) || is.null(sd2) || is.null(n2)) stop("Provide 'sd1', 'n1', 'sd2', and 'n2' for type 'difference_means'.")
    se <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
    explanation <- paste(
      "Standard error for difference in means = sqrt((sd1^2 / n1) + (sd2^2 / n2)):",
      "sd1 =", sd1, ", n1 =", n1, ", sd2 =", sd2, ", n2 =", n2, "=> SE =", round(se, 4)
    )
  } else if (type == "difference_proportions") {
    # Standard error for the difference in proportions
    if (is.null(p1) || is.null(n1) || is.null(p2) || is.null(n2)) stop("Provide 'p1', 'n1', 'p2', and 'n2' for type 'difference_proportions'.")
    se <- sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2))
    explanation <- paste(
      "Standard error for difference in proportions = sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2)):",
      "p1 =", p1, ", n1 =", n1, ", p2 =", p2, ", n2 =", n2, "=> SE =", round(se, 4)
    )
  } else {
    stop("Invalid 'type'. Choose from 'mean', 'proportion', 'difference_means', or 'difference_proportions'.")
  }

  return(list(
    standard_error = round(se, 4),
    explanation = explanation
  ))
}

#' @title Standard Error for a Single Mean
#' @description Calculates the standard error for a single mean.
#' @param sd Standard deviation of the sample.
#' @param n Sample size.
#' @return A list containing the standard error and detailed explanation.
#' @examples
#' standard_error_mean(sd = 10, n = 50)
#' @export
standard_error_mean <- function(sd, n) {
  if (is.null(sd) || is.null(n)) stop("Both 'sd' and 'n' are required.")
  se <- sd / sqrt(n)
  explanation <- paste(
    "Standard error for mean = sd / sqrt(n):",
    "sd =", sd, ", n =", n, "=> SE =", round(se, 4)
  )
  return(list(
    standard_error = round(se, 4),
    explanation = explanation
  ))
}

#' @title Standard Error for a Single Proportion
#' @description Calculates the standard error for a single proportion.
#' @param p Proportion.
#' @param n Sample size.
#' @return A list containing the standard error and detailed explanation.
#' @examples
#' standard_error_proportion(p = 0.5, n = 100)
#' @export
standard_error_proportion <- function(p, n) {
  if (is.null(p) || is.null(n)) stop("Both 'p' and 'n' are required.")
  se <- sqrt((p * (1 - p)) / n)
  explanation <- paste(
    "Standard error for proportion = sqrt(p * (1 - p) / n):",
    "p =", p, ", n =", n, "=> SE =", round(se, 4)
  )
  return(list(
    standard_error = round(se, 4),
    explanation = explanation
  ))
}

#' @title Standard Error for the Difference in Means
#' @description Calculates the standard error for the difference in two means.
#' @param sd1 Standard deviation of the first group.
#' @param n1 Sample size of the first group.
#' @param sd2 Standard deviation of the second group.
#' @param n2 Sample size of the second group.
#' @return A list containing the standard error and detailed explanation.
#' @examples
#' standard_error_difference_means(sd1 = 15, n1 = 30, sd2 = 20, n2 = 40)
#' @export
standard_error_difference_means <- function(sd1, n1, sd2, n2) {
  if (is.null(sd1) || is.null(n1) || is.null(sd2) || is.null(n2)) {
    stop("All of 'sd1', 'n1', 'sd2', and 'n2' are required.")
  }
  se <- sqrt((sd1^2 / n1) + (sd2^2 / n2))
  explanation <- paste(
    "Standard error for difference in means = sqrt((sd1^2 / n1) + (sd2^2 / n2)):",
    "sd1 =", sd1, ", n1 =", n1, ", sd2 =", sd2, ", n2 =", n2, "=> SE =", round(se, 4)
  )
  return(list(
    standard_error = round(se, 4),
    explanation = explanation
  ))
}

#' @title Standard Error for the Difference in Proportions
#' @description Calculates the standard error for the difference in two proportions.
#' @param p1 Proportion of the first group.
#' @param n1 Sample size of the first group.
#' @param p2 Proportion of the second group.
#' @param n2 Sample size of the second group.
#' @return A list containing the standard error and detailed explanation.
#' @examples
#' standard_error_difference_proportions(p1 = 0.6, n1 = 50, p2 = 0.5, n2 = 60)
#' @export
standard_error_difference_proportions <- function(p1, n1, p2, n2) {
  if (is.null(p1) || is.null(n1) || is.null(p2) || is.null(n2)) {
    stop("All of 'p1', 'n1', 'p2', and 'n2' are required.")
  }
  se <- sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2))
  explanation <- paste(
    "Standard error for difference in proportions = sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2)):",
    "p1 =", p1, ", n1 =", n1, ", p2 =", p2, ", n2 =", n2, "=> SE =", round(se, 4)
  )
  return(list(
    standard_error = round(se, 4),
    explanation = explanation
  ))
}

#' @title Conditional Probability Calculator
#' @description Calculates the conditional probability of an event given a condition from a contingency table.
#' @param table A contingency table (data frame or matrix) where rows represent events (e.g., age groups) and columns represent conditions (e.g., ice cream types).
#' @param event A character string or numeric index specifying the row (event) of interest.
#' @param condition A character string or numeric index specifying the column (condition) of interest.
#' @return A list containing the calculated probability and a detailed explanation.
#' @examples
#' # Define the table
#' data <- matrix(c(15, 8, 7, 20, 11, 15, 8, 7, 9), nrow = 3, byrow = TRUE,
#'                dimnames = list(c("Over 40", "20-40", "Under 20"),
#'                                c("Chocolate", "Vanilla", "Strawberry")))
#' # Calculate probability
#' conditional_probability(table = data, event = "Over 40", condition = "Chocolate")
#' @export
conditional_probability <- function(table, event, condition) {
  # Validate inputs
  if (!is.matrix(table) && !is.data.frame(table)) stop("Input table must be a matrix or data frame.")
  if (is.character(event)) event <- match(event, rownames(table))
  if (is.character(condition)) condition <- match(condition, colnames(table))
  if (is.na(event) || is.na(condition)) stop("Invalid event or condition specified.")

  # Extract relevant data
  condition_total <- sum(table[, condition])
  event_and_condition <- table[event, condition]

  # Calculate probability
  probability <- event_and_condition / condition_total

  # Explanation
  explanation <- paste(
    "Step 1: Identify the total for the condition (column).",
    "   Total for condition (", colnames(table)[condition], "):", condition_total,
    "",
    "Step 2: Identify the count for the event and condition (cell in table).",
    "   Count for event (", rownames(table)[event], ") and condition (", colnames(table)[condition], "):", event_and_condition,
    "",
    "Step 3: Calculate the conditional probability.",
    "   Probability = Count for event and condition / Total for condition.",
    "   Probability =", event_and_condition, "/", condition_total, "=", round(probability, 4)
  )

  return(list(
    probability = round(probability, 4),
    explanation = explanation
  ))
}

#' @title Cumulative Probability for a Discrete Random Variable
#' @description Calculates the cumulative probability \( P(X < x_0) \) for a random variable \( X \) with a given probability distribution.
#' @param x A numeric vector of possible values of the random variable.
#' @param p A numeric vector of probabilities corresponding to the values in `x`. These can be scaled values (e.g., involving k).
#' @param threshold A numeric value specifying the threshold \( x_0 \) for the calculation of \( P(X < x_0) \).
#' @return A list containing the cumulative probability and a detailed explanation.
#' @examples
#' # Example distribution: x = 0, 1, 2, 3 with P(X) = 4k, 5k, 8k, 3k
#' cumulative_probability(x = c(0, 1, 2, 3), p = c(4, 5, 8, 3), threshold = 2)
#' @export
cumulative_probability <- function(x, p, threshold) {
  if (length(x) != length(p)) stop("Vectors 'x' and 'p' must have the same length.")

  # Normalize probabilities
  total_p <- sum(p)
  normalized_p <- p / total_p

  # Filter values where X < threshold
  indices <- x < threshold
  cumulative_p <- sum(normalized_p[indices])

  # Explanation
  explanation <- paste(
    "Step 1: Normalize the probabilities.",
    "   Total sum of probabilities =", total_p,
    "   Normalized probabilities:", paste(round(normalized_p, 4), collapse = ", "),
    "",
    "Step 2: Identify values of X where X < threshold (", threshold, ").",
    "   Values of X considered:", paste(x[indices], collapse = ", "),
    "   Corresponding probabilities:", paste(round(normalized_p[indices], 4), collapse = ", "),
    "",
    "Step 3: Calculate cumulative probability.",
    "   Cumulative P(X < threshold) = Sum of probabilities where X < threshold.",
    "   Cumulative P(X < threshold) =", round(cumulative_p, 4)
  )

  return(list(
    cumulative_probability = round(cumulative_p, 4),
    explanation = explanation
  ))
}

#' @title Confidence Interval for Population Variance
#' @description Calculates the confidence interval for the population variance given the sample variance, sample size, and confidence level.
#' @param sd The population standard deviation.
#' @param n The sample size.
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#' @return A list containing the lower and upper bounds of the confidence interval for the variance.
#' @examples
#' variance_confidence_interval(sd = 5, n = 30, confidence_level = 0.95)
#' @export
variance_confidence_interval <- function(sd, n, confidence_level) {
  # Degrees of freedom
  df <- n - 1

  # Alpha value
  alpha <- 1 - confidence_level

  # Chi-squared critical values
  chi2_lower <- qchisq(alpha / 2, df)
  chi2_upper <- qchisq(1 - alpha / 2, df)

  # Sample variance
  sample_variance <- sd^2

  # Confidence interval bounds for variance
  lower_bound <- (df * sample_variance) / chi2_upper
  upper_bound <- (df * sample_variance) / chi2_lower

  # Explanation
  explanation <- paste(
    "Step 1: Degrees of freedom (df) = n - 1 =", df,
    "",
    "Step 2: Chi-squared critical values:",
    "   Lower critical value =", round(chi2_lower, 4),
    "   Upper critical value =", round(chi2_upper, 4),
    "",
    "Step 3: Calculate the bounds of the confidence interval for variance:",
    "   Lower bound = (df * sample_variance) / chi2_upper = (", df, "*", round(sample_variance, 4), ") /", round(chi2_upper, 4),
    "   Upper bound = (df * sample_variance) / chi2_lower = (", df, "*", round(sample_variance, 4), ") /", round(chi2_lower, 4),
    "",
    "Step 4: Confidence interval for variance = [", round(lower_bound, 4), ",", round(upper_bound, 4), "]"
  )

  # Return results
  return(list(
    lower_bound = round(lower_bound, 4),
    upper_bound = round(upper_bound, 4),
    explanation = explanation
  ))
}


#' @title Confidence Interval for Sample Variance
#' @description Calculates the confidence interval for the sample variance given the sample variance, sample size, and confidence level.
#' @param sample_variance The variance of the sample.
#' @param n The sample size.
#' @param confidence_level The confidence level as a decimal (e.g., 0.95 for 95 percent confidence).
#' @return A list containing the lower and upper bounds of the confidence interval for the sample variance.
#' @examples
#' sample_variance_confidence_interval(sample_variance = 25, n = 30, confidence_level = 0.95)
#' @export
sample_variance_confidence_interval <- function(sample_variance, n, confidence_level) {
  # Degrees of freedom
  df <- n - 1

  # Alpha value
  alpha <- 1 - confidence_level

  # Chi-squared critical values
  chi2_lower <- qchisq(alpha / 2, df)
  chi2_upper <- qchisq(1 - alpha / 2, df)

  # Confidence interval bounds for sample variance
  lower_bound <- (df * sample_variance) / chi2_upper
  upper_bound <- (df * sample_variance) / chi2_lower

  # Explanation
  explanation <- paste(
    "Step 1: Degrees of freedom (df) = n - 1 =", df,
    "",
    "Step 2: Chi-squared critical values:",
    "   Lower critical value =", round(chi2_lower, 4),
    "   Upper critical value =", round(chi2_upper, 4),
    "",
    "Step 3: Calculate the bounds of the confidence interval for sample variance:",
    "   Lower bound = (df * sample_variance) / chi2_upper = (", df, "*", round(sample_variance, 4), ") /", round(chi2_upper, 4),
    "   Upper bound = (df * sample_variance) / chi2_lower = (", df, "*", round(sample_variance, 4), ") /", round(chi2_lower, 4),
    "",
    "Step 4: Confidence interval for sample variance = [", round(lower_bound, 4), ",", round(upper_bound, 4), "]"
  )

  # Return results
  return(list(
    lower_bound = round(lower_bound, 4),
    upper_bound = round(upper_bound, 4),
    explanation = explanation
  ))
}

#' @title Probability for Sampling Distribution of the Sample Mean
#' @description Calculates the probability for the sample mean based on the sampling distribution.
#' @param sample_mean The hypothesized sample mean (e.g., total divided by sample size).
#' @param population_mean The population mean.
#' @param population_sd The population standard deviation.
#' @param n The sample size.
#' @param tail Specify "less", "greater", or "two-tailed" for the probability calculation.
#' @return A list containing the calculated probability and a detailed explanation.
#' @examples
#' # Example: Probability that the total grocery bill for one year is less than $7020
#' sampling_distribution_probability(
#'   sample_mean = 7020 / 52,
#'   population_mean = 140,
#'   population_sd = 10,
#'   n = 52,
#'   tail = "less"
#' )
#' @export
sampling_distribution_probability <- function(sample_mean, population_mean, population_sd, n, tail = "less") {
  # Standard error of the mean
  standard_error <- population_sd / sqrt(n)

  # Calculate the z-score
  z_score <- (sample_mean - population_mean) / standard_error

  # Calculate probability based on the tail type
  if (tail == "less") {
    probability <- pnorm(z_score)
    tail_explanation <- paste("P(X < sample_mean): Area to the left of z =", round(z_score, 4))
  } else if (tail == "greater") {
    probability <- 1 - pnorm(z_score)
    tail_explanation <- paste("P(X > sample_mean): Area to the right of z =", round(z_score, 4))
  } else if (tail == "two-tailed") {
    probability <- 2 * (1 - pnorm(abs(z_score)))
    tail_explanation <- paste("Two-tailed probability: 2 * Area beyond |z| =", round(abs(z_score), 4))
  } else {
    stop("Invalid tail type. Choose 'less', 'greater', or 'two-tailed'.")
  }

  # Explanation
  explanation <- paste(
    "Step 1: Calculate the standard error (SE).",
    "   SE = population_sd / sqrt(n) =", population_sd, "/ sqrt(", n, ") =", round(standard_error, 4),
    "",
    "Step 2: Compute the z-score.",
    "   z = (sample_mean - population_mean) / SE = (", sample_mean, "-", population_mean, ") /", round(standard_error, 4),
    "=", round(z_score, 4),
    "",
    "Step 3: Determine the probability based on the tail type.",
    tail_explanation,
    "",
    "Probability =", round(probability, 6)
  )

  # Return results
  return(list(
    probability = round(probability, 6),
    explanation = explanation
  ))
}
