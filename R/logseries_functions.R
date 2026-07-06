# Expected curve for the log-series distribution
curve_logSeries <- function(S, p, theta){
  (S - 1) * (1 - log(1 - theta * (1 - p))/log(1 - theta))
}

# M.l.e. for the log-series parameter
logser_mle <- function(data) {
  sample_mean <- mean(data)

  # Objective function
  objective <- function(p) {
    if (p <= 0 || p >= 1) {
      return(1e10)
    }
    return((p / (-(1 - p) * log(1 - p))) - sample_mean)
  }

  # Find the root
  root <- uniroot(objective, interval = c(1.e-6, 1.))$root
  return(root)
}

# Bootstrap for estimation
bootstrap_logser <- function(data, B = 1000) {
  n <- length(data)
  p_boot <- numeric(B)

  for (b in 1:B) {
    # resample with replacement
    sample_b <- sample(data, size = n, replace = TRUE)

    # Mle of the bootstrapped sample
    p_boot[b] <- logser_mle(sample_b)
  }

  return(p_boot)
}
