# Function to calculate the 2F1 hypergeometric, from python. 
hyp2f1_scipy <- function(a, b, c, z) {
  # Just in case z is complex, we force dtype complex of numpy
  if (is.complex(z)) {
    z <- np$array(z, dtype = "complex128")
  }

  sc$hyp2f1(a, b, c, z)
}

# Function to obtain the curve assuming species equivalence and beta
# distribution.
curve_beta <- function(S, p, N, alpha, beta) {
  # Guard: alpha, beta must be positive (optimizer can wander)
  if (any(alpha <= 0, beta <= 0)) return(rep(NaN, length(p)))

  hval <- hyp2f1_scipy(-N,alpha, alpha + beta,p)

  (S - 1) * (1 - hval)
}

# Function to estimate the beta distribution corresponding to the occupancies.
fit_beta_mle <- function(p_data) {
  # Remove boundary values (0 and 1 cause log(0) in Beta log-likelihood)
  p_clean <- p_data[p_data > 0 & p_data < 1]
  n_removed <- length(p_data) - length(p_clean)
  if (n_removed > 0)
    message(sprintf("Removed %d boundary value(s) (0 or 1) before fitting.", n_removed))

  fit <- fitdistr(p_clean, "beta",
                  start = list(shape1 = 1, shape2 = 1),
                  lower = 1e-6)

  list(
    alpha    = unname(fit$estimate["shape1"]),
    beta     = unname(fit$estimate["shape2"]),
    se_alpha = unname(fit$sd["shape1"]),
    se_beta  = unname(fit$sd["shape2"]),
    fit      = fit
  )
}

# Function to estimate the beta parameters by bootstrap. 
bootstrap_beta <- function(p_data, B = 2000, ci_level = 0.95, seed = 42) {
  set.seed(seed)

  n        <- length(p_data)
  p_clean  <- p_data[p_data > 0 & p_data < 1]
  n_clean  <- length(p_clean)

  alpha_boot <- numeric(B)
  beta_boot  <- numeric(B)
  failed     <- 0L

  for (b in seq_len(B)) {
    p_b <- sample(p_clean, size = n_clean, replace = TRUE)

    fit_b <- tryCatch(
      fitdistr(p_b, "beta",
               start = list(shape1 = 1, shape2 = 1),
               lower = 1e-6),
      error   = function(e) NULL,
      warning = function(w) NULL
    )

    if (!is.null(fit_b)) {
      alpha_boot[b] <- fit_b$estimate["shape1"]
      beta_boot[b]  <- fit_b$estimate["shape2"]
    } else {
      alpha_boot[b] <- NA_real_
      beta_boot[b]  <- NA_real_
      failed        <- failed + 1L
    }
  }

  if (failed > 0)
    message(sprintf("%d / %d bootstrap iterations failed and were excluded.", failed, B))

  converged <- !is.na(alpha_boot) & !is.na(beta_boot)
  alpha_boot <- alpha_boot[converged]
  beta_boot  <- beta_boot[converged]

  tail_prob <- (1 - ci_level) / 2
  probs     <- c(tail_prob, 0.5, 1 - tail_prob)
  ci_labels <- c("lower", "median", "upper")

  list(
    alpha      = alpha_boot,
    beta       = beta_boot,
    ci_alpha   = setNames(quantile(alpha_boot, probs), ci_labels),
    ci_beta    = setNames(quantile(beta_boot,  probs), ci_labels),
    n_success  = sum(converged),
    n_failed   = failed,
    ci_level   = ci_level
  )
}

