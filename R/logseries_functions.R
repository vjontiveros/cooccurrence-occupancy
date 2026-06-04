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

# Adjusting theta to the theoretical curve
bootstrap_nlsLM_logser <- function(y_data, S, p_data, B = 1000, start_theta) {
  
  n <- length(y_data)
  theta_boot <- numeric(B)
  
  for (b in 1:B) {
    
    # Re-muestreo con reemplazo
    idx <- sample(1:n, size = n, replace = TRUE)
    
    y_b <- y_data[idx]
    p_b <- p_data[idx]
    
    
    fit_b <- tryCatch( nlsLM(
      y_b ~ curve_logSeries(S, p_b, theta),
      start = list(theta = start_theta),
      lower = 1e-7,
      upper = 1,
      weights = p_b
    ),
    error = function(e) NULL
    )
    
    # Guardar resultado si converge
    if (!is.null(fit_b)) {
      theta_boot[b] <- coef(fit_b)["theta"]
    } else {
      theta_boot[b] <- NA
    }
  }
  
  # Eliminar fallos
  theta_boot <- theta_boot[!is.na(theta_boot)]
  
  return(theta_boot)
}
