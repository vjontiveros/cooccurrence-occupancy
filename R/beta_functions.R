library(reticulate)

py_require("scipy")
sc <- import("scipy.special")
np <- import("numpy")

hyp2f1_scipy <- function(a, b, c, z) {
  # Just in case z is complex, we force dtype complex of numpy
  if (is.complex(z)) {
    z <- np$array(z, dtype = "complex128")
  }
  
  sc$hyp2f1(a, b, c, z)
}

hyp1f1_scipy <- function(a, b, z) {
  # Just in case z is complex, we force dtype complex of numpy
  if (is.complex(z)) {
    z <- np$array(z, dtype = "complex128")
  }
  
  sc$hyp1f1(a, b, z)
}

curve_beta <- function(S, p, N, alpha, beta) {
  
  # Guard: alpha, beta must be positive (optimizer can wander)
  if (any(alpha <= 0, beta <= 0)) return(rep(NaN, length(p)))
  
  hval <- hyp2f1_scipy(-N,alpha, alpha + beta,p)
  
  (S - 1) * (1 - hval)
}

curve_beta_b <- function(S, p, N, alpha, beta) {
  # Guard: alpha, beta must be positive (optimizer can wander)
  if (any(alpha <= 0, beta <= 0)) return(rep(NaN, length(p)))
  
  hval <- hyp1f1_scipy(alpha, alpha + beta, N * log(1 - p))
  
  (S - 1) * (1 - hval)
}


bootstrap_nlsLM_beta <- function(y_data, S, p_data, N, B = 1000, start_alpha, 
                                 start_beta, p_sim){
  
  n <- length(y_data)
  values <- matrix(data = NA, ncol = length(p_sim), nrow = B)
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    y_b <- y_data[idx]
    p_b <- p_data[idx]
    
    fit_b <- tryCatch(
      nlsLM(
        y_b ~ curve_beta(S, p_b, N, alpha, beta),
        start   = list(alpha = start_alpha, beta = start_beta),
        lower   = c(0., 0.),
        upper   = c(1.e7, 1.e7),
        weights = p_b
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_b)) {
      alpha_boot <- coef(fit_b)["alpha"]
      beta_boot  <- coef(fit_b)["beta"]
    } else {
      next
    }
    values[b, ] <- curve_beta(S, p_sim, N, alpha_boot, beta_boot)
    
  }
  central_curve <- apply(values, MARGIN = 2, FUN = quantile, probs = 0.5, na.rm = TRUE)
  UCI95 <- apply(values, MARGIN = 2, quantile, probs = 0.975, na.rm = T)
  LCI95 <- apply(values, MARGIN = 2, quantile, probs = 0.025, na.rm = TRUE)
  
  data.frame(Freqs = p_sim, 
             Expected = central_curve,
             UCI95 = UCI95,
             LCI95 = LCI95, 
             Type = "Beta")
}

bootstrap_nlsLM_beta_b <- function(y_data, S, p_data, N, B = 1000, start_alpha, 
                                 start_beta, p_sim){
  
  n <- length(y_data)
  values <- matrix(data = NA, ncol = length(p_sim), nrow = B)
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    y_b <- y_data[idx]
    p_b <- p_data[idx]
    
    fit_b <- tryCatch(
      nlsLM(
        y_b ~ curve_beta_b(S, p_b, N, alpha, beta),
        start   = list(alpha = start_alpha, beta = start_beta),
        lower   = c(0., 0.),
        upper   = c(1.e7, 1.e7),
        weights = p_b
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_b)) {
      alpha_boot <- coef(fit_b)["alpha"]
      beta_boot  <- coef(fit_b)["beta"]
    } else {
      next
    }
    values[b, ] <- curve_beta_b(S, p_sim, N, alpha_boot, beta_boot)
    print(beta_boot)
  }
  central_curve <- apply(values, MARGIN = 2, FUN = quantile, probs = 0.5, na.rm = TRUE)
  UCI95 <- apply(values, MARGIN = 2, quantile, probs = 0.975, na.rm = T)
  LCI95 <- apply(values, MARGIN = 2, quantile, probs = 0.025, na.rm = TRUE)
  
  data.frame(Freqs = p_sim, 
             Expected = central_curve,
             UCI95 = UCI95,
             LCI95 = LCI95, 
             Type = "Beta")
}


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
bootstrap_beta_pars_curve <- function(S, p_data, y_data, N, B = 1000, start_alpha, 
                                      start_beta, ci_level = 0.95, seed = 42){
  set.seed(seed)
  
  n <- length(p_data)
  values <- matrix(data = NA, ncol = 2, nrow = B)
  
  for (b in 1:B) {
    
    idx <- sample(1:n, size = n, replace = TRUE)
    y_b <- y_data[idx]
    p_b <- p_data[idx]
    
    fit_b <- tryCatch(
      nlsLM(
        y_b ~ curve_beta_b(S, p_b, N, alpha, beta),
        start   = list(alpha = start_alpha, beta = start_beta),
        lower   = c(0., 0.),
        upper   = c(1.e7, 1.e7),
        weights = p_b
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_b)) {
      alpha_boot <- coef(fit_b)["alpha"]
      beta_boot  <- coef(fit_b)["beta"]
    }
    else {
      next
    }
    values[b, ] <- c(alpha_boot, beta_boot)
  }
  tail_prob <- (1 - ci_level) / 2
  probs     <- c(tail_prob, 0.5, 1 - tail_prob)
  ci_labels <- c("lower", "median", "upper")
  
  alpha_boot <- values[, 1]
  beta_boot <- values[, 2]
  list(
    alpha      = alpha_boot,
    beta       = beta_boot,
    ci_alpha   = setNames(quantile(alpha_boot, probs), ci_labels),
    ci_beta    = setNames(quantile(beta_boot,  probs), ci_labels),
    ci_level   = ci_level
  )
}
