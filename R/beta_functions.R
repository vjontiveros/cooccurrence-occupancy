library(reticulate)

py_require("scipy")
sc <- import("scipy.special")
np <- import("numpy")

hyp2f1_scipy <- function(a, b, c, z) {
  # Para argumentos complejos, forzamos dtype complejo de numpy
  if (is.complex(z)) {
    z <- np$array(z, dtype = "complex128")
  }
  
  sc$hyp2f1(a, b, c, z)
}

curve_beta <- function(S, p, N, alpha, beta) {
  
  # Guard: alpha, beta must be positive (optimizer can wander)
  if (any(alpha <= 0, beta <= 0)) return(rep(NaN, length(p)))
  
  hval <- hyp2f1_scipy(-N,alpha, alpha + beta,p)
  
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
