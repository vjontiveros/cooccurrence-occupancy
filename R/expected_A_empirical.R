expected_A_empirical <- function(p, q, S, N) {
  pq <- outer((1 - p), q * N, "^")       
  m  <- rowMeans(pq)   
  (S - 1) * (1 - m)
}


bootstrap_A_empirical <- function(n_boot, p, q, S, N){
  bootstrap_values <- matrix(NA, nrow = n_boot, ncol = length(p))
  for(b in 1:n_boot){
    q_boot <- sample(q, size = S, replace = TRUE)
    bootstrap_values[b, ] <- expected_A_empirical(p, q_boot, S, N)
  }
  bootstrap_values
}

