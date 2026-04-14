expected_A_empirical <- function(p, q, S, N) {
  pq <- outer(p, q, "*")       
  m  <- rowMeans((1 - pq)^N)   
  (S - 1) * (1 - m)
}

bootstrap_A_empirical <- function(n_boot, p, q, S, N){
  bootstrap_values <- matrix(NA, nrow = n_boot, ncol = length(p))
  for(b in 1:n_boot){
    q_boot <- sample(q, size = S, replace = TRUE)
    bootstrap_values[b, ] <- expected_A_empirical(p, q, S, N)
  }
  bootstrap_values
}

