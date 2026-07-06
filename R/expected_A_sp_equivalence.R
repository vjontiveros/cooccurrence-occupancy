expected_A_sp_equivalence <- function(p, q, S, N) {
  pq <- outer(p, q, "*")       
  m  <- rowMeans((1 - pq)^N)   
  (S - 1) * (1 - m)
}

bootstrap_A_sp_equivalence <- function(n_boot, p, q, S, N){
  bootstrap_values <- matrix(NA, nrow = n_boot, ncol = length(p))
  for(b in 1:n_boot){
    q_boot <- sample(q, size = S, replace = TRUE)
    bootstrap_values[b, ] <- expected_A_sp_equivalence(p, q_boot, S, N)
  }
  bootstrap_values
}
