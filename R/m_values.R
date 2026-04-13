

m_values <- function(sp.site, sp.names = NULL, mat.name = NULL){
   cols <- ncol(sp.site)
   rows <- nrow(sp.site)
   coocs <- tcrossprod(sp.site) # This gives a matrix for the cooccurrences of species i with species j.
   freqs <- diag(coocs)/cols

   Ns <- diag(coocs)
   diag(coocs) <- 0

   Sobs <- rowSums(coocs > 0)

   ExpMeans <- rep(NA, rows)
   ExpVar <- rep(NA, rows)
   for (i in 1:rows){
      temp <- 0
      temp3 <- 0
      for (j in 1:rows){
         if(i == j) next
         pij <- 1 - exp(lchoose(cols - Ns[i], Ns[j]) - lchoose(cols, Ns[j]))
         temp <- temp + pij
         temp2 <- 0
         for (k in (j + 1):rows){
            if(i == k) next

            temp2 <- sum(temp2, 1 - exp(lchoose(cols - Ns[k], Ns[i]) - lchoose(cols, Ns[i])), na.rm = T)
         }

         temp3 <- temp3 + pij * temp2

      }

      ExpMeans[i] <- temp
      ExpVar[i] <- temp * (1 - temp) + 2 * temp3
   }

   UL95 <- ExpMeans + 1.96 * sqrt(ExpVar)
   LL95 <- ExpMeans - 1.96 * sqrt(ExpVar)

   Index <- (Sobs - ExpMeans)/sqrt(ExpVar)

   output <- data.frame(Freqs = freqs,
                        Sobs = Sobs,
                        UCI95 = UL95,
                        LCI95 = LL95,
                        Expectation = ExpMeans,
                        Variance = ExpVar,
                        Index = Index)

   if(!is.null(sp.names)) output <- output %>% add_column(Species = sp.names)
   if(!is.null(mat.name)) output <- output %>% add_column(Habitat = mat.name)

   output
}
