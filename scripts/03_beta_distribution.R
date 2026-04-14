# Load packages -----------------------------------------------------------
library(tidyverse)
library(patchwork)

# U-shaped beta distribution ----------------------------------------------
p1 <- ggplot() + xlim(0, 1) + 
  geom_function(fun = dbeta, args = list(shape1 = .6, shape2 = .7)) +
  theme_bw() + xlab("Occupancy") + ylab("Probability density") +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618)

# Simulate it nrep times, for 100 sp and 100 sites. It takes quite some time, so
# I include it as an object in the results.

set.seed(7289495) #Set seed for reproducibility
nrep <- 1000
sims <- matrix(NA, nrow = 100 * nrep, ncol = 7)
for(rep in 1:nrep){
  occupancies <- rbeta(100, shape1 = .6, shape2 = .7) # Random occupancies
  r.numbers <- runif(10000) %>% matrix(nrow = 100, ncol = 100) 
  spxsite <- ((r.numbers) < occupancies) * 1.0 # We realize the occupancies
  c1 <- any(rowSums(spxsite) == 0) # Check no species goes without presence
  c2 <- any(colSums(spxsite) == 0) # Check no site is empty
  while(c1 | c2){
    occupancies <- rbeta(100, shape1 = .6, shape2 = .7)
    r.numbers <- runif(10000) %>% matrix(nrow = 100, ncol = 100)
    spxsite <- ((r.numbers) < occupancies) * 1.0
    print("Absent species or empty site")
    c1 <- any(rowSums(spxsite) == 0)
    c2 <- any(colSums(spxsite) == 0)
  }
  index <- c((1 + (rep - 1)*100):(rep*100))
  sims[index, 1:7] <- as.matrix(m_values(spxsite))
  print(rep)
}
colnames(sims) <- c("Freqs", "Sobs", "UCI95", "LCI95", "Expectation", 
                    "Variance", "Index")

# Save and load data in case that it takes too much time to run.
save(sims, file = "results/objects/sims_beta_U.RData")
# load("results/objects/sims_beta_U.RData") # uncomment if needed

df2 <- sims %>% as.data.frame() %>% group_by(Freqs) %>% 
  summarise(Q05 = quantile(Sobs, .05), Q95 = quantile(Sobs, .95), 
            mean = mean(Sobs))
df2 <- rbind(c(0, 0, 0, 0), df2) # Add 0 at the origin

p3 <- ggplot(df2, aes(x = Freqs, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = .5, fill = "mediumpurple3") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618) +
  xlab("Occupancy") + ylab("Cooccurrences") +
  scale_y_continuous(limits = c(0, 100))


# Another beta, decreading monotonically -----------------------------------------------

p2 <- ggplot() + xlim(0, 1) + 
  geom_function(fun = dbeta, args = list(shape1 = 1, shape2 = 3)) +
  theme_bw() + xlab("Occupancy") + ylab("Probability density") +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618)

# Simulate it nrep times, for 100 sp and 100 sites. It takes quite some time, so
# I include it as an object in the results. Same steps as above, just changing
# the beta parameters.
set.seed(611837)
nrep <- 1000
sims2 <- matrix(NA, nrow = 100 * nrep, ncol = 7)
for(rep in 1:nrep){
  occupancies <- rbeta(100, shape1 = 1, shape2 = 3)
  r.numbers <- runif(10000) %>% matrix(nrow = 100, ncol = 100)
  spxsite <- ((r.numbers) < occupancies) * 1.0
  c1 <- any(rowSums(spxsite) == 0)
  c2 <- any(colSums(spxsite) == 0)
  while(c1 | c2){
    occupancies <- rbeta(100, shape1 = 1, shape2 = 3)
    r.numbers <- runif(10000) %>% matrix(nrow = 100, ncol = 100)
    spxsite <- ((r.numbers) < occupancies) * 1.0
    print("Absent species or empty site")
    c1 <- any(rowSums(spxsite) == 0)
    c2 <- any(colSums(spxsite) == 0)
  }
  index <- c((1 + (rep - 1)*100):(rep*100))
  sims2[index, 1:7] <- as.matrix(m_values(spxsite))
  print(rep)
}
colnames(sims2) <- colnames(sims)

# Save and load data in case that it takes too much time to run.
save(sims2, file = "results/objects/sims_beta_decr.RData")
# load("results/objects/sims_beta_decr.RData") # uncomment if needed

df2 <- sims2 %>% as.data.frame() %>% group_by(Freqs) %>% 
  summarise(Q05 = quantile(Sobs, .05), Q95 = quantile(Sobs, .95), 
            mean = mean(Sobs))
df2 <- rbind(c(0, 0, 0, 0), df2)

# Plot
p4 <- ggplot(df2, aes(x = Freqs, y = mean)) +
  geom_line() +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = .5, fill = "mediumpurple3") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618) +
  xlab("Occupancy") + ylab("Cooccurrences") +
  scale_y_continuous(limits = c(0, 100))


# Merged plot -------------------------------------------------------------------
(p2 | p1)/(p4 | p3)

ggsave(filename = "results/figures/fig3_betas.png", width = 180, height = 143.6, 
       units = "mm", dpi = 300, )

