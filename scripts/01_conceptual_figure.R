
# Load packages and functions ---------------------------------------------

library(tidyverse)
source("R/m_values.R")

# Define the sp x site matrix in the conceptual figure --------------------

spxsite <- matrix(data = c(rep(1, 6), 
                           rep(1, 3), rep(0, 3),
                           1, 1, 0, 1, 0, 0,
                           1, 0, 0, 1, 0, 0,
                           0, 0, 1, 0, 1, 0,
                           rep(0, 5), 1), 
                  nrow = 6,
                  ncol = 6,
                  byrow = T)

# Calculate expectations

df <- m_values(spxsite)

# Adding 0 at the origin.
df2 <- df
df2[7, ] <- rep(0,7)
df2 <- df2 %>% mutate(UCI95 = pmin(UCI95, 5))

# Plot

ggplot(df, aes(x = Freqs, y = Sobs)) +
  geom_ribbon(data = df2, aes(ymin = LCI95, ymax = UCI95), alpha = .5, fill = "lightcoral") +
  geom_point() +
  geom_line(data = df2, aes(y = Expectation), linetype = 2) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618) +
  xlab("Occupancy") + ylab("Cooccurrences")

ggsave(
  filename = "results/figures/fig1_conceptual.png",
  width = 90, height = 56, units = "mm",
  dpi = 300,
  bg = "white"
)
