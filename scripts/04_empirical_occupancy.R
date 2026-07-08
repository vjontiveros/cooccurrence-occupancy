
# Load packages and functions ---------------------------------------------
library(tidyverse)
library(patchwork)

source("R/expected_A_empirical.R")
source("R/logseries_functions.R")


# Mediolittoral -----------------------------------------------------------
load("results/objects/resmedio.RData")

# Calculate the expected associations for the empirical distribution --------
S <- sum(res.medio$Freqs > 0) #Number of species present
N <- 148 #Number of sites
ps <- 5:1000 / 1000 #Occupancies to produce the expected associations
qs <- res.medio %>% filter(Freqs > 0) %>% pull(Freqs)

A_values <- expected_A_empirical(ps, qs, S, N)


# Bootstrapping -----------------------------------------------------------

# Bootstrapping the occupancies to obtain expected associations for the
# empirical distribution of occupancies.

bootstrap <- bootstrap_A_empirical(n_boot = 200, ps, qs, S, N)

df.A <- data.frame(Freqs = ps, Expected = A_values, 
                   UCI95 = apply(bootstrap, 2, quantile, probs = .95),
                   LCI95 = apply(bootstrap, 2, quantile, probs = .05), 
                   Type = "Empirical")


# Log-series --------------------------------------------------------------

# As we know the expression for cooccurrence-occupancy relationship given a
# log-series, we fit the theta to this expression.

p_data <- res.medio %>% filter(Freqs > 0) %>% pull(Freqs)
y_data <- res.medio %>% filter(Freqs > 0) %>% pull(Sobs)

S <- length(p_data)

# Getting the estimates for the logseries. First the bootstrap then the values
set.seed(42)
theta_boot <- bootstrap_logser(p_data * N)

logser_mle(p_data * N) 
quantile(theta_boot, probs = c(0.025, 0.975))

# Obtain the expected curve under a logSeries distribution
df.ls <- data.frame(Freqs = ps, 
                    Expected = curve_logSeries(S, ps, quantile(theta_boot, 0.5)),
                    UCI95 = curve_logSeries(S, ps, quantile(theta_boot, 0.975)),
                    LCI95 = curve_logSeries(S, ps, quantile(theta_boot, 0.025)), 
                    Type = "LogSeries")




# Plotting ----------------------------------------------------------------
cols <- c(
  # "Beta" = "#3366AA",  # muted blue
  "Empirical" = "#D55E00",  # neutral grey
  "LogSeries" = "#117733"   # muted purple
)
linetypes <- 1:2


p1 <- 
  ggplot(data = rbind(df.A, df.ls), 
       aes(x = Freqs, color = Type, linetype = Type)) + 
  geom_point(data = res.medio %>% filter(Freqs > 0), aes(x = Freqs, y = Sobs), inherit.aes = F) + 
  scale_x_log10() +
  geom_ribbon(aes(ymax = UCI95, ymin = LCI95, fill = Type), alpha = .35, color = NA) +
  geom_line(aes(y = Expected)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  scale_linetype_manual(values = linetypes) +
  theme(aspect.ratio = .618, legend.position = "none", 
        panel.grid.minor = element_blank()) +
  xlab("Occupancy") + ylab("Cooccurrences")

p2 <- ggplot(data = rbind(df.A, df.ls), 
       aes(x = Freqs, color = Type, linetype = Type)) + 
  geom_point(data = res.medio %>% filter(Freqs > 0), aes(x = Freqs, y = Sobs), inherit.aes = F) + 
  geom_ribbon(aes(ymax = UCI95, ymin = LCI95, fill = Type), alpha = .5, color = NA) +
  geom_line(aes(y = Expected)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  scale_linetype_manual(values = linetypes) +
  theme(aspect.ratio = .618, legend.position = "none", 
        panel.grid.minor = element_blank()) +
  xlab("Occupancy") + ylab("Cooccurrences")


# BCI ---------------------------------------------------------------------

load("results/objects/resbci.RData")

# Calculate the expected associations for the empirical distribution --------
S <- sum(res.bci$Freqs > 0) #Number of species present
N <- 1250 #Number of sites
ps <- 1:2000 / 2000 #Occupancies to produce the expected associations
qs <- res.bci %>% filter(Freqs > 0) %>% pull(Freqs)

A_values <- expected_A_empirical(ps, qs, S, N)


# Bootstrapping -----------------------------------------------------------

# Bootstrapping the expected associations for the empirical distribution of
# occupancies (occupancies now bootstrapped).

bootstrap <- bootstrap_A_empirical(n_boot = 200, ps, qs, S, N)

df.A <- data.frame(Freqs = ps, Expected = A_values, 
                   UCI95 = apply(bootstrap, 2, quantile, probs = .95),
                   LCI95 = apply(bootstrap, 2, quantile, probs = .05), 
                   Type = "Empirical")


# Log-series --------------------------------------------------------------

# As we know the expression for cooccurrence-occupancy relationship given a
# log-series, we fit the theta to this expression.

p_data <- res.bci %>% filter(Freqs > 0) %>% pull(Freqs)
y_data <- res.bci %>% filter(Freqs > 0) %>% pull(Sobs)

S <- length(p_data)

# Getting the estimates for the logseries. First the bootstrap then the values
set.seed(42)
theta_boot <- bootstrap_logser(p_data * N)

logser_mle(p_data * N) 
quantile(theta_boot, probs = c(0.025, 0.975))

# Obtain the expected curve under a logSeries distribution
df.ls <- data.frame(Freqs = ps, 
                    Expected = curve_logSeries(S, ps, quantile(theta_boot, 0.5)),
                    UCI95 = curve_logSeries(S, ps, quantile(theta_boot, 0.975)),
                    LCI95 = curve_logSeries(S, ps, quantile(theta_boot, 0.025)), 
                    Type = "LogSeries")



# Plots -------------------------------------------------------------------

p3 <- ggplot(data = rbind(df.A, df.ls), 
       aes(x = Freqs, color = Type, linetype = Type)) + 
  geom_point(data = res.bci %>% filter(Freqs > 0), aes(x = Freqs, y = Sobs), inherit.aes = F) + 
  scale_x_log10() +
  geom_ribbon(aes(ymax = UCI95, ymin = LCI95, fill = Type), alpha = .35, color = NA) +
  geom_line(aes(y = Expected)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  scale_linetype_manual(values = linetypes) +
  theme(aspect.ratio = .618, legend.position = "none", 
        panel.grid.minor = element_blank()) +
  xlab("Occupancy") + ylab("Cooccurrences")

p4 <- ggplot(data = rbind(df.A, df.ls), 
       aes(x = Freqs, color = Type, linetype = Type)) + 
  geom_point(data = res.bci %>% filter(Freqs > 0), aes(x = Freqs, y = Sobs), inherit.aes = F) + 
  geom_ribbon(aes(ymax = UCI95, ymin = LCI95, fill = Type), alpha = .5, color = NA) +
  geom_line(aes(y = Expected)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  scale_linetype_manual(values = linetypes) +
  theme(aspect.ratio = .618, legend.position = "none", 
        panel.grid.minor = element_blank()) +
  xlab("Occupancy") + ylab("Cooccurrences")


# Merged plot -------------------------------------------------------------

(p2 | p4)/(p1 | p3)

ggsave(filename = "results/figures/fig4_empirical.png", width = 180, height = 143.6, 
       units = "mm", dpi = 300, )

