# Load packages and sources -----------------------------------------------

library(MASS)
library(tidyverse)
library(patchwork)

# The next 4 lines are required to use scipy inside R
library(reticulate)
py_require("scipy")
sc <- import("scipy.special")
np <- import("numpy")

# Load functions
source("R/beta_functions.R")
source("R/expected_A_sp_equivalence.R")

# Load mediolitoral data --------------------------------------------------
load("results/objects/resmedio.RData")

p_data <- res.medio$Freqs
y_data <- res.medio$Sobs[p_data>0]
p_data <- p_data[p_data>0]
S <- length(p_data)
N <- 148

# Estimate beta from the occupancies, plot histogram, and perform bootstrap ----------------

beta_fit <- fit_beta_mle(p_data)

beta_fit


# Plot

bins <- 11
brks <- seq(0, max(p_data), length.out = bins + 1)

s3p1 <- ggplot(data.frame(x = p_data), aes(x = x)) + 
  geom_histogram(breaks = brks, boundary = 0, color = "white",
                 closed = "right", 
                 aes(y = after_stat(density)), linewidth = .2) +
  geom_function(fun = dbeta, 
                args = list(shape1 = beta_fit$alpha, shape2 = beta_fit$beta)) +
  theme_bw() + theme(aspect.ratio = .618, panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 10.5), oob = scales::oob_keep) +
  xlab("Occupancy") + ylab("Probability density")

#Log-scale

eps <- 0.01
s3p2 <- ggplot(data.frame(x = p_data), aes(x = x)) + 
  geom_rect(stat = "bin", breaks = brks, boundary = 0, color = "white",
    linewidth = .2, closed = "right", 
    aes(xmin = after_stat(xmin), xmax = after_stat(xmax),
        ymin = after_stat(ifelse(density > 0, eps, NA_real_)),
        ymax = after_stat(ifelse(density > 0, density, NA_real_)))) +
  geom_function(fun = dbeta,  
                args = list(shape1 = beta_fit$alpha, shape2 = beta_fit$beta)) +
  theme_bw() +
  theme(aspect.ratio = .618, panel.grid.minor = element_blank()) +
  scale_y_log10(limits = c(eps, 20)) +
  xlab("Occupancy") +
  ylab("Probability density")


# Parameter estimation ----------------------------------------------------

#Directly from the occupancies
p_bootstrap <- bootstrap_beta(p_data)
p_bootstrap

ps <- 1:396 / 396

temp <- matrix(NA, nrow = length(p_bootstrap$alpha), ncol = length(ps))
for(i in 1:length(p_bootstrap$alpha)){
  temp[i, ] <- curve_beta(S, ps, N, p_bootstrap$alpha[i], p_bootstrap$beta[i])
}


df.beta <- data.frame(Freqs = ps, 
                      Expected = apply(temp, 2, quantile, probs = .5), 
                      UCI95 = apply(temp, 2, quantile, probs = .975),
                      LCI95 = apply(temp, 2, quantile, probs = .025),
                      Type = "Beta")

# Bootstraping for the expectation under species equivalence --------------

A_values <- expected_A_sp_equivalence(ps, p_data, S, N)

bootstrap <- bootstrap_A_sp_equivalence(n_boot = 1000, ps, p_data, S, N)

df.A <- data.frame(Freqs = ps, Expected = A_values, 
                   UCI95 = apply(bootstrap, 2, quantile, probs = .975),
                   LCI95 = apply(bootstrap, 2, quantile, probs = .025), 
                   Type = "Empirical")


cols <- c(
  "Beta" = "#3366AA",  # muted blue
  "Empirical" = "#A23B72"#,  # muted plum
)
linetypes <- c(3, 1)

sf1 <- ggplot(data = rbind(df.A, df.beta), 
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

 sf3 <- ggplot(data = rbind(df.A, df.beta), 
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



# Load BCI data -----------------------------------------------------------
load("results/objects/resbci.RData")

p_data <- res.bci$Freqs
y_data <- res.bci$Sobs[p_data>0]
p_data <- p_data[p_data>0]
S <- length(p_data)
N <- 1250

# Estimate beta from the occupancies, plot histogram, and perform bootstrap ----------------

beta_fit <- fit_beta_mle(p_data)

beta_fit

# Plot
bins <- 17
brks <- seq(0, max(p_data), length.out = bins + 1)

s3p3 <- ggplot(data.frame(x = p_data), aes(x = x)) + 
  geom_histogram(breaks = brks, boundary = 0, color = "white",
                 closed = "right", 
                 aes(y = after_stat(density)), linewidth = .2) +
  geom_function(fun = dbeta, 
                args = list(shape1 = beta_fit$alpha, shape2 = beta_fit$beta)) +
  theme_bw() + theme(aspect.ratio = .618, panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0, 10.5), oob = scales::oob_keep) +
  xlab("Occupancy") + ylab("Probability density")

#Log-scale

eps <- 0.01
s3p4 <- ggplot(data.frame(x = p_data), aes(x = x)) + 
  geom_rect(stat = "bin", breaks = brks, boundary = 0, color = "white",
            linewidth = .2, closed = "right", 
            aes(xmin = after_stat(xmin), xmax = after_stat(xmax),
                ymin = after_stat(ifelse(density > 0, eps, NA_real_)),
                ymax = after_stat(ifelse(density > 0, density, NA_real_)))) +
  geom_function(fun = dbeta,  
                args = list(shape1 = beta_fit$alpha, shape2 = beta_fit$beta)) +
  theme_bw() +
  theme(aspect.ratio = .618, panel.grid.minor = element_blank()) +
  scale_y_log10(limits = c(eps, 10)) +
  xlab("Occupancy") +
  ylab("Probability density")

# Parameter estimation ----------------------------------------------------

#Directly from the occupancies
p_bootstrap <- bootstrap_beta(p_data)
p_bootstrap

ps <- 1:2500 / 2500

temp <- matrix(NA, nrow = length(p_bootstrap$alpha), ncol = length(ps))
for(i in 1:length(p_bootstrap$alpha)){
  temp[i, ] <- curve_beta(S, ps, N, p_bootstrap$alpha[i], p_bootstrap$beta[i])
}


df.beta <- data.frame(Freqs = ps, 
                      Expected = apply(temp, 2, quantile, probs = .5), 
                      UCI95 = apply(temp, 2, quantile, probs = .975),
                      LCI95 = apply(temp, 2, quantile, probs = .025),
                      Type = "Beta")

# Bootstraping for the expectation under species equivalence --------------

A_values <- expected_A_sp_equivalence(ps, p_data, S, N)

bootstrap <- bootstrap_A_sp_equivalence(n_boot = 1000, ps, p_data, S, N)

df.A <- data.frame(Freqs = ps, Expected = A_values, 
                   UCI95 = apply(bootstrap, 2, quantile, probs = .975),
                   LCI95 = apply(bootstrap, 2, quantile, probs = .025), 
                   Type = "Empirical")


sf2 <- ggplot(data = rbind(df.A, df.beta), 
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

sf4 <- ggplot(data = rbind(df.A, df.beta), 
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


# Saving plots ------------------------------------------------------------

(sf3 | sf4)/(sf1 | sf2)
ggsave(filename = "results/figures/fig_s2.png", width = 180, height = 143.6, 
       units = "mm", dpi = 300, )

(s3p1 | s3p3)/(s3p2 | s3p4)
ggsave(filename = "results/figures/fig_s3.png", width = 180, height = 143.6, 
       units = "mm", dpi = 300, )

