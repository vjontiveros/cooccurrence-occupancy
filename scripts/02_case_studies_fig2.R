# Load packages and functions -----------------------------------------------------------
library(readxl)
library(tidyverse)
library(patchwork)
source("R/m_values.R")

# Read data for Catalan mediolittoral habitat --------------------------------------------------------------- 

# The data consists in Catalan littoral habitat species (188) in the rows and 148
# samples for each of the 3 zones of the littoral habitat: infra-, medio-, and
# supralittoral. We keep just the mediolittoral for illustration purposes.

benthic <- read_excel("data/raw/supramedioinfra_labelled_rev.xlsx")

benthic[is.na(benthic)] <- 0 #NAs are actually 0s.

medio <- benthic[, 1:148*3] 
# supra <- benthic[, c(1:148*3 - 1)] #Supra
# infra <- benthic[, c(1:148*3 + 1)] #Infra

# medio <- ((medio + supra + infra) > 0) * 1.0 # Complete littoral

sp.site <- as.matrix(medio)
sp.site[sp.site > 0] <- 1 #transform to P/A

# Calculate the m-curve
res.medio <- m_values(sp.site, sp.names = unname(unlist(c(benthic[, 1]))), 
                      mat.name = "Medio")

# Save object for further analyses
save(res.medio, file = "results/objects/resmedio.RData")

# Plotting
p1 <- ggplot(data = res.medio, aes(x = Freqs, y = Sobs, text = Species)) +
  geom_ribbon(aes(x = Freqs, ymin = LCI95, ymax = UCI95), 
              fill = "lightcoral", alpha = 0.5, inherit.aes = F)  +
  geom_line(aes(x = Freqs, y = Expectation), inherit.aes = F, 
            lty = 2, col = "grey20", alpha = 0.9) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618) +
  xlab("Occupancy") + ylab("Cooccurrences")

p2 <- ggplot(data = res.medio, aes(x = Freqs, y = Index, text = Species)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -2,   ymax = 2,
           fill = "lightcoral", alpha = 0.5) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618) +
  xlab("Occupancy") + ylab("Association Index")


# Read BCI data -----------------------------------------------------------

# BCI data from version Jun 07, 2019 is available as a Dryad dataset located at
# https://datadryad.org/dataset/doi:10.15146/5xcp-0d46 . Specifically, we
# downloaded file bci.tree.zip and extracted file bci.tree8.rdata in folder
# data. This file contains a record for every tree ever recorded in the 2015 BCI
# census. Dead trees and other situations were also recorded, so for our
# purposes we kept just trees with a status that implies that the tree was alive
# (A, AR, AD, AM).

load("data/raw/bci.tree8.rdata")

bci <- 
  bci.tree8 %>% 
  filter(status %in% c("A", "AR", "AD", "AM")) %>% 
  select(quadrat, sp) %>% group_by(quadrat, sp) %>% 
  summarise(PA = (n() > 0) * 1.0, .groups = "drop_last") %>% 
  pivot_wider(names_from = quadrat, values_from = PA, values_fill = 0) %>% 
  ungroup() %>% as.data.frame()

bci.names <- bci %>% pull(sp)
bci <- bci %>% select(-1) %>% as.matrix()

res.bci <- m_values(bci, sp.names = bci.names, mat.name = "BCI")

# Save object for further analyses
save(res.bci, file = "results/objects/resbci.RData")

# Plots
p3 <-
  ggplot(data = res.bci, aes(x = Freqs, y = Sobs, text = Species)) +
  geom_ribbon(aes(x = Freqs, ymin = LCI95, ymax = UCI95), fill = "lightcoral", alpha = 0.5, inherit.aes = F)  +
  geom_line(aes(x = Freqs, y = Expectation), inherit.aes = F, lty = 2, col = "grey20", alpha = 0.9) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618) +
  xlab("Occupancy") + ylab("Cooccurrences")

p4 <-
  ggplot(data = res.bci, aes(x = Freqs, y = Index, text = Species)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = -2,   ymax = 2,
           fill = "lightcoral", alpha = 0.5) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = .618) +
  xlab("Occupancy") + ylab("Association Index")


# Plotting and saving all together ----------------------------------------

(p1 | p3)/(p2 | p4)

ggsave(filename = "results/figures/fig2_case_studies.png", width = 180, height = 143.6, units = "mm", dpi = 300)

