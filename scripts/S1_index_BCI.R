
# Load packages -----------------------------------------------------------

library(readxl)
library(tidyverse)


# Merging demography and association --------------------------------------

# Why BCI species have a negative or positive association index? Maybe
# demographic characteristics of the species can shed some light over this.
demography <- read_excel("data/raw/aaz4797_ruger_data_s1.xlsx", sheet = 2)
load("results/objects/resbci.RData")

bci.demo <- res.bci %>% select(Species, Index) %>% 
  inner_join(demography %>% mutate(Species = tolower(sp)))

# Explore models
lm(data = bci.demo, Index ~ PC1score + PC2score) %>% summary
lm(data = bci.demo, Index ~ PC1score + PC2score + PC2score:wd) %>% summary

# Compare different models
m0 <- lm(data = bci.demo, Index ~ PC1score + PC2score)
m1 <- lm(data = bci.demo, Index ~ PC1score + PC2score + PC2score:wd)
m2 <- lm(data = bci.demo, Index ~ PC1score)
m3 <- lm(data = bci.demo, Index ~ PC1score + PC2score + wd)
m4 <- lm(data = bci.demo, Index ~ PC2score)

anova(m2, m4) #Easiest models
anova(m2, m0) #Increasing complexity
anova(m0, m3) #Second model is not clearly better
anova(m0, m1) #m1 seems the best model

m1 %>% summary() #Reported in the main text.

# Plot
ggplot(bci.demo, aes(x = PC1score, y = PC2score, fill = Index)) + 
  geom_point(pch = 21, size = 2) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#D9D9D9", high = "#F59D0A") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1, 
        legend.position = "bottom")

ggsave(filename = "results/figures/Sfig1_tradeoffs.png", width = 90, height = 90, 
       units = "mm", dpi = 300, )
