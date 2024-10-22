library(tidyverse)
library(cmdstanr)
library(MCMCvis)

# now read in model data

# leps
leps_data <- readRDS("lepsEmergence-2024-07-17/lepEmergence-data-2024-07-17.rds")
leps_fit <- readRDS("lepsEmergence-2024-07-17/lepsEmergence-fit-2024-07-17.rds")


# supplemental response of leps by ows
gamma_df_leps <- data.frame(kappa = MCMCvis::MCMCpstr(leps_fit, params = 'gamma', 
                                                      func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                            kappa =MCMCvis::MCMCpstr(leps_fit, params = 'gamma', 
                                                     fun = median)) %>% 
  mutate(Pheno.Phase = c("Butterfly Emergence [Egg OWS]",
                         "Butterfly Emergence [Larvae OWS]",
                         "Butterfly Emergence [Pupae OWS]")) %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

theta_df_leps <- data.frame(kappa = MCMCvis::MCMCpstr(leps_fit, params = 'theta', 
                                                      func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                            kappa =MCMCvis::MCMCpstr(leps_fit, params = 'theta', 
                                                     fun = median)) %>% 
  mutate(Pheno.Phase = c("Butterfly Emergence [Egg OWS]",
                         "Butterfly Emergence [Larvae OWS]",
                         "Butterfly Emergence [Pupae OWS]")) %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

gddEffectPlot <- ggplot(gamma_df_leps) +
  geom_linerange(aes(x = median, xmin = CI_025, xmax = CI_975, y = Pheno.Phase)) +
  geom_linerange(aes(x = median, xmin = CI_25, xmax = CI_75, y = Pheno.Phase), linewidth = 1.5) +
  geom_point(aes(x = median, y = Pheno.Phase), shape = 3, size = 2.5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  labs(x = "Effect of GDD at median latitude", y = "") +
  theme_classic()
gddEffectPlot


LatEffectPlot <- 
  ggplot(theta_df_leps) +
  geom_linerange(aes(x = median, xmin = CI_025, xmax = CI_975, y = Pheno.Phase)) +
  geom_linerange(aes(x = median, xmin = CI_25, xmax = CI_75, y = Pheno.Phase), linewidth = 1.5) +
  geom_point(aes(x = median, y = Pheno.Phase), shape = 3, size = 2.5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  labs(x = "Effect of Latitude on GDD sensitivity", y = "") +
  theme_classic()

LatEffectPlot

# ggarrange
library(ggpubr)

ga <- ggarrange(gddEffectPlot, LatEffectPlot,
                labels = c("A",
                           "B"), 
                ncol = 1, nrow = 2)

ggsave(filename = "FigOutputs/LepOWSResponses_gddSens_LatSens.png", plot = ga, width = 6, height = 6, dpi = 400)
