library(tidyverse)
library(cmdstanr)
library(MCMCvis)

# now read in model data
# greenup
gu_data <- readRDS("greenup-2024-07-17/greenup-data-2024-07-17.rds")
gu_fit <- readRDS("greenup-2024-07-17/greenup-fit-2024-07-17.rds")
# leps
leps_data <- readRDS("lepsEmergence-2024-07-17/lepEmergence-data-2024-07-17.rds")
leps_fit <- readRDS("lepsEmergence-2024-07-17/lepsEmergence-fit-2024-07-17.rds")
# arrival
arr_data <- readRDS("arrive-2024-07-17/arrive-data-2024-07-17.rds")
arr_fit <- readRDS("arrive-2024-07-17/arrive-fit-2024-07-17.rds")
# fledge
fledge_data <- readRDS("fledge-2024-07-17/fledge-data-2024-07-17.rds")
fledge_fit <- readRDS("fledge-2024-07-17/fledge-fit-2024-07-17.rds")

###############################################################################
##########################GAMMAS###############################################
###############################################################################
# collect gammas -- effect of GDD on phenometric at median latitude
gamma_df_gu <- data.frame(kappa = MCMCvis::MCMCpstr(gu_fit, params = 'gamma', 
                                                    func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                          kappa =MCMCvis::MCMCpstr(gu_fit, params = 'gamma', 
                                                   fun = median)) %>% 
  mutate(Pheno.Phase = "Green-up") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

mu_gamma_df_leps <- data.frame(mu_gamma = MCMCvis::MCMCpstr(leps_fit, params = 'mu_gamma', 
                                                      func = function(x) quantile(x,  probs = c(0.025,0.25,0.75,0.975))),
                            mu_gamma = MCMCvis::MCMCpstr(leps_fit, params = 'mu_gamma', 
                                                     fun = median)) %>% 
  mutate(Pheno.Phase = "Butterfly Emergence") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)



kappa_df_arrival <- data.frame(kappa = MCMCvis::MCMCpstr(arr_fit, params = 'kappa', 
                                                        func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                              kappa =MCMCvis::MCMCpstr(arr_fit, params = 'kappa', 
                                                       fun = median)) %>% 
  mutate(Pheno.Phase = "Bird Arrival") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)


kappa_df_fledge <- data.frame(kappa = MCMCvis::MCMCpstr(fledge_fit, params = 'kappa', 
                                                    func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                              kappa =MCMCvis::MCMCpstr(fledge_fit, params = 'kappa', 
                                                       fun = median))%>% 
  mutate(Pheno.Phase = "Bird Breeding") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)


## combine together
comb_df <- bind_rows(gamma_df_gu, mu_gamma_df_leps,
                     kappa_df_arrival, kappa_df_fledge)

comb_df$Pheno.Phase <- factor(comb_df$Pheno.Phase, levels = c("Bird Breeding","Bird Arrival", "Butterfly Emergence","Green-up"))


#gamma is effect of GDD on fledge at median lat
gddEffectPlot <- ggplot(comb_df) +
  geom_linerange(aes(x = median, xmin = CI_025, xmax = CI_975, y = Pheno.Phase)) +
  geom_linerange(aes(x = median, xmin = CI_25, xmax = CI_75, y = Pheno.Phase), linewidth = 1.5) +
  geom_point(aes(x = median, y = Pheno.Phase), shape = 3, size = 2.5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  labs(x = "Effect of GDD at mean latitude", y = "") +
  scale_y_discrete(labels = scales::wrap_format(10))+
  theme_classic() +
  theme(axis.text.y = element_text(hjust=0.5))
gddEffectPlot

if (file.exists("FigOutputs/")){
  ggsave(filename = "FigOutputs/gddEffect.png", gddEffectPlot, width = 5, height = 4, dpi = 350)
} else{
  dir.create("FigOutputs")
  ggsave(filename = "FigOutputs/gddEffect.png", gddEffectPlot, width = 5, height = 4, dpi = 350)
}


###############################################################################
##########################THETAS###############################################
###############################################################################
# collect thetas -- effect of latitude on GDD
theta_df_gu <- data.frame(kappa = MCMCvis::MCMCpstr(gu_fit, params = 'theta', 
                                                    func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                          kappa =MCMCvis::MCMCpstr(gu_fit, params = 'theta', 
                                                   fun = median)) %>% 
  mutate(Pheno.Phase = "Green-up") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)
mu_theta_df_leps <- data.frame(kappa = MCMCvis::MCMCpstr(leps_fit, params = 'mu_theta', 
                                                      func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                            kappa =MCMCvis::MCMCpstr(leps_fit, params = 'mu_theta', 
                                                     fun = median)) %>% 
  mutate(Pheno.Phase = "Butterfly Emergence") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

theta_df_arrival <- data.frame(theta = MCMCvis::MCMCpstr(arr_fit, params = 'mu_theta', 
                                                         func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                               theta =MCMCvis::MCMCpstr(arr_fit, params = 'mu_theta', 
                                                        fun = median)) %>% 
  mutate(Pheno.Phase = "Bird Arrival") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)
theta_df_fledge <- data.frame(kappa = MCMCvis::MCMCpstr(fledge_fit, params = 'mu_theta', 
                                                        func = function(x) quantile(x, probs = c(0.025,0.25,0.75,0.975))),
                              kappa =MCMCvis::MCMCpstr(fledge_fit, params = 'mu_theta', 
                                                       fun = median))%>% 
  mutate(Pheno.Phase = "Bird Breeding") %>% 
  rename(CI_025 = 1,
         CI_25 = 2,
         CI_75 = 3,
         CI_975 = 4,
         median = 5)

## combine together
theta_comb_df <- bind_rows(theta_df_gu, mu_theta_df_leps,
                     theta_df_arrival, theta_df_fledge)

theta_comb_df$Pheno.Phase <- factor(theta_comb_df$Pheno.Phase, levels = c("Bird Breeding","Bird Arrival", "Butterfly Emergence","Green-up"))

#theta is effect of lat on GDD sensitivity
LatEffectPlot <- 
  ggplot(theta_comb_df) +
  geom_linerange(aes(x = median, xmin = CI_025, xmax = CI_975, y = Pheno.Phase)) +
  geom_linerange(aes(x = median, xmin = CI_25, xmax = CI_75, y = Pheno.Phase), linewidth = 1.5) +
  geom_point(aes(x = median, y = Pheno.Phase), shape = 3, size = 2.5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  labs(x = "Effect of Latitude on GDD sensitivity", y = "") +
  scale_y_discrete(labels = scales::wrap_format(10))+
  theme_classic() +
  theme(axis.text.y = element_text(hjust=0.5))

LatEffectPlot

ggsave(filename = "FigOutputs/LatitudeEffect.png", LatEffectPlot, width = 5, height = 4, dpi = 350)

# combine those two together 
# ggarrange
library(ggpubr)

ga <- ggarrange(gddEffectPlot, LatEffectPlot,
                labels = c("A",
                           "B"), 
                ncol = 1, nrow = 2)

ggsave(filename = "FigOutputs/gddSens_LatSens.png", plot = ga, width = 3, height = 4, dpi = 400)

