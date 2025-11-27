# running simulations for different data settings

library(here)
library(rjags)
library(rstan)
library(coda)
library(tidyverse)
library(ggplot2)
library(sf)
library(patchwork)
library(future.apply)
rstan_options(auto_write = TRUE) #saves compiled models to save re-doing each time
options(mc.cores = parallel::detectCores())

stan_model_compiled <- stan_model(file = here("Code", "code.stan"))

# Experiment settings
# L_values <- seq(10, 200, length.out = 5)

# Parallel setup
plan(multisession, workers = parallel::detectCores())



# data settings
n_s <- 1 # sites
t <- 4 # time points
# L <- 19 # locations per site
# n_t = 3 #number of temporal replicates at each site
# di <- rexp(n_s * L , rate = 1) # distance from source
# n_ecol <- n_s * t * L
S <- 3
S_star <- 0 #spike-in control species
M_n = 3 #temporal replicates per sample
# M <- rep(3, n_ecol) #temporal replicates per sample
N <- sum(M) #total temporal replicates
K_n = 4
# K <- rep(4, N) #PCR replicates per replicate sample
ncov_theta <- 2 #predictors for occupancy/detection part of model

# parameters
tau_true <- rep(.05, S) #this determines the variation or level of noise within species, between sites
phi_true <- rep(0.5, S)
sigma_true <- rep(0.02, S) #variation across samples
beta0_theta_true <- rep(0, S)
sigma_u_true <- 0.05
lambda_true <- rnorm(S + S_star, mean = 9.5, sd = 0.5) #this determines the average number of reads you get for each species - note this is exponential.
p_true <- c(rep(.95, S), rep(1, S_star))
q_true <- c(rep(.05, S), rep(1, S_star)) # prob of false positive
pi0_true <- .95
# lambda0 determines the number of reads you get, when there is a false positive.
# This can be played around with a bit, but depends on what our threshold is for removing samples with low read numbers.
lambda0_true <- 0.05


source(here("Code","function_pb.R"))

L_values = c(10, 50, 100, 150, 200)
iterations <- 20

# run for each value of L
results_all = map_df(L_values, function(L){
  res = run_simulation(iterations = iterations, cfvars, n_s, t, L, S, S_star,
                 M_n, N, K_n, ncov_theta,
                 tau_true, sigma_true, phi_true,
                 beta0_theta_true, lambda_true,
                 p_true, q_true, sigma_u_true,
                 pi0_true, lambda0_true,
  )
})

# plot
results_all$True <- rep(c(0,0,0,6,3,1,7,4,2,2,1,0,-1,-5,-10), times = length(L_values))
results_all$species <- rep(c("Brown Long-eared",
                             "Greater Horsehoe",
                             "Lesser Horseshoe"), times = (t+1)*length(L_values))
results_all$variable = rep(c("Visit: March", "Visit: May", "Visit: July", "Visit: September", "Distance"), each = 3, times = length(L_values))

{
  # beta_z_results$True <- as.vector(params$beta_z[L+2,])

  distanceef = ggplot(results_all[results_all$variable=='Distance',], aes(x = species,
                                          y = True,
                                          ymin = lower,
                                          ymax = upper)) +
    geom_errorbar(aes(colour = as.factor(L))) +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen",
                                              "10" = "#1b9e77",
                                              "50" = "#d95f02"#,
                                              # "100" = "#7570b3",
                                              # "150" = "#e7298a",
                                              # "200" = "#66a61e"
    )) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    xlab('Species') +
    ylab('Effect of distance from roost on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
      ggtitle('Distance')
  distanceef

}

{
  Visits = c("March", "May", "June", "September")
  results_all$visit = factor(rep(c(Visits, NA), each = 3, times = length(L_values)), levels = c("March", "May", "June", "September"))
  seasonef = ggplot(results_all[results_all$variable!='Distance',], aes(x = visit,
                                        y = True,
                                        ymin = lower,
                                        ymax = upper)) +
    geom_errorbar(aes(colour = as.factor(L))) +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen",
                                              "10" = "#1b9e77",
                                              "50" = "#d95f02"#,
                                              # "100" = "#7570b3",
                                              # "150" = "#e7298a",
                                              # "200" = "#66a61e"
    )) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    # scale_x_continuous(breaks = 1:t, labels = Visits) +
    xlab('Visit') +
    ylab('Effect of season on DNA biomass') +
    facet_wrap(~species, scales = 'free') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14),
          strip.background = element_blank()) +
    ylim(-1,8) +
    ggtitle("Season")
  seasonef
}
