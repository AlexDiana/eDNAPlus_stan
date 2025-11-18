library(here)
library(rjags)
library(rstan)
library(coda)
library(tidyverse)
library(ggplot2)
library(sf)
library(patchwork)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# setwd("C:/Users/Alex/Google Drive/R Folder/Lecturership/ednaplus")

source(here("Code","function_pb.R"))


cfvars = read.csv(here("Data","allsitevars_CF.csv"))
#remove bat roost & visit that hasn't happpened yet
# cfvars = cfvars %>% filter(Date < Sys.time()) %>% filter(Name != 'Bat roost')
cfvars = cfvars %>% filter(Visit != 5) %>% filter(Name != 'Bat roost')
# add technical replicates to dataset
# cfvars_rep <- cfvars %>%
#   uncount(weights = 4, .id = "TechRep")

# scale the covariates.
# scale distance
cfvars$dist_m_scale = scale(cfvars$dist_m)#[1]

# data settings
n_s <- 1 # sites
t <- 4 # time points
L <- 19 # locations per site
# n_t = 3 #number of temporal replicates at each site
# di <- rexp(n_s * L , rate = 1) # distance from source
n_ecol <- n_s * t * L

S <- 3
S_star <- 0 #spike-in control species
M <- rep(3, n_ecol) #temporal replicates per sample
N <- sum(M) #total temporal replicates
K <- rep(4, N) #PCR replicates per replicate sample
ncov_theta <- 2 #predictors for occupancy/detection part of model

# parameters
tau_true <- rep(.5, S) #this determines the variation or level of noise within species, between sites
phi_true <- rep(1, S) # response of the occupancy probability to the biomass
sigma_true <- rep(0.02, S) #variation across samples
sigma_y_true <- rep(0.5, S) #variation across samples
beta0_theta_true <- rep(3, S) # baseline detection rate
sigma_u_true <- 0.05 # PCR noise
lambda_true <- rnorm(S + S_star, mean = 20, sd = 1)
p_true <- c(rep(.95, S), rep(1, S_star))

q_true <- c(rep(0, S), rep(1, S_star)) # prob of false positive
# pi0_true <- .95
# lambda0 determines the number of reads you get, when there is a false positive.
# This can be played around with a bit, but depends on what our threshold is for removing samples with low read numbers.
# lambda0_true <- 0.05

mu0_true <- 1
sd0_true <- 1

list_simdata <- simulateData(cfvars, n_s, t, L,
                             S, S_star, M, N, K,
                             ncov_theta,
                             tau_true, sigma_true, sigma_y_true,
                             phi_true, beta0_theta_true,
                             lambda_true, p_true, q_true,
                             sigma_u_true, mu0_true, sd0_true)

stan_data <- list_simdata$stan_data
params <- list_simdata$params

stan_data$mu_mu0 <- 1
stan_data$sd_mu0 <- 1

# MODEL --------

stan_model_compiled <- stan_model(file = here("Code","code_new.stan"))

sampling <- T

# init_fun <- function() {
#   list(
#     lambda = stan_data$lambda_prior
#   )
# }

if(sampling){
  results_stan <-
    rstan::sampling(
      stan_model_compiled,
      data = stan_data,
      pars = c("beta_z", "beta_theta",
               # "logl","v_im",
               "lambda","beta0_theta",
               "tau","sigma","sigma_y",
               "p","q",
               # "logit_p","logit_q",
               "phi","mu0","sigma0"
      ),
      # init = init_fun,
      chains = 1,
      iter = 3000)
} else {
  results_stan <-
    rstan::vb(
      stan_model_compiled,
      data = stan_data,
      pars = c("beta_z", "beta_theta","logl","v_im","lambda","beta0_theta",
               "tau","sigma","sigma_y","u_imk",
               "p","q",
               # "logit_p","logit_q",
               "phi","mu0","sigma0"
      ),
      # init = init_fun,
      algorithm = "meanfield",
      iter = 55000,
      # elbo_samples = 500,
      # adapt_engaged = F,
      tol_rel_obj = 0.00001,
      output_samples = 500,
      eta = 2)
}

# RESULTS --------

matrix_results <- as.matrix(results_stan)

matrix_results2 <- matrix_results %>%
  as.data.frame

texts <- colnames(matrix_results2)

# effect of distance
{
  beta_z_results <- matrix_results %>%
    as.data.frame

  # pattern <- sprintf("\\[%d,", L + 2)
  pattern <- sprintf("\\[%d,", 2)
  selected <- grepl(pattern, texts) & grepl("beta_z", texts)
  # on each column, calculate the 95% credible intervals
  beta_z_results <- beta_z_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  # beta_z_results$True <- as.vector(params$beta_z[L+2,])
  beta_z_results$True <- as.vector(params$beta_z[2,])

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  distanceef = ggplot(beta_z_results, aes(x = 1:nrow(beta_z_results),
                             y = True,
                             ymin = `2.5%`,
                             ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
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

# effect of season
{
  beta_z_results <- matrix_results %>%
    as.data.frame
  pattern <- sprintf("\\[%d,", 1)
  selected <- grepl(pattern, texts) & grepl("beta_z", texts)
  # on each column, calculate the 95% credible intervals
  beta_z_results <- beta_z_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  beta_z_results$True <- as.vector(params$beta_z[1,])

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  seasonef = ggplot(beta_z_results, aes(x = 1:nrow(beta_z_results),
                                          y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("Season")

  seasonef
}

# phi
{
  phi_results <- matrix_results %>%
    as.data.frame
  selected <- grepl("phi", texts)
  # on each column, calculate the 95% credible intervals
  phi_results <- phi_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  phi_results$True <- as.vector(params$phi)

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")

  phi_plot = ggplot(phi_results, aes(x = 1:nrow(phi_results),
                                          y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("Season")

  phi_plot
}

# sigma_y
{
  sigma_y_results <- matrix_results %>%
    as.data.frame
  selected <- grepl("sigma_y", texts)
  # on each column, calculate the 95% credible intervals
  sigma_y_results <- sigma_y_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  sigma_y_results$True <- as.vector(params$sigma_y)

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  sigmay_plot = ggplot(sigma_y_results, aes(x = 1:nrow(sigma_y_results),
                                          y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("Sigma y")

  sigmay_plot
}

# tau
{
  tau_results <- matrix_results %>%
    as.data.frame
  selected <- grepl("tau", texts)
  # on each column, calculate the 95% credible intervals
  tau_results <- tau_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  tau_results$True <- as.vector(params$tau)

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  tau_plot = ggplot(tau_results, aes(x = 1:nrow(tau_results),
                                          y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("Sigma y")

  tau_plot
}

# lambda
{
  lambda_results <- matrix_results %>%
    as.data.frame
  selected <- grepl("lambda", texts)
  # on each column, calculate the 95% credible intervals
  lambda_results <- lambda_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  lambda_results$True <- as.vector(params$lambda)

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  lambda_plot = ggplot(lambda_results, aes(x = 1:nrow(lambda_results),
                                          y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("Sigma y")

  lambda_plot
}

# p
{
  p_results <- matrix_results %>%
    as.data.frame
  selected <- grepl("p", texts) & !grepl("phi", texts) & !grepl("lp", texts)
  # on each column, calculate the 95% credible intervals
  p_results <- p_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  p_results$True <- as.vector(params$p)

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  p_plot = ggplot(p_results, aes(x = 1:nrow(p_results),
                                          y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("p")

  p_plot
}

# q
{
  q_results <- matrix_results %>%
    as.data.frame
  selected <- grepl("q", texts)
  # on each column, calculate the 95% credible intervals
  q_results <- q_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  q_results$True <- as.vector(params$q)

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  q_plot = ggplot(q_results, aes(x = 1:nrow(q_results),
                                          y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("Sigma y")

  q_plot
}

# u_imk
{
  uimk_results <- matrix_results %>%
    as.data.frame
  selected <- grepl("u", texts)
  # on each column, calculate the 95% credible intervals
  uimk_results <- uimk_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  uimk_results$True <- as.vector(params$u)

  species_names <- c("Brown Long-eared",
                     "Greater Horsehoe",
                     "Lesser Horseshoe")
  q_plot = ggplot(q_results, aes(x = 1:nrow(q_results),
                                          # y = True,
                                          ymin = `2.5%`,
                                          ymax = `97.5%`)) +
    geom_errorbar() +
    # geom_point(aes(colour = "Latent (true) effect"), size = 3) +
    scale_colour_manual(name = "", values = c("Latent (true) effect" = "darkgreen")) + # legend text + colour
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S, labels = species_names) +
    xlab('Species') +
    ylab('Effect of season on DNA biomass') +
    theme_classic() +
    theme(axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title =element_text(size = 14),
          legend.text = element_text(size = 14)) +
    ggtitle("Sigma y")

  q_plot
}

View(t(apply(matrix_results, 2, function(x) quantile(x, probs = c(0.025, 0.975)))))

distanceef+seasonef + plot_layout(guides = 'collect')

# logl is the predicted log DNA biomass of a site. there is one prediction per site for each time interval.
# map of DNA biomass

logl_results = matrix_results %>%
    as.data.frame()

texts <- colnames(logl_results)
loglsel <- grepl('logl', texts)
# for each species, there are 76 estimates of logl
# this will be in the order of 1-19 for visit 1, visit 2 etc, or visit 1-4 for site 1, site 2 etc.
# it is the same format as the original data (cfvars)

logl_results = logl_results[,loglsel] %>%
  apply(., 2, function(x) c(
    mean = mean(x),
    quantile(x, probs = c(0.025, 0.975)))) %>% t %>%
    as.data.frame()

logl_results$species = rep(1:S, each = n_s * t * L)
cfvars_sub = cfvars %>% group_by(Name,Visit) %>% filter(row_number() == 1)
#merge cfvarssub with loglresults (each needs to repeat)
cfvars_species = do.call(rbind, replicate(3,cfvars_sub, simplify = FALSE))
logl_results = cbind(logl_results, cfvars_species)

logl_sf = st_as_sf(logl_results, coords = c("Longitude", "Latitude"), crs = 4326)
ggplot(data = logl_sf) +
    geom_sf(aes(size = mean), alpha = 0.6, colour = 'grey') +
    theme_minimal() + facet_grid(species ~ Visit)


# plot DNA biomass over distance and time
species_labels <- c('1' = "Brown Long-eared",
                   "2" = "Greater Horsehoe",
                    "3" = "Lesser Horseshoe")

# relabel visit
logl_sf$Visit <- factor(
  logl_sf$Visit,
  levels = 1:4,
  labels = c("March", "May", "June", "September")
)
library(ggh4x)
ggplot(data = logl_sf, aes(x = dist_m, y = mean, ymax = `97.5%`, ymin = `2.5%`,
                           fill = Visit,
                           colour = Visit)) +
  geom_smooth() +
  ggh4x::facet_grid2(species~Visit, labeller = labeller(species = species_labels), scales = 'free_x', independent = "x") +
  labs(fill = 'Visit', colour = 'Visit', x = 'Distance from roost (m)',
       y = 'DNA Biomass') +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
        )

# logl is a proxy for how much biomass is at the site. but we don't know amplification rates, we don't
# know shedding rate. so its on an abstract scale because we cannot actually measure how much DNA is in the area.

# ------

# oim_CI <- matrix_results %>%
#   as.data.frame %>%
#   select(starts_with("o_im")) %>%
#   # apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
#   apply(., 2, function(x) quantile(x, probs = .5)) %>%
#   as.data.frame
#
# logl_CI <- matrix_results %>%
#   as.data.frame %>%
#   select(starts_with("logl")) %>%
#   apply(., 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) %>% t %>%
#   # apply(., 2, function(x) quantile(x, probs = .5)) %>%
#   as.data.frame
#
# logl_CI_mat <- array(NA, dim = c(3, n, S))
# for (i in 1:n) {
#   for (s in 1:S) {
#     logl_CI_mat[1,i,s] <- logl_CI$`2.5%`[(s-1)*n + i]
#     logl_CI_mat[2,i,s] <- logl_CI$`50%`[(s-1)*n + i]
#     logl_CI_mat[3,i,s] <- logl_CI$`97.5%`[(s-1)*n + i]
#   }
# }
#
# s <- 4
#
# colnames(OTU)[s]
#
# data.frame(
#   log_CI_min = logl_CI_mat[1,,s],
#   log_CI_median = logl_CI_mat[2,,s],
#   log_CI_max = logl_CI_mat[3,,s],
#   Dist = sitesInfo$DistFromPAedge) %>%
#   filter(!is.na(Dist)) %>%
#   ggplot(aes(x = Dist,
#              y = log_CI_median,
#              ymin = log_CI_min,
#              ymax = log_CI_max)) + geom_line() + geom_point(size = 2)+ geom_errorbar() + theme_bw()
# # logl_CI_mat[,1] <-
#
# sitesInfo <- data %>% group_by(site) %>% summarise(DistFromPAedge = DistFromPAedge[1])
#
# data$DistFromPAedge
#
#
#
# qplot(oim_CI$., o_true)
#
# # ggplot(data = oim_CI, aes(x = 1:N,
# #                           y = o_true,
# #                           ymin = `2.5%`,
# #                           ymax = `97.5%`)) + geom_errorbar() + geom_point()
#
# # traceplot(results_stan, pars = c("beta_z[1,1]","beta_theta[1,1]"))
#
#
#
#
#
# #
#
#
# qplot(rep(o_true, K), logy1[,S+1])
#
# matrix_results %>%
#   as.data.frame %>%
#   select(starts_with("phi")) %>%
#   apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975)))
