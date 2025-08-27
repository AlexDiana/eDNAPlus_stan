library(here)
library(rjags)
library(rstan)
library(coda)
library(tidyverse)
library(ggplot2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# setwd("C:/Users/Alex/Google Drive/R Folder/Lecturership/ednaplus")

source(here("Code","function_pb.R"))


cfvars = read.csv(here("Data","allsitevars_CF.csv"))
#remove bat roost & visit that hasn't happpened yet
cfvars = cfvars %>% filter(Date < Sys.time()) %>% filter(Name != 'Bat roost')
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
tau_true <- rep(.2, S) #this determines the variation or level of noise within species, between sites
phi_true <- rep(1, S)
sigma_true <- rep(0.2, S)
beta0_theta_true <- rep(0, S)
sigma_u_true <- 1
lambda_true <- rnorm(S + S_star, mean = 9, sd = 1)
p_true <- c(rep(.95, S), rep(1, S_star))

q_true <- c(rep(.05, S), rep(1, S_star)) # prob of false positive
pi0_true <- .95
# lambda0 determines the number of reads you get, when there is a false positive.
# This can be played around with a bit, but depends on what our threshold is for removing samples with low read numbers.
lambda0_true <- 3

list_simdata <- simulateData(cfvars, n_s, t, L,
                             S, S_star, M, N, K,
                             ncov_theta,
                             tau_true, sigma_true, phi_true, beta0_theta_true,
                             lambda_true, p_true, q_true,
                             sigma_u_true, pi0_true, lambda0_true)
stan_data <- list_simdata$stan_data
params <- list_simdata$params

# MODEL --------

stan_model_compiled <- stan_model(file = here("Code","code.stan"))


sampling <- F

if(sampling){
  results_stan <-
    rstan::sampling(
      stan_model_compiled,
      data = stan_data,
      pars = c("beta_z", "beta_theta","logl","v_im","lambda",
               "sigma0","sigma","sigma1","p","phi"
      ),
      # init = init_fun,
      chains = 1,
      iter = 15000)
} else {
  results_stan <-
    rstan::vb(
      stan_model_compiled,
      data = stan_data,
      pars = c("beta_z", "beta_theta","logl","v_im","lambda",
               "sigma0","sigma","sigma1","p","phi"
      ),
      # init = init_fun,
      algorithm = "meanfield",
      iter = 15000,
      # elbo_samples = 500,
      tol_rel_obj = 0.00001,
      output_samples = 500)
}

# RESULTS --------

matrix_results <- as.matrix(results_stan)

# effect of distance
{
  beta_z_results <- matrix_results %>%
    as.data.frame

  texts <- colnames(beta_z_results)
  pattern <- sprintf("\\[%d,", L + 2)
  selected <- grepl(pattern, texts) & grepl("beta_z", texts)
  # on each column, calculate the 95% credible intervals
  beta_z_results <- beta_z_results[,selected] %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  beta_z_results$True <- as.vector(params$beta_z[L+2,])

  ggplot(beta_z_results, aes(x = 1:nrow(beta_z_results),
                             y = True,
                             ymin = `2.5%`,
                             ymax = `97.5%`)) + geom_errorbar() + geom_point() +
    # ylim(-1, 0) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    scale_x_continuous(breaks = 1:S)

}


# map of DNA biomass
{
  logl_results = matrix_results %>%
    as.data.frame()

  texts <- colnames(logl_results)
  loglsel <- grepl('logl', texts)
  # for each species, there are 76 estimates of logl
  # this will be in the order of 1-19 for visit 1, visit 2 etc, or visit 1-4 for site 1, site 2 etc.

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


  library(sf)

  logl_sf = st_as_sf(logl_results, coords = c("Longitude", "Latitude"), crs = 4326)

  ggplot(data = logl_sf) +
    geom_sf(aes(size = mean), alpha = 0.6, colour = 'grey') +
    theme_minimal() + facet_grid(species ~ Visit)
}

# plot effect of distance over time
ggplot(data = logl_sf, aes(x = dist_m, y = mean, ymax = `97.5%`, ymin = `2.5%`,
                           group = as.factor(Visit), fill = as.factor(Visit),
                           colour = as.factor(Visit))) +
  geom_smooth() +
  facet_wrap(species ~ Visit)



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
