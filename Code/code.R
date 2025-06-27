library(here)
library(rjags)
library(rstan)
library(coda)
library(tidyverse)
library(ggplot2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("C:/Users/Alex/Google Drive/R Folder/Lecturership/ednaplus")

source(here("Code","function.R"))

# data settings
n <- 100
S <- 20
S_star <- 1
M <- rep(2, n)
N <- sum(M)
K <- rep(2, N)
ncov_z <- 1
ncov_theta <- 2

# parameters
tau_true <- rep(1, S)
phi_true <- rep(1, S)
sigma_true <- rep(1, S)
beta0_theta_true <- rep(0, S)
sigma_u_true <- 1
lambda_true <- rnorm(S + S_star, mean = 12, sd = 1)
p_true <- c(rep(.95, S), rep(1, S_star))
q_true <- c(rep(.05, S), rep(1, S_star))
pi0_true <- .95
lambda0_true <- 7

stan_data <- simulateData(n, S, S_star, M, N, K,
                          ncov_z, ncov_theta,
                          tau_true, sigma_true, phi_true, beta0_theta_true,
                          lambda_true, p_true, q_true,
                          sigma_u_true, pi0_true, lambda0_true)

# MODEL --------

stan_model_compiled <- stan_model(file = here("Code","code.stan"))

results_stan <- rstan::vb(
  stan_model_compiled,
  data = stan_data,
  algorithm = "meanfield",
  pars = c("beta_z", "beta_theta","logl","o_im","v_im","lambda",
           "sigma0","sigma","sigma1","p","phi"
  ),
  # init = init_fun,
  iter = 10000,
  # elbo_samples = 500,
  tol_rel_obj = 0.001,
  output_samples = 500)

save(results_stan, file = "results_stan.rda")

# results_stan <- sampling(
#   stan_model_compiled,
#   data = stan_data,
#   chains = 1,
#   iter = 4000,   # Total iterations: 7000 (burn-in) + 20000 (sampling)
#   warmup = 2000,
#   thin = 1,       # No thinning by default, Stan often doesn't need it
#   seed = 123,
#   # init = list(inits_stan),
#   # init = 0,
#   pars = c("beta_z", "beta_theta","logl") # Parameters to extract, similar to JAGS `variable.names`
# )

matrix_results <- as.matrix(results_stan)

oim_CI <- matrix_results %>%
  as.data.frame %>%
  select(starts_with("o_im")) %>%
  # apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
  apply(., 2, function(x) quantile(x, probs = .5)) %>%
  as.data.frame

logl_CI <- matrix_results %>%
  as.data.frame %>%
  select(starts_with("logl")) %>%
  apply(., 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))) %>% t %>%
  # apply(., 2, function(x) quantile(x, probs = .5)) %>%
  as.data.frame

logl_CI_mat <- array(NA, dim = c(3, n, S))
for (i in 1:n) {
  for (s in 1:S) {
    logl_CI_mat[1,i,s] <- logl_CI$`2.5%`[(s-1)*n + i]
    logl_CI_mat[2,i,s] <- logl_CI$`50%`[(s-1)*n + i]
    logl_CI_mat[3,i,s] <- logl_CI$`97.5%`[(s-1)*n + i]
  }
}

s <- 4

colnames(OTU)[s]

data.frame(
  log_CI_min = logl_CI_mat[1,,s],
  log_CI_median = logl_CI_mat[2,,s],
  log_CI_max = logl_CI_mat[3,,s],
  Dist = sitesInfo$DistFromPAedge) %>%
  filter(!is.na(Dist)) %>%
  ggplot(aes(x = Dist,
             y = log_CI_median,
             ymin = log_CI_min,
             ymax = log_CI_max)) + geom_line() + geom_point(size = 2)+ geom_errorbar() + theme_bw()
# logl_CI_mat[,1] <-

sitesInfo <- data %>% group_by(site) %>% summarise(DistFromPAedge = DistFromPAedge[1])

data$DistFromPAedge



qplot(oim_CI$., o_true)

# ggplot(data = oim_CI, aes(x = 1:N,
#                           y = o_true,
#                           ymin = `2.5%`,
#                           ymax = `97.5%`)) + geom_errorbar() + geom_point()

# traceplot(results_stan, pars = c("beta_z[1,1]","beta_theta[1,1]"))


# beta_z
{
  beta_z_results <- matrix_results %>%
    as.data.frame %>%
    select(starts_with("beta_z")) %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  beta_z_results$True <- as.vector(beta_z_true)


  ggplot(beta_z_results, aes(x = 1:nrow(beta_z_results),
                             y = True,
                             ymin = `2.5%`,
                             ymax = `97.5%`)) + geom_errorbar() + geom_point()
}

# beta_theta
{
  beta_theta_results <- matrix_results %>%
    as.data.frame %>%
    select(starts_with("beta_theta")) %>%
    apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975))) %>% t %>%
    as.data.frame

  beta_theta_results$True <- as.vector(beta_theta_true)


  ggplot(beta_z_results, aes(x = 1:nrow(beta_z_results),
                             y = True,
                             ymin = `2.5%`,
                             ymax = `97.5%`)) + geom_errorbar() + geom_point()
}


#


qplot(rep(o_true, K), logy1[,S+1])

matrix_results %>%
  as.data.frame %>%
  select(starts_with("phi")) %>%
  apply(., 2, function(x) quantile(x, probs = c(0.025, 0.975)))
