# runmodel

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

# data
sd = read.csv("output/spatial_data_longform_withquant_speciesID.csv")
siteinfo = read.csv("output/all_site_vars2.csv")

# new sampleID
sd$SampleIDrep = paste0(sd$SampleID, "_", sd$pcr_replicate)


# for now, just do canadafarm
sd = sd %>% filter(Location == "Canada Farm")
siteinfo = siteinfo %>% filter(Location == "Canada Farm") %>% droplevels()
# need to get data into correct format for model
# this is a list of all data

# to get number of ecological samples, we need to condense the siteinfo df.
sites_ecol = siteinfo %>%
  select(Point_ID, Visit, dist_m) %>%
  distinct()


# for canada farm, there are 19 locations
L = n_distinct(siteinfo$Point_ID)
# 5 time points
t = n_distinct(siteinfo$Visit)
# 3 temporal replicates
# so n_ecol is not exactly L*T because some locations are missing replicates - AS has no data for visit 5
n_ecol = nrow(sites_ecol)

# there are 28 components of the stan data
# [1] "n" the number of samples
# "ncov_z"  number of variables?
# "X_z" matrix of samples x variables
# "N" total number of temporal replicates (n_ecol x n_days)
# "M" for each ecological sample, how many temporal reps?
# "N2" total number of pcr replicates (N x pcr replicate)
# "sumM" cumulative sum of M
# "ncov_theta" number of secondary variables
# "X_theta" dataframe of secondary variables - length of M
# "o_im" i don't think this is needed, to do with estimating variation in simulated data
# "K" number of pcr replicates per temporal replicate
# "im_idx" how to map the field sample to technical replicate
# "sumK" cumulative sum of K
# "S" number of species
# "S_star"  number of spike-ins
#  "logy1" read counts for each sample - one column per species, one row per sample.
#  "biomassInSample"
#   a_sigma0"        "b_sigma0"        "a_sigma1"        "b_sigma1"
# [22] "a_p"             "b_p"             "a_q"             "b_q"             "a_theta0"        "b_theta0"        "lambda_prior"


n = n_ecol #94
# n covariates = number of visits (minus reference value) + 1 for distance
ncov_z = (t-1) + 1
# order sites ecol by site then visit
sites_ecol = sites_ecol %>%
  arrange(Point_ID, Visit)
X_z <- model.matrix(~ factor(Visit) + dist_m - 1, sites_ecol)
# remove the first column (visit 1) as this is our reference value
X_z <- X_z[,-1]
# this is the number of samples in siteinfo
# its not quite 285 because some samples are missing.
N = nrow(siteinfo) # 273
# for each ecological sample, how many temporal reps did we take?
sites_M = siteinfo %>%
  group_by(Point_ID, Visit) %>%
  summarise(M = n()) %>%
  arrange(Point_ID, Visit)
M = sites_M$M
length(M)
# total number of PCR replicates
sites_K = sd %>%
  group_by(Point_ID, Visit, Night) %>%
  summarise(K = n_distinct(pcr_replicate)) %>%
  arrange(Point_ID, Visit, Night)
K <- sites_K$K
length(K)
N
# length K and N should be the same length.

#PCR replicates per replicate sample
N2 = n_distinct(sd$SampleIDrep) # 817
sumM <- c(0, cumsum(M)[-n])
sumK <- c(0, cumsum(K)[-N])

# needs to ordered the same as other dataframes
x_theta_df = siteinfo %>%
  select(Point_ID, Visit, Night,
         t2m_c_mean,
         total_precipitation_mm,
         wind10m_ms_mean) %>%
  arrange(Point_ID, Visit, Night)
X_theta <- as.matrix(x_theta_df[, c("t2m_c_mean",
                              "total_precipitation_mm", "wind10m_ms_mean")])
X_theta = scale(X_theta)
ncov_theta = ncol(X_theta)

# i'm not sure what this is, but it needs to be length n x N_covariates
o_im = rep(0, N)

im_idx <- rep(1:n_ecol, M)
length(im_idx)# map field sample to technical replicate

# number of species
S = n_distinct(sd$SpeciesID)
S_star = 0

# data frame of read count data.
# each column is a species, each row is a sample
# first, lets sum read numbers per sample by SpeciesID
sd_summ = sd %>%
  group_by(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  summarise(reads = sum(read_count)) %>%
  arrange(Point_ID, Visit, Night, pcr_replicate, SpeciesID) %>%
  ungroup()
sd_summ2 = pivot_wider(sd_summ, names_from = SpeciesID, values_from = reads, values_fill = 0)
logy1 = as.matrix(sd_summ2[, -(1:4)]) # remove first four columns which are not species data
# sample concentration
sq_avg = sd %>%
  group_by(Point_ID, Visit, Night, SampleID) %>%
  summarise(biomassInSample = mean(PreIndexConcentration.ng.ul., na.rm = TRUE)) %>%
  arrange(Point_ID, Visit, Night)
biomassInSample = sq_avg$biomassInSample
# the model does not support NA
biomassInSample[is.na(biomassInSample)] = 0
biomassInSample = ifelse(biomassInSample > 0, 1,0)

lambda_start <- apply(logy1, 2, function(x) mean(x[x > 1]))
lambda_start[is.na(lambda_start)] <- 1

# --- Output for Stan ---
stan_data <- list(n = n,
                  ncov_z = ncov_z,
                  X_z = X_z,
                  N = N,
                  M = M,
                  N2 = N2,
                  sumM = sumM,
                  ncov_theta = ncov_theta,
                  X_theta = X_theta[im_idx, ],
                  o_im = o_im,
                  K = K,
                  im_idx = im_idx,
                  sumK = sumK,
                  S = S,
                  S_star = S_star,
                  # y = yimk_true,
                  logy1 = logy1,
                  biomassInSample = biomassInSample,
                  a_sigma0 = 1,
                  b_sigma0 = 1,
                  a_sigma1 = 0.1,
                  b_sigma1 = 0.1,
                  a_p = 10,
                  b_p = 1,
                  a_q = 1,
                  b_q = 10,
                  a_theta0 = 1,
                  b_theta0 = 10,
                  lambda_prior = lambda_start
)
# run mcmc?
sampling = F
stan_model_compiled <- stan_model(file = here("Code","code_new.stan"))
stan_data$mu_mu0 <- 1
stan_data$sd_mu0 <- 1

if(sampling){
  results_stan <-
    rstan::sampling(
      stan_model_compiled,
      data = stan_data,
      pars = c("beta_z", "beta_theta",
               "logl","v_im",
               "lambda","beta0_theta",
               "tau","sigma","sigma_y",
               "p","q",
               # "logit_p","logit_q",
               "phi","mu0","sigma0"
      ),
      # init = init_fun,
      chains = 4,
      iter = 5000,
      cores = 4,
      control = list(max_treedepth = 15))
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

