
logistic <- function(x) 1 / (1 + exp(-x))


# modifying the simulate data function to incorportate my existing data

# this function takes in a list of parameters:
# N_s = number of sites
#t = number of time points, l = number of locations per site, di = distance covariates,
# S = number of species, S_star = number of spike-ims, M = X technical replicates per site, N = total number of replicate samples (sum(M)),
# K = sequencing replicates per replicate sample.
# ncov_theta,
# tau, sigma, phi, beta0_theta,
# lambda, p, q, sigma_u, pi0, lambda0
simulateData <- function(data, n_s, t, L, n_t,
                         S, S_star, M, N, K,
                         ncov_theta,
                         tau, sigma, phi, beta0_theta,
                         lambda, p, q, sigma_u, pi0, lambda0){


  n_ecol <- n_s * t * L     # total ecological field samples (site × visit × location × temporal replicate)
  N <- sum(M)                # total number of temporal replicates
  N2 <- sum(K)               # total number of PCR replicates


  # --- Index helpers ---
  # site_id     <- rep(1:L, each = t * n_t)
  # visit_id    <- rep(rep(1:t, each = n_t), times = L)
  # temporal_id <- rep(1:n_t, times = L * t)

  im_idx <- rep(1:n_ecol, M)           # map field sample to technical replicate
  sumM <- c(0, cumsum(M)[-n_ecol])
  sumK <- c(0, cumsum(K)[-N])

  # create dataframe of each visit to each site
  # Build design matrix from real cfvars
  # Do i need to include temporal replicates in this matrix?
  data_sub = data %>% group_by(Name,Visit) %>% filter(row_number() == 1)
  data_sub$Visit <- scale(data_sub$Visit)
  # X_z <- model.matrix(~ factor(Name) + factor(Visit) + dist_m_scale - 1, data = data_sub)
  X_z <- model.matrix(~ factor(Name) + Visit + dist_m_scale - 1, data = data_sub)
  ncov_z = ncol(X_z)

  # Choose which covariates to include in detection
  # Once this has one row per technical replicate it will need to be subsetted again.
  X_theta <- as.matrix(data[, c("t2m_c_mean",
                                  "total_precipitation_mm", "wind10m_ms_mean")])
  X_theta = scale(X_theta)
  ncov_theta <- ncol(X_theta)

  # --- Species-level coefficients ---
  # Here we are generating coefficients i.e how each variable affects abundance of each species

  #generate coefficients of each variable on each species
  # this was previously random, but i'm adding in directions
  beta_z_true = matrix(0, ncol(X_z), S)
  # no effect of site
  beta_z_true[grep("Name", colnames(X_z)), ] <- 0
  # Positive effect of visit
  beta_z_true[grep("Visit", colnames(X_z)), ] <- 1   # or any positive value you like
  # Negative effect of distance
  beta_z_true[grep("dist_m", colnames(X_z)), ] <- -1 # adjust magnitude if you want


  # randomly generate the effects of environmental variables (e.g. temp and precipitation) on OCCUPANCY
  #i.e is the species present or not
  # At the moment this is random, but you could assign directions/intensity across species.
  beta_theta_true <- matrix(sample(c(-1,1), ncov_theta * S, replace = TRUE), ncov_theta, S)

  # generate the effects of environmental variables (e.g. temp and precipitation) on DETECTION PROBABILITY
  #i.e how strong is the signal of species presence?
  beta_w_true <- matrix(sample(c(-1,1), ncov_theta * S, replace = TRUE), ncov_theta, S)

  # --- Latent abundance ---
  # Latent state is the 'true' state - before noise from other processess are added.
  # we are randomly generating this based on our random coefficients and other parameters
  logl_true <- X_z %*% beta_z_true + sapply(1:S, function(s) rnorm(n_ecol, sd = tau[s]))
  # so now we have a latent abundance prediction for each species at each ecological sample

  # at the moment this generates negative abundances.

  # --- Occupancy probabilities and presence/absence ---

  # Copy ecological values for each technical replicate
  theta_true <- logistic(matrix(beta0_theta_true, N, S, byrow = TRUE) + #create a matrix of S species and N (tech sample) rows
                           X_theta[] %*% beta_theta_true + #matrix of covariate effects on each species in each sample
                           logl_true[im_idx,] * matrix(phi, N, S, byrow = TRUE))
  # theta_true represents occupancy probability of each species in each sample, before technical-replicate noise is added.
  delta_true <- t(sapply(1:N, function(i) {
    sapply(1:S, function(s) rbinom(1, 1, theta_true[i, s]))
  }))
  # delta_true represents actual presence absence of each species in each sample

  # --- Latent expression values ---
  # X_theta[] %*% beta_w_true is estimating how much the covariates impact DETECTION
  v_true <- X_theta %*% beta_w_true + logl_true[im_idx, ] + sapply(1:S, function(s) rnorm(N, sd = sigma[s]))
  # This adds species specific noise to each sample to simulate observation noise.
  v_true[delta_true == 0] <- NA
  # remove noise from species that are absent

  # Spike-in species
  v_spike <- matrix(0, N, S_star)
  v_all <- cbind(v_true, v_spike)

  # --- Observation offset ---
  o_true <- -log(apply(v_all, 1, function(x) sum(exp(x[!is.na(x)]))))

  # --- Technical + sequencing replicates ---
  u_true <- rnorm(N2, sd = sigma_u)
  yimk_true <- matrix(NA, N2, S + S_star)
  # yimk will be our observed data i.e read numbers, after noise from lab process has been added.

  for (s in 1:(S + S_star)) { #for each species
    for (i in 1:n_ecol) { #for each ecological sample
      for (m in 1:M[i]) { #for each technical sample of an ecological sample
        for (k in 1:K[sumM[i] + m]) { #for each sequence replicate of a technical replicate..
          if (s <= S) {
            # real species
            cimk <- ifelse(delta_true[i, s] == 1, rbinom(1, 1, p[s]), 0) #use detection probability to determine if species is present or not
            if (cimk == 1) {
              lambda_simk <- exp(lambda[s] + v_true[i, s] + o_true[i] + u_true[sumK[sumM[i] + m] + k])
              yimk_true[sumK[sumM[i] + m] + k, s] <- rpois(1, lambda_simk) #generate read count

            } else {
              yimk_true[sumK[sumM[i] + m] + k, s] <- rbinom(1, 1, 1 - pi0) * rpois(1, lambda0) #simulate occasional false positives
            }

          } else {
            # spike-in species
            cimk <- rbinom(1, 1, p[s])
            if (cimk == 1) {
              lambda_simk <- exp(lambda[s] + o_true[i] + u_true[sumK[sumM[i] + m] + k])
              yimk_true[sumK[sumM[i] + m] + k, s] <- rpois(1, lambda_simk)
            } else {
              yimk_true[sumK[sumM[i] + m] + k, s] <- rbinom(1, 1, 1 - pi0) * rpois(1, lambda0)
            }
          }
        }
      }
    }
  }

  # --- Log-transformed counts ---
  logy1 <- log(yimk_true + 1)
  lambda_start <- log(apply(yimk_true, 2, function(x) mean(x[x > 10])))

  # --- Output for Stan ---
  stan_data <- list(n = n_ecol,
                    ncov_z = ncov_z,
                    X_z = X_z,
                    N = N,
                    N2 = N2,
                    M = M,
                    sumM = sumM,
                    ncov_theta = ncov_theta,
                    X_theta = X_theta[im_idx, ],
                    K = K,
                    im_idx = im_idx,
                    sumK = sumK,
                    S = S,
                    S_star = S_star,
                    y = yimk_true,
                    logy1 = logy1,
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

  params <- list(beta_z = beta_z_true,
                 beta_theta = beta_theta_true,
                 beta_w = beta_w_true,
                 logl = logl_true)

  list("stan_data" = stan_data, "params" = params)
}




