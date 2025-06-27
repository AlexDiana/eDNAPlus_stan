
logistic <- function(x) 1 / (1 + exp(-x))


simulateData <- function(n, S, S_star, M, N, K,
                         ncov_z, ncov_theta,
                         tau, sigma, phi, beta0_theta,
                         lambda, p, q, sigma_u, pi0, lambda0){


  N2 <- sum(K)
  im_idx <- rep(1:n, M)

  sumM <- c(0, cumsum(M)[-n])
  sumK <- c(0, cumsum(K)[-N])

  X_z <- matrix(rnorm(n * ncov_z), n, ncov_z)
  X_theta <- matrix(rnorm(N * ncov_theta), N, ncov_theta)

  beta_z_true <- matrix(sample(c(-1,1), ncov_z * S, replace = T), ncov_z, S)
  beta_theta0_true <- rep(3, S)
  beta_theta_true <- matrix(sample(c(-1,1), ncov_theta * S, replace = T), ncov_theta, S)
  beta_w_true <- matrix(sample(c(-1,1), ncov_theta * S, replace = T), ncov_theta, S)

  lambda0_true <- 7

  # simulate logl
  if(S > 0){
    logl_true <- X_z %*% beta_z_true + sapply(seq_len(S), function(s) rnorm(n, sd = tau[s]))
  } else {
    logl_true <- matrix(0, n, S)
  }

  theta_true <- logistic(
    matrix(beta0_theta, N, S, byrow = T) + X_theta %*% beta_theta_true + logl_true[im_idx,] * matrix(phi, N, S, T)
  )

  delta_true <- t(sapply(1:N, function(i){
    sapply(1:S, function(s) {
      rbinom(1, 1, theta_true[i,s])
    })
  }))

  v_true <- X_theta %*% beta_w_true + logl_true[im_idx,] + sapply(1:S, function(s) rnorm(N, sd = sigma[s]))
  v_true[delta_true == 0] <- NA

  v_spike <- matrix(0, N, S_star)

  v_all <- cbind(v_true, v_spike)

  o_true <- - log(apply(v_all, 1, function(x){
    sum(exp(x[!is.na(x)]))
  }))

  # observation process 1
  if(T) {

    u_true <- rnorm(N2, sd = sigma_u)

    yimk_true <- matrix(NA, N2, S + S_star)
    for (s in 1:S) {
      for (i in 1:n) {
        for (m in 1:M[i]) {
          for (k in 1:K[sumM[i] + m]) {

            if(delta_true[sumM[i] + m, s] == 1){
              cimk <- rbinom(1, 1, p[s])
            } else {
              cimk <- 0
            }

            if(cimk == 1){
              lambda_simk <- exp(lambda[s] + v_true[sumM[i] + m, s] + o_true[sumM[i] + m] + u_true[sumK[sumM[i] + m] + k])
              yimk_true[sumK[sumM[i] + m] + k, s] <- rpois(1, lambda_simk)
            } else {
              yimk_true[sumK[sumM[i] + m] + k, s] <-
                rbinom(1, 1, 1 - pi0) * rpois(1, lambda0)
            }

          }
        }
      }
    }

    for (s in 1:S_star) {
      for (i in 1:n) {
        for (m in 1:M[i]) {
          for (k in 1:K[sumM[i] + m]) {

            cimk <- rbinom(1, 1, p[S + s])

            if(cimk == 1){
              lambda_simk <- exp(lambda[S + s] + o_true[sumM[i] + m] + u_true[sumK[sumM[i] + m] + k])
              yimk_true[sumK[sumM[i] + m] + k, S + s] <- rpois(1, lambda_simk)
            } else {
              yimk_true[sumK[sumM[i] + m] + k, S + s] <-
                rbinom(1, 1, 1 - pi0) * rpois(1, lambda0)
            }



          }
        }
      }
    }
  }

  # multinomial observation process
  if(F) {
    sum_yimk <- rpois(N2, lambda = 10^5)

    yimk_true <- matrix(NA, N2, S + S_star)
    for (i in 1:n) {
      for (m in 1:M[i]) {
        for (k in 1:K[sumM[i] + m]) {

          cimk_all <- rep(NA, S + S_star)
          for (s in 1:(S + S_star)) {
            cimk_all[s] <- rbinom(1, 1, p[s])
          }

          v_current <- c(v_true[sumM[i] + m, ], v_spike[sumM[i] + m, ])
          v_current[is.na(v_current)] <- -exp(10)
          probs <- lambda * exp(v_current)
          probs <- probs / sum(probs)

          yimk_true[sumK[sumM[i] + m] + k,] <- rmultinom(1, sum_yimk[sumK[sumM[i] + m] + k], probs)

          # if(cimk == 1){
          #   lambda_simk <- exp(lambda_true[s] + v_true[sumM[i] + m, s] + o_true[sumM[i] + m] + u_true[sumK[sumM[i] + m] + k])
          #   yimk_true[sumK[sumM[i] + m] + k, s] <- rpois(1, lambda_simk)
          # } else {
          #   yimk_true[sumK[sumM[i] + m] + k, s] <-
          #     rbinom(1, 1, 1 - pi0_true) * rpois(1, lambda0_true)
          # }

        }
      }
    }

  }

  logy1 <- log(yimk_true + 1)

  lambda_start <- log(apply(yimk_true, 2, function(x){
    mean(x[x > 10])
  }))

  stan_data <- list(n = n,
                    ncov_z = ncov_z,
                    X_z = X_z,
                    N = N,
                    N2 = N2,
                    M = M,
                    sumM = sumM,
                    ncov_theta = ncov_theta,
                    X_theta = X_theta,
                    K = K,
                    im_idx = rep(1:N, K),
                    sumK = sumK,
                    S = S,
                    S_star = S_star,
                    logy1 = logy1,
                    a_sigma0 = 1,
                    b_sigma0 = 1,
                    a_sigma1 = .1,
                    b_sigma1 = .1,
                    a_p = 10,
                    b_p = 1,
                    a_q = 1,
                    b_q = 10,
                    a_theta0 = 1,
                    b_theta0 = 10,
                    lambda_prior = lambda_start
  )

}
