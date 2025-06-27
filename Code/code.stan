data {
  // number of sites
  int<lower = 1> n;
  int<lower = 0> ncov_z;
  matrix[n, ncov_z] X_z;

  // total number of samples
  int<lower = 1> N;
  // total number of technical replicates
  int<lower = 1> N2;

  // samples per site
  int<lower = 0> M[n];
  int<lower = 0> sumM[n];
  int<lower = 0> im_idx[N2 ];
  // sample-level detection covariates
  int<lower = 0> ncov_theta;
  matrix[N, ncov_theta] X_theta;

  // PCR per marker
  int<lower = 1> K[N];
  int<lower = 0> sumK[N];

  // number of species
  int<lower = 0> S;
  // number of spike ins
  int<lower = 0> S_star;

  // survey level information
  matrix[N2, S + S_star] logy1;
  // int logy_na[N2, S];

  // int delta[N, S];

  // priors
  // real prior_beta_psi;
  // real<lower = 0> prior_beta_psi_sd;
  // real prior_beta_theta;
  // real<lower = 0> prior_beta_theta_sd;

  real a_sigma0;
  real b_sigma0;
  real a_sigma1;
  real b_sigma1;
  real a_p;
  real b_p;
  real a_q;
  real b_q;
  real a_theta0;
  real b_theta0;
  real<lower = 0> lambda_prior[S + S_star];

}

parameters {

  matrix[n, S] logl;
  matrix[N, S] v_im;
  vector[N2] u_imk;
  vector[N] o_im; // offset

  real<lower = 0> tau[S];
  real<lower = 0> sigma[S];

  vector[S + S_star] lambda;
  matrix[ncov_z, S] beta_z;
  matrix[ncov_theta, S] beta_w;
  matrix[ncov_theta, S] beta_theta;
  // vector[S] beta0_z;
  vector[S] beta0_theta;
  vector<lower=0>[S] phi;

  real<lower = 0, upper = 1> pi0;
  real<lower = 0> sigma0;

  // false positive parameter
  real<lower = 0, upper = 1> p[S + S_star];
  real<lower = 0, upper = 1> q[S + S_star];
  real<lower = 0, upper = 1> theta0[S];
  real<lower = 0> sigma1[S + S_star];

}

transformed parameters {

  // matrix[n, S] Xzbetaz = rep_matrix(beta0_z', n) + X_z * beta_z;
  matrix[n, S] Xzbetaz = X_z * beta_z;
  matrix[N, S] Xthetabetaw = X_theta * beta_w;
  matrix[N, S] Xthetabetatheta = X_theta * beta_theta;

  matrix[N, S] logit_theta;
  for (i in 1:n) {
    for(m in 1:M[i]){
      for (s in 1:S) {
        logit_theta[sumM[i] + m, s] = beta0_theta[s] + Xthetabetatheta[sumM[i] + m,s] + phi[s] * logl[i,s];
      }
    }
  }

  real log_theta0[S] = log(theta0);
  real log1m_theta0[S]  = log1m(theta0);

  matrix[N, S] log_theta = log_inv_logit(logit_theta);
  matrix[N, S] log1m_theta = log1m_inv_logit(logit_theta);

}

model {

  real log_p_ydelta1;
  real log_p_ydelta0;

  sigma0 ~ gamma(a_sigma0, b_sigma0);

  to_vector(beta_z) ~ normal(0, 10);
  to_vector(beta_w) ~ normal(0, 10);
  to_vector(beta_theta) ~ normal(0, 10);

  phi ~ normal(0, 10);
  tau ~ gamma(1, 1);
  sigma ~ gamma(1, 1);

  lambda[1:(S + S_star)] ~ normal(lambda_prior[1:(S + S_star)], 0.1);
  beta0_theta ~ normal(0, .1);
  theta0 ~ beta(a_theta0, b_theta0);
  p ~ beta(a_p, b_p);
  q ~ beta(a_q, b_q);
  sigma1 ~ gamma(a_sigma1, b_sigma1);

  o_im ~ normal(0, 10);
  u_imk ~ normal(0, 10);

  // likelihood of logl

  for(s in 1:S){

    for (i in 1:n) {
      target += normal_lpdf(logl[i, s] | Xzbetaz[i, s], tau[s]);
    }

  }

  // likelihood of v

  for(s in 1:S){

    // target += normal_lpdf(v_im[, s] | rep_matrix(logl[im_idx, s], N) + Xthetabetaw[, s], sigma[s]);

    for (i in 1:n) {

      for(m in 1:M[i]){

        target += normal_lpdf(v_im[sumM[i] + m,s] | logl[i,s] + Xthetabetaw[sumM[i] + m,s], sigma[s]);

      }

    }

  }


  ////////// precomputation of lambda_s /////////////////////////////////////////////////////

  // matrix[N2, S] lambda_simk = rep_matrix(lambda[1:S]', N2) +
    //   rep_matrix(o_im[im_idx], S) +
    //   v_im[im_idx, 1:S];

  // matrix[N2, S] lambda_simk;
  //
    // for(s in 1:S){
      //
        //   for (i in 1:n) {
          //
            //     for(m in 1:M[i]){
              //
                //       for(k in 1:K[sumM[i] + m]){
                  //
                    //         lambda_simk[sumK[sumM[i] + m] + k,s] = lambda[s] + o_im[sumM[i] + m] + v_im[sumM[i] + m,s];// + u_imk[sumK[sumM[i] + m] + k];
                    //         // real lambda_simk = lambda[s] + o_im[sumM[i] + m] + v_im[sumM[i] + m,s] + u_imk[sumK[sumM[i] + m] + k];
                    //
                      //       }
              //     }
          //   }
      // }

  // precomputation of individual likelihoods

  // matrix[N2, S] logpyc1_precomp;
  // for (i in 1:N2) {
    //   for (s in 1:S) {
      //     logpyc1_precomp[i, s] = normal_lpdf(logy1[i,s] | lambda_simk[i, s], sigma1[s]);
      //   }
    // }
  //
    // matrix[N2, S + S_star] logpyc0_precomp;
  // for (i in 1:N2) {
    //   for (s in 1:(S + S_star)) {
      //     logpyc0_precomp[i, s] =
        //     log_sum_exp(
          //       log(pi0) + normal_lpdf(logy1[i,s] | log(1), .00001),
          //       log(1 - pi0) + normal_lpdf(logy1[i,s] | 0, sigma0)
          //       );
      //   }
    // }


  /////////// loglikelihood ////////////////////////

    for(s in 1:S){

      for (i in 1:n) {

        for(m in 1:M[i]){

          // log probability of y given delta = 0
          log_p_ydelta0 = 0;

          for(k in 1:K[sumM[i] + m]){

            log_p_ydelta0 +=
              // logpyc0_precomp[sumK[sumM[i] + m] + k, s];
            log_sum_exp(
              log(pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | log(1), .00001),
              // log(pi0) + logy1_lpdf_0[sumK[sumM[i] + m] + k,s],
              // log(1 - pi0) + logy1_lpdf_1[sumM[i] + m] + k,s]
log(1 - pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | 0, sigma0)
            );

          }

          // log_p_ydelta0 = sum(logpyc0_precomp[(sumK[sumM[i] + m] + 1):(sumK[sumM[i] + m] + K[sumM[i] + m]), s]);

          // log probability of y given delta = 1
          log_p_ydelta1 = 0;

          for(k in 1:K[sumM[i] + m]){

            // real lambda_simk = lambda[s];// + o_im[sumM[i] + m] + v_im[sumM[i] + m,s];// + u_imk[sumK[sumM[i] + m] + k];
            real lambda_simk = lambda[s] + o_im[sumM[i] + m] + v_im[sumM[i] + m,s];// + u_imk[sumK[sumM[i] + m] + k];
            // real lambda_simk = lambda[s] + o_im[sumM[i] + m] + v_im[sumM[i] + m,s] + u_imk[sumK[sumM[i] + m] + k];

            real logpyc1 = normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | lambda_simk, sigma1[s]);

            real logpyc0 =
              log_sum_exp(
                log(pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | log(1), .00001),
                log(1 - pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | 0, sigma0)
              );

            log_p_ydelta1 +=
              log_sum_exp(
                logpyc1 + log(p[s]),
                // logpyc1_precomp[sumK[sumM[i] + m] + k, s] + log(p[s]),
                logpyc0 + log(1 - p[s])
                // logpyc0_precomp[sumK[sumM[i] + m] + k, s] + log(1 - p[s])
              );

          }

          target += log_sum_exp(
            log_theta[sumM[i] + m,s] + log_p_ydelta1,
            log1m_theta[sumM[i] + m,s] + log_p_ydelta0);

        }

      }

    }

  ///////////// loglikelihood of spikeins ///////////////////

    for(s in 1:S_star){

      for (i in 1:n) {

        for(m in 1:M[i]){

          for(k in 1:K[sumM[i] + m]){

            real lambda_simk_spike = lambda[S + s] + o_im[sumM[i] + m];// + u_imk[sumK[sumM[i] + m] + k];
            // real lambda_simk = lambda[S + s] + o_im[sumM[i] + m] + u_imk[sumK[sumM[i] + m] + k];

            real logpyc1 = normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | lambda_simk_spike, sigma1[S + s]);

            real logpyc0 =
              // logpyc0_precomp[sumK[sumM[i] + m] + k, S + s];
            log_sum_exp(
              log(pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | log(1), .00001),
              log(1 - pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | 0, sigma0)
            );

            target +=
              // normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | lambda_simk_spike, sigma1[S + s]);
            // normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | lambda_simk_spike, .1);
            // logpyc1;
            log_sum_exp(
              logpyc1 + log(p[S + s]),
              logpyc0 + log(1 - p[S + s])
            );

          }

        }
      }
    }




  // for (s in 1:S_star) {
    //   // Precompute lambda_simk for all relevant indices
    //   vector[N2] lambda_simk;
    //   int pos = 1;
    //
      //   for (i in 1:n) {
        //     for (m in 1:M[i]) {
          //       int im_index = sumM[i] + m;
          //       for (k in 1:K[im_index]) {
            //         lambda_simk[pos] = lambda[S + s] + o_im[im_index]; // + u_imk[sumK[im_index] + k] if needed
            //         pos += 1;
            //       }
          //     }
        //   }
    //
      //   // Vectorized computations
    //   vector[N] logpyc1 = normal_lpdf(logy1[1:sumK[sumM[n] + M[n]], S + s] | lambda_simk, sigma1[S + s]);
    //   vector[sumK[sumM[n] + M[n]]] logpyc0 = logpyc0_precomp[1:sumK[sumM[n] + M[n]], S + s];
    //
      //   // Vectorized log_sum_exp
    //   vector[sumK[sumM[n] + M[n]]] log_terms = log_sum_exp(
      //     logpyc1 + log(p[S + s]),
      //     logpyc0 + log(1 - p[S + s])
      //     );
    //
      //     // Accumulate the sum
    //     target += sum(log_terms);
    // }

}
