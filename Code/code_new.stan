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
  // int<lower = 0> im_idx[N2 ];
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
  vector[N] o_im; // offset
  // int logy_na[N2, S];

   int<lower = 0> biomassInSample[N];

  // priors
  // real prior_beta_psi;
  // real<lower = 0> prior_beta_psi_sd;
  // real prior_beta_theta;
  // real<lower = 0> prior_beta_theta_sd;

  real mu_mu0;
  real sd_mu0;
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

  real<lower = 0> tau[S];
  real<lower = 0> sigma[S];
  real<lower = 0> sigma_y[S + S_star];

  vector[S + S_star] lambda;
  matrix[ncov_z, S] beta_z;
  matrix[ncov_theta, S] beta_w;
  matrix[ncov_theta, S] beta_theta;
  // vector[S] beta0_z;
  vector[S] beta0_theta;
  vector<lower=0>[S] phi;

  // real<lower = 0, upper = 1> pi0;
  real mu0;
  real<lower = 0> sigma0;

  // false positive parameter
  // real<lower = 0, upper = 1> p[S + S_star];
  // real<lower = 0, upper = .2> q[S + S_star];
  real logit_p[S + S_star];
  real logit_q[S + S_star];
  real<lower = 0, upper = 1> theta0[S];

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


  mu0 ~ normal(mu_mu0, sd_mu0);
  phi ~ normal(0, 10);
  tau ~ gamma(1, 1);
  sigma ~ gamma(1, 1);

  lambda[1:(S + S_star)] ~ normal(lambda_prior[1:(S + S_star)], 0.1);
  beta0_theta ~ normal(0, .1);
  theta0 ~ beta(a_theta0, b_theta0);
  // p ~ beta(a_p, b_p);
  // q ~ beta(a_q, b_q);
  logit_p ~ normal(2.2, 1);
  logit_q ~ normal(-2.2, 1);
  sigma_y ~ gamma(a_sigma1, b_sigma1);

  u_imk ~ normal(0, 10);

  to_vector(v_im) ~ normal(0, 10);

  // likelihood of logl

  for(s in 1:S){

    for (i in 1:n) {
      target += normal_lpdf(logl[i, s] | Xzbetaz[i, s], tau[s]);
    }

  }

   for(s in 1:S){

      for (i in 1:n) {

        for(m in 1:M[i]){

          target += normal_lpdf(v_im[sumM[i] + m,s] | logl[i,s] + Xthetabetaw[sumM[i] + m,s], sigma[s]);

        }
      }
   }

    /////////// loglikelihood ////////////////////////

    for(s in 1:S){

      for (i in 1:n) {

        for(m in 1:M[i]){

          // log probability of y given delta = 0
          log_p_ydelta0 = 0;

          for(k in 1:K[sumM[i] + m]){

            if(logy1[sumK[sumM[i] + m] + k,s] > 0){

              // log_p_ydelta0 += log(q[s]) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | mu0, sigma0);
              log_p_ydelta0 += log_inv_logit(logit_q[s]) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | mu0, sigma0);

            } else {

              log_p_ydelta0 += log1m_inv_logit(logit_q[s]);
              // log_p_ydelta0 += log(1 - q[s]);

            }

            // log_p_ydelta0 += normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | log(1), .00001);

          }

          // log probability of y given delta = 1
          log_p_ydelta1 = 0;

          // log_p_ydelta1 += normal_lpdf(v_im[sumM[i] + m,s] | logl[i,s] + Xthetabetaw[sumM[i] + m,s], sigma[s]);

          if(biomassInSample[sumM[i] + m]){

            for(k in 1:K[sumM[i] + m]){

              // real lambda_simk = lambda[s];// + o_im[sumM[i] + m] + v_im[sumM[i] + m,s];// + u_imk[sumK[sumM[i] + m] + k];
              real lambda_simk = lambda[s] + o_im[sumM[i] + m] + v_im[sumM[i] + m,s] + u_imk[sumK[sumM[i] + m] + k];
              // real lambda_simk = lambda[s] + o_im[sumM[i] + m] + v_im[sumM[i] + m,s] + u_imk[sumK[sumM[i] + m] + k];

              // real logpyc1 = normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | lambda_simk, sigma_y[s]);

              if(logy1[sumK[sumM[i] + m] + k,s] > 0){

                // log_p_ydelta1 += log(p[s]) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | lambda_simk, sigma_y[s]);
                log_p_ydelta1 += log_inv_logit(logit_p[s]) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | lambda_simk, sigma_y[s]);

              } else {

                log_p_ydelta1 += log1m_inv_logit(logit_p[s]);;
                // log_p_ydelta1 += log(1 - p[s]);

              }

              // real logpyc0 =
              // log_sum_exp(
                //   log(pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | log(1), .00001),
                //   log(1 - pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,s] | 0, sigma0)
                //   );

                // log_p_ydelta1 +=
                // log_sum_exp(
                  //   logpyc1 + log(p[s]),
                  //   // logpyc1_precomp[sumK[sumM[i] + m] + k, s] + log(p[s]),
                  //   logpyc0 + log(1 - p[s])
                  //   // logpyc0_precomp[sumK[sumM[i] + m] + k, s] + log(1 - p[s])
                  //   );

            }

          } else {

            log_p_ydelta1 = -1000000;

          }

          target += log_sum_exp(
            log_theta[sumM[i] + m,s] + log_p_ydelta1,
            log1m_theta[sumM[i] + m,s] + log_p_ydelta0);

        }

      }

    }

    ///////////// loglikelihood of spikeins ///////////////////

    // for(s in 1:S_star){
      //
      //   for (i in 1:n) {
        //
        //     for(m in 1:M[i]){
          //
          //       for(k in 1:K[sumM[i] + m]){
            //
            //         real lambda_simk_spike = lambda[S + s] + o_im[sumM[i] + m];// + u_imk[sumK[sumM[i] + m] + k];
            //         // real lambda_simk = lambda[S + s] + o_im[sumM[i] + m] + u_imk[sumK[sumM[i] + m] + k];
            //
            //         real logpyc1 = normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | lambda_simk_spike, sigma_y[S + s]);
            //
            //         real logpyc0 =
            //         // logpyc0_precomp[sumK[sumM[i] + m] + k, S + s];
            //         log_sum_exp(
              //           log(pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | log(1), .00001),
              //           log(1 - pi0) + normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | 0, sigma0)
              //           );
              //
              //           target +=
              //           // normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | lambda_simk_spike, sigma1[S + s]);
              //           // normal_lpdf(logy1[sumK[sumM[i] + m] + k,S + s] | lambda_simk_spike, .1);
              //           // logpyc1;
              //           log_sum_exp(
                //             logpyc1 + log(p[S + s]),
                //             logpyc0 + log(1 - p[S + s])
                //             );
                //
                //       }
                //
                //     }
                //   }
                // }


}
