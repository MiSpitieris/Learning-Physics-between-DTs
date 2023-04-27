functions {
  vector mu_fn(real mu_wk2, 
               real R, 
               int nP, 
               int nI){
                 vector[nP]  mP = rep_vector(mu_wk2,nP); 
                 vector[nI]  mI = rep_vector((1/R)*mu_wk2,nI);
                 vector[nP + nI] mu;
                 mu = append_row(mP, mI);
                 return(mu);
  }
  // Physics informed prior kernel of the WK2 model
  matrix K_wk2(vector tP, 
               vector tI,
               real rho,
               real alpha,
               real sigmaP,
               real sigmaI,
               real rho_d,
               real alpha_d,
               real R,
               real C) {
    int nP = rows(tP);
    int nI = rows(tI);
    matrix[nP + nI, nP + nI] K;
    matrix[nP, nP] KB;

// KP
for (i in 1:(nP-1)){
  K[i,i] = pow(alpha, 0.2e1);
  for (j in (i+1):nP){
    K[i,j] = exp(-pow(tP[i] - tP[j], 0.2e1) * pow(rho, -0.2e1));
    K[i,j] = pow(alpha, 0.2e1) * K[i,j];
    K[j,i] = K[i,j];
  }
  K[nP,nP] = pow(alpha, 0.2e1);
}
K[1:nP, 1:nP] = K[1:nP, 1:nP] + diag_matrix(rep_vector(pow(sigmaP, 0.2e1), nP));  

// K_delta
// press_Bias
    for (i in 1:(nP-1)){
      KB[i,i] = pow(alpha_d, 0.2e1);
      for (j in (i+1):nP){
        KB[i,j] = exp(-pow(tP[i] - tP[j], 0.2e1) * pow(rho_d, -0.2e1));
        KB[i,j] = pow(alpha_d, 0.2e1) * KB[i,j];
        KB[j,i] = KB[i,j];
      }
      KB[nP,nP] = pow(alpha_d, 0.2e1);
    }
    K[1:nP, 1:nP] = K[1:nP, 1:nP] + KB[1:nP, 1:nP];

// KPI
for (i in 1:nP){
  for (j in 1:nI){
    K[i, nP + j] = 0.1e1 / R * exp(-pow(tP[i] - tI[j], 0.2e1) * pow(rho, -0.2e1)) 
    + 0.2e1 * C * (tP[i] - tI[j])* pow(rho, -0.2e1) 
    * exp(-pow(tP[i] - tI[j], 0.2e1) * pow(rho, -0.2e1));
    K[i, nP + j] = pow(alpha, 0.2e1) * K[i, nP + j];              
  }
}

        // KIP (KIP = KPI')
        K[(nP + 1):(nP + nI), 1:nP] = K[1:nP, (nP + 1):(nP + nI)]';
     
        // KI
        for (i in 1:(nI-1)){
          K[nP + i, nP +i] = 
            pow(R, -0.2e1) * exp(-pow(tI[i] - tI[i], 0.2e1) * pow(rho, -0.2e1)) 
          + C * C * (0.2e1 * pow(rho, -0.2e1) * exp(-pow(tI[i] - tI[i], 0.2e1) 
                                                        * pow(rho, -0.2e1)) - 0.4e1 * pow(tI[i] - tI[i], 0.2e1) 
                         * pow(rho, -0.4e1) * exp(-pow(tI[i] - tI[i], 0.2e1) * pow(rho, -0.2e1)));
          K[nP + i, nP +i] = pow(alpha, 0.2e1) * K[nP + i, nP +i];
          for (j in (i+1):nI){
            K[nP + i, nP +j] = 
              pow(R, -0.2e1) * exp(-pow(tI[i] - tI[j], 0.2e1) * pow(rho, -0.2e1)) 
            + C * C * (0.2e1 * pow(rho, -0.2e1) * exp(-pow(tI[i] - tI[j], 0.2e1) 
                                                          * pow(rho, -0.2e1)) - 0.4e1 * pow(tI[i] - tI[j], 0.2e1) 
                           * pow(rho, -0.4e1) * exp(-pow(tI[i] - tI[j], 0.2e1) * pow(rho, -0.2e1)));
            K[nP + i, nP +j] = pow(alpha, 0.2e1) * K[nP + i, nP +j];
            K[nP + j, nP +i] = K[nP + i, nP +j];
          }
          K[nP + nI, nP +nI] = 
            pow(R, -0.2e1) * exp(-pow(tI[nI] - tI[nI], 0.2e1) * pow(rho, -0.2e1)) 
          + C * C * (0.2e1 * pow(rho, -0.2e1) * exp(-pow(tI[nI] - tI[nI], 0.2e1) 
                                                        * pow(rho, -0.2e1)) - 0.4e1 * pow(tI[nI] - tI[nI], 0.2e1) 
                         * pow(rho, -0.4e1) * exp(-pow(tI[nI] - tI[nI], 0.2e1) * pow(rho, -0.2e1)));
          K[nP + nI, nP +nI] = pow(alpha, 0.2e1) * K[nP + nI, nP +nI];
        }
        K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)] = K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)]
        + diag_matrix(rep_vector(pow(sigmaI, 0.2e1), nI));
        return cholesky_decompose(K);
  }
}
data { 
  int<lower=1> nP;
  int<lower=1> nI;
  int<lower=1> Ns;
  int<lower=1,upper=Ns> id[Ns];
  vector[nP] tP[Ns]; 
  vector[nI] tI[Ns];
  matrix[Ns, nP] yP;
  matrix[Ns, nI] yI;
}

transformed data {
  int<lower=1> N_tot=nP+nI;
  matrix[Ns,N_tot] y;
  y=append_col(yP, yI);
}
parameters {
  real<lower=0> R_tilde[Ns];
  real<lower=0> C_tilde[Ns];
  
  real<lower=0> mu_wk2_tilde[Ns];
  real rho_tilde[Ns];   //non-centered sd of rho
  real alpha_tilde[Ns]; //non-centered sd of alpha
  // delta
  real rho_d_tilde;   //non-centered sd of rho (delta)
  real alpha_d_tilde; //non-centered sd of marginal sd (delta)
  // noise sds
  real<lower=0,upper=10> sigmaP;
  real<lower=0,upper=20> sigmaI;
  // global parameters
  real<lower=0> rho_m;   // median of individual prior on rho 
  real<lower=0> rho_s;   //sd of individual prior on rho
  real<lower=0> alpha_m; // median of individual prior on alpha
  real<lower=0> alpha_s; //sd of individual prior on alpha
  // delta
  real<lower=0.05,upper=1> rho_d;  // length scale for delta
  real<lower=5, upper=40> alpha_d; // marginal standard deviation
  // global parameters
  real<lower=0.5,upper=2> mu_R;
  real<lower=1,upper=2> tau_R;
  real<lower=0.5,upper=2> mu_C;
  real<lower=1,upper=2> tau_C;
  real<lower=60,upper=100> mu_muWK2;
  real<lower=20> tau_muWK2;
}
transformed parameters {
  real<lower=0> mu_wk2[Ns]; // PI prior mean
  real<lower=0> rho[Ns];    // length scale
  real<lower=0> alpha[Ns]; //marginal standard deviation
  real<lower=0> R[Ns]; // physical parameters
  real<lower=0> C[Ns]; // physical parameters
  // Non-centered parameterization of individual parameters
  for (s in 1:Ns) {
    rho[s] = exp(log(rho_m) + rho_s * rho_tilde[s]);
    alpha[s] = exp(log(alpha_m) + alpha_s * alpha_tilde[s]);
    R[s] = mu_R + tau_R * R_tilde[s];
    C[s] = mu_C + tau_C * C_tilde[s];
    mu_wk2[s] = mu_muWK2 + tau_muWK2 * mu_wk2_tilde[s];
  }
}
model {
  matrix[nP + nI, nP + nI] L_K[Ns];
  vector[nP + nI] mu[Ns];
  for (s in 1:Ns) {
    L_K[s] = K_wk2(tP[s], tI[s], rho[s], alpha[s], sigmaP, sigmaI, rho_d, alpha_d, R[s], C[s]);
    mu[s] = mu_fn(mu_wk2[s], R[s], nP, nI);
  }
  // Global level
  rho_m ~ inv_gamma(2, 2);
  alpha_m ~ normal(0,20); 
  rho_s ~ normal(0, 0.5);
  alpha_s ~ normal(0, 10);
  // delta
  rho_d~ inv_gamma(2, 0.5);
  alpha_d ~ normal(0,20);

  // non centered parameterization
  rho_tilde ~ normal(0, 1); 
  alpha_tilde ~ normal(0, 1); 
  // delta
  rho_d_tilde ~ normal(0, 1);
  alpha_d_tilde ~ normal(0, 1);
  R_tilde ~ normal(0, 1);  
  C_tilde ~ normal(0, 1);
  mu_wk2_tilde ~ normal(0, 1);

  // likelihood
  for (i in 1:Ns){
    y[i] ~ multi_normal_cholesky(mu[id[i]], L_K[id[i]]);
  }
}
