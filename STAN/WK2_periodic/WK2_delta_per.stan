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
  matrix K_wk2(matrix tP, 
               matrix tI,
               real rho,
               real alpha,
               real p,
               real sigmaP,
               real sigmaQ,
               real rho_d,
               real alpha_d,
               real R,
               real C,
               vector yI) {
    int nP = rows(tP);
    int nI = rows(tI);
    matrix[nP + nI, nP + nI] K;
    matrix[nP, nP] KB;
    real a;
    real b;
    // KP
    for (i in 1:(nP-1)){
      K[i,i] = pow(alpha, 0.2e1);
      for (j in (i+1):nP){
        a=tP[i,1]; 
        b=tP[j,1];
        K[i,j] = exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1));
        K[i,j] = pow(alpha, 0.2e1) * K[i,j];
        K[j,i] = K[i,j];
      }
      K[nP,nP] = pow(alpha, 0.2e1);
    }
    K[1:nP, 1:nP] = K[1:nP, 1:nP] + diag_matrix(rep_vector(pow(sigmaP, 0.2e1), nP));
    
    // press_Bias
    for (i in 1:(nP-1)){
      KB[i,i] = pow(alpha_d, 0.2e1);
      for (j in (i+1):nP){
        a=tP[i,1]; 
        b=tP[j,1];
        KB[i,j] = exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho_d, -0.2e1));
        KB[i,j] = pow(alpha_d, 0.2e1) * KB[i,j];
        KB[j,i] = KB[i,j];
      }
      KB[nP,nP] = pow(alpha_d, 0.2e1);
    }
    K[1:nP, 1:nP] = K[1:nP, 1:nP] + KB[1:nP, 1:nP];
    
    // KPI
    for (i in 1:nP){
      for (j in 1:nI){
        a=tP[i,1];  
        b=tI[j,1]; 
        K[i, nP + j] = 
        0.1e1 / R * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) + 
        0.4e1 * C * sin(pi() * (a - b) * p) * pow(rho, -0.2e1) * pi() * p * 
        cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1));
        K[i, nP + j] = pow(alpha, 0.2e1) * K[i, nP + j];              
      }
    }
    
    // KIP (KIP = KPI')
    K[(nP + 1):(nP + nI), 1:nP] = K[1:nP, (nP + 1):(nP + nI)]';
   //  for (i in 1:nI){
   //    for (j in 1:nP){
   //     a=tI[i,1];  
   //     b=tP[j,1]; 
   //     K[nP + i, j] = 
   //     0.1e1 / R * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) - 
   //     0.4e1 * C * sin(pi() * (a - b) * p) * pow(rho, -0.2e1) * pi() * p * 
   //     cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1));
   //     K[nP + i, j] = pow(alpha, 0.2e1) * K[nP + i, j];
   //   }
   // }        
   // KI
    for (i in 1:(nI-1)){
     K[nP + i, nP +i] = 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(rho, -0.2e1)) * pow(R, -0.2e1) + 
     C * C * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * 0 * p), 0.2e1) * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(rho, -0.2e1)) * pow(rho, -0.2e1) - 
     0.4e1 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(rho, -0.2e1)) * pow(rho, -0.2e1) - 
     0.16e2 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     pow(cos(pi() * 0 * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * 
     pow(rho, -0.2e1)) * pow(rho, -0.4e1));
      K[nP + i, nP +i] = pow(alpha, 0.2e1) * K[nP + i, nP +i];
     for (j in (i+1):nI){
      a=tI[i,1];
      b=tI[j,1];
      K[nP + i, nP +j] = 
      exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) * pow(R, -0.2e1) + 
      C * C * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * (a - b) * p), 0.2e1) * 
      exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) * pow(rho, -0.2e1) - 
      0.4e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
      exp(-0.2e1 * pow(sin(pi() * (a - b) *       p), 0.2e1) * pow(rho, -0.2e1)) * 
      pow(rho, -0.2e1) - 0.16e2 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
      pow(cos(pi() * (a - b) * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * 
      pow(rho, -0.2e1)) * pow(rho, -0.4e1));

      K[nP + i, nP +j] = pow(alpha, 0.2e1) * K[nP + i, nP +j];
      K[nP + j, nP +i] = K[nP + i, nP +j];
     }
     K[nP + nI, nP +nI] = 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(rho, -0.2e1)) * pow(R, -0.2e1) + 
     C * C * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * 0 * p), 0.2e1) * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(rho, -0.2e1)) * pow(rho, -0.2e1) - 
     0.4e1 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(rho, -0.2e1)) * pow(rho, -0.2e1) - 
     0.16e2 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     pow(cos(pi() * 0 * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * 
     pow(rho, -0.2e1)) * pow(rho, -0.4e1));
     K[nP + nI, nP +nI] = pow(alpha, 0.2e1) * K[nP + nI, nP +nI];
     }
     K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)] = K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)]
            + diag_matrix(rep_vector(pow(sigmaQ, 0.2e1), nI));
     return cholesky_decompose(K);
  }
  
}

data {
  int<lower=1> nP;
  int<lower=1> nI;
  matrix[nP,1] tP;
  matrix[nI,1] tI;
  vector[nP] yP;
  vector[nI] yI;
  real p;
}

transformed data {
  vector[nP + nI] y = append_row(yP, yI);
}

parameters {
  real<lower=0> rho;
  real<lower=0,upper=50> alpha;
  real<lower=0> rho_d;
  real<lower=0,upper=30> alpha_d;
  real<lower=0, upper=10> sigmaP;
  real<lower=3, upper=40> sigmaQ;
  real<lower=0> mu_wk2;
  real<lower=0.5, upper=3> R;
  real<lower=0.5, upper=6> C;
}

model {
  // Chol. of PI kernel
  matrix[nP + nI, nP + nI] L_K = K_wk2(tP, tI, rho, alpha, p, sigmaP, sigmaQ,rho_d, alpha_d, R, C, yI);
  // mean vector
  vector[nP + nI] mu = mu_fn(mu_wk2, R, nP, nI);
  // priors
  rho~normal(0,1);
  alpha ~ normal(0,20);
  rho_d~normal(0,1.0/3);
  alpha_d ~ normal(0,10);
  mu_wk2 ~ normal(mean(yP), 10);
  
  y ~ multi_normal_cholesky(mu, L_K);
}

