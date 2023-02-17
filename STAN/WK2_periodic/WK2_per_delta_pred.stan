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
               real rho_d, 
               real alpha_d,
               real p,
               real sigmaP,
               real sigmaI,
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
        a=tP[i]; 
        b=tP[j];
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
        a=tP[i]; 
        b=tP[j];
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
        a=tP[i];  
        b=tI[j]; 
        K[i, nP + j] = 
        0.1e1 / R * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) + 
        0.4e1 * C * sin(pi() * (a - b) * p) * pow(rho, -0.2e1) * pi() * p * 
        cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1));
        K[i, nP + j] = pow(alpha, 0.2e1) * K[i, nP + j];              
      }
    }
    
    // KIP (KIP = KPI')
    K[(nP + 1):(nP + nI), 1:nP] = K[1:nP, (nP + 1):(nP + nI)]';
     
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
      a=tI[i];
      b=tI[j];
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
            + diag_matrix(rep_vector(pow(sigmaI, 0.2e1), nI));
     return cholesky_decompose(K);
  }
  
  matrix KP_pred(vector tP,
                 vector tP_pred,
                 real rho,
                 real alpha,
                 real p){
    int nP = rows(tP);
    int nP_pred = rows(tP_pred);
    real a;
    real b;
    matrix[nP_pred, nP] KP;
    
    for (i in 1:nP_pred){
      for (j in 1:nP){
        a=tP_pred[i];
        b=tP[j];
        KP[i,j] = exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1));
        KP[i,j] = pow(alpha, 0.2e1) * KP[i,j];
      }
    }
    return KP;
  }
  
  matrix KP_Bias(vector tP,
                 vector tP_pred,
                 real rho_d,
                 real alpha_d,
                 real p
     ){
       int nP = rows(tP);
       int nP_pred = rows(tP_pred);
       matrix[nP_pred, nP] KB;
       real a;
       real b;
       for (i in 1:nP_pred){
         for (j in 1:nP){
          a=tP_pred[i];
          b=tP[j];
          KB[i,j] = exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho_d, -0.2e1));
          KB[i,j] = pow(alpha_d, 0.2e1) * KB[i,j];
         }
       }
       return KB;
  }
  
  // KPI_pred
  matrix KPI_pred(vector tI,
                  vector tP_pred,
                  real rho,
                  real alpha,
                  real p,
                  real R,
                  real C
  ){
    int nP = rows(tP_pred);
    vector[nP] tP = tP_pred;
    int nI = rows(tI);
    matrix[nP, nI] KPI;
    real a;
    real b;

    for (i in 1:nP){
      for (j in 1:nI){
        a=tP[i]; 
        b=tI[j];
        KPI[i, j] = 
        0.1e1 / R * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) + 
        0.4e1 * C * sin(pi() * (a - b) * p) * pow(rho, -0.2e1) * pi() * p * 
        cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1));
        KPI[i, j] = pow(alpha, 0.2e1) * KPI[i, j];
      }
    }
    return KPI; 
  }
  
  // KIP_pred
  matrix KIP_pred(vector tI_pred,
                  vector tP,
                  real rho,
                  real alpha,
                  real p,
                  real R,
                  real C){
    int nI = rows(tI_pred);
    vector[nI] tI = tI_pred;
    int nP = rows(tP);
    matrix[nI, nP] KIP;
    real a;
    real b;
    for (i in 1:nI){
      for (j in 1:nP){
        a=tI[i]; 
        b=tP[j];
        KIP[i, j] = 
        0.1e1 / R * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) - 
        0.4e1 * C * sin(pi() * (a - b) * p) * pow(rho, -0.2e1) * pi() * p * 
        cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1));
        KIP[i, j] = pow(alpha, 0.2e1) * KIP[i, j];
      }
    }
    return KIP;
  }
  
  matrix KI_pred(vector tI_pred,
                 vector tI,
                 real rho,
                 real alpha,
                 real p,
                 real R,
                 real C){
    int nI = rows(tI);
    int nI_pred = rows(tI_pred);
    matrix[nI_pred, nI] KI;
    real a;
    real b;
    // KI
    for (i in 1:nI_pred){
      for (j in 1:nI){
        a=tI_pred[i]; 
        b=tI[j];
        KI[i, j] = 
        exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) * pow(R, -0.2e1) + 
        C * C * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * (a - b) * p), 0.2e1) * 
        exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(rho, -0.2e1)) * pow(rho, -0.2e1) - 
        0.4e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
        exp(-0.2e1 * pow(sin(pi() * (a - b) *       p), 0.2e1) * pow(rho, -0.2e1)) * 
        pow(rho, -0.2e1) - 0.16e2 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
        pow(cos(pi() * (a - b) * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * 
        pow(rho, -0.2e1)) * pow(rho, -0.4e1));
        KI[i, j] = pow(alpha, 0.2e1) * KI[i, j];
      }
    }
    return KI;
  }
  
  vector Pgp_pred_rng(vector tP,
                      vector tI,
                      vector yP,
                      vector yI,
                      vector tP_pred,
                      real rho,
                      real alpha,
                      real rho_d, 
                      real alpha_d,
                      real p,
                      real R,
                      real C,
                      real sigmaP,
                      real sigmaI) {
    int nP = rows(tP);
    int nI = rows(tI);
    int N = nP + nI;
    int nP_pred = rows(tP_pred);
    vector[nP_pred] f2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[nP_pred, N] qT_p;
      matrix[N, nP_pred] v_pred;
      vector[nP_pred] f2_mu;
      matrix[nP_pred, nP_pred] cov_f2;
      matrix[nP_pred, nP_pred] diag_delta;
      vector[N] y;
      y[1:nP] = yP[1:nP];
      y[(nP+1):N] = yI[1:nI];
       
      L_K =  K_wk2(tP, tI, rho, alpha, rho_d, alpha_d, p, sigmaP, sigmaI, R, C, yI);
      K_div_y = mdivide_left_tri_low(L_K, y);
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';

      qT_p[1:nP_pred, 1:nP] = KP_pred(tP, tP_pred, rho, alpha, p)+KP_Bias(tP, tP_pred, rho_d, alpha_d,p);
      qT_p[1:nP_pred, (nP+1):N] = KPI_pred(tI, tP_pred, rho, alpha, p, R, C);
      f2_mu = (qT_p * K_div_y);
      v_pred = mdivide_left_tri_low(L_K, qT_p');
      cov_f2 = KP_pred(tP_pred, tP_pred, rho, alpha, p)+KP_Bias(tP_pred, tP_pred, rho_d, alpha_d,p) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(1e-6, nP_pred));
                                                                      
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
  
  vector Igp_pred_rng(vector tP,
                      vector tI,
                      vector yP,
                      vector yI,
                      vector tI_pred,
                      real rho,
                      real alpha,
                      real rho_d, 
                      real alpha_d,
                      real p,
                      real R,
                      real C,
                      real sigmaP,
                      real sigmaI) {
    int nP = rows(tP);
    int nI = rows(tI);
    int N = nP + nI;
    int nI_pred = rows(tI_pred);
    vector[nI_pred] f2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[nI_pred, N] qT_I;
      matrix[N, nI_pred] v_pred;
      vector[nI_pred] f2_mu;
      matrix[nI_pred, nI_pred] cov_f2;
      matrix[nI_pred, nI_pred] diag_delta;
      vector[N] y;

      y[1:nP] = yP[1:nP];
      y[(nP+1):N] = yI[1:nI];
      
      L_K =  K_wk2(tP, tI, rho, alpha, rho_d, alpha_d, p, sigmaP, sigmaI, R, C, yI);
      K_div_y = mdivide_left_tri_low(L_K, y);
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';
      qT_I[1:nI_pred, 1:nP] = KIP_pred(tI_pred, tP, rho, alpha, p, R, C);
      qT_I[1:nI_pred, (nP+1):N] = KI_pred(tI_pred, tI, rho, alpha, p, R, C);
      f2_mu = (qT_I * K_div_y);
      v_pred = mdivide_left_tri_low(L_K, qT_I');
      cov_f2 = KI_pred(tI_pred, tI_pred, rho, alpha, p, R, C) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(1e-6, nI_pred));
                                                                    
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}

data {
  int<lower=1> nP;
  int<lower=1> nI;
  int<lower=1> nP_pred;
  int<lower=1> nI_pred;
  vector[nP] tP;
  vector[nI] tI;
  vector[nP_pred] tP_pred;
  vector[nI_pred] tI_pred;
  vector[nP] yP;
  vector[nI] yI;
  real<lower=0> p;
  
  // posterior samples
  int N_samples;
  vector[N_samples] rho;
  vector[N_samples] alpha;
  vector[N_samples] sigmaP;
  vector[N_samples] sigmaI;
  vector[N_samples] R;
  vector[N_samples] C;
  
  vector[N_samples] rho_d;
  vector[N_samples] alpha_d;
  
}
parameters {
}
model {
}
generated quantities {
  vector[nP_pred] f_P;
  matrix[N_samples, nP_pred] y_P;
  vector[nI_pred] f_I;
  matrix[N_samples, nI_pred] y_I;

  for(n in 1:N_samples) {
    f_P = Pgp_pred_rng(tP, tI, yP, yI, tP_pred[1:nP_pred], rho[n], alpha[n], rho_d[n], alpha_d[n], p, R[n], C[n], sigmaP[n], sigmaI[n]);
    for(np in 1:nP_pred) {
      y_P[n, np] = normal_rng(f_P[np], sigmaP[n]);
    }  
  }
  
  for(n in 1:N_samples) {
    f_I = Igp_pred_rng(tP, tI, yP, yI, tP_pred[1:nI_pred], rho[n], alpha[n], rho_d[n], alpha_d[n], p, R[n], C[n], sigmaP[n], sigmaI[n]);
    for(np in 1:nI_pred) {
      y_I[n, np] = normal_rng(f_I[np], sigmaI[n]);
    }  
  }
}



