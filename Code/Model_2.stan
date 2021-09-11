//
// This Stan program corresponds to the Second Model (see master thesis) that is used for Model 2 and Model 2bis.
//

//  The data block specifies the data that is conditioned upon in Bayes Rule.
data {
  int<lower=0> N;    // number of  individuals
  int<lower=1> I;    // number of individual predictors
  int<lower=1> S;    // number of site
  int<lower=1> J;    // number of site predictors
  int<lower=1,upper=S> site[N];  // site of individual i
  matrix[N, I] X;    // individual-level predictors
  matrix[J, S] tZ;   // site-level predictors transposed
  vector[N] y;       // individual-level outcomes
}

// The parameters block declares the parameters which posterior distribution is sought.
parameters {
  matrix[I, S] u;
  cholesky_factor_corr[I] L_Omega;  // Cholesky factor for Omega
  vector<lower=0,upper=pi()/2>[I] theta_unif;  // uniform distribution used to build the Cauchy priors for theta
  matrix[I, J] gamma;   // site-level coefficients
  real<lower=0> sigma_s[S];  
}

// We use a "non-centered parameterization" for theta and
// a "Cholesky parametrization" for Omega:
transformed parameters {
  vector<lower=0>[I] theta = 2.5 * tan(theta_unif); //Cauchy priors for theta
  matrix[I, S] beta = gamma * tZ + diag_pre_multiply(theta, L_Omega) * u;
}

// We finally write our hierarchical model, with the priors:
model {
  vector[N] mu;
  vector[N] sigma; // variable used to merge the sigma_s of different sites
  sigma_s ~ uniform(0,100000);  // prior for sigma_s
  for(n in 1:N) {
    mu[n] = X[n, ] * beta[, site[n]];
    sigma[n] = sigma_s[site[n]];
  }
  to_vector(u) ~ std_normal(); 
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, sqrt(5)); // N(0,5)
  y ~ normal(mu, sigma);
}

