//
// This Stan program corresponds to the First Model (see master thesis) that is used for Model 1.
//

//  The data block specifies the data that is conditioned upon in Bayes Rule.
data {       
  int<lower=0> S;  // number of sites: must be positive
  real hat_tau_s[S];  // vector of estimated treatment effects
  real<lower=0> sigma_s[S];  // vector of standard errors of estimated treatment effects: components must all be positive
}

// The parameters block declares the parameters which posterior distribution is sought.
parameters {         
  real tau;  // mean of the treatment effect (over all sites)
  real<lower=0> sigma;   // standard deviation of the treatment effect (over all sites)
  vector[S] eta;    // tau_s = mu + tau * eta is the mean of the treatment effect (for each site)
}

// We use a "non-centered parameterization":
transformed parameters { 
  vector[S] tau_s;
  tau_s = tau + sigma * eta;
}

// We finally write our hierarchical model, with the priors:
model {  
  tau ~ normal(0, sqrt(5)); // hyperprior 1 : N(0,5)
  sigma ~ cauchy(0, 5); // hyperprior 2 : half-Cauchy(0,5)
  target += normal_lpdf(eta | 0, 1);    
  target += normal_lpdf(hat_tau_s | tau_s, sigma_s);
}



