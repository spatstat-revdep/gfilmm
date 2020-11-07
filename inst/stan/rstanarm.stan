
data {
  int<lower=1> J; // number of between variances
  ...
}

parameters {
  real<lower=0> tau;
  vector<lower=0>[J] theta;
  corr_matrix[J] Omega;
  ...
}

transformed parameters {
  real<lower=0> sigma_between_total;
  vector<lower=0>[J] sigmas_between;
  cov_matrix[J] Sigma; 
  sigma_between_total = sqrt(J) * tau; 
  Sigma = quad_form_diag(Omega, sigma_between_total * sqrt(theta));
  for(j in 1:J){
    sigmas_between[j] = Sigma[j,j];
  }
}

model {
  tau ~ gamma(shape, scale);
  theta ~ dirichlet_rng(rep_vector(concentration, J));
  Omega ~ lkj_corr(regularization);
  ...
}

