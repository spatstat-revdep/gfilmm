// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=1> N;
  real y[N];
  int<lower=1> I;
  int<lower=1> groupID[N];
  real<lower=0> shape;
  real<lower=0> rate;
}

parameters {
  real mu;
  real alpha[I];
  real<lower=0> sigma_error;
  real<lower=0> sigma_group;
}

model {
  mu ~ normal(0, 100);
  sigma_error ~ cauchy(0, 5);
  sigma_group ~ gamma(shape, rate);
  for(i in 1:I){
    alpha[i] ~ normal(0, sigma_group);
  }
  for(k in 1:N){
    y[k] ~ normal(mu + alpha[groupID[k]], sigma_error);
  }
}

generated quantities {
  real<lower=0> sigma_total;
  sigma_total = sqrt(sigma_group^2 + sigma_error^2);
}
