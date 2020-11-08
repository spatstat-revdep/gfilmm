data {
  int<lower=1> N;
  real y[N];
  int<lower=1> I;
  int<lower=1> J;
  int<lower=1> PartID[N];
  int<lower=1> OperatorID[N];
  int<lower=1> InteractionID[N];
  real<lower=0> shape;
  real<lower=0> scale;
  real<lower=0> concentration;
}

parameters {
  real<lower=0> tau;
  simplex[2] theta;
  vector[I] PartA;
  vector[J] Op;
  vector[I*J] OpPart;
  real<lower=0> sigmaE;
}

transformed parameters {
  real<lower=0> sigma_between_total;
  vector<lower=0>[2] sigmas_between;
  real<lower=0> sigmaO;
  real<lower=0> sigmaOP;
  sigma_between_total = sqrt(2) * tau; 
  sigmas_between = sigma_between_total * sqrt(theta);
  sigmaO = sigmas_between[1];
  sigmaOP = sigmas_between[2];
}

model {
  tau ~ gamma(shape, scale);
  theta ~ dirichlet(rep_vector(concentration, 2));
  Op ~ normal(0, sigmaO);
  OpPart ~ normal(0, sigmaOP);
  sigmaE ~ cauchy(0, 5);
  PartA ~ normal(0, 100); 
  for(k in 1:N){
    y[k] ~ normal(
      PartA[PartID[k]] + 
        Op[OperatorID[k]] + 
          OpPart[InteractionID[k]], 
      sigmaE
    );
  }
}

generated quantities {
  real<lower=0> sigmaTotal;
  sigmaTotal = sqrt(sigmaE^2 + sigmaO^2 + sigmaOP^2);
}
