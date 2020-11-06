// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

data {
  int<lower=1> N;
  real y[N];
  int<lower=1> I;
  int<lower=1> J;
//  int<lower=1> Kij[I,J];
  int<lower=1> PartIndicator[N];
  int<lower=1> OperatorIndicator[N];
  int<lower=1> InteractionIndicator[N];
}

parameters {
  real PartA[I];
  real<lower=0> sigmaE;
  real<lower=0> sigmaO;
  real<lower=0> sigmaOP;
  real Op[J];
  real OpPart[I*J];
  
}

model {
  PartA ~ normal(0, 100);
  sigmaE ~ cauchy(0, 5);
  sigmaO ~ cauchy(0, 5);
  sigmaOP ~ cauchy(0, 5);
  Op ~ normal(0, sigmaO);
  OpPart ~ normal(0, sigmaOP);
  for(k in 1:N){
    y[k] ~ normal(
      PartA[PartIndicator[k]] + 
        Op[OperatorIndicator[k]] + 
          OpPart[InteractionIndicator[k]], 
      sigmaE
    );
  }
}

generated quantities {
  real<lower=0> sigmaTotal;
  sigmaTotal = sqrt(sigmaE^2 + sigmaO^2 + sigmaOP^2);
}
