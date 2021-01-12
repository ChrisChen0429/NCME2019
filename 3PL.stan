data {
  int<lower=0> I;
  int<lower=0> J; 
  int<lower=0,upper=1> obs_3pl[I,J];
}

parameters {
  vector[I] theta;
  vector<lower=0>[J] a;
  vector[J] b; 
  vector<lower=0,upper=1>[J] c;
}
transformed parameters{
  real<lower=0,upper=1> p[I,J]; 
  for (i in 1:I){
    for (j in 1:J){
      p[i,j] =inv_logit(a[j]*(theta[i] - b[j]));
    }
  }
}

model {
  theta ~ normal(0,1);
  b ~ normal(0,2);
  a ~ lognormal(0,0.5);
  c ~ uniform(0,1);
  for (i in 1:I){
    for (j in 1:J){
      obs_3pl[i,j] ~ bernoulli(c[j] + (1 - c[j])*p[i,j]);
    }
  }
}
