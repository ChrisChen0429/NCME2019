data {
  int<lower=0> I; //number of exaimee
  int<lower=0> J; //number of item
  int<lower=0> N; //number of alternatives
  int y[I,J,N];
  int z[J,N];
}

parameters {
  real c[J,N];
  real<lower=0> d[J];
  real<lower=0,upper=1> lambda[I,J];
}
model {
  d ~ lognormal(0,1);
  for (j in 1:J){
    for (n in 1:N){
      c[j,n] ~ normal((1/2)*d[j]*z[j,n],1);
     for (i in 1:I){
            y[i,j,n] ~ bernoulli( lambda[i,j] * (1 - normal_cdf(c[j,n] - d[j] * z[j,n],0,1)) + (1 - lambda[i,j]) * (1 - normal_cdf(c[j,n],0,1)));
          }
        }
    }
}
generated quantities{
  real beta[J];
  real theta[I];
  real middle[J,N];
  
  for (j in 1:J){
    for (n in 1:N){
     middle[j,n] =  c[j,n]*z[j,n];
    }
    beta[j] = mean(middle[j,]);  
  }
  for (i in 1:I){
    theta[i] = mean(lambda[i,]);
  }
}
