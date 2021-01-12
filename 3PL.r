library(R2jags)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rjags)
library(mcmcplots)
library(truncnorm)
library(Rlab)

## simulate the fake data
set.seed(88)
J <- 20           # number of item
N <- 4            # number of alternatives
I <- 500          # number of examees

### SDT part

z = array(NA, c(J,N))    # true answer 
for (j in 1:J){
  for (n in 1:N){
    this <- rbern(1,0.5)
    if (this > 0.5){
      z[j,n] <- 1
    }else{
      z[j,n] <- -1
    }
  }
}

lambda = array(NA, c(I,J))
y = array(NA, c(I,J,N))
c = array(NA, c(J,N))
delta = array(NA, c(I,J,N))
p = array(NA, c(I,J,N))

d = rlnorm(J,0,0.5)
theta <- rnorm(I,0,0.5)
for (j in 1:J){
  for (n in 1:N){
    c[j,n] <- rnorm(n=1,(1/2)*d[j]*z[j,n],1)    
  }
}

b <- apply(c*z,1,mean)
for (i in 1:I){
  for (j in 1:J){
    lambda[i,j] <-1 / (1 + exp( b[j] - theta[i]))
    for (n in 1:N){
      p[i,j,n] <- lambda[i,j] * (1 -pnorm(c[j,n] - d[j] * z[j,n],0,1)) + (1 - lambda[i,j]) * (1- pnorm(c[j,n],0,1))
      y[i,j,n] <- rbern(1,p[i,j,n])        
    }
  }
}


obs_3pl= matrix(NA, nrow = I, ncol = J)

for (i in 1:I){
  for (j in 1:J){
    obs_3pl[i,j] =  as.numeric(all(as.numeric(z[j,] == 1) == y[i,j,]))
  }
}
PL.dat <- list("I"=I,"J"=J,"obs_3pl"=obs_3pl)
PL_fit <- stan("3PL.stan")


estimated_beta <- summary(PL_fit, pars = c("b"))$summary[,c("mean")]
estimated_theta <- summary(PL_fit, pars = c("theta"))$summary[,c("mean")]
estimated_a <- summary(PL_fit, pars = c("a"))$summary[,c("mean")]

m <- cbind(rank(estimated_theta), rank(theta))
cor(m) 

m <- cbind(rank(estimated_beta), rank(b))
cor(m) 

m <- cbind(rank(estimated_a), rank(d))
cor(m) 
