library(R2jags)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rjags)
library(mcmcplots)
library(truncnorm)
library(Rlab)

## simulate the fake data
set.seed(93)
J <- 40           # number of multiple choice could take
N <- 5            # number of items in the exam
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
theta = rnorm(I,0,1)
b = rnorm(J,0,2)

for (i in 1:I){
  for (j in 1:J){
    lambda[i,j] <- 1 / (1 + exp(b[j] -  theta[i]))
    for (n in 1:N){
      c[j,n] <- rnorm(n=1,(1/2)*d[j]*z[j,n],1)
      p[i,j,n] <- lambda[i,j] * (1 -pnorm(c[j,n] - d[j] * z[j,n],0,1)) + (1 - lambda[i,j]) * (1- pnorm(c[j,n],0,1))
      y[i,j,n] <- rbern(1,p[i,j,n])        
    }
  }
}


HMM.SDT.dat <- list("N"=N,"I"=I,"J"=J,"y"=y,'z'=z)
HMM.SDT.fit <- stan("HMM-SDT.stan", data=HMM.SDT.dat)
saveRDS(HMM.SDT.fit,'HMM_SDT_fit.rds')