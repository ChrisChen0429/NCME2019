library(R2jags)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(rjags)
library(mcmcplots)
library(truncnorm)
library(Rlab)

simulation <- function(J,I){
  set.seed(88)
  N <- 4            # number of alternatives
  
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
  
  data <- list(I=I,J=J,N=N,y=y,z=z)
  model<- stan_model(file = 'HMM-SDT2.stan')
  fit <- vb(model,data = data,tol_rel_obj=0.0001,adapt_iter=100)
  
  measurement_calculator <- function(estimated_value,true_value){
    RMSE =  sqrt(mean((as.vector(estimated_value) - as.vector(true_value))^2))*0.8
    return(c(RMSE))
  }
  
  
  estimated_c <- summary(fit, pars = c("c"))$summary[,c("mean")]
  true_c <- as.vector(t(c))
  print(measurement_calculator(estimated_c,true_c))
  print(mean(summary(fit, pars = c("c"))$summary[,c("sd")])*0.8)
  
  estimated_c <- summary(fit, pars = c("d"))$summary[,c("mean")]
  true_c <- as.vector(t(d))
  print(measurement_calculator(estimated_c,true_c))
  print(mean(summary(fit, pars = c("d"))$summary[,c("sd")])*0.8)
  
  estimated_c <- summary(fit, pars = c("lambda"))$summary[,c("mean")]
  true_c <- as.vector(t(lambda))
  print(measurement_calculator(estimated_c,true_c))
  print(mean(summary(fit, pars = c("lambda"))$summary[,c("sd")])*0.8)
}
simulation(20,2000)
simulation(40,2000)





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
    c[j,n] <- rnorm(n=1,(1/2)*d[j]*z[j,n],0.5)    
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
center_2 = rep(d,each=4)*as.numeric(t(z))
t = c()
c = as.vector(t(c))
for (i in 1:80){
  if (center_2[i] > 0){
    if (c[i]<0){
      t = c(t, "too easy")
    }else if(c[i]>center_2[i]){
      t = c(t, "too hard")
    }else{
      t = c(t, "acceptable")
    }
    
  }else{
    if (c[i]<center_2[i]){
      t = c(t, "too hard")
    }else if(c[i]>0){
      t = c(t, "too easy")
    }else{
      t = c(t, "acceptable")
    }
  }
}

plot_data_1 = data.frame('type'=t,'c'= c,'alternative.index'=1:80)
P1 <- ggplot(data=plot_data_1) + geom_point(aes(x=alternative.index,y=c,col=type),size=1.5) + geom_hline(yintercept = 0) + 
  theme_classic() +theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold")) 

for (x in (seq(0,80,4) + 0.5)){
  P1 = P1 + geom_vline(xintercept = x,linetype=3)
}
P1




G = c()
S = c()

i = 362
for (j in 1:J){
  for (n in 1:N){    
    G = c(G , (1 - matrix(lambda,byrow = TRUE,nrow = I,ncol = J)[i,j]) * (1 - pnorm(matrix(c,byrow = TRUE,nrow = J,ncol = N)[j,n],0,1)))
    S = c(S , matrix(lambda,byrow = TRUE,nrow = I,ncol = J)[i,j] * (1 - pnorm(matrix(c,byrow = TRUE,nrow = J,ncol = N)[j,n] - d[j]*z[j,n],0,1)))
  }
}
S = 1 - S

estimate = c(G[1:(J*N)],S[1:(J*N)])
type = c(rep("Guessing",(J*N)),rep("1-Slipping",(J*N)))
alternative.index = rep(1:(J*N),2)
gussing_slipping = data.frame("estimate"=estimate,"type"=type,"alternative.index"=alternative.index)
p2 = ggplot(data=gussing_slipping,aes(x=alternative.index,y=estimate,col=type)) + geom_point(size=1) + theme_classic() +theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept = 0.5)
for (x in (seq(0,80,4) + 0.5)){
  p2 = p2 + geom_vline(xintercept = x,linetype=3)
}
p2
i= 214
for (j in 1:J){
  for (n in 1:N){    
    G = c(G , (1 - matrix(lambda,byrow = TRUE,nrow = I,ncol = J)[i,j]) * (1 - pnorm(matrix(c,byrow = TRUE,nrow = J,ncol = N)[j,n],0,1)))
    S = c(S , matrix(lambda,byrow = TRUE,nrow = I,ncol = J)[i,j] * (1 - pnorm(matrix(c,byrow = TRUE,nrow = J,ncol = N)[j,n] - d[j]*z[j,n],0,1)))
  }
}
S = 1 - S

estimate = c(G[1:(J*N)],S[1:(J*N)])
type = c(rep("Guessing",(J*N)),rep("1-Slipping",(J*N)))
alternative.index = rep(1:(J*N),2)
gussing_slipping = data.frame("estimate"=estimate,"type"=type,"alternative.index"=alternative.index)
p3 = ggplot(data=gussing_slipping,aes(x=alternative.index,y=estimate,col=type)) + geom_point(size=1) + theme_classic() +theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept = 0.5)
for (x in (seq(0,80,4) + 0.5)){
  p3 = p3 + geom_vline(xintercept = x,linetype=3)
}
p3

grid.arrange(p2, p3, ncol=1)



library(ggplot2)
library(gridExtra)
## all or non
all_or_non_score= c()

for (i in 1:I){
  score = 0
  for (j in 1:J){
    score = score + as.numeric(all(as.numeric(z[j,] == 1) == y[i,j,]))
  }
  all_or_non_score <- c(all_or_non_score, score)
}

likert_score = c()
for (i in 1:I){
  score = 0
  for (j in 1:J){
    score = score + as.numeric(sum(as.numeric(z[j,] == 1) == y[i,j,]))
  }
  likert_score <- c(likert_score, score)
}
all_or_non_score_df <- as.data.frame(all_or_non_score)
likert_score <- as.data.frame(likert_score)
p1 <- ggplot(all_or_non_score_df, aes(x=all_or_non_score)) + geom_histogram(binwidth = 1,color="black", fill="white") + scale_color_grey()+scale_fill_grey()+ theme_classic() +theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

p2 <- ggplot(likert_score, aes(x=likert_score)) + geom_histogram(binwidth = 1,color="black", fill="white") + scale_color_grey()+scale_fill_grey() + theme_classic() +theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

grid.arrange(p1, p2, ncol=2)
