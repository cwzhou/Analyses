gdata <- function(N=10,lambda_0D=0.1,lambda_0R=4,beta_R=log(2),beta_D=log(3),
                  ztype=0,zparam=0.5,ctype=1,cparam=2,
                  gaptype=0,gapparam1=0.2,gapparam2=0.25,#rho;
                  num_A=2,tau=10){
  
  # ztype indicates the distribution for covariate z: 0=normal(0,1),1=binary(zparam),2=uniform(0,1)
  # ctype indicates the distribution for censoring (0=exponential,1=uniform,
  #                                                 9=no censoring)
  
  # generating covariates
  if (ztype == 0){z <- rnorm(N)} 
  if (ztype == 1){z <- rbinom(N, 1, zparam)}
  if (ztype == 2){z <- runif(N)}
  
  # generating treatment
  A = rbinom(N, num_A-1, 0.5)
  
  # generating censoring time
  if (ctype == 0){cc <- rexp(N,cparam)}
  if (ctype == 1){cc <- runif(N,min=0,max=tau)}
  
  # generating gap time 1 using Gumbel bivariate exponential model
  #if alpha=0 then independent and don't need to do this.
  alpha1 = 4*gapparam1
  alpha2 = 4*gapparam2
  lambda_D = lambda_0D*exp(beta_D*z)
  lambda_R = lambda_0R*exp(beta_R*z)
  gaptime1 = rexp(N,lambda_R)
  if (alpha1 == 0){
    tt = rexp(N,lambda_D)
  } else{
    tt <- gumbelbiexp(N=N,z=z,alpha=alpha1,lambda_D=lambda_D,lambda_R=lambda_R,beta=beta_R,y=gaptime1)$tt
  }
}

Surv(tt,rep(1, 100))
s1 <- survfit(Surv(tt, rep(1,100)) ~ 1)
str(s1)
s1
# function for generating failure time from Gumbel bivariate exponential distribution (1960)
gumbelbiexp <- function(u, alpha, y) {
  cc = 1-2*exp(-y); 
  print("1-2exp(-y)")
  # print(cc)
  rootsqd = (u-1)/(cc*alpha) + ((1+cc*alpha)/(2*cc*alpha))^2; #print(rootsqd)
  root = sqrt(rootsqd);
  # print("root")
  # print(root)
  expnegx_plus = (1+cc*alpha)/(2*cc*alpha) + root;
  expnegx_minus = (1+cc*alpha)/(2*cc*alpha) - root;
  print((expnegx_plus > 1))
  print(expnegx_minus)
  print((expnegx_plus > 1))
  expnegx = expnegx_minus*(expnegx_plus>1)+expnegx_plus*(expnegx_plus<=1)
  print(sprintf("exp(-x): %s", expnegx))
  print("tt")
  tt = -log(expnegx)
  return(tt)
}

N=1000
set.seed(2023)
u_N = runif(N)
y_N = rexp(N)
alpha_seq = 4*c(0.0625,0.125, 0.25)
alpha = alpha_seq[2]; print(alpha)
print(sprintf("#### alpha: %s #### ", alpha))
tt = gumbelbiexp(u_N, alpha, y_N)

set.seed(1)
alpha_seq = seq(-1, 1, 0.25); alpha <- alpha_seq[-which(alpha_seq == 0)]
z = rnorm(1) #covariates
lambda_0D = abs(rnorm(5));
sorted_lam0D <- lambda_0D[order(lambda_0D)]
lambda_0R = abs(rnorm(5)); 
sorted_lam0R <- lambda_0R[order(lambda_0R)]
beta_D = c(-1,0.2,3)
beta_R = c(-0.2,0.1,3)
df1 = data.frame()
length(z)*length(alpha)*length(lambda_0D)*length(lambda_0R)*length(beta_D)*length(beta_R)
for (alpha_val in alpha){
  for (z_val in z){
    for (lambda_0D_val in sorted_lam0D){
      for (lambda_0R_val in sorted_lam0R){
        for (beta_D_val in beta_D){
          for (beta_R_val in beta_R){
            lambda_D_val = lambda_0D_val*exp(beta_D_val*z_val)
            lambda_R_val = lambda_0R_val*exp(beta_R_val*z_val)
            y_val = rexp(N, lambda_R_val)
            x = gumbelbiexp_x_given_y(N=1,z=z_val,alpha=alpha_val,lambda_D=lambda_D_val,lambda_R=lambda_R_val,y=y_val)
            row = data.frame(N = N, z = z_val, alpha = alpha_val, 
                             lambda_0D = lambda_0D_val, beta_D = beta_D_val, lambda_D = lambda_D_val, 
                             lambda_0R = lambda_0R_val, beta_R = beta_R_val, lambda_R = lambda_R_val, 
                             gaptime = y_val, generated_failure_time = x)
            # print(row)
            df1 <- rbind(df1, row) 
            # print(df1)
          }
        }
      }
    }
  }
}
View(df1)
dim(df1)

df1$g_ft = as.numeric(df1$generated_failure_time)
df1 %>% filter(g_ft>0) %>% nrow()

library(dplyr)
# Calculate proportions 
df1 %>%
  summarize(
    pos_prop = sum(g_ft > 0, na.rm = TRUE)/nrow(df1),
    neg_prop = sum(g_ft < 0, na.rm = TRUE)/nrow(df1),
    miss_prop = sum(is.na(g_ft))/nrow(df1)
  )


  
  tt = -(1/lambda_D)*
  #suma <- matrix(rep(0,ngrp),ngrp,1) 
  tt[,1] <- -log(1-u[,1])*exp(-beta*z[,1])
  suma <- suma+exp(tt[,i-1]*exp(beta*z[,i-1])/theta)
  bl <- log((i-1)-suma+(suma-(i-2))*(1-u[,i])^(-(theta+i-1)^(-1)))
  tt[,i] <- theta*bl*exp(-beta*z[,i])
  list(tt=tt) 
}



gaptime = function(gaptype=0,N,z,lambda_0R,rho,alpha1,beta_R){
  
}

gdata <- function(ngrp=100,nmember=2,theta=0.8,beta=log(2), 
                  ztype=0,zparam=0.4,ctype=0,cparam=2,stoptime=5) {
 
  # generating failure time from clayton model
  
  tt <- clayton(ngrp,nmember,z,theta,beta)$tt
  
  # generating the observed time and indicator variable
  
  if (ctype == 9) {
    x <- tt
    d <- matrix(rep(1,ngrp*nmember),ngrp,nmember) 
  } else if (ctype != 9)  {
    cc <- (cc <= stoptime) * cc + (cc > stoptime) * stoptime  
    x <- (tt <= cc) * tt + (tt > cc) * cc
    d <- (tt <= cc)
  }
  list(x=x,d=d,z=z) 
}

gdatap <- function(ngrp=100,nmember=2,theta=0.8,beta=log(2), 
                   ztype=0,zparam=0.4,ctype=0,cparam=2,stoptime=5,ttype=1,v=0,lambda=1) {
  # ztype indicates the distribution for z 
  #    (0-2 are for cluster varying covariates: 
  #               0=normal(0,1),1=binary,2=uniform(0,1))
  #     3-5 are for cluster constant covariates:
  #               3=normal(0,1),4=binary,5=uniform(0,1))
  #     6 is for assigned two treatment groups (half each).
  # ctype indicates the distribution for censoring (0=exponential,1=uniform,
  #                                                 9=no censoring)
  # ttype indicates the hazard function for failure (1=exponential,
  #                                                  2=piecewise exponential) 
  
  # generating covariates
  
  if (ztype == 0) z <- matrix(rnorm(ngrp*nmember),ngrp,nmember)
  if (ztype == 1) z <- matrix(rbinom(ngrp*nmember, 1, zparam),ngrp,nmember)
  if (ztype == 2) z <- matrix(runif(ngrp*nmember),ngrp,nmember)
  if (ztype == 3) z <- matrix(rnorm(ngrp),ngrp,nmember)
  if (ztype == 4) z <- matrix(rbinom(ngrp, 1, zparam),ngrp,nmember)
  if (ztype == 5) z <- matrix(runif(ngrp),ngrp,nmember)
  if (ztype == 6) z <- matrix(c(rep(0,ngrp/2),rep(1,ngrp/2)),ngrp,nmember)
  
  # generating censoring time
  
  if (ctype == 0) cc <- matrix(rexp(ngrp*nmember,cparam),ngrp,nmember)
  if (ctype == 1) cc <- matrix(runif(ngrp*nmember,max=cparam),ngrp,nmember)
  
  # generating failure time from clayton model
  
  if (ttype == 1) tt <- clayton(ngrp,nmember,z,theta,beta)$tt
  if (ttype == 2) tt <- claytonp(ngrp,nmember,z,theta,beta,v,lambda)$tt
  
  # generating the observed time and indicator variable
  
  if (ctype == 9) {
    x <- tt
    d <- matrix(rep(1,ngrp*nmember),ngrp,nmember) 
  } else if (ctype != 9)  {
    cc <- (cc <= stoptime) * cc + (cc > stoptime) * stoptime  
    x <- (tt <= cc) * tt + (tt > cc) * cc
    d <- (tt <= cc)
  }
  list(x=x,d=d,z=z) 
}

# function for generating multivariate failure time from Clayton model

clayton <- function(ngrp=100,nmember=2,z,theta=0.8,beta) { 
  u <- matrix(runif(ngrp*nmember),ngrp,nmember)
  tt <- matrix(rep(0,ngrp*nmember),ngrp,nmember)
  suma <- matrix(rep(0,ngrp),ngrp,1) 
  tt[,1] <- -log(1-u[,1])*exp(-beta*z[,1])
  if (nmember > 1) for (i in 2:nmember) {
    suma <- suma+exp(tt[,i-1]*exp(beta*z[,i-1])/theta)
    bl <- log((i-1)-suma+(suma-(i-2))*(1-u[,i])^(-(theta+i-1)^(-1)))
    tt[,i] <- theta*bl*exp(-beta*z[,i])
  }
  list(tt=tt) 
}

# function for generating multivariate failure time from Clayton model
# with piecewise exponential marginals

claytonp <- function(ngrp=100,nmember=2,z,theta=0.8,beta=0,v=0,lambda=1) { 
  nv <- length(v)
  u <- matrix(runif(ngrp*nmember),ngrp,nmember)
  tt <- matrix(rep(0,ngrp*nmember),ngrp,nmember)
  
  sumv <- 0
  bl <- matrix(rep(0,ngrp*nv),ngrp,nv)
  bl1 <- -log(1-u[,1])*exp(-beta[1]*z[,1])/lambda[1]
  if (nv==1) tt[,1] <- bl1 else if (nv > 1) {
    bl[,1] <- bl1*(bl1 >= v[1])*(bl1 < v[2])
    for (m in 2:nv) {
      sumv <- sumv + (v[m]-v[m-1])*exp(beta[m-1]*z[,1])*lambda[m-1]
      bl1 <- -(log(1-u[,1])+sumv)*exp(-beta[m]*z[,1])/lambda[m]+v[m]
      if (m==nv) bl[,m] <- bl1*(bl1 >= v[m]) else if (m < nv) 
        bl[,m] <- bl1*(bl1 >= v[m])*(bl1 < v[m+1])
    }
    tt[,1] <- apply(bl,1,sum)
  }
  suma <- (1-u[,1])^(-1/theta)
  if (nmember > 1) for (i in 2:nmember) {
    bl2 <- (i-1)-suma+(suma-(i-2))*(1-u[,i])^(-(theta+i-1)^(-1))
    
    sumv <- 0 
    bl <- matrix(rep(0,ngrp*nv),ngrp,nv)
    bl1 <- theta*log(bl2)*exp(-beta[1]*z[,i])/lambda[1]
    if (nv==1) tt[,i] <- bl1 else if (nv > 1) { 
      bl[,1] <- bl1*(bl1 >= v[1])*(bl1 < v[2])
      for (m in 2:nv) {
        sumv <- sumv + (v[m]-v[m-1])*exp(beta[m-1]*z[,i])*lambda[m-1]
        bl1 <- (theta*log(bl2)-sumv)*exp(-beta[m]*z[,i])/lambda[m]+v[m]
        if (m==nv) bl[,m] <- bl1*(bl1 >= v[m]) else if (m < nv) 
          bl[,m] <- bl1*(bl1 >= v[m])*(bl1 < v[m+1])
      }
      tt[,i] <- apply(bl,1,sum)
    }
    suma <- suma+bl2
  }
  list(tt=tt) 
}

logrank <- function(x,d,z,w=1,w1=1,w2=1,w0=1)  {
  failure <- sort(unique(x[d==1]))
  yi <- ifelse(outer(x,failure,">="),1,0)
  di <- ifelse(outer(x,failure,"=="),1,0)
  ybar <- apply(yi,2,sum)
  dbar <- apply(di,2,sum)
  yi1 <- apply(yi[z==1,],2,sum)
  di1 <- apply(di[z==1,],2,sum)
  q0 <- (sum(w0*(di1-dbar*yi1/ybar)))^2 # log-rank type
  qs0 <- max((cumsum(w0*(di1-dbar*yi1/ybar)))^2) # sup log-rank type
  q1 <- (sum(w1*(di1-dbar*yi1/ybar)))^2 # log-rank type
  qs1 <- max((cumsum(w1*(di1-dbar*yi1/ybar)))^2) # sup log-rank type
  q2 <- (sum(w2*(di1-dbar*yi1/ybar)))^2 # log-rank type
  qs2 <- max((cumsum(w2*(di1-dbar*yi1/ybar)))^2) # sup log-rank type
  q <- (sum(w*(di1-dbar*yi1/ybar)))^2 # log-rank type
  qs <- max((cumsum(w*(di1-dbar*yi1/ybar)))^2) # sup log-rank type
  list(q=q,qs=qs,q1=q1,qs1=qs1,q2=q2,qs2=qs2,q0=q0,qs0=qs0)
}    

logrankB <- function(x,d,z,w=1,w1=1,w2=1,w0=1,B=100)  {
  failure <- sort(unique(x[d==1]))
  yi <- ifelse(outer(x,failure,">="),1,0)
  di <- ifelse(outer(x,failure,"=="),1,0)
  ybar <- apply(yi,2,sum)
  dbar <- apply(di,2,sum)
  yi1 <- NULL
  di1 <- NULL
  for (i in 1:B) {
    yi1 <- cbind(yi1,apply(yi[z[,i]==1,],2,sum))
    di1 <- cbind(di1,apply(di[z[,i]==1,],2,sum)) }
  q0 <- (apply(w0*(di1-dbar*yi1/ybar),2,sum))^2 # log-rank type
  qs0 <- apply((apply(w0*(di1-dbar*yi1/ybar),2,cumsum))^2,2,max)
  # sup log-rank type
  q1 <- (apply(w1*(di1-dbar*yi1/ybar),2,sum))^2 # log-rank type
  qs1 <- apply((apply(w1*(di1-dbar*yi1/ybar),2,cumsum))^2,2,max)
  # sup log-rank type
  q2 <- (apply(w2*(di1-dbar*yi1/ybar),2,sum))^2 # log-rank type
  qs2 <- apply((apply(w2*(di1-dbar*yi1/ybar),2,cumsum))^2,2,max)
  # sup log-rank type
  q <- (apply(w*(di1-dbar*yi1/ybar),2,sum))^2 # log-rank type
  qs <- apply((apply(w*(di1-dbar*yi1/ybar),2,cumsum))^2,2,max)
  # sup log-rank type
  list(q=q,qs=qs,q1=q1,qs1=qs1,q2=q2,qs2=qs2,q0=q0,qs0=qs0)
}    

