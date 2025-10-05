# Title: Functions for Comparator Methods Script
# Date: 01.21.2024
##############################################################################################
##############################################################################################
##############################################################################################
### PMCR: function for Precision Medicine with Competing Risks (2 Med options, 2 risks)

# pmcr can't handle factors with only two levels aka indicators
convert_factors_to_numeric <- function(x) {
  if (is.factor(x) && length(levels(x)) == 2) {
    as.numeric(x) - 1
  } else {
    x
  }
}


# k=21;beta = beta0[k,]; mydata=mydata; pix.hat=mydata[["pix.hat"]]; Xp=Xp; tt0=t0
# function calculating the CIFs
calc.F<-function(beta,mydata,pix.hat,Xp){
  # print('calculating CIFs')
  # View(Xp)
  # View(beta)
  xb<-Xp %*% beta
  n<-nrow(mydata)
  # smoothed version of weight
  h<-sd(xb)*(n/4)^(-1/3)
  w.numer<-pnorm(xb/h)*mydata$A+(1-pnorm(xb/h))*(1-mydata$A)
  w.denom<-pix.hat*mydata$A+(1-pix.hat)*(1-mydata$A)
  wi<-w.numer/w.denom
  # print("cc1")

  ord.t0<-sort(mydata$time)
  ord.d0<-as.numeric(mydata$status>0)[order(mydata$time)]
  ord.ebs<-mydata$status[order(mydata$time)]
  ord.wi<-wi[order(mydata$time)]
  # print("cc2")

  tt<-ord.t0[ord.d0==1]
  # tt_tmp <<-tt
  SIt<-numeric()
  SIt[1]<-1
  for(this_k in 1:length(tt)){
    SIt[this_k+1]<-SIt[this_k]*(1-sum(ord.wi[ord.t0==tt[this_k]])/sum(ord.wi[ord.t0>=tt[this_k]]))
    # cat("\n",this_k,":\n", SIt, "\n")
  }
  # print("cc3")
  # sit <<-SIt

  tt1<-ord.t0[ord.ebs==1]
  tt2<-ord.t0[ord.ebs==2]
  lamb1<-lamb2<-numeric()
  for(k in 1:length(tt1)) lamb1[k]<-sum(ord.wi[ord.t0==tt1[k]])/sum(ord.wi[ord.t0>=tt1[k]])
  for(k in 1:length(tt2)) lamb2[k]<-sum(ord.wi[ord.t0==tt2[k]])/sum(ord.wi[ord.t0>=tt2[k]])
  # print("cc4")
  # tmp_lamb1 <<- lamb1
  # tmp_lamb2 <<- lamb2
  
  # print(tt1)
  if(length(tt1)>0){ F1.tt0<-stepfun(tt1,cumsum(c(0,stepfun(tt,SIt)(tt1)*lamb1)))}else{F1.tt0<-0}
  # print("cc4.5")
  if(length(tt2)>0){ F2.tt0<-stepfun(tt2,cumsum(c(0,stepfun(tt,SIt)(tt2)*lamb2)))}else{F2.tt0<-0}
  # print("cc5")
return(list(F1t0=F1.tt0,F2t0=F2.tt0,St0=stepfun(tt,SIt)))
}

# beta = beta0[k,]; mydata=mydata; pix.hat=mydata[["pix.hat"]]; Xp=Xp; tt0=t0
# optim function
opt.beta1<-function(beta,mydata,pix.hat,Xp,tt0=5){
  # message("opt.beta1")
   cif.fit<<-calc.F(beta,mydata,pix.hat,Xp)
   # print(cif.fit)
   if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
return(F1.tt0)
}
opt.beta2<-function(beta,mydata,pix.hat,Xp,tt0=5){
  # message('opt.beta2')
   cif.fit<-calc.F(beta,mydata,pix.hat,Xp)
   if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
return(F2.tt0)
}
opt.beta<-function(beta,mydata,pix.hat,Xp,tt0=5,alp=0.1,M=1000){
   # message('opt.beta')
   cif.fit<-calc.F(beta,mydata,pix.hat,Xp)
   if(class(cif.fit$F1t0)[1]!="numeric"){F1.tt0<-cif.fit$F1t0(tt0)}else{F1.tt0<-0}
   if(class(cif.fit$F2t0)[1]!="numeric"){F2.tt0<-cif.fit$F2t0(tt0)}else{F2.tt0<-0}
return(F1.tt0+M*(F2.tt0-alp)*as.numeric(F2.tt0>alp))
}

####################################################################################################
# Arguments in function PMCR                                                                       #
# Time: observed event or censoring time                                                           #
# Event: indicator variable of events or censoring(0)                                              #
# formula: treatment group indicator ~ predictors                                                  #
# data: data frame with all variables needed                                                       #
# rgenoud: if TRUE, the rgenoud function will be applied to search beta                            #
# Restrict: if TRUE, restricted regime is estimated. otherwise, unrestricted regimes are returned. #
# propscore: method for calculating the propensity scores. "logistic" or "tree".                   #
#            Note, the option "true" can be specified when the data is generated using gen.data,   #
#            and the true propensity score will be used.                                           #
# t0: time for CIF to be evaluated at.                                                             #
# alp: alpha value as the upper limit of risk2 CIF                                                 #
# M: a large number used as penalty for the restriction being violated                             #
####################################################################################################

# Time="obs_time"; Event="status";rgenoud=FALSE; Restrict=FALSE; propscore="logistic"; 
# formula=modelPr_PMCR; data=train_pmcr;
# t0=t0_pmcr
PMCR<-function(Time,
               Event,
               formula=formula(data),
               data=parent.frame(),
               rgenoud=TRUE,
               Restrict=TRUE,
               propscore="logistic",
               t0=5,alp=0.1,M=1000){
  # message("VERY FIRST ARGUMENT t0:",t0)
  mydata<-data
  call <- match.call()
  avars<-all.vars(formula)
  mydata$time<-data[,Time]
  mydata$status<-data[,Event]
  mydata$A<-data[,avars[1]]
  if(propscore=="true")     mydata$pix.hat<-mydata$piz
  if(propscore=="logistic") mydata$pix.hat<-predict(glm(formula,mydata,family=binomial("logit")),type="response")
  if(propscore=="tree")     mydata$pix.hat<-predict(rpart(formula,mydata,method="class"))[,"1"]
  temp.x<-terms(formula, data = mydata)
  Xp <- model.matrix(temp.x,mydata)
  npar<-ncol(Xp)
  beta0<-matrix(0,2*npar,npar)
  for(i in 1:npar){
    beta0[2*i-1,i]<-1
    beta0[2*i,i]<--1
  }
  outest<-list(data=mydata, cov = as.data.frame(Xp))
  if(!Restrict){
    # message("NOT restrict t0:",t0)
    f1val<-1
    # View(mydata)
    # beta0 <<- beta0
    # opt.beta1 <<- opt.beta1
    # mydata <<- mydata
    # pix.hat <<-mydata$pix.hat
    # Xp<<-Xp
    # tt0<<-t0
    # View(Xp)
    # print(t0)
    # print(beta0)
    # print(nrow(beta0))
    for(k in 1:nrow(beta0)){
      # message("k:",k)
      # message("166 t0:",t0)
      fit1<-try(optim(par=beta0[k,],
                      fn=opt.beta1,
                      mydata=mydata,
                      pix.hat=mydata[["pix.hat"]],
                      Xp=Xp,
                      tt0=t0),
                silent=FALSE)
      # print(opt.beta1)
      # print(fit1)
      # print(t0)
      # print(cif.fit$F1t0)
      # print(cif.fit$F1t0(t0))
      if(fit1$value<f1val){
        # message("true")
        outest$beta1<-fit1$par/sqrt(sum(fit1$par^2))
        f1val<-fit1$value
        # message("bottom k=",k,"\n")
        # print(c(outest$beta1,f1val))
      }
      # message("end of k")
      # message("179 t0:",t0)
    }
    # print(5)
    f2val<-1
    for(k in 1:nrow(beta0)){
      # message("k:",k)
      # t0 = 1
      # print(t0)
      # message("188 t0:",t0)
      # print(cif.fit$F2t0)
      fit2<-try(optim(par=beta0[k,],fn=opt.beta2,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0),silent=TRUE)
      # print(fit2)
      if(fit2$value<f2val){
        outest$beta2<-fit2$par/sqrt(sum(fit2$par^2))
        f2val<-fit2$value
        # cat(paste("k=",k,"\n"))
        #print(c(outest$beta2,f2val))
      }
    }
    # message("198 t0:",t0)
    # print(6)
  if(rgenoud){
    fit1<-try(genoud(opt.beta1,nvars=npar,max=FALSE,starting.values=fit1$par,max.generations=30,print.level=0,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0),silent=TRUE)
    if(!is.character(fit1)){
      if(fit1$value<f1val){
        outest$beta1<-fit1$par/sqrt(sum(fit1$par^2))
        f1val<-fit1$value
        #cat("rgenoud replace\n")
        #print(c(outest$beta1,f1val))
      }
    }
    # print(7)
    fit2<-try(genoud(opt.beta2,nvars=npar,max=FALSE,starting.values=fit2$par,max.generations=30,print.level=0,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0),silent=TRUE)
    if(!is.character(fit2)){
      if(fit2$value<f2val){
        outest$beta2<-fit2$par/sqrt(sum(fit2$par^2))
        f2val<-fit2$value
        #cat("rgenoud replace\n")
        #print(c(outest$beta2,f2val))
      }
    }
  }
  # message('FF 179')
  FF<-calc.F(beta=outest$beta1,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp)
  outest$Fbeta1<-c(FF$F1t0(t0),FF$F2t0(t0))
  outest$f1val<-f1val
  # message('FF 183')
  FF<-calc.F(beta=outest$beta2,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp)
  outest$Fbeta2<-c(FF$F1t0(t0),FF$F2t0(t0))
  outest$f2val<-f2val
  }

  # print(8)
  if(Restrict){
    # message("Restrict t0:",t0)
    f3val<-1
    # print(opt.beta)
    for(k in 1:nrow(beta0)){
      fit<-try(optim(par=beta0[k,],fn=opt.beta,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0,alp=alp,M=M),silent=TRUE)
      # print(fit)
      if(fit$value<f3val){
        outest$beta3<-fit$par/sqrt(sum(fit$par^2))
        f3val<-fit$value
        #cat(paste("k=",k,"\n"))
        #print(c(outest$beta3,f3val))
      }
    }
    # print(9)
    if(rgenoud){
      fit<-try(genoud(opt.beta,nvars=npar,max=FALSE,starting.values=fit$par,max.generations=30,print.level=0,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp,tt0=t0,alp=alp,M=M),silent=TRUE)
      if(!is.character(fit)){
        if(fit$value<f3val){
          outest$beta3<-fit$par/sqrt(sum(fit$par^2))
          f3val<-fit$value
          #cat("rgenoud replace\n")
          #print(c(outest$beta3,f3val))
        }
      }
    }
    # message('FF 211')
    FF<-calc.F(beta=outest$beta3,mydata=mydata,pix.hat=mydata$pix.hat,Xp=Xp)
    outest$f3val<-f3val
    outest$Fbeta3<-c(FF$F1t0(t0),FF$F2t0(t0))
  }
return(outest)
}
##############################################################################################
##############################################################################################
##############################################################################################



#### AIPWE:::::
#data.df = train
aipwe_data_format = function(data.df){
  tempdata = NULL
  # View(data.df)
  A <- data.df %>%
    dplyr::select("Trt") %>%
    mutate(Trt = ifelse(Trt == 1, 1, 0)) %>% # Trt has to be 0/1 # LIMITATION TO AIPWE: can only consider binary treatments
    dplyr::pull(Trt)
  colnames(A) = NULL
  xs <- data.df %>%
    dplyr::select(-c(obs_time, Trt, status, D.0, D.1, D.2)) %>%
    # model.matrix(~ . -1 , data = .) # %>%
    as.matrix()
  ixn = matrix(NA, nrow = nrow(xs), ncol = ncol(xs))
  # View(xs)
  # View(ncol(xs))
  for (c in 1:ncol(xs)){
    print(c)
    cc = as.numeric(xs[,c])
    cc1 = cbind(cc,A) %>%
      as.data.frame() %>%
      mutate(ix = ifelse(A == 1, cc, 0)) %>%
      dplyr::select(ix)
    ixn[,c] = cc1$ix
  }
  #covariates, treatment assignment (if of interest: and interactions)
  Z1 = cbind(xs, A, ixn)
  # View(Z1)
  Z1_Vnames <- c(paste0("V", 1:ncol(xs)),
                 "A",
                 paste0("V", (ncol(xs)+1+1):(ncol(xs)+1+ncol(ixn))))
  colnames(Z1) = Z1_Vnames
  # View(Z1)
  X <- data.df %>% dplyr::select("obs_time") %>% as.vector()
  cause <- data.df %>% dplyr::select("status")
  tempdata <- data.frame(time=X,cause=cause, Z1)
  colnames(tempdata) = c("time", "cause", colnames(Z1))
  tempdata<-tempdata[order(tempdata$time),]
  return(list(data = tempdata, Z1=Z1))
}

#' returns weight calculated by bi-level weight lasso selection method
weightAGB <- function(cj,gamma,tune,beta,j,nu=1.5){
  message("weightAGB: returns weight calculated by bi-level weight lasso selection method")
  #j is a vector with group index
  #cj is a vector with c-value for each group
  #tune equals to ((1-gamma)/tau*gamma)^(gamma-1)
  ngroup <- length(unique(j))
  theta <- sapply(1:ngroup, function(i){
    normbeta <- (sum(abs(beta[which(j==i)]))^(gamma-1))
    normbeta <- normbeta/abs(beta[which(j==i)])^nu
    return(cj[i]*tune*normbeta)
  })
  return(unlist(theta))
}

wsvmL1.group <- function(x,y,w,gamma){
  # View(x)
  #dimension of covariate
  pp <- ncol(x)
  # message("pp_wsvmL1.group:", pp)
  #sample size
  nn <- nrow(x)
  ##objective vector
  obj <- c(c(w),gamma,rep(0,1+pp))
  ##inequality matrix G
  G <- matrix(0,nrow = nn+pp+pp, ncol = nn+pp+1+pp)
  G[1:nn,] <- cbind(diag(nn),matrix(0,nn,pp),y,x*y)
  G[(nn+1):(nn+pp),] <- cbind(matrix(0,nrow=pp,ncol=nn),
                              cbind(
                                diag(pp),
                                matrix(0,nrow = pp,ncol=1),
                                -diag(pp)))
  G[(nn+pp+1):(nn+2*pp),] <- cbind(matrix(0,nrow=pp,ncol=nn),
                                   cbind(diag(pp),
                                         matrix(0,nrow = pp,ncol = 1),
                                         diag(pp)))

  h <- c(rep(1,nn),rep(0,2*pp))
  dir <- c(rep('>=',nn+2*pp))
  lwr <- list(ind = c((nn+pp+1):(nn+pp+1+pp)),
              val=rep(-Inf,1+pp))
  bds = list(lower = lwr)
  soln <- Rglpk::Rglpk_solve_LP(
    obj,G,dir,h,bounds = bds,max=F,types = rep('C',nn+1)
  )
  #soln <- lpSolve::lp('min',obj,G,dir,h)
  b0 <- soln$solution[(nn+pp+1)]
  b <- soln$solution[(nn+pp+1+1):(nn+pp+1+pp)]
  M <- matrix(c(b0,b),pp+1,1)
  rownames(M) <- c('intercept',paste0('X',c(1:pp),seperate = ' '))
  return(list(sv=M))
}


#' function generate event time for cause 1
#' event time follows Fine and Gray model
Tcause1.cll <- function(p, beta, Z, omega=1){
  n <- nrow(Z)
  U <- runif(n)
  prob1 <- (1-(1-p)^(exp(Z%*%beta)*omega))
  W <- (1-U*prob1)^(1/exp(Z%*%beta)/omega)
  T <- -log(1-(1-W)/p)
  return(T)
}
#' function generate event time for cause 2
Tcause2.cll <- function(beta, Z,omega=1)
{
  n <- nrow(Z)
  U <- runif(n)
  T <- -log(U)/exp(Z%*%beta)
  return(T)
}
CensTime <- function(lambdaC, betaC, ZC, omega=1){
  n <- nrow(ZC)
  U <- runif(n)
  C <- -log(U)/exp(ZC%*%betaC)/lambdaC/omega
  return(C)
}


# function to generate data for the simulation studies
dataGen <- function(m,
                    mu.m,
                    sigma.m,
                    gammaA,
                    p,
                    beta11){
  tempdata <- NULL
  #cont variables
  xs.cont <- mvrnorm(m,mu.m,sigma.m)
  #generate categorical variables
  xs.cat <- sapply(1:m,function(i){
    #start with all te categories
    xs.cats1 <- rep(0,2*3)
    xs.cats2 <- rep(0,2*4)
    xs.cats3 <- rep(0,2*6)
    cat.temp1 <- sample(c(1:3),2,replace = TRUE);
    cat.temp2 <- sample(c(1:4),2,replace=TRUE)
    cat.temp3 <- sample(c(1:6),2,replace=TRUE)
    xs.cats1[c(cat.temp1[1],cat.temp1[2]+3)] <- 1
    xs.cats2[c(cat.temp2[1],cat.temp2[2]+4)] <- 1
    xs.cats3[c(cat.temp3[1],cat.temp3[2]+6)] <- 1
    #delete dummies
    xs.cats1 <- xs.cats1[-c(3,6)]
    xs.cats2 <- xs.cats2[-c(4,8)]
    xs.cats3 <- xs.cats3[-c(6,12)]
    return(c(xs.cats1,xs.cats2,xs.cats3))
  })
  xs.cat <- t(xs.cat)
  xs.cat <- xs.cat[,-c(16:20)]
  #combine cont and cat
  xs <- cbind(xs.cont,xs.cat)
  ##with intercept
  xs.1 <- cbind(1,xs)

  logitA<-xs.1%*%gammaA
  pA<-exp(logitA)/(1+exp(logitA))

  A<-rbinom(m,1,prob=pA)
  #Z1 has 65 col x1-x32 A A*x1-A*x32
  #covariates, treatment assignment and interactions
  # View(A)
  # View(xs)
  # View(A*xs)
  Z1 <- cbind(xs,A,A*xs)
  #probabiliy of cause 1
  prob1 <- c(1-(1-p)^(exp(Z1%*%beta11)))
  cause <- sapply(1:m,function(i){
    rbinom(1,1,1-prob1[i])+1
  })
  T <- numeric(m)
  T[cause==1] <- #T1.tmp[cause==1]
    Tcause1.cll(p, beta11, Z=Z1, omega=1)[cause==1]
  T[cause==2] <- #T2.tmp[cause==2]
    Tcause2.cll(beta11, Z=Z1, omega=1)[cause==2]
  C <- runif(m)*6
  X <- pmin(T,C)
  cause <- cause*(T <= C)
  tempdata <- data.frame(time=X,cause=cause, Z1)
  tempdata<-tempdata[order(tempdata$time),]
  # print(head(tempdata))
  return(list(data = tempdata,Z1=Z1))
}
# all0<-AIPWEIF.CIFnew(
#   tempdata = tempdata,
#   a=0,
#   beta.c=fitcensor,
#   tau0=tau1,
#   proA=proA,
#   censurv.u=censurv0,
#   pred1=pred01,
#   pred2=pred02,
#   fitcrr1=fitcrr1
# )
AIPWEIF.CIFnew<-function(tempdata,
                         a,
                         beta.c,
                         tau0,
                         proA,
                         censurv.u,
                         pred1,
                         pred2,
                         fitcrr1){
  # message("AIPWEIF.CIFnew")
  A<-tempdata$A; nn<-nrow(tempdata)
  pp <- which(names(tempdata)=='A')-1-2
  restriction<-(tempdata$time<=tau0)
  fcause1<-1*(tempdata$cause==1)*1*restriction
  censurv <- diag(censurv.u)
  beta <- matrix(beta.c,ncol=1)
  p <- nrow(beta)
  z <- as.matrix(cbind(tempdata[,c(3:(3+pp-1))],a,a*tempdata[,c(3:(3+pp-1))]))
  # print("cifnew0")
  ############### first term##################################
  IF1<-A^a*(1-A)^(1-a)*fcause1/proA^a/(1-proA)^(1-a)/censurv
  IF1[censurv==0]<-0
  ############### second term##################################
  # message("tau0:",tau0)
  tau0loc<-length(tempdata$time[tempdata$time<=tau0])
  # message("tau0loc:", tau0loc)
  EF1<-pred1[tau0loc,]
  IF2<--(2*a-1)*(A-proA)/proA^a/(1-proA)^(1-a)*EF1 # nn*1
  # print("cifnew1")
  ############### third term##################################
  delta_c <- (tempdata$cause == 0)*1 # nn*1
  #Calculate Lambda_c
  Y <-(do.call(cbind, replicate(nn, tempdata$time, simplify = FALSE)) >=
         do.call(rbind, replicate(nn, tempdata$time, simplify = FALSE))) *1
  #print(dim(Y))
  expz <- exp(z%*%beta)
  expzp <- do.call(cbind, replicate(p, expz, simplify = FALSE))
  zexpz <- expzp*z
  temp0 <- t(expz)%*%Y
  S0hat <- temp0 + (temp0==0)
  # print("cifnew2")
  #survival probability
  surv.prob.t<-t(1-pred1-pred2)#nn*nn_t
  # print("cifnew3")
  EF.T<-t((
    do.call(rbind,replicate(nn,EF1,simplify = FALSE))
  )- pred1)/(surv.prob.t)#nn*nn_t
  # print("cifnew4")
  # print(tau0loc+1)
  # print(nn)
  # View(EF.T)
  EF.T[,(tau0loc+1):nn]<-0
  # print("cifnew5")
  #martingale for censoring, dMc
  dNc = (do.call(cbind,replicate(nn,tempdata$time,simplify = FALSE))
         == do.call(rbind,replicate(nn,tempdata$time,simplify = FALSE)))*1
  # print("cifnew6")
  dNc[tempdata$cause != 0,] = 0
  # print("cifnew7")
  dlamb10c = colSums(dNc)/S0hat ##1*n
  dMc = dNc - Y*(do.call(rbind,replicate(nn,dlamb10c,simplify=FALSE)))*
    (do.call(cbind,replicate(nn,expz,simplify=FALSE))) #nn*nn_t
  censurv.G<-t(censurv.u)
  dMc_G = dMc/censurv.G*EF.T #nn*nn_t
  dMc_G[censurv.G==0]<-0
  dMc_G[surv.prob.t==0]<-0
  #change NA to 0 (small sample NC issue)
  dMc_G[is.na(dMc_G)] <-0
  dMc_G[dMc_G==-Inf] <- -1e100
  dMc_G[dMc_G==Inf] <- 1e1000

  Mc_G_EF= rowSums(dMc_G)
  IF3_temp = A^a*(1-A)^(1-a)/proA^a/(1-proA)^(1-a)*Mc_G_EF

  # print("cifnew5")
  IF3 = IF3_temp[1:nn]
  return(cbind(IF1,IF2,IF3))
}

betaUpdateBLWL <- function(xs,y.aipwe,w.aipwe,tune,pp,gamma,cj,j){
  # message("betaUpdateBLWL")
  # message("ncol(xs):", ncol(xs), " part2")
  new.gamma <- rep(0.001,pp)
  beta0 <- wsvmL1.group(x = xs, y = y.aipwe, w = w.aipwe, new.gamma)$sv
  # print('deep1')
  beta0 <- sapply(1:length(tune),function(t){
    if(any(abs(beta0[-1])>=10^(-5))){
      new.gamma <- t(weightAGB(cj,gamma,tune[t],beta0[-1],j))
      # new.gamma <- t(weightAGBCpp(cj,gamma,tune[t],beta0[-1],j))
      #weight.AGB(cj,gamma,tune[t],beta0[-1],j)
      if(all(new.gamma==Inf |abs(beta0[-1])<=10^(-5)
      )){
        # print('deep2')
        beta0.new <- matrix(c(1,rep(0,pp)),ncol=1)
      }else if(any(new.gamma == Inf | abs(beta0[-1]) <=10^(-5)
      )){
        keep.ind <- which(new.gamma != Inf & abs(beta0[-1])>10^(-5)
        )
        # print('deep3')
        beta0.tmp <- wsvmL1.group(x = as.matrix(xs[,keep.ind]), y = y.aipwe,
                                  w = w.aipwe, new.gamma[keep.ind])$sv
        # print('deep4')
        beta0.new <- matrix(0,ncol=1,nrow=pp+1)
        beta0.new[c(1,keep.ind+1)] <- beta0.tmp
      }else{
        beta0.new <- wsvmL1.group(x = xs, y = y.aipwe,
                                  w = w.aipwe, new.gamma)$sv
        # print('deep5')
      }
      beta0 <- beta0.new
    }else{
      beta0 <- c(1,beta0[-1])
    }
    # print('deep7')
    return(beta0)
  })
  # print("deep10")
  return(beta0)
}
#
# data_list = data.df.aipwe # this is for simulations
# data_list = train_aipwe; # data_list = uh # this is for RDA
# # pp.v = 30;#c(5,5,5,5,5,5)#ncov
# tau1 = as.numeric(t0_aipwe)
# tune = c(0.001,0.01,0.5,1,seq(0.1,400,length.out=16))
aipwe.fit <- function(data_list,
                      # beta1A, beta11,
                      pp.v, tau1, tune){
  print("starting aipwe.fit")
  m = dim(data_list[["data"]])[1]
  tempdata <- data_list$data
  Z1 <- data_list$Z1
  #group index
  j <- rep(1:length(pp.v),times=pp.v)
  # number of cov and group
  pp<- sum(pp.v); ngroup<- length(unique(j))

  # print("correct PS")
  # Correct PS
  proA<-glm(paste0(paste0('A~'),
                   paste0('V',c(1:pp), seperate='',collapse = '+')),
            family=binomial(link = "logit"),data=tempdata)$fitted.values
  # print("wrong PS")
  # # Wrong PS
  # proAwrong<- 0.1
  #Obtain censoring survival probability
  objC <- with(tempdata,Surv(time, (cause==0)))
  fit.cox<-coxph(as.formula(
    paste0(
      'objC~',
      paste0('V', c(1:((ncol(Z1)-1)/2)), seperate = '', collapse = '+'),
      '+A+',
      paste0('V', c((((ncol(Z1)-1)/2)+2):ncol(Z1)), seperate = '', collapse = '+')
    )
  ),data=tempdata)
  censoring <-survfit(fit.cox)
  fitcensor<-fit.cox$coef
  print("a=0")
  #a=0
  data0 <- tempdata
  #set trt = 0 and interactions all zero
  data0[,c((pp+3):(pp*2+3))] <- 0
  #deleting time and cause
  data0[,c(1,2)] <- NULL
  print("a=1")
  #a=1
  data1 <- tempdata
  data1[,(pp+3)] <- 1
  #trt = 1 and interactions are the same as the covariates
  data1[,c((pp+3+1):(pp*2+3))] <- tempdata[,c(3:(pp+2))]
  #deleting time and cause
  data1[,c(1,2)] <- NULL
  censoring0 <-survfit(fit.cox,newdata=data0)
  censoring1 <-survfit(fit.cox,newdata=data1)
  kmest <- stepfun(censoring$time, c(1, censoring$surv))
  censurv<-kmest(tempdata$time)
  censurv <- censurv0 <- censurv1 <- matrix(0,m,m)
  for(i in 1:m){
    kmest0 <- stepfun(censoring0$time, c(1, censoring0$surv[,i]))
    censurv0[,i]<-kmest0(tempdata$time)
    kmest1 <- stepfun(censoring1$time, c(1, censoring1$surv[,i]))
    censurv1[,i]<-kmest1(tempdata$time)
  }
  print("fit fine and gray model")
  # fit Fine and Gray model
  # View(tempdata$time)
  # View(tempdata$cause)
  # View(tempdata[,3:(ncol(Z1)+2)])
  # message("dim(Z1):", dim(Z1))
  # message("dim(tempdata):", dim(tempdata))
  # print('fg1')

  # # Create a model.matrix for factors:
  # cov1 <- model.matrix(~ .,
  #                      data = tempdata[,3:(ncol(Z1)+2)])[, -1]
  # print(tempdata$time)
  # print(tempdata$cause)
  # print(tempdata[,3:(ncol(Z1)+2)])
  fitcrr1<-crr(tempdata$time,
               tempdata$cause,
               cov1 = tempdata[,3:(ncol(Z1)+2)]) #,cov1) 
  # View(fitcrr1)
  # print('fg2')
  fitcrr2<-crr(
    tempdata$time,
    tempdata$cause,
    tempdata[,3:(ncol(Z1)+2)],
    failcode=2)
  # View(fitcrr2)
  # print('preds')
  pred01.temp<-predict(fitcrr1,as.matrix(data0)) #trt0
  pred02.temp<-predict(fitcrr2,as.matrix(data0)) #trt0
  pred11.temp<-predict(fitcrr1,as.matrix(data1)) #trt1
  pred12.temp<-predict(fitcrr2,as.matrix(data1)) #trt1
  #CIF
  #make nn_t *nn
  true.pred02 <- true.pred12 <- pred01<-pred11<-pred02<-pred12<-matrix(0,m,m)
  for (i in 1:m){
    pred01.temp.step <- stepfun(pred01.temp[,1], c(0, pred01.temp[,(1+i)]))
    pred01[,i]<-pred01.temp.step(tempdata$time)
    pred11.temp.step <- stepfun(pred11.temp[,1], c(0, pred11.temp[,(1+i)]))
    pred11[,i]<-pred11.temp.step(tempdata$time)
    pred02.temp.step <- stepfun(pred02.temp[,1], c(0, pred02.temp[,(1+i)]))
    pred02[,i]<-pred02.temp.step(tempdata$time)
    pred12.temp.step <- stepfun(pred12.temp[,1], c(0, pred12.temp[,(1+i)]))
    pred12[,i]<-pred12.temp.step(tempdata$time)
    # true.pred02[i,] <- (1-p)^exp(as.matrix(data0)%*%beta11)*
    #   (1-exp(-tempdata$time[i]*exp(as.matrix(data0)%*%beta11)))
    # true.pred12[i,] <- (1-p)^exp(as.matrix(data1)%*%beta11)*
    #   (1-exp(-tempdata$time[i]*exp(as.matrix(data1)%*%beta11)))
  }
  # print("correct propensity score")
  #correct propensity score
  #trt0
  all0<-AIPWEIF.CIFnew(
    tempdata,0,fitcensor,tau1,proA,censurv0,pred01,pred02,fitcrr1
  )
  #trt1
  all1<-AIPWEIF.CIFnew(
    tempdata,1,fitcensor,tau1,proA,censurv1,pred11,pred12,fitcrr1
  )
  print("t6")
  IPWIF0<-all0[,1]
  IPWIF1<-all1[,1]
  AIPWIF02<-IPWIF0+all0[,2]
  AIPWIF12<-IPWIF1+all1[,2]
  AIPWIF0<-AIPWIF02+all0[,3]
  AIPWIF1<-AIPWIF12+all1[,3]
  # message("pp:", pp)
  xs <- as.matrix(tempdata[,c(3:(pp+2))])
  xs.1 <- cbind(1,xs)
  # trueA <- as.numeric(xs.1%*%beta1A>0)
  CF.aipwe<-AIPWIF1-AIPWIF0
  w.aipwe<-abs(CF.aipwe); y.aipwe<-2*(CF.aipwe>0)-1
  cj <- rep(1,ngroup);
  gamma <- 0.5;
  # message("ncol(xs):",ncol(xs))
  beta0 <- betaUpdateBLWL(xs,y.aipwe,w.aipwe,tune,pp,gamma,cj,j)
  estA <- matrix(as.numeric(xs.1%*%beta0>0),m,length(tune))
  select <- (abs(beta0) > 10^(-5))
  select.noint <- select[-1,]
  group.select <- sapply(unique(j),function(i){
    sapply(1:length(tune),function(l){
      as.numeric(any(abs(select.noint[which(i==j),l])>10^(-5)))
    }
    )})
  print("emat")
  emat <- apply(abs(w.aipwe)*
                  (sign(xs.1%*%beta0)-y.aipwe)^2,2,sum)
  cn <- log(pp)
  bic <- log(emat/m) + cn*log(m)*apply(abs(beta0)>10^(-5),2,sum)/m
  #refit------------------
  min.bic <- which.min(bic)
  refit.group <- group.select[min.bic,]
  group.num <- which(refit.group==1)
  print("storing beta estimates into vector")
  #store the beta estimates in a vector of length pp+1
  beta.temp <- rep(0,pp+1)
  # print("t9")
  if(length(group.num)==0){
    print("length(group.num) == 0")
    beta.refit <- c(1,rep(0,pp))
    estA.refit <- as.numeric(xs.1%*%beta.refit>0)
    beta.temp <- beta.refit
  }else{
    # print("not 0")
    x.ind <- sapply(1:length(group.num),function(i){
      ind <- which(j==group.num[i])
      return(ind)
    })
    x.ind <- unlist(x.ind)
    beta.refit <- wsvmL1.group(
      x = xs[,x.ind],
      y = y.aipwe,
      w = w.aipwe,
      rep(0,length(x.ind))
    )$sv
    #beta.refit.norm <- beta.refit/sqrt(sum(beta.refit^2))
    beta.temp[c(1,x.ind+1)] <- beta.refit
  }
  # print(estA.refit)
  message("end of aipwe")
  return(beta.temp)
}

message("End of F02.ComparatorMethod_Functions.R")


# End of script -------------------------------------------------------------

