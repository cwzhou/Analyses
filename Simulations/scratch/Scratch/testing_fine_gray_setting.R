q = 0.3
truemod = TRUE
eta = c(0,-1,-1.4, 0.5)

# gen.data<-function(n,q=0.8,eta=c(0,1,-1,0.5),c0=8,truemod=TRUE){
  Z1<-runif(n)#runif(n,-2,2)
  Z2<-runif(n)#runif(n,-2,2)
  if(truemod) eze<-exp(cbind(1,Z1,Z2)%*%eta[-length(eta)])
  if(!truemod) eze<-exp(cbind(1,Z1,Z2,Z1^2)%*%eta)
  A<-rbinom(n,1,eze/(1+eze))
  
  eza1<-exp(-Z1+A*(Z1-Z2))      # exponential of linear predictors in PH model of risk 1
  eza2<-exp(-A*(2+Z1+2*Z2))     # rate of expoential dist. for risk2
  
  p1<-1-(1-q)^eza1
  rsk1<-rbinom(n,1,p1)
  
  U<-runif(n)
  T1<--log(1-(1-(1-p1*U)^(1/eza1))/q)
  T2<--log(1-U)/eza2
  T<-rsk1*T1+(1-rsk1)*T2
  ebs<-rsk1+2*as.numeric(rsk1==0)
  
  C<-runif(n,0,c0)
  obst<-apply(cbind(T,C),1,min)
  delta<-as.numeric(T<=C)
  status<-delta*ebs
  return(data.frame(id=1:n,time=obst,status=status,Z1=Z1,Z2=Z2,A=A,piz=eze/(1+eze)))
# }

#######################################
ZZ<-expand.grid(rep(list(seq(-2,2,.01)),2))
#Z1<-ZZ[,1]
#Z2<-ZZ[,2]
#######################################
resltabl<-function(q,alps,prefix,mydata,truth,Z1=ZZ[,1],Z2=ZZ[,2]){
  Beta<-read.table(paste(prefix,"_beta.txt",sep=""))
  Fval<-read.table(paste(prefix,"_fval.txt",sep=""))
  F1est<-read.table(paste(prefix,"_F1est.txt",sep=""))
  F2est<-read.table(paste(prefix,"_F2est.txt",sep=""))
  attain<-read.table(paste(prefix,"_attain.txt",sep=""))
  prop<-read.table(paste(prefix,"_prop.txt",sep=""))
  rectm<-read.table(paste(prefix,"rectm.txt",sep=""))
  
  subdata1<-mydata[attain[,1]==1]
  subdata2<-mydata[attain[,2]==1]
  subdata3<-mydata[attain[,3]==1]
  
  Beta1<-Beta[attain[,1]==1,1:3]
  Beta2<-Beta[attain[,2]==1,4:6]
  Beta3<-Beta[attain[,3]==1,7:9]
  
  F1est1<-F1est[attain[,1]==1,1]
  F1est2<-F1est[attain[,2]==1,2]
  F1est3<-F1est[attain[,3]==1,3]
  
  est1<-c(apply(Beta1,2,mean),mean(F1est1))
  est2<-c(apply(Beta2,2,mean),mean(F1est2))
  est3<-c(apply(Beta3,2,mean),mean(F1est3))
  
  bias1<-est1-truth[1,]
  bias2<-est2-truth[2,]
  bias3<-est3-truth[3,]
  
  stder1<-c(apply(Beta1,2,sd),sd(F1est1))
  stder2<-c(apply(Beta2,2,sd),sd(F1est2))
  stder3<-c(apply(Beta3,2,sd),sd(F1est3))
  
  blk1<-cbind(truth[1,],est1,stder1)
  blk2<-cbind(truth[2,],est2,stder2)
  blk3<-cbind(truth[3,],est3,stder3)
  ##########
  est2.1<-matrix(nrow=nrow(Beta1),ncol=2)
  pcc1<-numeric()
  for(i in 1:nrow(Beta1)){
    est2.1[i,]<-search.beta(beta=as.numeric(Beta1[i,]),Z1=Z1,Z2=Z2,n=length(Z1),q=q,t0=2)
    Zmati<-as.matrix(cbind(1,subdata1[[i]][,c("Z1","Z2")]))
    estregm<-(Zmati%*%as.numeric(Beta1[i,]))>0
    tregm<-(Zmati%*%as.numeric(truth[1,-4]))>0
    pcc1[i]<-mean(estregm==tregm)
  }
  est2.2<-matrix(nrow=nrow(Beta2),ncol=2)
  pcc2<-numeric()
  for(i in 1:nrow(Beta2)){
    est2.2[i,]<-search.beta(beta=as.numeric(Beta2[i,]),Z1=Z1,Z2=Z2,n=length(Z1),q=q,t0=2)
    Zmati<-as.matrix(cbind(1,subdata2[[i]][,c("Z1","Z2")]))
    estregm<-(Zmati%*%as.numeric(Beta2[i,]))>0
    tregm<-(Zmati%*%as.numeric(truth[2,-4]))>0
    pcc2[i]<-mean(estregm==tregm)
  }
  
  est2.3<-matrix(nrow=nrow(Beta3),ncol=2)
  pcc3<-numeric()
  for(i in 1:nrow(Beta3)){
    est2.3[i,]<-search.beta(beta=as.numeric(Beta3[i,]),Z1=Z1,Z2=Z2,n=length(Z1),q=q,t0=2)
    Zmati<-as.matrix(cbind(1,subdata3[[i]][,c("Z1","Z2")]))
    estregm<-(Zmati%*%as.numeric(Beta3[i,]))>0
    tregm<-(Zmati%*%as.numeric(truth[3,-4]))>0
    pcc3[i]<-mean(estregm==tregm)
  }
  
  est5<-c(mean(est2.1[est2.1[,2]<=alps[1],1]),mean(est2.2[est2.2[,2]<=alps[2],1]),mean(est2.3[est2.3[,2]<=alps[3],1]))
  stder5<-c(sd(est2.1[est2.1[,2]<=alps[1],1]),sd(est2.2[est2.2[,2]<=alps[2],1]),sd(est2.3[est2.3[,2]<=alps[3],1]))
  
  est6<-c(mean(pcc1[est2.1[,2]<=alps[1]]),mean(pcc2[est2.2[,2]<=alps[2]]),mean(pcc3[est2.3[,2]<=alps[3]]))
  stder6<-c(sd(pcc1[est2.1[,2]<=alps[1]]),sd(pcc2[est2.2[,2]<=alps[2]]),sd(pcc3[est2.3[,2]<=alps[3]]))
  
  blk4<-round(rbind(est5,est6,apply(attain,2,mean)),3)
  blk5<-round(rbind(stder5,stder6),3)
  
  out1<-round(rbind(blk1,cbind(rep(NA,3),blk4[,1],c(blk5[,1],NA))),3)
  out2<-round(rbind(blk2,cbind(rep(NA,3),blk4[,2],c(blk5[,2],NA))),3)
  out3<-round(rbind(blk3,cbind(rep(NA,3),blk4[,3],c(blk5[,3],NA))),3)
  
  rownames(out1)<-rownames(out2)<-rownames(out3)<-c("beta0","beta1","beta2","hatF1","F1","PCD","F2<alp")
  colnames(out1)<-colnames(out2)<-colnames(out3)<-c("Truth","Est","Std")
  return(list(alps=alps,alpha1tab=out1,alpha2tab=out2,alpha3tab=out3))
}
