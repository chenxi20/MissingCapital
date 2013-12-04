#Simulation Experiments -- Simulated NLS for inverse production function 
#Xi Chen, STATEC Luxembourg
#Version 1.1, 08/09/2013 - 28/11/2013
#This code generate a fully simulated data set.

rm(list = ls())
library(plyr)
library(plm)
library(foreign)
library(xtable)
library(devtools)
library(optimx)
#Work directoty:
#setwd("D:/Rcenter/MyRoutines")
#source("D:/Rcenter/MyRoutines/Prod_DGP_v1.1.r")

set.seed(123456)
source_url("https://raw.github.com/chenxi20/MissingCapital/master/Prod_DGP_v1.1.r")

Data=DGP(N=500,Tt=5,a0=0.1,a1=0.95,sexi=0.1,w0bar=1,sew0=1,K0bar=5,seK0=1,delta=0.08,b0=0.5,b1=0.9,seI=0.01,c0=5,c1=0.5,c2=0.8,seL=0.9,d0=10,d1=0.5,d2=0.5,seM=0.9,e0=0.2,e1=0.3,e2=0.5,se=1,f0=0.5,f1=0.5,stata=F,stoTC=F,ltt=F)
names(Data)

Qtrue=function(theta,N,TT,data){
  alpha=theta[1]
  beta=theta[2]
  gamma=theta[3]
  d=data
  delta=0.08
  t=c(1:TT)
  res2_NT=0
  for (n in 1:N){
    dn=subset(d,id==n)
    res2_sum=matrix(0,ncol=1,nrow=TT-1)
    Kn=(dn$Y*dn$u^(-1)*dn$L^(-beta)*dn$M^(-gamma))^(1/alpha)
    res2_sum=res2_sum+(dn$I[2:TT]-Kn[2:TT]+(1-delta)*Kn[1:(TT-1)])*(dn$I[2:TT] -Kn[2:TT]+(1-delta)*Kn[1:(TT-1)])
    res2_NT=rbind(res2_NT,res2_sum)
  }
  res2=sum(res2_NT)/(2*N*(TT-1))
  return(res2)
}

Qnaive=function(theta,N,TT,data){
  alpha=theta[1]
  beta=theta[2]
  gamma=theta[3]
  d=data
  delta=0.08
  t=c(1:TT)
  res2_NT=0
  for (n in 1:N){
    dn=subset(d,id==n)
    res2_sum=matrix(0,ncol=1,nrow=TT-1)
    Kn=(dn$Y*dn$L^(-beta)*dn$M^(-gamma))^(1/alpha)
    res2_sum=res2_sum+(dn$I[2:TT]-Kn[2:TT]+(1-delta)*Kn[1:(TT-1)])*(dn$I[2:TT] -Kn[2:TT]+(1-delta)*Kn[1:(TT-1)])
    res2_NT=rbind(res2_NT,res2_sum)
  }
  res2=sum(res2_NT)/(2*N*(TT-1))
  return(res2)
}




#Code testing infrastructure:

#alpha=0.2
#beta=0.3
#gamma=0.5
#see=1
#data=Data
#delta=0.08
#N=500
#TT=5
#t=c(1:TT)
#S=2



#simulated objective function over S replications VERSION 2
Qsim2=function(theta,S,N,TT,data){
  alpha=theta[1]
  beta=theta[2]
  gamma=theta[3]
  see=1
  #d=data
  delta=0.08
  t=c(1:TT)
  h_S=matrix(0,ncol=1,nrow=N*(TT-1))
  for (s in 1:S){
    set.seed(1234+s)
    v_s=rnorm(N*TT,0,1)
    data_s=data.frame(cbind(data,v_s))
    for (n in 1:N){
       dn_s=subset(data_s,id==n)
       Kn_s=(dn_s$Y*exp(dn_s$v_s)^(-see)*dn_s$L^(-beta)*dn_s$M^(-gamma))^(1/alpha)
       hn_s_tmp = Kn_s[2:TT]-(1-delta)*Kn_s[1:(TT-1)]
       if(n==1){hn_s=hn_s_tmp}
       if(n>1){hn_s=c(hn_s,hn_s_tmp)}
      }
    h_S=h_S+as.matrix(hn_s)
   }
  H_S=h_S/S
  res2=sum((data$I[2:TT]-H_S)*(data$I[2:TT]-H_S))
  return(res2)
}


v0=c(0.2,0.3,0.5)
Qtrue(theta=v0,N=500,TT=5,data=Data)
Qnaive(theta=v0,N=500,TT=5,data=Data)
#Qsim(theta=v0,S=50,N=500,TT=5,data=Data)
Qsim2(theta=v0,S=100,N=500,TT=5,data=Data)

est1 = optimx(v0,Q,N=500,TT=5,data=Data, method ="BFGS",control=list(abstol=1e-10,reltol=1e-10)) 
est2 = optimx(v0,Qsim2,S=5,N=500,TT=5,data=Data, method ="BFGS",control=list(abstol=1e-6,reltol=1e-6)) 

write.table(est2, file = "est2.csv") 

