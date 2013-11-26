#It is working but slow!!!

rm(list = ls())
library(plyr)
library(plm)
library(foreign)
library(xtable)
library(mgcv)
#Work directoty:
getwd()
setwd("/Users/chen/Desktop/Works_25112013/MyRoutines_25112013")
set.seed(123456)

source("/Users/chen/Desktop/Works_25112013/MyRoutines_25112013/MCE_DGP_v1.0.r")

Data=DGP(N=500,Tt=5,a0=0.1,a1=0.95,sexi=0.1,w0bar=1,sew0=1,K0bar=5,seK0=1,delta=0.08,b0=0.5,b1=0.9,seI=0.01,c0=5,c1=0.5,c2=0.8,seL=0.9,d0=10,d1=0.5,d2=0.5,seM=0.9,e0=0.2,e1=0.3,e2=0.5,se=0.1,f0=0.5,f1=0.5,stata=F,stoTC=F,ltt=F) 

YY=(Data$K^0.2)*(Data$L^0.3)*(Data$M^0.5)
Data=data.frame(cbind(Data,YY))


#simulated objective function over S replications
Q=function(theta,S,N,TT,data){
alpha=theta[1]
beta=theta[2]
gamma=theta[3]
#a0=theta[4]
#a1=theta[5]  
d=data
delta=0.08
t=c(1:TT)
res2_NT=0
  for (n in 1:N){
    dn=subset(d,id==n)
    res2_sum=matrix(0,ncol=1,nrow=TT-1)
#    for (s in 1:S){
#    set.seed(123+s)  
#    U_s=rnorm(T,0,0)    
#    K=(dn$YY/exp(a0+a1*t+U_s)*dn$L^(-beta)*dn$M^(-gamma))^(1/alpha)  
     #Kn=(dn$YY*dn$L^(-beta)*dn$M^(-gamma))^(1/alpha)  
     #res2_sum=res2_sum+(dn$I[2:TT]-Kn[2:TT]+(1-delta)*Kn[1:(TT-1)])*(dn$I[2:TT] -Kn[2:TT]+(1-delta)*Kn[1:(TT-1)])  
     lnK=(log(dn$YY)-beta*log(dn$L)-gamma*log(dn$M))/alpha
     res2_sum=res2_sum+(dn$I[2:TT]-exp(lnK[2:TT])+(1-delta)*exp(lnK[1:(TT-1)]))*(dn$I[2:TT] -exp(lnK[2:TT])+(1-delta)*exp(lnK[1:(TT-1)]))  
#    }
#    res2_mean=res2_sum/S
     res2_NT=rbind(res2_NT,res2_sum)    
  }
res2=sum(res2_NT)/(N*(TT-1))
return(res2)
}
v0=c(0.1,0.2,0.3)
Q(theta=v0,S=1,N=500,TT=5,data=Data)

bnlm = optim(v0,Q,S=1,N=500,TT=5,data=Data, method ="L-BFGS-B",control = list(factr=1e10,pgtol=0))

