#Monte Carlo Experiments -- Data Generating Process with function DGP(), old name: MCE_DGP_v1.1.R
#Xi Chen, STATEC Luxembourg
#Version 1.1, 08/09/2013 - 27/11/2013
#This code generate a fully simulated production data set.


#Code testing infrastructure
#----------------------------------------------------------------------------
#rm(list = ls())
#library(plyr)
#library(plm)
#library(foreign)
#library(xtable)
#set.seed(12345)
#N=3  	    #Number of individual
#Tt=2		    #Number of periods
#J=8		    #Number of variables
#a0=0.1	    #Productivity Markov process
#a1=0.9
#sexi=0.1
#w0bar=1	  #Initial productivity level in log
#sew0=1
#K0bar=5	  #Initial capital stock in log-normal distribution
#seK0=1
#delta=0.08	#Fixed depreciation rate
#b0=0.4		  #Investment rule
#b1=0.6
#c0=50		  #Labor demande
#c1=0.5
#c2=0.8
#seI=0.01
#seL=0.1
#d0=1		    #Material demande
#d1=0.5
#d2=0.5
#seM=0.5
#e0=0.2     #Production function
#e1=0.4
#e2=0.7
#se=1
#stoTC=T    #Stochastic technical change
#ltt=F      #Linear time trend
#----------------------------------------------------------------------------


DGP = function(N,Tt,a0,a1,sexi,w0bar,sew0,K0bar,seK0,delta,b0,b1,seI,c0,c1,c2,seL,d0,d1,d2,seM,e0,e1,e2,se,f0,f1,stata=F,stoTC=T,ltt=F){
  require(plyr)
  require(plm)
  require(foreign)
  
  ind=1:N           #individual indicator                        
  year=1:Tt         #time indicator        
  J=8               #number of variables    
  
  #Generating productivity
  #----------------------------------------------------------------------------
  #intial value of prodcutivity (normal distributed)
  wi0=as.matrix(rnorm(N,w0bar,sew0))
  #productivity 1-order Markov process
  wit=wi0
  for (t in 1 : (Tt-1)){
    wit_tmp=a0+a1*wit[,t]+rnorm(N, 0, sexi)
    wit=cbind(wit, wit_tmp)}
  
  
  #Generating investment and capital stock
  #----------------------------------------------------------------------------
  #intial value of capital stock (normal distributed)
  Ki0=as.matrix(exp(rnorm(N,K0bar,seK0)))
  Ii0= (Ki0^b0)*(exp(wi0)^b1)*exp(rnorm(N,0,seI))
  #Investment rule and Capital stock:
  Kit=Ki0
  Iit=Ii0
  for (t in 1 : (Tt-1)){
    Kit_tmp= (1-delta)*Kit[,t]+Iit[,t]
    Iit_tmp= (Kit_tmp^b0)*(exp(wit[,t+1])^b1)*exp(rnorm(N,0,seI))
    Kit=cbind(Kit, Kit_tmp)
    Iit=cbind(Iit, Iit_tmp)
  }
  #Iit=(Kit^b0)*(exp(wit)^b1)*exp(rnorm(T*N,0,0.01))
  
  
  #Generating Labor and Materail demande
  #----------------------------------------------------------------------------
  Lit0=c0*(Kit[,1]^c1)*(exp(wit[,1])^c2)*exp(rnorm(N,0,seL))
  Mit0=d0*(Kit[,1]^d1)*(exp(wit[,1])^d2)*exp(rnorm(N,0,seM))
  Lit=Lit0
  Mit=Mit0
  for (t in 2 : Tt){
    Lit_tmp=c0*(Kit[,t]^c1)*(exp(wit[,t])^c2)*exp(rnorm(N,0,seL))
    Mit_tmp=d0*(Kit[,t]^d1)*(exp(wit[,t])^d2)*exp(rnorm(N,0,seM))
    Lit=cbind(Lit, Lit_tmp)
    Mit=cbind(Mit, Mit_tmp)
  }
  
  
  #Generating output
  #----------------------------------------------------------------------------
  #without time trend without stochastic technical change
  if (stoTC==F & ltt==F){ Yit=(Kit^e0)*(Lit^e1)*(Mit^e2)}
  
  #without time trend with stochastic technical change
  if (stoTC==T & ltt==F){ Yit=(Kit^e0)*(Lit^e1)*(Mit^e2)*exp(wit)}
  
  #without time trend with stochastic technical change
  if (stoTC==F & ltt==T){ Yit=exp(f0+f1*year)*(Kit^e0)*(Lit^e1)*(Mit^e2)}
  
  #with time trend with stochastic technical change
  if (stoTC==T & ltt==T){ Yit=exp(f0+f1*year)*(Kit^e0)*(Lit^e1)*(Mit^e2)*exp(wit)}
  
  #Transform the matrix to panel data structure
  #----------------------------------------------------------------------------
  Panel=matrix (NA, nrow=Tt,ncol=J)
  Panel[,1]=t(rep(1,Tt))
  Panel[,2]=t(c(1:Tt))
  Panel[,3]=t(Yit[1,])
  Panel[,4]=t(Lit[1,])
  Panel[,5]=t(Mit[1,])
  Panel[,6]=t(Kit[1,])
  Panel[,7]=t(Iit[1,])
  Panel[,8]=t(wit[1,])
  for (i in 2:N){
    Panel_tmp=matrix (NA, nrow=Tt,ncol=J)
    Panel_tmp[,1]=t(rep(i,Tt))
    Panel_tmp[,2]=t(c(1:Tt))
    Panel_tmp[,3]=t(Yit[i,])
    Panel_tmp[,4]=t(Lit[i,])
    Panel_tmp[,5]=t(Mit[i,])
    Panel_tmp[,6]=t(Kit[i,])
    Panel_tmp[,7]=t(Iit[i,])
    Panel_tmp[,8]=t(wit[i,])
    Panel=rbind(Panel,Panel_tmp)
  }
  u=exp(rnorm(N*Tt,0,se))
  PanelData=data.frame(cbind(Panel, Panel[,3]*u,log(Panel[,3]), log(Panel[,3]*u), log(Panel[,4]), log(Panel[,5]), log(Panel[,6]), log(Panel[,7]), u))
  PanelData=rename(PanelData, c("V1"="id","V2"="year","V3"="Ystar","V4"="L","V5"="M","V6"="K","V7"="I","V8"="lnw","V9"="Y","V10"="lnYstar","V11"="lnY","V12"="lnL","V13"="lnM","V14"="lnK","V15"="lnI"))
  
  #the option "stata" allows us to save the data in STATA form: 
  if (stata==TRUE) {write.dta(PanelData, file = "PanelData_Prod_DGP.dta")}
  return(PanelData)
}
