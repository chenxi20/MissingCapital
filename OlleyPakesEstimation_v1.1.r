#Production function Estimation (by Olley and Pakes, 1996) with panel bootstraping
#Xi Chen,  STATEC Luxembourg 
#Version 1.1, 15/09/2013-04/12/2013
#Functions for estimating the production function, the package "plyr", "mgcv", "plm" is required


Estimat = function(panel, GAM=F){
  #Estimation with Olley and Pakes method without bootstraping
  #"panel" should be an objet of class "data.frame" and variables are in LOG and named as "ln__" 
  #"GAM" use of mgcv for generalized addditive model estimation (slower)
  require(plyr)
  require(mgcv)
  require(plm)
  d<- panel
  #The first-stage estimation
  if (GAM==F){
  lnK2=d$lnK*d$lnK
  lnK3=d$lnK*d$lnK*d$lnK
  lnI2=d$lnI*d$lnI
  lnI3=d$lnI*d$lnI*d$lnI
  reg1= lm (lnY ~ lnL + lnM +lnK +lnK2+ lnK3 +lnI +lnI2+ lnI3, data=d )}
  else{reg1= gam(lnY ~ lnL + lnM + s(lnK,lnI),data=d)}
  #summary(reg1)
  bl=reg1$coefficients[2]
  bm=reg1$coefficients[3]
  phi=reg1$fitted.values-bl*d$lnL-bm*d$lnM
  panel2=pdata.frame(cbind(d,phi), c("id","year"))
  
  #The second-stage estimation
  sres2 = function (betak,data){
    #the current value
    w = data$phi-betak*data$lnK
    #the lag value
    wl = lag(w)
    wl2 = wl*wl
    wl3 = wl*wl*wl
    balanced=subset(data.frame(cbind(data,w,wl,wl2,wl3)),year!=1)
    reg = lm(w~wl+wl2+wl3,data=balanced)
    summary(reg)$coefficients
    
    res2=reg$residuals*reg$residuals
    Sres2=sum(res2)
    return(Sres2) } 
#plot for visualizing the optimization
#  test=matrix(NA, ncol=2,nrow=150)
#  for (i in -50:100){
#    test[i+50,1]=i/100+0.025
#    test[i+50,2]=sres2(bk=(i/100+0.025),data=panel2)
#  }
#  plot(test[,1],test[,2])  
  b0=0
  bnlm = optim(b0,sres2,data=panel2, method ="L-BFGS-B")
  bk=c("lnK"=bnlm$par) 
  return(c(bl,bm,bk))
} 


EstBoot = function(NBoot,RBoot,panelBoot){
  #Panel Bootraping for estimating the standard errors of estimator
  #"NBoot" number of individuals
  #"RBoot" number of replications
  #"panelBoot" the data set in question
N=NBoot
R=RBoot
d=panelBoot
ResBoot=Estimat(panel=d, GAM=F)
for (j in 1:R){ 
  set.seed(1234+j)
  random=as.integer(runif(N, min=1, max=N))
  resampled=subset(d,id==random[1])
  for (i in 2:N){
    Resampled=subset(d,id==random[i])
    resampled=rbind(resampled,Resampled)
  }   
  ResBoot=cbind(ResBoot,Estimat(panel=resampled, GAM=F))
} 
#write.csv(ResBoot, file = "ResBoot.csv")
SdBoot=matrix(NA,ncol=1,nrow=dim(ResBoot)[1])
for (i in 1:dim(ResBoot)[1]){
  SdBoot[i]=sd(ResBoot[i,])
}  
BootRes=cbind(ResBoot[,1],SdBoot,ResBoot[,1]/SdBoot)
return(BootRes)
}

EstProd = function (panel,GAM=F){
  #productivity estimation (include the exogeneous shocks)
  #"panel" should be an objet of class "data.frame" and variables are in LOG and named as "ln__" 
  #"GAM" use of mgcv for generalized addditive model estimation (slower)
  d=panel
  gam=GAM
  est=Estimat(panel=d, GAM=gam)
  productiv=d$lnY-est[1]*d$lnL-est[2]*d$lnM-est[3]*d$lnK
  return(productiv)
}

