#A collection of panel data descriptive functions 
#Xi Chen,  STATEC Luxembourg 
#Version 1.0, 15/09/2013


#Code testing infrastructure
#----------------------------------------------------------------------------
#rm(list = ls())
#library(plyr)
#library(plm)
#source("D:/Rcenter/MyRoutines/Prod_DGP_v1.1.r")
#set.seed(123456)
#Data=DGP(N=1000,Tt=10,a0=0.1,a1=0.95,sexi=0.1,w0bar=1,sew0=1,K0bar=5,seK0=1,delta=0.08,b0=0.5,b1=0.9,seI=0.01,c0=5,c1=0.5,c2=0.8,seL=0.9,d0=10,d1=0.5,d2=0.5,seM=0.9,e0=0.2,e1=0.3,e2=0.5,se=1,f0=0.5,f1=0.5,stata=T,stoTC=T,ltt=F)


FirstDes= function(data){
  #Some essential descriptive statistics 
  #Arguments: "data" is a data.frame
  v1=names(data)
  v2=dim(data)  
  v3=summary(data)
  v4=head(data)
  newList <- list("Dimension of dataset" = v2, "Names of variables" = v1, "Descriptive Statistics" = v3, "First few lines" = v4)
  return(newList)  
}



PanelAverage= function(data,N,T,over="id"){
  require(plyr)
  require(plm)
  #Averaging the panel data over ind/year 
  #Arguments: "data" is a data.frame; 
  #           "N" number of individuals
  #           "T" number of periods
  #           "over" the option on averaging function [id: Over individual averaging/year: Over time averaging]  
#Over time averaging
if (over=="year"){
  AData=apply(subset(data,id==1),2, mean)
  for (i in 2:N){
  AData_tmp=apply(subset(data, id==i), 2, mean)  
  AData=rbind(AData,AData_tmp)
  }
}
#Over individual averaging
if (over=="id"){  
  AData=apply(subset(data,year==1),2, mean)
  for (i in 2:T){
  AData_tmp=apply(subset(data, year==i), 2, mean)  
  AData=rbind(AData,AData_tmp)
  }
}
return(data.frame(AData, row.names = NULL))
}



PanelPlot= function(data,N,T,over="id",var=1){
  #Averaging the panel data over ind/year and plot
  #Arguments: "data" is a data.frame; 
  #           "N" number of individuals
  #           "T" number of periods
  #           "over" a string, the option on ploting [id: Cross-sectional distribution by averaging over time/year: Evoluation by averaging over individuals]
  #           "var" a numeric, select variable of interest
  AData=PanelAverage(data,N=N,T=T,over=over)
  varname=names(data)[var]
  #plot individual distribution
  if (over=="year"){
    plot1=plot(density(AData[,var]),main="Cross-sectional distribution by averaging over time",xlab="",ylab="",col=2,cex=0.1)
    newList <- list("the varaible of interest" = varname, plot1)
    return(newList)  
  }
  #plot time evoluation
  if (over=="id"){  
    plot2=plot(x=AData$year, y=AData[,var],main="Evoluation by averaging over individual",xlab="",ylab="",type="l",col=2,cex=0.1)
    newList <- list("the varaible of interest" = varname, plot2)
    return(newList)  
    }
  }
