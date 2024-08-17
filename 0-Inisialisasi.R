inisial=function(data,alfa0)
{
  library(pracma)
  library(MASS)
  data=data.frame(data)
  n=nrow(data)
  y1=as.matrix(data[,1])
  y2=as.matrix(data[,2])
  x=as.matrix(data[,-c(1,2)])
  
  #Inisialisasi Parameter dari Poisson Regression 
  f1=glm(formula=y1~x,family=quasipoisson(link=log))
  f2=glm(formula=y2~x,family=quasipoisson(link=log))
  beta10=f1$coefficients
  beta20=f2$coefficients
  x=as.matrix(cbind(rep(1,n),x))
  p=ncol(x)
  miu10=exp(x%*%beta10)
  miu20=exp(x%*%beta20)
  alfa1=summary(f1)$dispersion
  alfa2=summary(f2)$dispersion  
  alfa=c(alfa0,alfa1,alfa2)
  alfa=as.matrix(alfa)
  lamda0=cov(scale(y1,center=T,scale=T),scale(y2,center=T,scale=T))
  rownames(alfa)<-c('alfa0', 'alfa1','alfa2')
  start=as.matrix(c(beta10,beta20,lamda0,alfa))
  return(start)
}
