#########################################################################################################
BGPR=function(data,start,epsilon,maximum.iteration)
{
  eps=epsilon
  maxit=maximum.iteration
  
  library(MASS)
  data=data.frame(data)
  n=nrow(data)
  y1=as.matrix(data[,1])
  y2=as.matrix(data[,2])
  x=as.matrix(data[,-c(1,2)])
  x=cbind(1,x)
  p=ncol(x)
  
  beta10=start[1:p]
  beta20=start[(p+1):(2*p)]
  miu10=exp(x%*%beta10)
  miu20=exp(x%*%beta20)
  
  
  #Hasil
  Parameter=matrix(0,ncol=2*p+4,nrow=1)
  Std.Error=matrix(0,ncol=2*p+4,nrow=1)
  Z.Value=matrix(0,ncol=2*p+4,nrow=1)
  P.Value=matrix(0,ncol=2*p+4,nrow=1)
  Info=data.frame(matrix(0,ncol=1,nrow=8))
  colnames(Parameter)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  colnames(Std.Error)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  colnames(Z.Value)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  colnames(P.Value)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  y1hat=c()
  y2hat=c()
  
    fit=popBGPR(y1,y2,x,start,maxit,eps)
    parameter=as.matrix(fit$param)
    hes=fit$Hess
    inv=diag(pinv(-hes))
    se=round(as.matrix(sqrt(abs(inv))),5)
    z=round(parameter/se,5)
    con=fit$converged
    pv=round(2*pnorm(abs(z),lower.tail=FALSE),5)
    
    Parameter=parameter
    Std.Error=se
    Z.Value=z
    P.Value=pv
    
    p=ncol(x)
    p0=1
    x0=as.matrix(rep(1,n))
    par0=parameter[-c((2:(p)),((p+2):(p*2)))]
    
    ln.H1=round(Q_BGPR(y1,y2,x,parameter),3)
    ln.H0=round(Q_BGPR(y1,y2,x0,par0),3)
    
    G2=round(-2*(ln.H0-ln.H1),5)
    v=2*(ncol(data)-2)
    pvalF=round(pchisq((G2),v,lower.tail=FALSE),5)
    aic=round(n*(ln.H1)+(2*length(parameter)),3)
    
    Info=rbind(fit$iter,fit$converged,round(fit$norm.iter,3),ln.H1,ln.H0,G2,pvalF,aic)
    rownames(Info)=c("Number of Iteration","Converged/Not","Norm of Last Iteration","ln.H1","ln.H0","G^2","P.Value of F","AIC")
    colnames(Info)="Value"
    Info=noquote(Info)
    
    y1hat=exp(x%*%as.matrix(parameter[1:p]))
    y2hat=exp(x%*%as.matrix(parameter[(p+1):(2*p)]))
  
  
  Hasil=cbind(y1,round(miu10),round(y1hat),y2,round(miu20),round(y2hat))
  colnames(Hasil)=c("Y1","Y1.Pois-Reg","Y1.BGPR","Y2","Y2.Pois-Reg","Y2.BGPR")
  
  error1=(log(y1)-log(y1hat))
  error2=(log(y2)-log(y2hat))
  
  Tampil=data.frame(cbind(Parameter,Std.Error,Z.Value,P.Value))
  colnames(Tampil)=c('Parameter','Std.Error','Z.Value','P.Value')
  cat('   ','\n')
  cat('   ','\n')
  cat('******** Bivariate Generalized Poisson Regression ********','\n')
  cat('   ','\n')
  cat('-------------------------------------------','\n')
  cat('       Hasil Penghitungan Y.hat BGPR       ','\n')
  cat('-------------------------------------------','\n')
  print(Hasil)
  cat('   ','\n')
  cat('------------------------------------','\n')
  cat('       Hasil Uji Parsial BGPR       ','\n')
  cat('------------------------------------','\n')
  print(Tampil)
  cat('   ','\n')
  cat('----------------------------------------------------------','\n')
  cat('       Informsasi Iterasi & Hasil Uji Serentak BGPR       ','\n')
  cat('----------------------------------------------------------','\n')
  print(Info)
  list(Y1.hat=y1hat,Y2.hat=y2hat,Hasil=Hasil,Parameter=Parameter,Std.Error=Std.Error,Z.Value=Z.Value,P.Value=P.Value,Info=Info,AIC=aic,Error1=error1,Error2=error2)
}
