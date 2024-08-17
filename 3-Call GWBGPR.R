#########################################################################################################
GWBGPR=function(data,coord,w,start,epsilon,maximum.iteration,Folder_simpan)
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
  w=as.matrix(w)
  
  beta10=start[1:p]
  beta20=start[(p+1):(2*p)]
  miu10=exp(x%*%beta10)
  miu20=exp(x%*%beta20)
  
  
  #Hasil
  Parameter=matrix(0,ncol=2*p+4,nrow=n)
  Std.Error=matrix(0,ncol=2*p+4,nrow=n)
  Z.Value=matrix(0,ncol=2*p+4,nrow=n)
  P.Value=matrix(0,ncol=2*p+4,nrow=n)
  Info=data.frame(matrix(0,ncol=8,nrow=n))
  colnames(Parameter)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  colnames(Std.Error)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  colnames(Z.Value)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  colnames(P.Value)=c(paste("Beta1",c(0:(p-1)),sep=""),paste("Beta2",c(0:(p-1)),sep=""),"Lamda0",paste("Alfa",c(0:2),sep=""))
  colnames(Info)=c("Number of Iteration","Converged/Not","Norm of Last Iteration","ln.H1","ln.H0","G^2","P.Value of F","AIC")
  y1hat=c()
  y2hat=c()
  
  for (l in 1:n)
  {
    print(noquote(paste("Analyzing location number : ",l)))
    fit=pop(y1,y2,x,w,start,l,maxit,eps)
    parameter=as.matrix(fit$param)
    hes=fit$Hess
    inv=diag(pinv(-hes))
    se=round(as.matrix(sqrt(abs(inv))),5)
    z=round(parameter/se,5)
    con=fit$converged
    pv=round(2*pnorm(abs(z),lower.tail=FALSE),5)
    
    Parameter[l,]=parameter
    Std.Error[l,]=se
    Z.Value[l,]=z
    P.Value[l,]=pv
    
    p=ncol(x)
    p0=1
    x0=as.matrix(rep(1,n))
    par0=parameter[-c((2:(p)),((p+2):(p*2)))]
    
    ln.H1=round(Q(y1,y2,x,w,l,parameter),3)
    ln.H0=round(Q(y1,y2,x0,w,l,par0),3)
    
    G2=round(-2*(ln.H0-ln.H1),5)
    v=2*(ncol(data)-2)
    pvalF=round(pchisq((G2),v,lower.tail=FALSE),5)
    aic=round((n*ln.H1)+(2*length(parameter)),3)
    
    Info[l,]=c(fit$iter,fit$converged,round(fit$norm.iter,3),ln.H1,ln.H0,G2,pvalF,aic)
    
    y1hat[l]=exp(x[l,]%*%as.matrix(parameter[1:p]))
    y2hat[l]=exp(x[l,]%*%as.matrix(parameter[(p+1):(2*p)]))
  }
  
  Lokasi=1:n
  Parameter=data.frame(Lokasi,Parameter)
  Std.Error=data.frame(Lokasi,Std.Error)
  Z.Value=data.frame(Lokasi,Z.Value)
  P.Value=data.frame(Lokasi,P.Value)
  Info=data.frame(Lokasi,Info)
  
  if (exists("Hasil_BGPR$Y1.hat")==TRUE && exists("Hasil_BGPR$Y2.hat")==TRUE)
  {
    Hasil=cbind(y1,round(miu10),round(Hasil_BGPR$Y1.hat),round(y1hat),y2,round(miu20),round(Hasil_BGPR$Y2.hat),round(y2hat))
    colnames(Hasil)=c("Y1","Y1.Pois-Reg","Y1.BGPR","Y1.GWBGPR","Y2","Y2.Pois-Reg","Y2.BGPR","Y2.GWBGPR")
  } else {
    Hasil=cbind(y1,round(miu10),round(y1hat),y2,round(miu20),round(y2hat))
    colnames(Hasil)=c("Y1","Y1.Pois-Reg","Y1.GWBGPR","Y2","Y2.Pois-Reg","Y2.GWBGPR")
  }
  
  error1=(log(y1)-log(y1hat))
  error2=(log(y2)-log(y2hat))
  
  write.table(Parameter,paste(Folder_simpan,"Parameter_GWBGPR.csv",sep=""),sep=";",row.names = F)
  write.table(Std.Error,paste(Folder_simpan,"StdError_GWBGPR.csv",sep=""),sep=";",row.names = F)
  write.table(Z.Value,paste(Folder_simpan,"Z-Value_GWBGPR.csv",sep=""),sep=";",row.names = F)
  write.table(P.Value,paste(Folder_simpan,"P-Value_GWBGPR.csv",sep=""),sep=";",row.names = F)
  write.table(Info,paste(Folder_simpan,"Info_GWBGPR.csv",sep=""),sep=";",row.names = F)
  
  Tampil=data.frame(cbind(Parameter,Std.Error,Z.Value,P.Value))
  colnames(Tampil)=c('Parameter','Std.Error','Z.Value','P.Value')
  cat('   ','\n')
  cat('   ','\n')
  cat('******** Geographically Weighted Bivariate Generalized Poisson Regression ********','\n')
  cat('   ','\n')
  cat('-------------------------------------------','\n')
  cat('       Hasil Penghitungan Y.hat GWBGPR       ','\n')
  cat('-------------------------------------------','\n')
  print(Hasil)
  cat('   ','\n')
  cat('------------------------------------','\n')
  cat('       Hasil Uji Parsial GWBGPR       ','\n')
  cat('------------------------------------','\n')
  print(Tampil)
  cat('   ','\n')
  cat('----------------------------------------------------------','\n')
  cat('       Informsasi Iterasi & Hasil Uji Serentak GWBGPR       ','\n')
  cat('----------------------------------------------------------','\n')
  print(Info)
  list(Y1.hat=y1hat,Y2.hat=y2hat,AIC=sum(as.numeric(Info[,9])),Info=Info,Error1=error1,Error2=error2)
}