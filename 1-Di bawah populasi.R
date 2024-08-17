#Fungsi dibawah H1 (Populasi)
pop=function(y1,y2,x,w,param,num.location,maxit,eps)
{
  l=num.location
  norm.iter=100000
  iter=1
  Hess=matrix(0,length(param),length(param))
  
  y1=as.matrix(y1)
  y2=as.matrix(y2)
  x=as.matrix(x)
  w=as.matrix(w)
  p=ncol(x)
  n=nrow(x)
  
  param0=param
  be1=as.matrix(param[1:p])
  be2=as.matrix(param[(p+1):(2*p)])
  lamda0=param[(2*p+1)]
  alfa0=param[(2*p+2)]
  alfa1=param[(2*p+3)]
  alfa2=param[(2*p+4)]
  start=param
  
  while (iter<=maxit && ((alfa1<0 || alfa2<0) || norm.iter>eps))
  {
    param0=param
    be1=as.matrix(param[1:p])
    be2=as.matrix(param[(p+1):(2*p)])
    lamda0=param[(2*p+1)]
    alfa0=param[(2*p+2)]
    alfa1=param[(2*p+3)]
    alfa2=param[(2*p+4)]
    
    while (length(which((exp(x%*%be1)==Inf)==TRUE))!=0)
      be1=0.1*be1
    while (length(which((exp(x%*%be2)==Inf)==TRUE))!=0)
      be2=0.1*be2
    while (length(which((1/exp(x%*%be1)==Inf)==TRUE))!=0)
      be1=0.1*be1
    while (length(which((1/exp(x%*%be2)==Inf)==TRUE))!=0)
      be2=0.1*be2
    while (length(which((exp(x%*%be1)==0)==TRUE))!=0)
      be1=0.1*be1
    while (length(which((exp(x%*%be2)==0)==TRUE))!=0)
      be2=0.1*be2
    
    param0=as.matrix(c(be1,be2,lamda0,alfa0,alfa1,alfa2))
    
    miu1=exp(x%*%be1)
    miu2=exp(x%*%be2)
  
    
    ############################## TURUNAN PERTAMA #########################
    
    #Turunan Pertama Q terhadap Lamda0
    fungsi1=0
    fungsi2=0
    
    for (i in 1:n)
    {
      fungsi1=fungsi1+(((-n/lamda0)-(3*n))*w[i,l])
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi2=fungsi2+((((-y1[i]+k+1)/((miu1[i]-lamda0)+((y1[i]-k)*alfa1)))+((-y2[i]+k+1)/((miu2[i]-lamda0)+((y2[i]-k)*alfa2))))*(k-1)/(lamda0+(k*alfa0))*w[i,l])
        }
      }
    }
    diff1.Q.lamda0=fungsi1+fungsi2
    
    #Turunan Pertama Q terhadap Beta1
    fungsi1=0
    fungsi2=0
    fungsi3=0
    
    for (i in 1:n)
    {
      fungsi1=fungsi1+(1/miu1[i]*(miu1[i]*x[i,])*w[i,l])
      fungsi2=fungsi2+(miu1[i]*x[i,]*w[i,l])
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi3=fungsi3+(((y1[i]-k-1)*miu1[i]*x[i,]/((miu1[i]-lamda0)+((y1[i]-k)*alfa1)))*w[i,l])
        }
      }
    }
    diff1.Q.beta1=as.matrix(fungsi1+fungsi2+fungsi3)
    
    #Turunan Pertama Q terhadap Beta2
    fungsi1=0
    fungsi2=0
    fungsi3=0
    
    for (i in 1:n)
    {
      fungsi1=fungsi1+(1/miu2[i]*(miu2[i]*x[i,])*w[i,l])
      #print(i)
      #print((1/miu2[i]*(miu2[i]*x[i,])*w[i,l]))
      fungsi2=fungsi2+(miu2[i]*x[i,]*w[i,l])
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi3=fungsi3+(((y2[i]-k-1)*miu2[i]*x[i,]/((miu2[i]-lamda0)+((y2[i]-k)*alfa2)))*w[i,l])
        }
      }
    }
    diff1.Q.beta2=as.matrix(fungsi1+fungsi2+fungsi3)
    
    #Turunan Pertama Q terhadap Alfa1
    fungsi1=0
    fungsi2=0
    
    for (i in 1:n)
    {
      fungsi1=fungsi1+(-y1[i]*w[i,l])
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi2=fungsi2+((((y1[i]-k-1)*(y1[i]-k)/((miu1[i]-lamda0)+((y1[i]-k)*alfa1)))+1)*w[i,l])
        }
      }
    }
    diff1.Q.alfa1=fungsi1+fungsi2
    
    #Turunan Pertama Q terhadap Alfa2
    fungsi1=0
    fungsi2=0
    
    for (i in 1:n)
    {
      fungsi1=fungsi1+(-y2[i]*w[i,l])
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi2=fungsi2+(((y2[i]-k-1)*(y2[i]-k)/((miu2[i]-lamda0)+((y2[i]-k)*alfa2)))*w[i,l])
        }
      }
    }
    diff1.Q.alfa2=fungsi1+fungsi2
    
    #Turunan Pertama Q terhadap Alfa0
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+((-1+(k*(k-1)/(lamda0+(k*alfa0))))*w[i,l])
        }
      }
    }
    diff1.Q.alfa0=fungsi1
    
    gtheta=rbind(diff1.Q.beta1,diff1.Q.beta2,diff1.Q.lamda0,diff1.Q.alfa0,diff1.Q.alfa1,diff1.Q.alfa2)
    
    ############################## TURUNAN KEDUA #########################   
    #Turunan Kedua Q terhadap Lamda0
    fungsi1=0
    fungsi2=0
    
    for (i in 1:n)
    {
      fungsi1=fungsi1+(n/lamda0*w[i,l])
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi2=fungsi2+((((-y1[i]+k+1)/(((miu1[i]-lamda0)+((y1[i]-k)*alfa1))^2))+(((-y2[i]+k+1)/(((miu2[i]-lamda0)+((y2[i]-k)*alfa2))^2))+((k-1)/((lamda0+(k*alfa0))^2))))*w[i,l])
        }
      }
    }
    diff2.Q.lamda0=fungsi1+fungsi2
    
    #Turunan Kedua Q terhadap Lamda0 - Beta1
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(((y1[i]-k-1)*miu1[i]*x[i,]/(((miu1[i]-lamda0)+((y1[i]-k)*alfa1))^2))*w[i,l])
        }
      }
    }
    diff2.Q.lamda0.beta1=as.matrix(fungsi1)
    
    #Turunan Kedua Q terhadap Lamda0 - Beta2
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(((y2[i]-k-1)*miu2[i]*x[i,]/(((miu2[i]-lamda0)+((y2[i]-k)*alfa2))^2))*w[i,l])
        }
      }
    }
    diff2.Q.lamda0.beta2=as.matrix(fungsi1)
    
    #Turunan Kedua Q terhadap Lamda0 - Alfa1
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(((y1[i]-k)*(y1[i]-k-1)/(((miu1[i]-lamda0)+((y1[i]-k)*alfa1))^2))*w[i,l])
        }
      }
    }
    diff2.Q.lamda0.alfa1=fungsi1
    
    #Turunan Kedua Q terhadap Lamda0 - Alfa2
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(((y2[i]-k)*(y2[i]-k-1)/(((miu2[i]-lamda0)+((y2[i]-k)*alfa2))^2))*w[i,l])
        }
      }
    }
    diff2.Q.lamda0.alfa2=fungsi1
    
    #Turunan Kedua Q terhadap Lamda0 - Alfa0
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(((-k*(k-1))/((lamda0+(k*alfa0))^2))*w[i,l])
        }
      }
    }
    diff2.Q.lamda0.alfa0=fungsi1
    
    #Turunan Kedua Q terhadap Beta1
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(x[i,]%*%t(x[i,])*miu1[i]*(1+((-y1[i]+k+1)*miu1[i]/(((miu1[i]-lamda0)+((y1[i]-k)*alfa1))^2))+((y1[i]-k-1)/((miu1[i]-lamda0)+((y1[i]-k)*alfa1))))*w[i,l])
        }
      }
    }
    diff2.Q.beta1=as.matrix(fungsi1)
    
    #Turunan Kedua Q terhadap Beta2
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(x[i,]%*%t(x[i,])*miu2[i]*(1+((-y2[i]+k+1)*miu2[i]/(((miu2[i]-lamda0)+((y2[i]-k)*alfa2))^2))+((y2[i]-k-1)/((miu2[i]-lamda0)+((y2[i]-k)*alfa2))))*w[i,l])
        }
      }
    }
    diff2.Q.beta2=as.matrix(fungsi1)
    
    #Turunan Kedua Q terhadap Beta1 - Alfa1
    fungsi1=0
    
    for (i in 1:n)
    { kk=min(y1[i],y2[i])
    if (kk!=0)
    {for (k in  0:kk)
    {fungsi1=fungsi1+(((-y1[i]+k)*(y1[i]-k-1)*miu1[i]*x[i,]/(((miu1[i]-lamda0)+((y1[i]-k)*alfa1))^2))*w[i,l])}}
    }
    diff2.Q.beta1.alfa1=as.matrix(fungsi1)
    
    #Turunan Kedua Q terhadap Beta2 - Alfa2
    fungsi1=0
    
    for (i in 1:n)
    {kk=min(y1[i],y2[i])
    if (kk!=0)
    {for (k in  0:kk)
    {fungsi1=fungsi1+(((-y2[i]+k)*(y2[i]-k-1)*miu2[i]*x[i,]/(((miu2[i]-lamda0)+((y2[i]-k)*alfa2))^2))*w[i,l])}}
    }
    diff2.Q.beta2.alfa2=as.matrix(fungsi1)
    
    #Turunan Kedua Q terhadap Beta1 - Beta2
    diff2.Q.beta1.beta2=as.matrix(matrix(0,p,p))
    
    #Turunan Kedua Q terhadap Beta1 - Alfa2
    diff2.Q.beta1.alfa2=as.matrix(rep(0,p))
    
    #Turunan Kedua Q terhadap Beta2 - Alfa1
    diff2.Q.beta2.alfa1=as.matrix(rep(0,p))
    
    #Turunan Kedua Q terhadap Beta1 - Alfa0
    diff2.Q.beta1.alfa0=as.matrix(rep(0,p))
    
    #Turunan Kedua Q terhadap Beta2 - Alfa0
    diff2.Q.beta2.alfa0=as.matrix(rep(0,p))
    
    #Turunan Kedua Q terhadap Alfa1 - Alfa2
    diff2.Q.alfa1.alfa2=0
    
    #Turunan Kedua Q terhadap Alfa1 - Alfa0
    diff2.Q.alfa1.alfa0=0
    
    #Turunan Kedua Q terhadap Alfa2 - Alfa0
    diff2.Q.alfa2.alfa0=0
    
    #Turunan Kedua Q terhadap Alfa1
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+((((y1[i]-k)^2)*(-y1[i]+k+1)/(((miu1[i]-lamda0)+((y1[i]-k)*alfa1))^2))*w[i,l])
        }
      }
    }
    diff2.Q.alfa1=fungsi1
    
    #Turunan Kedua Q terhadap Alfa2
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+((((y2[i]-k)^2)*(-y2[i]+k+1)/(((miu2[i]-lamda0)+((y2[i]-k)*alfa2))^2))*w[i,l])
        }
      }
    }
    diff2.Q.alfa2=fungsi1
    
    #Turunan Kedua Q terhadap Alfa0
    fungsi1=0
    
    for (i in 1:n)
    {
      kk=min(y1[i],y2[i])
      if (kk!=0)
      {
        for (k in  0:kk)
        {
          fungsi1=fungsi1+(((-(k^2))*(k-1)/(lamda0+(k*alfa0)))*w[i,l])
        }
      }
    }
    diff2.Q.alfa0=fungsi1
    
    Hess.iter=matrix(0,length(param),length(param))
    
    Hess.iter[1:p,1:p]=diff2.Q.beta1
    Hess.iter[1:p,(p+1):(2*p)]=diff2.Q.beta1.beta2
    Hess.iter[(2*p+1),1:p]=diff2.Q.lamda0.beta1
    Hess.iter[1:p,(2*p+1)]=diff2.Q.lamda0.beta1
    Hess.iter[1:p,(2*p+2)]=diff2.Q.beta1.alfa0
    Hess.iter[(2*p+2),1:p]=diff2.Q.beta1.alfa0
    Hess.iter[1:p,(2*p+3)]=diff2.Q.beta1.alfa1
    Hess.iter[(2*p+3),1:p]=diff2.Q.beta1.alfa1
    Hess.iter[1:p,(2*p+4)]=diff2.Q.beta1.alfa2
    Hess.iter[(2*p+4),1:p]=diff2.Q.beta1.alfa2
    
    Hess.iter[(p+1):(2*p),(p+1):(2*p)]=diff2.Q.beta2
    Hess.iter[(p+1):(2*p),1:p]=diff2.Q.beta1.beta2
    Hess.iter[(2*p+1),(p+1):(2*p)]=diff2.Q.lamda0.beta2
    Hess.iter[(p+1):(2*p),(2*p+1)]=diff2.Q.lamda0.beta2
    Hess.iter[(p+1):(2*p),(2*p+2)]=diff2.Q.beta2.alfa0
    Hess.iter[(2*p+2),(p+1):(2*p)]=diff2.Q.beta2.alfa0
    Hess.iter[(p+1):(2*p),(2*p+3)]=diff2.Q.beta2.alfa1
    Hess.iter[(2*p+3),(p+1):(2*p)]=diff2.Q.beta2.alfa1
    Hess.iter[(p+1):(2*p),(2*p+4)]=diff2.Q.beta2.alfa2
    Hess.iter[(2*p+4),(p+1):(2*p)]=diff2.Q.beta2.alfa2
    
    Hess.iter[(2*p+1),(2*p+1)]=diff2.Q.lamda0
    Hess.iter[(2*p+1),(2*p+2)]=diff2.Q.lamda0.alfa0
    Hess.iter[(2*p+2),(2*p+1)]=diff2.Q.lamda0.alfa0
    Hess.iter[(2*p+1),(2*p+3)]=diff2.Q.lamda0.alfa1
    Hess.iter[(2*p+3),(2*p+1)]=diff2.Q.lamda0.alfa1
    Hess.iter[(2*p+1),(2*p+4)]=diff2.Q.lamda0.alfa2
    Hess.iter[(2*p+4),(2*p+1)]=diff2.Q.lamda0.alfa2
    
    Hess.iter[(2*p+2),(2*p+2)]=diff2.Q.alfa0
    Hess.iter[(2*p+2),(2*p+3)]=diff2.Q.alfa1.alfa0
    Hess.iter[(2*p+3),(2*p+2)]=diff2.Q.alfa1.alfa0
    Hess.iter[(2*p+2),(2*p+4)]=diff2.Q.alfa2.alfa0
    Hess.iter[(2*p+4),(2*p+2)]=diff2.Q.alfa2.alfa0
    
    Hess.iter[(2*p+3),(2*p+3)]=diff2.Q.alfa1
    Hess.iter[(2*p+3),(2*p+4)]=diff2.Q.alfa1.alfa2
    Hess.iter[(2*p+4),(2*p+3)]=diff2.Q.alfa1.alfa2
    
    Hess.iter[(2*p+4),(2*p+4)]=diff2.Q.alfa2
    
    param=param0-(ginv(Hess.iter)%*%as.matrix(gtheta))
    param[(2*p+2):(2*p+4)]=abs(param[(2*p+2):(2*p+4)])
    norm.iter=sqrt(sum((param-param0)^2))
    
    if (iter==1) {
      param.iter=param
    } else {
      param.iter=cbind(param.iter,param)
    }
    
    if (iter>1) {
      adj=ifelse(sum(param0-param.iter[,(iter-1)])!=0,"Adjusted","Non-Adjusted")
      adj=adj[1]
    } else {
      adj=ifelse(sum(param0-start)!=0,"Adjusted","Non-Adjusted")
      adj=adj[1]
    }
    
    
    if (adj=="Adjusted") {
      param.iter[,(iter-1)]=param0
    } else {
      param.iter[,(iter-1)]=param.iter[,(iter-1)]
    }
    
    iter=iter+1
    }
  converged=ifelse((iter-1)<maxit,"Converged","Not-Converged")
  list(iter=iter-1,converged=converged,param=param,param.iter=param.iter,norm.iter=norm.iter,Hess=Hess.iter)
}

Q=function(y1,y2,x,w,l,param)
{
  p=ncol(x)
  n=nrow(x)
  be1=as.matrix(param[1:p])
  miu1=exp((x)%*%be1)
  be2=as.matrix(param[(p+1):(2*p)])
  miu2=exp(x%*%be2)
  miu0=param[(2*p+1)]
  alfa0=param[(2*p+2)]
  alfa1=param[(2*p+3)]
  alfa2=param[(2*p+4)]
  A1=0; A2=0; A3=0; A4=0; A5=0; A6=0; A7=0; A8=0; B=0
  for (i in 1:n)
  {
    A1=A1+(log(miu0)*w[i,l])
    A2=A2+(log(abs(miu1[i]-miu0))*w[i,l])
    A3=A3+(log(abs(miu2[i]-miu0))*w[i,l])
    A4=A4+(miu0*w[i,l])
    A5=A5+((miu1[i]-miu0)*w[i,l])
    A6=A6+((miu2[i]-miu0)*w[i,l])
    A7=A7+(y1[i]*alfa1*w[i,l])
    A8=A8+(y2[i]*alfa2*w[i,l])
    kk=min(y1[i],y2[i])
    Bi=matrix(0,ncol=1,nrow=kk+1)
    for (k in 0:kk)
      {
      B1=((y1[i]-k-1)*log(abs(miu1[i]-miu0)+((y1[i]-k)*alfa1)))-lfactorial(y1[i]-k)+(k*(alfa1+alfa2-alfa0))
      B2=((y2[i]-k-1)*log(abs(miu2[i]-miu0)+((y2[i]-k)*alfa2)))-lfactorial(y2[i]-k)+((k-1)*(miu0+(k*miu0)))-lfactorial(k)
      Bi[k+1]=B1+B2
    }
    B=B+(sum(Bi)*w[i,l])
  }
  Q=A1+A2+A3-A4+A5+A6-A7-A8+B
  return(Q)
}
