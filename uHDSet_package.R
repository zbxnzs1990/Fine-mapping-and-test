library('MASS')
library("foreach")
library("doParallel")

single_marker<-function(y,x){
n=length(y);
y=(y-mean(y))/sd(y);
# x=(x-mean(x))/sd(x)
ro2=((t(y)%*%x)/n)^2;
stat=(n-2)*ro2/(1-ro2);
# stat(is.nan(stat))=0;
p=pchisq(stat,df=1,lower.tail = FALSE);
# p(is.nan(p))=1;
return(p)
}

Si_dev_v<-function(x,p,mu){
y=x;
id=(x>mu);
y[id]=p*sign(x[id])*abs(x[id])^(p-1);
y[!id]=p*(x[!id]^2/(2*mu)+mu/2)^(p-1)*x[!id]/mu;
return(y)
}

Si_v<-function(x,p,mu){
y=x;
id=(x>mu);
y[id]=abs(x[id])^p;
y[!id]=(x[!id]^2/(2*mu)+mu/2)^p;
return(y)
}

smooth_fun<-function(x,A,Kin,y,lambda,p,mu,weight){
f=t(y-A%*%x)%*%Kin%*%(y-A%*%x)+lambda*sum(Si_v(x*weight,p,mu));
return(f)
}

smooth_fun_o<-function(x,A,Kin,y,lambda,p,mu){
f=t(y-A%*%x)%*%Kin%*%(y-A%*%x)+lambda*sum(Si_v(x,p,mu));
return(f)
}

smooth_fun_dev<-function(x,A,Kin,y,lambda,p,mu,weight){
df=2*t(A)%*%Kin%*%(A%*%x-y)+lambda*sum(Si_dev_v(x*weight,p,mu));
return(df)
}

smooth_fun_dev_o<-function(x,A,Kin,y,lambda,p,mu){
df=2*t(A)%*%Kin%*%(A%*%x-y)+lambda*sum(Si_dev_v(x,p,mu));
return(df)
}

USR_unweighted<-function(A,y,lambda,Kin,p){
  # % Function for USR unweighted method, which is a special case of USR and
  # % computed faster than USR
  # 
  # % Input
  # % A: independent matrix (genetic score matrix)
  # % y: response variable (phenotype data)
  # % lambda: tuning parameter that control the sparsity level
  # % Kin: inverse of Kinship matrix
  # % p: Lp norm based regularization
  # % weight: weight coefficient for weighted method
  # 
  # % Output
  # % solution: the sparse solution of unweighted USR
  
  n=nrow(A);
  n_snp=ncol(A);
  
  esp0=1e-3;
  r=1;
  ro=0.8;
  delta=0.3;
  mu=0.1;
  
  x0=L_half_F(A,y,lambda,Kin);
  support=which(x0!=0);
  if(length(support)>0){
  B=as.matrix(A[,support]);
  xb=x0[support];
  mu_loop=0;
  x_mu=xb;
  d_mu=as.matrix(1);
  
#   print(xb)
#   print(class(B))
#   print(support)
#   print(B[1])
  
  # Lower bound
  L1=matrix(0,length(support),1);
  L2=(lambda*p*sqrt(length(support))/(2*max(svd(B)$d)*max(svd(Kin)$d)*sqrt(t(y)%*%Kin%*%y)))^(1/(1-p));
  for (j in 1:length(support)){
    L1[j]=(lambda*p*(1-p)/(2*t(B[,j])%*%Kin%*%B[,j]))^(1/(2-p));
  }
  
  while (norm(d_mu,"F")>1e-5*n_snp){
    x_mu_old=x_mu;
    g0=smooth_fun_dev_o(x_mu,B,Kin,y,lambda,p,mu);
    d0=-g0;
    k=0;
    x=x_mu;
    d=d0;
    g=g0;
    sk=as.matrix(1);
    while (norm(sk,"F")>1e-5*n_snp){
      d_old=d;
      g_old=g;
      x_old=x;
      i=0;
      while (smooth_fun_o(x+(ro^i)*d,B,Kin,y,lambda,p,mu)>smooth_fun_o(x,B,Kin,y,lambda,p,mu)+delta*(ro^i)*t(g)%*%d){
        i=i+1;
      }
      x=x+(ro^i)*d;
      g=smooth_fun_dev_o(x,B,Kin,y,lambda,p,mu);
      yk=g-g_old;
      sk=x-x_old;
      tk=esp0*norm(g,"F")^r+max(0,t(sk)%*%yk/(t(sk)%*%sk));
      zk=yk-tk*sk;
      beta=t(g)%*%zk/(t(d_old)%*%zk)-2*(t(zk)%*%zk)*(t(g)%*%d_old)/(t(d_old)%*%zk)^2;
      theta=t(g)%*%d_old/(t(d_old)%*%zk);
      d=-g+beta[1]*d_old+theta[1]*zk;
      k=k+1;
    }
    x_mu=x;
    for (i in 1:length(support)){
      if (abs(x_mu[i])<max(L1[i],L2)){
        x_mu[i]=0;
      }
    }
    d_mu=x_mu-x_mu_old;
    mu=mu/5;
    mu_loop=mu_loop+1;
  }
  solution=matrix(0,n_snp,1);
  solution[support]=x_mu;
}else{solution=x0;}
  return(solution)
}

USR<-function(A,y,lambda,Kin,p,weight){
  # % Function for USR unweighted method, which is a special case of USR and
  # % computed faster than USR
  # 
  # % Input
  # % A: independent matrix (genetic score matrix)
  # % y: response variable (phenotype data)
  # % lambda: tuning parameter that control the sparsity level
  # % Kin: inverse of Kinship matrix
  # % p: Lp norm based regularization
  # % weight: weight coefficient for weighted method
  # 
  # % Output
  # % solution: the sparse solution of unweighted USR
  
  n=nrow(A);
  n_snp=ncol(A);
  
  esp0=1e-3;
  r=1;
  ro=0.8;
  delta=0.3;
  mu=0.1;
  
  x0=L_half_Fw(A,y,lambda,Kin,weight);
  support=which(x0!=0);
  B=as.matrix(A[,support]);
  xb=x0[support];
  weight_supp=weight[support];
  mu_loop=0;
  x_mu=xb;
  d_mu=as.matrix(1);
  
  # Lower bound
  L1=matrix(0,length(support),1);
  L2=(lambda*p*sqrt(length(support))/(2*max(svd(B)$d)*max(svd(Kin)$d)*sqrt(t(y)%*%Kin%*%y)))^(1/(1-p));
  for (i in 1:length(support)){
    L1[i]=(lambda*weight_supp[i]*p*(1-p)/(2*t(B[,i])%*%Kin%*%B[,i]))^(1/(2-p));
  }
  
  while (norm(d_mu,"F")>1e-5*n_snp){
    x_mu_old=x_mu;
    g0=smooth_fun_dev(x_mu,B,Kin,y,lambda,p,mu,weight_supp);
    d0=-g0;
    k=0;
    x=x_mu;
    d=d0;
    g=g0;
    sk=as.matrix(1);
    while (norm(sk,"F")>1e-5*n_snp){
      d_old=d;
      g_old=g;
      x_old=x;
      i=0;
      while (smooth_fun(x+(ro^i)*d,B,Kin,y,lambda,p,mu,weight_supp)>smooth_fun(x,B,Kin,y,lambda,p,mu,weight_supp)+delta*(ro^i)*t(g)%*%d){
        i=i+1;
      }
      x=x+(ro^i)*d;
      g=smooth_fun_dev(x,B,Kin,y,lambda,p,mu,weight_supp);
      yk=g-g_old;
      sk=x-x_old;
      tk=esp0*norm(g,"F")^r+max(0,t(sk)%*%yk/(t(sk)%*%sk));
      zk=yk-tk*sk;
      beta=t(g)%*%zk/(t(d_old)%*%zk)-2*(t(zk)%*%zk)*(t(g)%*%d_old)/(t(d_old)%*%zk)^2;
      theta=t(g)%*%d_old/(t(d_old)%*%zk);
      d=-g+beta[1]*d_old+theta[1]*zk;
      k=k+1;
    }
    x_mu=x;
    for (i in 1:length(support)){
      if (abs(x_mu[i])<max(L1[i],L2)){
        x_mu[i]=0;
      }
    }
    d_mu=x_mu-x_mu_old;
    mu=mu/5;
    mu_loop=mu_loop+1;
  }
  solution=matrix(0,n_snp,1);
  solution[support]=x_mu;
  return(solution)
}

Scaled_USR<-function(A,y,lambda,Kin,p){
  # % Function for USR unweighted method, which is a special case of USR and
  # % computed faster than USR
  # 
  # % Input
  # % A: independent matrix (genetic score matrix)
  # % y: response variable (phenotype data)
  # % lambda: tuning parameter that control the sparsity level
  # % Kin: inverse of Kinship matrix
  # % p: Lp norm based regularization
  # % weight: weight coefficient for weighted method
  # 
  # % Output
  # % solution: the sparse solution of unweighted USR
  
  n=nrow(A);
  n_snp=ncol(A);
  
  esp0=1e-3;
  r=1;
  ro=0.8;
  delta=0.3;
  mu=0.1;
  
  x0=L_half_F(A,y,lambda,Kin);
  support=which(x0!=0);
  if(length(support)>0){
    B=as.matrix(A[,support]);
    xb=x0[support];
    mu_loop=0;
    x_mu=xb;
    d_mu=as.matrix(1);
    
    #   print(xb)
    #   print(class(B))
    #   print(support)
    #   print(B[1])
    
    # Lower bound
    L1=matrix(0,length(support),1);
    L2=(lambda*p*sqrt(length(support))/(2*max(svd(B)$d)*max(svd(Kin)$d)*sqrt(t(y)%*%Kin%*%y)))^(1/(1-p));
    for (j in 1:length(support)){
      L1[j]=(lambda*p*(1-p)/(2*t(B[,j])%*%Kin%*%B[,j]))^(1/(2-p));
    }
    
    while (norm(d_mu,"F")>1e-5*n_snp){
      x_mu_old=x_mu;
      g0=smooth_fun_dev_o(x_mu,B,Kin,y,lambda,p,mu);
      d0=-g0;
      k=0;
      x=x_mu;
      d=d0;
      g=g0;
      sk=as.matrix(1);
      while (norm(sk,"F")>1e-5*n_snp){
        d_old=d;
        g_old=g;
        x_old=x;
        i=0;
        while (smooth_fun_o(x+(ro^i)*d,B,Kin,y,lambda,p,mu)>smooth_fun_o(x,B,Kin,y,lambda,p,mu)+delta*(ro^i)*t(g)%*%d){
          i=i+1;
        }
        x=x+(ro^i)*d;
        g=smooth_fun_dev_o(x,B,Kin,y,lambda,p,mu);
        yk=g-g_old;
        sk=x-x_old;
        tk=esp0*norm(g,"F")^r+max(0,t(sk)%*%yk/(t(sk)%*%sk));
        zk=yk-tk*sk;
        beta=t(g)%*%zk/(t(d_old)%*%zk)-2*(t(zk)%*%zk)*(t(g)%*%d_old)/(t(d_old)%*%zk)^2;
        theta=t(g)%*%d_old/(t(d_old)%*%zk);
        d=-g+beta[1]*d_old+theta[1]*zk;
        k=k+1;
      }
      x_mu=x;
      for (i in 1:length(support)){
        if (abs(x_mu[i])<max(L1[i],L2)){
          x_mu[i]=0;
        }
      }
      d_mu=x_mu-x_mu_old;
      mu=mu/5;
      mu_loop=mu_loop+1;
    }
    solution=matrix(0,n_snp,1);
    solution[support]=x_mu;
  }else{solution=x0;}
  return(solution)
}

half_threshold<-function(x,lambda,mu){
  n=length(x);
  y=matrix(0,nrow=n,ncol=1);
  t=54^(1/3)/4*(lambda*mu)^(2/3);
  for (i in 1:n){
    if (abs(x[i])>t){
      y[i]=f_l_mu(x[i],lambda,mu);
    }
  }
  return(y)
}

half_threshold_w<-function(x,lambda,mu){
  #   x=matrix(c(0.01,2,3),nrow=3,ncol=1);
  #   lambda=1;
  #   mu=1;
  n=length(x);
  y=matrix(0,nrow=n,ncol=1);
  t=54^(1/3)/4*(lambda*mu)^(2/3);
  for (i in 1:n){
    if (abs(x[i])>t[i]){
      y[i]=f_l_mu(x[i],lambda[i],mu);
    }
  }
  return(y)
}

soft_threshold<-function(x,y){
  if (x>0 & y<x){
    z=x-y
  }else if (x<0 & y<abs(x)){
    z=x+y;
  }else{
    z=0;
  }
  return(z)
}

f_l_mu<-function(x,lambda,mu){
  y=2/3*x*(1+cos(2*pi/3-2/3*acos(lambda*mu/8/(abs(x)/3)^(3/2))));
  return(y)
}

L_half_F<-function(A,y,lambda,Kin,option="typical",k=1){
  # % Function for unified L1/2 thresholding algorithm
  # 
  # % Input
  # % A: independent matrix (genetic score matrix)
  # % y: response variable (phenotype data)
  # % lambda: tuning parameter that control the sparsity level
  # % Kin: inverse of Kinship matrix
  # % option: if option='typical', solve the unified L1/2 thresholding
  # % algorithm; if option~='typical', solve the L1/2 thresholding under
  # % predetermined sparsity level
  # % k: predetermined sparsity level (only works when option~='typical')
  # 
  # % Output
  # % solution: the sparse solution of L1/2 norm based regularization
  
  n=nrow(A);
  n_snp=ncol(A);
  NormA=max(eigen(A%*%t(A))$values);
  mu=1/NormA
  x=matrix(0,n_snp,1);
  esp=1;step=0;
  if (option=="typical"){
    print("typical")
    while (esp>0.001 && step<100){
      x_old=x;
      B=x_old+mu*t(A)%*%Kin%*%(y-A%*%x_old);
      x=half_threshold(B,lambda,mu);
      step=step+1;
      esp=norm(x-x_old,"F");
    }
  }else{  
    while (esp>0.001 && step<100){
      x_old=x;
      B=x_old+mu*t(A)%*%Kin%*%(y-A%*%x_old);
      Babs=abs(B);
      SB=sort(Babs,decreasing=TRUE);
      lambda=sqrt(96)/9*NormA*SB[k]^1.5;
      x=half_threshold(B,lambda,mu);
      step=step+1;
      esp=norm(x-x_old,"F");
    }
  }
  return(x)
}

Scaled_L_half<-function(A,y,lambda0,option="typical",k=2){
# % Function for Scaled L1/2 thresholding algorithm
# 
# % Input
# % A: independent matrix (genetic score matrix)
# % y: response variable (phenotype data)
# % lambda0: initial tuning parameter that control the sparsity level
# 
# % Output
# % solution: the sparse solution of Scaled L1/2 norm based regularization
# % sigma: estimated variance
# % lambda: final lambda

n=nrow(A);
n_snp=ncol(A);
NormA=max(eigen(A%*%t(A))$values);
mu=1/NormA;
sigma=norm(y,"F")/sqrt(n);
lambda=lambda0;
x=matrix(0,n_snp,1);
esp=1;step=0;

if (option=="typical"){
#   print("typical")
  while (esp>0.001 & step<100){
    lambda_old=lambda;
    sigma_old=sigma;
    x_old=x;
    B=x_old+mu*t(A)%*%(y-A%*%x_old)/(2*n*sigma);
    x=half_threshold(B,lambda,mu/(2*n*sigma));
    step=step+1;
    
    sigma=norm(y-A%*%x,"F")/sqrt(n);
    lambda=lambda0*sigma;
    esp=norm(x-x_old)+abs(sigma-sigma_old)+abs(lambda-lambda_old)
  }
}else{
  print("sparse level")
  while (esp>0.001 && step<100){
    lambda_old=lambda;
    sigma_old=sigma;
    x_old=x;
    B=x_old+mu*t(A)%*%(y-A%*%x_old)/(2*n*sigma);
    Babs=abs(B);
    SB=sort(Babs,decreasing=TRUE);
    lambda=sqrt(96)/9*NormA*SB[k]^1.5;
    x=half_threshold(B,lambda*(2*n*sigma),mu/(2*n*sigma));
    step=step+1;
    sigma=norm(y-A%*%x,"F")/sqrt(n);
    lambda=lambda0*sigma
    esp=norm(x-x_old)+abs(sigma-sigma_old)+abs(lambda-lambda_old);
  }
}

solution=list(x=x,sigma=sigma,lambda=lambda);
return(solution)
}

Scaled_L_half_F<-function(A,y,lambda0,Kin,option="typical",k=1){
  # % Function for unified L1/2 thresholding algorithm
  # 
  # % Input
  # % A: independent matrix (genetic score matrix)
  # % y: response variable (phenotype data)
  # % lambda: tuning parameter that control the sparsity level
  # % Kin: inverse of Kinship matrix
  # % option: if option='typical', solve the unified L1/2 thresholding
  # % algorithm; if option~='typical', solve the L1/2 thresholding under
  # % predetermined sparsity level
  # % k: predetermined sparsity level (only works when option~='typical')
  # 
  # % Output
  # % solution: the sparse solution of L1/2 norm based regularization
  
  n=nrow(A);
  n_snp=ncol(A);
  NormA=max(eigen(A%*%t(A))$values);
  mu=1/NormA
  sigma=norm(y,"F")/sqrt(n)[1]
  lambda=lambda0
  x=matrix(0,n_snp,1);
  esp=1;step=0;
  
  if (option=="typical"){
    print("typical")
    while (esp>0.001 && step<100){
      lambda_old=lambda
      sigma_old=sigma
      x_old=x;
      B=x_old+mu*t(A)%*%Kin%*%(y-A%*%x_old)/(2*n*sigma);
      x=half_threshold(B,lambda,mu/(2*n*sigma));
      step=step+1;
      
      sigma=sqrt(t(y-A%*%x)%*%Kin%*%(y-A%*%x))[1]/sqrt(n)
      lambda=lambda0*sigma
      esp=norm(x-x_old)+abs(sigma-sigma_old)+abs(lambda-lambda_old)
    }
  }else{  
    print("sparse level")
    B=mu*t(A)%*%Kin%*%y/(2*n*sigma);
    Babs=abs(B);
    SB=sort(Babs,decreasing=TRUE);
    lambda0=sqrt(96)/9*NormA*SB[k]^1.5;
    while (esp>0.001 && step<100){
      lambda_old=lambda;
      sigma_old=sigma;
      x_old=x;
      B=x_old+mu*t(A)%*%Kin%*%(y-A%*%x_old)/(2*n*sigma);
      Babs=abs(B);
      SB=sort(Babs,decreasing=TRUE);
      lambda=sqrt(96)/9*NormA*SB[k]^1.5;
      x=half_threshold(B,lambda*(2*n*sigma),mu/(2*n*sigma));
      step=step+1;
      sigma=sqrt(t(y-A%*%x)%*%Kin%*%(y-A%*%x))[1]/sqrt(n)
      lambda=lambda0*sigma
      esp=norm(x-x_old)+abs(sigma-sigma_old)+abs(lambda-lambda_old)
    }
  }
  return(list(x=x,sigma=sigma,lambda=lambda))
}

p_value_L_half<-function(y,X){
  # % Algorithm for L_half p-value
  # 
  # % Input
  # % X: independent matrix
  # % y: response variable
  # 
  # % Output
  # % P: the p-value vector for each selected independent variables
  # % clc;clear;
  # % n=10;
  # % p=50;
  # % X=randn(n,p);
  # % y=randn(n,1);
  
  n=nrow(X);
  p=ncol(X);
  lambda0=2*sqrt(log(p)/n);
  solution=Scaled_L_half(X,y,lambda0);
  gama=solution$lambda/solution$sigma;
  
  Z=t(X)%*%X/n;
  
  M=diag(p);
  if (gama<1){
    fr<-function(x){
      y=t(x)%*%Z%*%x;
      return(y)
    }
    gr<-function(x){
      y=2*Z%*%x;
      return(y)
    }
    
    for (i in 1:p){
      x0=matrix(0,p,1);
      x0[i]=1/Z[i,i];
      ei=matrix(0,p,1);
      ei[i]=1;
      ui=rbind(-Z,Z);
      ci=rbind(-ei-gama*matrix(1,p,1),ei-gama*matrix(1,p,1));
      if (length(which(ui%*%x0-ci<=0))>0){
#         print(which(ui%*%x0-ci<=0))
#           print(i)
      }else{
        opt=constrOptim(x0,fr,gr,ui,ci);
        M[,i]=as.matrix(opt$par);
      }
    }
    if (sum(sum(is.nan(M)))>0){
      print('M is infeasible')
      M=eye(p,p);
    }
  }
  theta=solution$x+M%*%t(X)%*%(y-X%*%solution$x)/n;
  P=theta;
  MZM=M%*%Z%*%t(M);
  for (i in 1:p){
    P[i]=2*(1-pnorm(sqrt(n)*abs(theta[i])/(solution$sigma*sqrt(MZM[i,i])+1e-6)));
  }
  return(P);
}

Scaled_Lp<-function(A,y,p_norm,lambda0,option="typical",k=2){
# % Function of Scaled Lp norm regularization
# 
# % Input
# % A: independent matrix
# % y: response variable
# % lambda0: initial tuning parameter that control the sparsity level lambda0=2*sqrt(log(p)/n);
# % p_norm: Lp(0<p<1) norm
# 
# % Output
# % solution: the sparse solution of Scaled Lp
# % sigma: estimated variance
# % lambda: final lambda
  L_half=Scaled_L_half(A,y,lambda0,option,k);
  n=nrow(A);
  n_snp=ncol(A);
  
  esp0=1e-3;
  r=1;
  ro=0.8;
  delta=0.3;
  mu=0.1;
  Kin=diag(n);
  x0=L_half$x;
  sigma=L_half$sigma;
  lambda=L_half$lambda;
  support=which(x0!=0);
  B=as.matrix(A[,support]);
  xb=x0[support];
  mu_loop=0;
  x_mu=xb;
  d_mu=as.matrix(n_snp);
  
  
  # % Lower bound
  L1=matrix(0,length(support),1);
  L2=L1;
  if(ncol(B)>0){
    L2=(lambda*p_norm*sqrt(length(support))/(2*norm(B,"F")*sqrt(sum(y*y))))^(1/(1-p_norm));
    for (i in 1:length(support)){
      L1[i]=(lambda*p_norm*(1-p_norm)/(2*sum(B[,i]*B[,i])))^(1/(2-p_norm));
    }
  
  
  while (d_mu>1e-5*n_snp){
    lambda_old=lambda;
    sigma_old=sigma;
    x_mu_old=x_mu;
    x=x_mu;
    g0=smooth_fun_dev_Scaled(x,B,y,lambda,sigma,n,p_norm,mu);
    
    d0=-g0;
    k=0;    
    d=d0;
    g=g0;
    sk=matrix(1,1,1);
    while (norm(sk,"F")>1e-5*n_snp){
      d_old=d;
      g_old=g;
      x_old=x;
      i=0;
      while (smooth_fun_Scaled(x+(ro^i)*d,B,y,lambda,sigma,n,p_norm,mu)>smooth_fun_Scaled(x,B,y,lambda,sigma,n,p_norm,mu)+delta*(ro^i)*sum(g*d)){
        i=i+1;
      }
      x=x+(ro^i)*d;
      g=smooth_fun_dev_Scaled(x,B,y,lambda,sigma,n,p_norm,mu);
      yk=g-g_old;
      sk=x-x_old;
      tk=esp0*norm(g,"F")^r+max(0,-sum(sk*yk)/sum(sk*sk));
      zk=yk-tk*sk;
      beta=sum(g*zk)/sum(d_old*zk)-2*sum(zk*zk)*sum(g*d_old)/sum(d_old*zk)^2;
      theta=sum(g*d_old)/sum(d_old*zk);
      d=-g+beta*d_old+theta*zk;
      k=k+1;
    }
    x_mu=x;
    for (i in 1:length(support)){
      if (abs(x_mu[i])<max(L1[i],L2)){
        x_mu[i]=0;
      }
    }
    sigma=norm(y-B%*%x,"F")/sqrt(n);
    lambda=lambda0*sigma;
    d_mu=norm(x_mu-x_mu_old,"F")+abs(sigma-sigma_old)+abs(lambda-lambda_old);
    
    mu=mu/5;
    
    mu_loop=mu_loop+1;
  }
  sol=matrix(0,n_snp,1);
  sol[support]=x_mu;
  }else{
    sol=matrix(0,n_snp,1);
  }
  solution=list(x=sol,sigma=sigma,lambda=lambda);
  return(solution)
}

uHDSet<-function(y,X,Kin,p_norm=0.3,option="typical",k=2,permute=100000){
  # Algorithm for uHDSet test
  # 
  # Input
  # X: independent matrix (genotype matrix)
  # y: response variable
  # Kinship: inverse of Kinship matrix (diagnal element should be 0.5)
  # 
  # Output
  # P.uFineMap: the p-value vector for the uFineMap test
  # P.uHDSet: p-value for the uHDSet test
  
  Kin=solve(2*Kinship)
  n=nrow(X);
  p=ncol(X);
  lambda0=2*sqrt(log(p)/n)
  solution1=Scaled_L_half_F(X,y,lambda0,Kin,option="t",k=(floor(0.03*p)+1))
  solution2=Scaled_L_half_F(X,y,lambda0,Kin,option="typical",k=(floor(0.02*p)+1))
  x=(solution1$x+solution2$x)/2
  sigma=(solution1$sigma+solution2$sigma)/2
  Z=t(X)%*%Kin%*%X/n;
  
  M=InverseLinfty(Z, n, resol=1.3, mu=NULL, maxiter=50, threshold=1e-2, verbose=TRUE)
  
  #sigma=(sigma+sum(y^2)/n)/2
  theta=x+M%*%t(X)%*%Kin%*%(y-X%*%x)/n
  P=theta
  MZM=M%*%Z%*%t(M)
  denom=sigma^2*diag(MZM)
  Stat1=sqrt(n)*abs(theta)/(sqrt(denom)+1e-8)
  P=2*(1-pnorm(Stat1));
  
  eKin=eigen(Kin)
  eKinl=diag(sqrt(eKin$values))
  eKinQ=eKin$vectors
  Kin_half=eKinQ%*%eKinl%*%t(eKinQ)
  
  A=sigma*M%*%t(X)%*%Kin_half
  
  ZA=A%*%t(A)
  eZA=eigen(ZA)
  rankZ=sum(eZA$values>1E-5)
  ZA.d=eZA$values[1:rankZ]
  
  eZAl=diag(sqrt(ZA.d))
  U1=eZA$vectors[,1:rankZ]
  B=U1%*%eZAl
  
  Z=n*solve(t(B)%*%B)%*%t(B)%*%theta
  maxZ=max(Z^2)
  sumZ=sum(Z^2)
  P_sum=pchisq(sumZ,rankZ,lower.tail=FALSE)
  P_multiple=1-pchisq(maxZ,1)^rankZ
  P_adjust=1-pchisq(max(Stat1^2),1)^rankZ
  P_adjust_sum=pchisq(sum(Stat1^2),rankZ,lower.tail=FALSE)
  
  Simu=permute;
  max_distribution=matrix(0,Simu,1)
  
  denom=sigma^2*diag(MZM)
  multi=mvrnorm(Simu,matrix(0,p,1),sigma^2*MZM)
  for (i in 1:Simu){
    max_distribution[i]=max(multi[i,]^2/denom)
  }
  Stat=max(n*theta^2/denom)
  P.uHDSet=sum(Stat<max_distribution)/Simu
  
  return(list(P.uFineMap = P, P.uHDSet = P.uHDSet))
}

uHDSet_parallel<-function(y,X,Kin,p_norm=0.3,option="typical",k=2,permute=100000){
  # Algorithm for uHDSet test
  # 
  # Input
  # X: independent matrix (genotype matrix)
  # y: response variable
  # Kinship: inverse of Kinship matrix (diagnal element should be 0.5)
  # 
  # Output
  # P.uFineMap: the p-value vector for the uFineMap test
  # P.uHDSet: p-value for the uHDSet test
  
  Kin=solve(2*Kinship)
  n=nrow(X);
  p=ncol(X);
  lambda0=2*sqrt(log(p)/n)
  solution1=Scaled_L_half_F(X,y,lambda0,Kin,option="t",k=(floor(0.03*p)+1))
  solution2=Scaled_L_half_F(X,y,lambda0,Kin,option="typical",k=(floor(0.02*p)+1))
  x=(solution1$x+solution2$x)/2
  sigma=(solution1$sigma+solution2$sigma)/2
  Z=t(X)%*%Kin%*%X/n;
  
  M=InverseLinfty_parallel(Z, n, resol=1.3, mu=NULL, maxiter=50, threshold=1e-2, verbose=TRUE)
  
  #sigma=(sigma+sum(y^2)/n)/2
  theta=x+M%*%t(X)%*%Kin%*%(y-X%*%x)/n
  P=theta
  MZM=M%*%Z%*%t(M)
  denom=sigma^2*diag(MZM)
  Stat1=sqrt(n)*abs(theta)/(sqrt(denom)+1e-8)
  P=2*(1-pnorm(Stat1));
  
  eKin=eigen(Kin)
  eKinl=diag(sqrt(eKin$values))
  eKinQ=eKin$vectors
  Kin_half=eKinQ%*%eKinl%*%t(eKinQ)
  
  A=sigma*M%*%t(X)%*%Kin_half
  
  ZA=A%*%t(A)
  eZA=eigen(ZA)
  rankZ=sum(eZA$values>1E-5)
  ZA.d=eZA$values[1:rankZ]
  
  eZAl=diag(sqrt(ZA.d))
  U1=eZA$vectors[,1:rankZ]
  B=U1%*%eZAl
  
  Z=n*solve(t(B)%*%B)%*%t(B)%*%theta
  maxZ=max(Z^2)
  sumZ=sum(Z^2)
  P_sum=pchisq(sumZ,rankZ,lower.tail=FALSE)
  P_multiple=1-pchisq(maxZ,1)^rankZ
  P_adjust=1-pchisq(max(Stat1^2),1)^rankZ
  P_adjust_sum=pchisq(sum(Stat1^2),rankZ,lower.tail=FALSE)
  
  Simu=permute;
  max_distribution=matrix(0,Simu,1)
  
  denom=sigma^2*diag(MZM)
  multi=mvrnorm(Simu,matrix(0,p,1),sigma^2*MZM)
  for (i in 1:Simu){
    max_distribution[i]=max(multi[i,]^2/denom)
  }
  Stat=max(n*theta^2/denom)
  P.uHDSet=sum(Stat<max_distribution)/Simu
  
  return(list(P.uFineMap = P, P.uHDSet = P.uHDSet))
}

smooth_fun_Scaled<-function(x,A,y,lambda,sigma,n,p,mu){
  f=t(y-A%*%x)%*%(y-A%*%x)/(2*n*sigma)+lambda*sum(Si_v(x,p,mu));
  return(f)
}

smooth_fun_Scaled_fam<-function(x,A,Kin,y,lambda,sigma,n,p,mu){
  f=t(y-A%*%x)%*%Kin%*%(y-A%*%x)/(2*n*sigma)+lambda*sum(Si_v(x,p,mu));
  return(f)
}

smooth_fun_dev_Scaled<-function(x,A,y,lambda,sigma,n,p,mu){
  y=t(A)%*%(A%*%x-y)/(n*sigma)+lambda*sum(Si_dev_v(x,p,mu));
  return(y)
}

smooth_fun_dev_Scaled_fam<-function(x,A,Kin,y,lambda,sigma,n,p,mu){
  y=t(A)%*%Kin%*%(A%*%x-y)/(n*sigma)+lambda*sum(Si_dev_v(x,p,mu));
  return(y)
}

SoftThreshold <- function( x, lambda ) {
  #
  # Standard soft thresholding
  #
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i){
        v <- v+1;
      }
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval!=beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j]
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"%",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }                        
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

InverseLinfty_parallel <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl, cores = detectCores() - 1)
  data <- foreach (i = 1:p, .combine = rbind) %dopar% {
    InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
      p <- nrow(sigma);
      rho <- max(abs(sigma[i,-i])) / sigma[i,i];
      mu0 <- rho/(1+rho);
      beta <- rep(0,p);
      
      if (mu >= mu0){
        beta[i] <- (1-mu0)/sigma[i,i];
        returnlist <- list("optsol" = beta, "iter" = 0);
        return(returnlist);
      }
      
      diff.norm2 <- 1;
      last.norm2 <- 1;
      iter <- 1;
      iter.old <- 1;
      beta[i] <- (1-mu0)/sigma[i,i];
      beta.old <- beta;
      sigma.tilde <- sigma;
      diag(sigma.tilde) <- 0;
      vs <- -sigma.tilde%*%beta;
      
      while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
        
        for (j in 1:p){
          oldval <- beta[j];
          v <- vs[j];
          if (j==i){
            v <- v+1;
          }
          beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
          if (oldval!=beta[j]){
            vs <- vs + (oldval-beta[j])*sigma.tilde[,j]
          }
        }
        
        iter <- iter + 1;
        if (iter==2*iter.old){
          d <- beta - beta.old;
          diff.norm2 <- sqrt(sum(d*d));
          last.norm2 <-sqrt(sum(beta*beta));
          iter.old <- iter;
          beta.old <- beta;
          if (iter>10)
            vs <- -sigma.tilde%*%beta;
        }
      }
      
      returnlist <- list("optsol" = beta, "iter" = iter)
      return(returnlist)
    }
    SoftThreshold <- function( x, lambda ) {
      #
      # Standard soft thresholding
      #
      if (x>lambda){
        return (x-lambda);}
      else {
        if (x< (-lambda)){
          return (x+lambda);}
        else {
          return (0); }
      }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }                        
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  M=as.matrix(data)
  stopCluster(cl)
  return(M)
}

NoiseSd <- function( yh, A, n ){
  ynorm <- sqrt(n)*(yh/sqrt(diag(A)));
  sd.hat0 <- mad(ynorm);
  
  zeros <- (abs(ynorm)<3*sd.hat0);
  y2norm <- sum(yh[zeros]^2);
  Atrace <- sum(diag(A)[zeros]);
  sd.hat1 <- sqrt(n*y2norm/Atrace);
  
  ratio <- sd.hat0/sd.hat1;
  if (max(ratio,1/ratio)>2)
    print("Warning: Noise estimate problematic");
  
  s0 <- sum(zeros==FALSE);
  return (list( "sd" = sd.hat1, "nz" = s0));
}
