library("foreach")
library("doParallel")
library("ggplot2")
library("expm")
library("fda")
set.seed(97);

obtainbeta=function(X,t,l,tanchor,K,allPhi,allPhip,Phi,lambda){
  minbeta=-2;#becareful with this setting
  maxbeta=0;
  while((maxbeta-minbeta)>0.001){
    temp=searchbeta(X,t,l,tanchor,K,allPhi,allPhip,Phi,lambda,minbeta,maxbeta);
    minbeta=temp[[1]];
    maxbeta=temp[[2]];
  }
  return((maxbeta+minbeta)/2);
};
searchbeta=function(X,t,l,tanchor,K,allPhi,allPhip,Phi,lambda,minbeta,maxbeta){
  candidatebetalength=5;
  candidatebeta=seq(minbeta,maxbeta,length.out=candidatebetalength);
  Lerror=array(0,candidatebetalength);
  for(betait in 1:candidatebetalength){
    beta=candidatebeta[betait];
    Lindierror=array(0,n);
    for(i in 1:n){
      Lindierror[i]=sum((X[i,]-Phi[i,,]%*%  solve(t(Phi[i,,])%*%Phi[i,,]+lambda*colSums(allPhip[i,,]+beta*allPhi[i,,])%*%t(colSums(allPhip[i,,]+beta*allPhi[i,,]))/(l*l))  %*%t(Phi[i,,])%*%X[i,])**2);
    }
    Lerror[betait]=sum(Lindierror);
  }
  inde=which.min(Lerror);
  if(inde==1){inde=inde+1;}
  if(inde==candidatebetalength){inde=inde-1;}
  return(list(candidatebeta[inde-1],candidatebeta[inde+1]));
}

#sample size
for(n in c(200)){
for(scorenumber in c(1,3,5)){
for(m in c(3,5,7,9,11)){
      
#parallel
parallelobtainhatbeta=function(it){
library("fda")
#generate data
  if(scorenumber==1){zeta=runif(n,1,2);}
  if(scorenumber==3){zeta=rnorm(n,1,0.1);for(i in 1:n){zeta[i]=max(zeta[i],0);}}
  if(scorenumber==5){zeta=rlnorm(n,0.5,0.1);}
t=array(0,c(n,m));
X=array(0,c(n,m));
for(i in 1:n){
  t[i,]=sort(runif(m,0,1));
  for(j in 1:m){
    X[i,j]=zeta[i]*exp(t[i,j])+rnorm(1,0,0.1);
  }
}
#B spline Basis
l=20;
tanchor=c(1:l)/l;
norder=2;
K=norder+1;if(m<=5){K=norder;}
allPhi=array(0,c(n,l,K));
allPhip=array(0,c(n,l,K));
for(i in 1:n){
  allPhi[i,,]=bsplineS(x=tanchor,breaks=seq(0,1,length.out=K-norder+2),norder=norder );
  allPhip[i,,]=bsplineS(x=tanchor,breaks=seq(0,1,length.out=K-norder+2),norder=norder,nderiv=1);
}
Phi=array(0,c(n,m,K));
for(i in 1:n){for(j in 1:m){for(k in 1:K){
    Phi[i,j,k]=allPhi[i,ceiling(t[i,j]*l),k];
}}}
# estimation
candidatelambda=c(0.1,1,10,100,1000);
candidatelambdalength=length(candidatelambda);
GCV=array(0,candidatelambdalength);
hatbetalist=array(0,candidatelambdalength);
for(lambdait in 1:candidatelambdalength){
  lambda=candidatelambda[lambdait];
  err=try(obtainbeta(X,t,l,tanchor,K,allPhi,allPhip,Phi,lambda));
  if("try-error" %in% class(err)){  
    hatbetalist[lambdait]=0;
    GCV[lambdait]=10**8;
  }
  else{
    hatbetalist[lambdait]=obtainbeta(X,t,l,tanchor,K,allPhi,allPhip,Phi,lambda);
    #GCV
    Lindierror=array(0,n);
    traceindi=array(0,n);
    for(i in 1:n){
      Lindierror[i]=sum((X[i,]-Phi[i,,]%*% solve(t(Phi[i,,])%*%Phi[i,,]+lambda*colSums(allPhip[i,,]+hatbetalist[lambdait]*allPhi[i,,])%*%t(colSums(allPhip[i,,]+hatbetalist[lambdait]*allPhi[i,,]))/(l*l)) %*%t(Phi[i,,])%*%X[i,])**2);
      traceindi[i]=sum(diag( Phi[i,,]%*% solve(t(Phi[i,,])%*%Phi[i,,]+lambda*colSums(allPhip[i,,]+hatbetalist[lambdait]*allPhi[i,,])%*%t(colSums(allPhip[i,,]+hatbetalist[lambdait]*allPhi[i,,]))/(l*l)) %*%t(Phi[i,,]) ));
    }
    GCV[lambdait]=sum(Lindierror)/((n*m-sum(traceindi))**2);
  }
}
hatbeta=hatbetalist[which.min(GCV)];
#MSE
profileMSE=array(0,n);
for(i in 1:n){
  profileMSE[i]=sum( ( zeta[i]*exp(tanchor) - allPhi[i,,]%*% solve(t(Phi[i,,])%*%Phi[i,,]+lambda*colSums(allPhip[i,,]+hatbeta*allPhi[i,,])%*%t(colSums(allPhip[i,,]+hatbeta*allPhi[i,,]))/(l*l))%*%t(Phi[i,,])%*%X[i,] )**2 )/l;
}
return(mean(profileMSE));
}

itermax=100;
closeAllConnections();closeAllConnections();
cl <- makeCluster(100);
registerDoParallel(cl);
temp <- foreach(it=c(1:itermax), .combine='c') %dopar% parallelobtainhatbeta(it);
stopCluster(cl);
      
result_recoveryMSE=temp[which(temp>0)];
cat(n,scorenumber,m,round(mean(result_recoveryMSE),4),round(sd(result_recoveryMSE)/sqrt(length(result_recoveryMSE)),4),"\n");
}}}#end for n,m, scorenumber


