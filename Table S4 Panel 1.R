library("foreach")
library("doParallel")
library("expm")
Ker=function(s,h){
  temp=ifelse((s-h)>0,0,ifelse((s+h)<0,0,((1-(abs(s/h))**3)**3)*70/(h*81)))
  return(temp);
}
V=array(0,10);
B=array(0,10);
V0=175/247;
V[1]=420/187;
V[2]=7840/243;
B0=35/243;
B[1]=-35/81;
B[2]=140/81;

#d=1, q=1
select_h_one_dimen_one_poly=function(n,m,t,X,h_length,h_candidate,order){
  indexorder=sample(c(1:n),n)
  trainindex=indexorder[1:(n/2)];
  testindex=indexorder[(n/2+1):n]; 
  testerror=array(0,h_length);
  for(hit in 1:h_length){
    h=h_candidate[hit];
    for(i in testindex){
      for(j in 1:m){
        if(t[i,j]==(-100)){next;}
        testerror[hit]=testerror[hit]+(X[i,j]-one_dimen_one_poly(t[i,j],length(trainindex),m,t[trainindex,],X[trainindex,],h,0))**2;
      }
    }
  }
  factor=1;
  if(order>0){
    factor=(((2*order+1)*(order+2)*V[order]*B0*B0)/(2*V0*B[order]*B[order]))**(1/(2*order+5));
  }
  return(h_candidate[which.min(testerror)]*factor);
}
one_dimen_one_poly=function(tau,n,m,t,X,h,order){
  u0=sum(Ker(t-tau,h));
  u1=sum(Ker(t-tau,h)*(t-tau));
  u2=sum(Ker(t-tau,h)*(t-tau)*(t-tau));
  R0=sum(Ker(t-tau,h)*X);
  R1=sum(Ker(t-tau,h)*X*(t-tau));
  if(order==0){
    ini_order=(R0*u2-R1*u1)/(u0*u2-u1*u1)
    return(ini_order);
  }
  if(order==1){
    fir_order=(-R0*u1+R1*u0)/(u0*u2-u1*u1)
    return(fir_order);
  }
}

#d=1, q=2
select_h_one_dimen_sec_poly=function(n,m,t,X,h_length,h_candidate,order){
  indexorder=sample(c(1:n),n)
  trainindex=indexorder[1:(n/2)];
  testindex=indexorder[(n/2+1):n]; 
  testerror=array(0,h_length);
  for(hit in 1:h_length){
    h=h_candidate[hit];
    for(i in testindex){
      for(j in 1:m){
        if(t[i,j]==(-100)){next;}
        testerror[hit]=testerror[hit]+(X[i,j]-one_dimen_sec_poly(t[i,j],length(trainindex),m,t[trainindex,],X[trainindex,],h,0))**2;
      }
    }
  }
  factor=1;
  if(order>0){
    factor=(((2*order+1)*(order+2)*V[order]*B0*B0)/(2*V0*B[order]*B[order]))**(1/(2*order+5));
  }
  return(h_candidate[which.min(testerror)]*factor);
}
one_dimen_sec_poly=function(tau,n,m,t,X,h,order){
  u0=sum(Ker(t-tau,h));
  u1=sum(Ker(t-tau,h)*(t-tau));
  u2=sum(Ker(t-tau,h)*(t-tau)*(t-tau));
  u3=sum(Ker(t-tau,h)*(t-tau)*(t-tau)*(t-tau));
  u4=sum(Ker(t-tau,h)*(t-tau)*(t-tau)*(t-tau)*(t-tau));
  R0=sum(Ker(t-tau,h)*X);
  R1=sum(Ker(t-tau,h)*X*(t-tau));
  R2=sum(Ker(t-tau,h)*X*(t-tau)*(t-tau));
  if(order==0){
    ini_order=(R0*(u2*u4-u3*u3)-R1*(u1*u4-u2*u3)+R2*(u1*u3-u2*u2))/(u0*(u2*u4-u3*u3)-u1*(u1*u4-u2*u3)+u2*(u1*u3-u2*u2))
    return(ini_order);
  }
  if(order==1){
    fir_order=(-R0*(u1*u4-u2*u3)+R1*(u0*u4-u2*u2)-R2*(u0*u3-u1*u2))/(u0*(u2*u4-u3*u3)-u1*(u1*u4-u2*u3)+u2*(u1*u3-u2*u2))
    return(fir_order);
  }
  if(order==2){
    sec_order=2*(R0*(u1*u3-u2*u2)-R1*(u0*u3-u1*u2)+R2*(u0*u2-u1*u1))/(u0*(u2*u4-u3*u3)-u1*(u1*u4-u2*u3)+u2*(u1*u3-u2*u2))
    return(sec_order);
  }
}

#d=2, q=1
select_h_two_dimen_one_poly=function(n,m,t,X,h_length,h_candidate,order){
  indexorder=sample(c(1:n),n)
  trainindex=indexorder[1:(n/2)];
  testindex=indexorder[(n/2+1):n]; 
  testerror=array(0,h_length);
  for(hit in 1:h_length){
    h=h_candidate[hit];
    for(i in testindex){
      for(j1 in 1:m){
        for(j2 in 1:m){
          if(t[i,j1]==(-100)){next;}
          if(t[i,j2]==(-100)){next;}
          testerror[hit]=testerror[hit]+(X[i,j1]*X[i,j2]-two_dimen_one_poly(t[i,j1],t[i,j2],length(trainindex),m,t[trainindex,],X[trainindex,],h,0))**2;
        }}
    }
  }
  factor=1;
  if(order>0){
    factor=(((2*order+1)*(order+2)*V[order]*B0*B0)/(2*V0*B[order]*B[order]))**(1/(2*order+5));
  }
  return(h_candidate[which.min(testerror)]*factor);
}
two_dimen_one_poly=function(tau1,tau2,n,m,t,X,h,order){
  S00=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)));
  S10=sum(rowSums(Ker(t-tau1,h)*(t-tau1))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*(t-tau1)*Ker(t-tau2,h)));
  S01=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h)*(t-tau2))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)*(t-tau2)));
  S20=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*Ker(t-tau2,h)));
  S11=sum(rowSums(Ker(t-tau1,h)*(t-tau1))*rowSums(Ker(t-tau2,h)*(t-tau2))-rowSums(Ker(t-tau1,h)*(t-tau1)*Ker(t-tau2,h)*(t-tau2)));
  S02=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)*(t-tau2)*(t-tau2)));
  R00=sum(rowSums(Ker(t-tau1,h)*X)*rowSums(Ker(t-tau2,h)*X)-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)*X*X));
  R10=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*X)*rowSums(Ker(t-tau2,h)*X)-rowSums(Ker(t-tau1,h)*X*(t-tau1)*Ker(t-tau2,h)*X));
  R01=sum(rowSums(Ker(t-tau1,h)*X)*rowSums(Ker(t-tau2,h)*X*(t-tau2))-rowSums(Ker(t-tau1,h)*X*Ker(t-tau2,h)*X*(t-tau2)));
  if(order==0){
    ini_order=(R00*(S20*S02-S11*S11)-R10*(S10*S02-S01*S11)+R01*(S10*S11-S01*S20))/(S00*(S20*S02-S11*S11)-S10*(S10*S02-S01*S11)+S01*(S10*S11-S01*S20))
    return(ini_order);
  }
  if(order==1){
    firini_order=(-R00*(S10*S02-S11*S01)+R10*(S00*S02-S01*S01)-R01*(S00*S11-S10*S01))/(S00*(S20*S02-S11*S11)-S10*(S10*S02-S01*S11)+S01*(S10*S11-S01*S20))
    return(firini_order);
  }
}


parallelobtainhatbeta=function(it){
  write.csv(it,"process.csv");
  # generate data
  if(scorenumber==1){zeta=runif(n,1,2);}
  if(scorenumber==3){zeta=rnorm(n,1,0.1);for(i in 1:n){zeta[i]=max(zeta[i],0);}}
  if(scorenumber==5){zeta=rlnorm(n,0.5,0.1);}
  t=array(0,c(n,m));
  X=array(0,c(n,m));
  for(i in 1:n){
    t[i,]=runif(m,0,1);
    for(j in 1:m){
      X[i,j]=zeta[i]*exp(t[i,j])+rnorm(1,0,0.1);
    }
  }
  # varying m
  for(i in 1:n){
    if(varyingm==1){remainingm=floor(runif(1,0,1)*6)+4;}
    if(varyingm==2){remainingm=floor(runif(1,0,1)*6)+9;}
    if(varyingm==3){remainingm=rpois(1,4);
    remainingm=min(remainingm,15);
    remainingm=max(remainingm,3);}
    if(varyingm==4){remainingm=rpois(1,8);
    remainingm=min(remainingm,15);
    remainingm=max(remainingm,3);}
    for(j in c((remainingm+1):m)){
      t[i,j]=-100
    }
  }
  
#estimation procedure
recovery_procedure=function(n,m,testt,testX,valit,valiX){
  # estimate beta
  X=testX;
  t=testt;
  
  h_length=5;
  h_candidate=exp(seq(log(0.3),log(0.6),length.out=h_length));
  h_EX=select_h_one_dimen_one_poly(n,m,t,X,h_length,h_candidate,0);  
  h_EXp=select_h_one_dimen_one_poly(n,m,t,X,h_length,h_candidate,1);
  l=20;
  H=min(max(h_EXp,h_EX),0.45);
  betatanchor=seq(H,1-H,length.out=l);
  # EX
  EX=array(0,l);
  for(o in 1:l){
    EX[o]=one_dimen_one_poly(betatanchor[o],n,m,t,X,h_EX,0);
  }
  # EXp
  EXp=array(0,l);
  for(o in 1:l){
    EXp[o]=one_dimen_one_poly(betatanchor[o],n,m,t,X,h_EXp,1);
  }
  # hatbeta
  hatbeta=-sum(EX*EXp)/sum(EX*EX);
  
  #recovery, new discretion parameter l
  X=valiX;
  t=valit;
  l=20;
  tanchor=c(1:l)/l;
  centerX=array(0,c(n,m));
  fitEX=array(0,c(n,m));
  for(i in 1:n){
    for(j in 1:m){
      if(t[i,j]==(-100)){next;}
      fitEX[i,j]=one_dimen_one_poly(t[i,j],n,m,t,X,h_EX,0);
      centerX[i,j]=X[i,j]-fitEX[i,j];
    }
  }
  EX=array(0,l);
  for(o in 1:l){
    EX[o]=one_dimen_one_poly(tanchor[o],n,m,t,X,h_EX,0);
  }
  EXp=array(0,l);
  for(o in 1:l){
    EXp[o]=one_dimen_one_poly(tanchor[o],n,m,t,X,h_EXp,1);
  }
  C=array(0,c(l,l));
  h_C=select_h_two_dimen_one_poly(n,m,t,centerX,h_length,h_candidate,0);
  for(o1 in 1:l){
    for(o2 in 1:l){
      C[o1,o2]=two_dimen_one_poly(tanchor[o1],tanchor[o2],n,m,t,centerX,h_C,0);
    }
  }
  temp=eigen(C,symmetric=TRUE)
  eigenvalues=temp$values; #eigenvalues[1]/sum(eigenvalues[which(eigenvalues>0)])>0.85
  eigenvectors=temp$vectors; #eigenvectors[,k] is the kth eigenvector 
  if(eigenvalues[1]/sum(eigenvalues[which(eigenvalues>0)])<0.85){
    return(-1);
  }
  
  Cp=array(0,c(l,l));
  h_Cp=select_h_two_dimen_one_poly(n,m,t,centerX,h_length,h_candidate,1);
  for(o1 in 1:l){
    for(o2 in 1:l){
      Cp[o1,o2]=two_dimen_one_poly(tanchor[o1],tanchor[o2],n,m,t,centerX,h_Cp,1);
    }
  }
  eigenvectorsp=array(0,c(l,10)); #the derivative of eigenvector
  for(k in 1:10){
    eigenvectorsp[,k]=(Cp%*%eigenvectors[,k])/eigenvalues[k];
  }
  
  # calculate MSE for different kappa
  kappalength=5;
  kappalist=exp(seq(log(0.001),log(0.1),length.out=kappalength));
  GCV=array(0,kappalength);
  for(kappait in 1:kappalength){
    kappa=kappalist[kappait];
    #estimate xi_{i,1}
    hatxi1=array(0,n);
    for(i in 1:n){
      numeritor=0;
      dominant=0;
      for(j in 1:m){
        if(t[i,j]==(-100)){next;}
        numeritor=numeritor+(fitEX[i,j]-X[i,j])*eigenvectors[ceiling(t[i,j]*l),1];
        dominant=dominant+eigenvectors[ceiling(t[i,j]*l),1]*eigenvectors[ceiling(t[i,j]*l),1];
      }
      numeritor=numeritor*l/(m*kappa);
      dominant=dominant*l/(m*kappa);
      for(o in 1:l){
        numeritor=numeritor+(EXp[o]+hatbeta*EX[o])*(eigenvectorsp[o,1]+hatbeta*eigenvectors[o,1]);
        dominant=dominant+(eigenvectorsp[o,1]+hatbeta*eigenvectors[o,1])*(eigenvectorsp[o,1]+hatbeta*eigenvectors[o,1]);
      }
      hatxi1[i]=-numeritor/dominant;
    }
    #fit X on t
    fitX=array(0,c(n,m));
    for(i in 1:n){
      for(j in 1:m){
        if(t[i,j]==(-100)){next;}
        fitX[i,j]=fitEX[i,j]+hatxi1[i]*eigenvectors[ceiling(t[i,j]*l),1];
      }
    }
    GCV[kappait]=sum((fitX[,]-valiX[,])**2)/(n*m);
  }#end of kappait
  return(kappalist[which.min(GCV)]);
}
#select kappa
kappa=recovery_procedure(n/2,m,t[1:(n/2),],X[1:(n/2),],t[(n/2+1):n,],X[(n/2+1):n,])
if(kappa==-1){return(-1);}

#recovery and mse
# estimate beta
h_length=5;
h_candidate=exp(seq(log(0.3),log(0.6),length.out=h_length));
h_EX=select_h_one_dimen_one_poly(n,m,t,X,h_length,h_candidate,0);  
h_EXp=select_h_one_dimen_one_poly(n,m,t,X,h_length,h_candidate,1);
l=20;
H=min(max(h_EXp,h_EX),0.45);
betatanchor=seq(H,1-H,length.out=l);
# EX
EX=array(0,l);
for(o in 1:l){
  EX[o]=one_dimen_one_poly(betatanchor[o],n,m,t,X,h_EX,0);
}
# EXp
EXp=array(0,l);
for(o in 1:l){
  EXp[o]=one_dimen_one_poly(betatanchor[o],n,m,t,X,h_EXp,1);
}
# hatbeta
hatbeta=-sum(EX*EXp)/sum(EX*EX);

#recovery, new discretion parameter l
l=20;
tanchor=c(1:l)/l;
realX=array(0,c(n,l));
for(i in 1:n){
  for(o in 1:l){
    realX[i,o]=zeta[i]*exp(tanchor[o]);
  }
}
centerX=array(0,c(n,m));
fitEX=array(0,c(n,m));
for(i in 1:n){
  for(j in 1:m){
    if(t[i,j]==(-100)){next;}
    fitEX[i,j]=one_dimen_one_poly(t[i,j],n,m,t,X,h_EX,0);
    centerX[i,j]=X[i,j]-fitEX[i,j];
  }
}
EX=array(0,l);
for(o in 1:l){
  EX[o]=one_dimen_one_poly(tanchor[o],n,m,t,X,h_EX,0);
}
EXp=array(0,l);
for(o in 1:l){
  EXp[o]=one_dimen_one_poly(tanchor[o],n,m,t,X,h_EXp,1);
}
C=array(0,c(l,l));
h_C=select_h_two_dimen_one_poly(n,m,t,centerX,h_length,h_candidate,0);
for(o1 in 1:l){
  for(o2 in 1:l){
    C[o1,o2]=two_dimen_one_poly(tanchor[o1],tanchor[o2],n,m,t,centerX,h_C,0);
  }
}
temp=eigen(C,symmetric=TRUE)
eigenvalues=temp$values; #eigenvalues[1]/sum(eigenvalues[which(eigenvalues>0)])>0.85
eigenvectors=temp$vectors; #eigenvectors[,k] is the kth eigenvector 
if(eigenvalues[1]/sum(eigenvalues[which(eigenvalues>0)])<0.85){
  return(-1);
}

Cp=array(0,c(l,l));
h_Cp=select_h_two_dimen_one_poly(n,m,t,centerX,h_length,h_candidate,1);
for(o1 in 1:l){
  for(o2 in 1:l){
    Cp[o1,o2]=two_dimen_one_poly(tanchor[o1],tanchor[o2],n,m,t,centerX,h_Cp,1);
  }
}
eigenvectorsp=array(0,c(l,10)); #the derivative of eigenvector
for(k in 1:10){
  eigenvectorsp[,k]=(Cp%*%eigenvectors[,k])/eigenvalues[k];
}

#estimate xi_{i,1}
hatxi1=array(0,n);
for(i in 1:n){
  numeritor=0;
  dominant=0;
  for(j in 1:m){
    if(t[i,j]==(-100)){next;}
    numeritor=numeritor+(fitEX[i,j]-X[i,j])*eigenvectors[ceiling(t[i,j]*l),1];
    dominant=dominant+eigenvectors[ceiling(t[i,j]*l),1]*eigenvectors[ceiling(t[i,j]*l),1];
  }
  numeritor=numeritor*l/(m*kappa);
  dominant=dominant*l/(m*kappa);
  for(o in 1:l){
    numeritor=numeritor+(EXp[o]+hatbeta*EX[o])*(eigenvectorsp[o,1]+hatbeta*eigenvectors[o,1]);
    dominant=dominant+(eigenvectorsp[o,1]+hatbeta*eigenvectors[o,1])*(eigenvectorsp[o,1]+hatbeta*eigenvectors[o,1]);
  }
  hatxi1[i]=-numeritor/dominant;
}
#recover X on tanchor
recoveryX=array(0,c(n,l));
for(i in 1:n){
  for(o in 1:l){
    recoveryX[i,o]=EX[o]+hatxi1[i]*eigenvectors[o,1];
  }
}
return(sum((recoveryX[,]-realX[,])**2)/(n*l));
}



#sample size
for(scorenumber in c(1,3,5)){
  for(n in c(100,200)){
    m=16;
    for(varyingm in c(1:4)){
      
  set.seed(97);
  itermax=100;
  #parallel
  closeAllConnections();
  closeAllConnections();
  cl <- makeCluster(100);
  registerDoParallel(cl);
  result <- foreach(it=c(1:itermax), .combine='c') %dopar% parallelobtainhatbeta(it);
  stopCluster(cl);
  closeAllConnections();
  closeAllConnections();
  
  result_recoveryMSE=result[which(result>0)];
  cat(scorenumber,n,varyingm,mean(result_recoveryMSE),sd(result_recoveryMSE)/sqrt(length(result_recoveryMSE)),"\n");
}}}

