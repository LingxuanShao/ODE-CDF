library("foreach")
library("doParallel")
require("ggplot2")
library("expm")
set.seed(97);

Ker=function(s,h){
  temp=ifelse((s-h)>0,0,ifelse((s+h)<0,0,((1-(abs(s/h))**3)**3)*70/(h*81)))
  return(temp);
}
#V,B for bandwidth selection
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
  if(order==1){#derivative of the first component
    firini_order=(-R00*(S10*S02-S11*S01)+R10*(S00*S02-S01*S01)-R01*(S00*S11-S10*S01))/(S00*(S20*S02-S11*S11)-S10*(S10*S02-S01*S11)+S01*(S10*S11-S01*S20))
    return(firini_order);
  }
}

#d=2, q=2
select_h_two_dimen_sec_poly=function(n,m,t,X,h_length,h_candidate,order){#target on the diagonal (tau,tau)
  indexorder=sample(c(1:n),n)
  trainindex=indexorder[1:(n/2)];
  testindex=indexorder[(n/2+1):n]; 
  testerror=array(0,h_length);
  for(hit in 1:h_length){
    h=h_candidate[hit];
    for(i in testindex){
      for(j in 1:m){
        if(t[i,j]==(-100)){next;}
        testerror[hit]=testerror[hit]+(X[i,j]*X[i,j]-two_dimen_sec_poly(t[i,j],length(trainindex),m,t[trainindex,],X[trainindex,],h,0))**2;
      }
    }
  }
  factor=1;
  if(order>0){
    factor=(((2*order+1)*(order+2)*V[order]*B0*B0)/(2*V0*B[order]*B[order]))**(1/(2*order+5));
  }
  return(h_candidate[which.min(testerror)]*factor);
}
two_dimen_sec_poly=function(tau,n,m,t,X,h,order){#target on the diagonal (tau,tau)
  tau1=tau;
  tau2=tau;
  S00=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)));
  S10=sum(rowSums(Ker(t-tau1,h)*(t-tau1))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*(t-tau1)*Ker(t-tau2,h)));
  S01=S10;
  #S01=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h)*(t-tau2))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)*(t-tau2)));
  S20=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*Ker(t-tau2,h)));
  S11=sum(rowSums(Ker(t-tau1,h)*(t-tau1))*rowSums(Ker(t-tau2,h)*(t-tau2))-rowSums(Ker(t-tau1,h)*(t-tau1)*Ker(t-tau2,h)*(t-tau2)));
  S02=S20;
  #S02=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)*(t-tau2)*(t-tau2)));
  S30=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*(t-tau1)*Ker(t-tau2,h)));
  S21=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h)*(t-tau2))-rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*Ker(t-tau2,h)*(t-tau2)));
  S12=S21;
  S03=S30;
  #S12=sum(rowSums(Ker(t-tau1,h)*(t-tau1))*rowSums(Ker(t-tau2,h)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*(t-tau1)*Ker(t-tau2,h)*(t-tau2)*(t-tau2)));
  #S03=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h)*(t-tau2)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)*(t-tau2)*(t-tau2)*(t-tau2)));
  S40=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h))-rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*(t-tau1)*(t-tau1)*Ker(t-tau2,h)));
  S31=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h)*(t-tau2))-rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*(t-tau1)*Ker(t-tau2,h)*(t-tau2)));
  S22=sum(rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*(t-tau1)*(t-tau1)*Ker(t-tau2,h)*(t-tau2)*(t-tau2)));
  S13=S31;
  S04=S40;
  #S13=sum(rowSums(Ker(t-tau1,h)*(t-tau1))*rowSums(Ker(t-tau2,h)*(t-tau2)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*(t-tau1)*Ker(t-tau2,h)*(t-tau2)*(t-tau2)*(t-tau2)));
  #S04=sum(rowSums(Ker(t-tau1,h))*rowSums(Ker(t-tau2,h)*(t-tau2)*(t-tau2)*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*Ker(t-tau2,h)*(t-tau2)*(t-tau2)*(t-tau2)*(t-tau2)));
  R00=sum(rowSums(Ker(t-tau1,h)*X)*rowSums(Ker(t-tau2,h)*X)-rowSums(Ker(t-tau1,h)*X*Ker(t-tau2,h)*X));
  R10=sum(rowSums(Ker(t-tau1,h)*X*(t-tau1))*rowSums(Ker(t-tau2,h)*X)-rowSums(Ker(t-tau1,h)*X*(t-tau1)*Ker(t-tau2,h)*X));
  R01=R10;
  #R01=sum(rowSums(Ker(t-tau1,h)*X)*rowSums(Ker(t-tau2,h)*X*(t-tau2))-rowSums(Ker(t-tau1,h)*X*Ker(t-tau2,h)*X*(t-tau2)));
  R20=sum(rowSums(Ker(t-tau1,h)*X*(t-tau1)*(t-tau1))*rowSums(Ker(t-tau2,h)*X)-rowSums(Ker(t-tau1,h)*X*(t-tau1)*(t-tau1)*Ker(t-tau2,h)*X));
  R11=sum(rowSums(Ker(t-tau1,h)*X*(t-tau1))*rowSums(Ker(t-tau2,h)*X*(t-tau2))-rowSums(Ker(t-tau1,h)*X*(t-tau1)*Ker(t-tau2,h)*X*(t-tau2)));
  R02=R20;
  #R02=sum(rowSums(Ker(t-tau1,h)*X)*rowSums(Ker(t-tau2,h)*X*(t-tau2)*(t-tau2))-rowSums(Ker(t-tau1,h)*X*Ker(t-tau2,h)*X*(t-tau2)*(t-tau2)));
  hatH=array(0,c(4,4));
  hatV=array(0,4);
  hatH[1,1]=S00;
  hatH[1,2]=S10+S01;
  hatH[1,3]=S20+S02;
  hatH[1,4]=S11;
  hatH[2,1]=S10+S01;
  hatH[2,2]=S20+S02+2*S11;
  hatH[2,3]=S30+S21+S12+S03;
  hatH[2,4]=S21+S12;
  hatH[3,1]=S20+S02;
  hatH[3,2]=S30+S21+S12+S03;
  hatH[3,3]=S40+S04+2*S22;
  hatH[3,4]=S31+S13;
  hatH[4,1]=S11;
  hatH[4,2]=S21+S12;
  hatH[4,3]=S31+S13;
  hatH[4,4]=S22;
  hatV[1]=R00;
  hatV[2]=R10+R01;
  hatV[3]=R20+R02;
  hatV[4]=R11;
  hattheta=solve(hatH)%*%hatV;
  if(order==0){return(hattheta[1]);}
  if(order==1){return(hattheta[2]);}#hattheta[2] for \partial_{s_{1}}C=\partial_{s_{2}}C
  if(order==2){return(hattheta[4]);}#hattheta[4] for \partial_{s_{1}s_{2}}C 
  if(order==3){return(2*hattheta[3]);}#hattheta[3] for \partial_{s_{1}s_{1}}C/2=\partial_{s_{2}s_{2}}C/2
}

#sample size
for(scorenumber in c(1,3,5)){
  for(n in c(100,200)){
    m=16;
    for(varyingm in c(1:4)){
      
itermax=100;
#parallel
parallelobtainhatbeta=function(it){
  
#generate data
if(scorenumber==1){zeta=runif(n,1,2);}
if(scorenumber==3){zeta=rnorm(n,1,0.1);for(i in 1:n){zeta[i]=max(zeta[i],0);}}
if(scorenumber==5){zeta=rlnorm(n,0.5,0.1);}
t=array(0,c(n,m));
X=array(0,c(n,m));
for(i in 1:n){
  for(j in 1:m){
    t[i,j]=runif(1,0,1);
    X[i,j]=t[i,j]*t[i,j]+2*t[i,j]*sqrt(zeta[i]+1)+zeta[i]+3+rnorm(1,0,0.1);
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

#estimation
#h selection
h_length=5;
h_candidate1=exp(seq(log(0.4),log(0.6),length.out=h_length));
h_EXpXp=select_h_two_dimen_sec_poly(n,m,t,X,h_length,h_candidate1,2);
h_candidate2=exp(seq(log(0.25),log(0.35),length.out=h_length));
h_EX=select_h_one_dimen_one_poly(n,m,t,X,h_length,h_candidate2,0);
#tanchor
l=10;
H=min(max(h_EXpXp,h_EX),0.45);
tanchor=seq(H,1-H,length.out=l);
#
EXpXp=array(0,l);
for(o in 1:l){
  EXpXp[o]=two_dimen_sec_poly(tanchor[o],n,m,t,X,h_EXpXp,2);
}
#
EX=array(0,l);
for(o in 1:l){
  EX[o]=one_dimen_one_poly(tanchor[o],n,m,t,X,h_EX,0);
}
#compute hatbeta
hatbeta=-sum((EXpXp+8)*EX)/sum(EX*EX);
#relative error
beta=-4;
relaerror=((hatbeta-beta)*(hatbeta-beta))/(beta*beta);
write.csv(it,"process.csv");
return(relaerror);
}
closeAllConnections();closeAllConnections();
cl <- makeCluster(100);
registerDoParallel(cl);
temp <- foreach(it=c(1:itermax), .combine='c') %dopar% parallelobtainhatbeta(it);
stopCluster(cl);
closeAllConnections();closeAllConnections();

cat(scorenumber, n, varyingm, mean(temp), sd(temp)/sqrt(itermax), "\n");

}}}#end for n,m, scorenumber


