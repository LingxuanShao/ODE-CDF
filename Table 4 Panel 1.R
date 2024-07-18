library("readxl");
library("foreach")
library("doParallel")
require("ggplot2")
library("expm")
set.seed(1);

data <- read_excel("/home/shaolx/Documents/Project 4 ODE-CDF/data.xlsx")
# 4 for cured; 6 for current confirmed; 9 for dead; 14 for name
n=33;
m=40;
#wash data
index=which(as.matrix(data[,8])==20200210);
name=as.matrix(data[index,14]);
nohubeiindex=index[-18];
nohubeiname=as.matrix(data[nohubeiindex,14]);
Xcurrent=array(0,c(n,m));
Xcured=array(0,c(n,m));
Xdead=array(0,c(n,m));
Xcheck=array(0,c(n,m));
for(i in 1:n){
  Xcurrent[i,]=as.matrix(data[c(nohubeiindex[i]:(nohubeiindex[i]+m-1)),6]);
  Xcured[i,]=as.matrix(data[c(nohubeiindex[i]:(nohubeiindex[i]+m-1)),4]);
  Xdead[i,]=as.matrix(data[c(nohubeiindex[i]:(nohubeiindex[i]+m-1)),9]);
  Xcheck[i,]=as.matrix(data[c(nohubeiindex[i]:(nohubeiindex[i]+m-1)),14]);
}
t=array(0,c(n,m));
for(i in 1:n){
  for(j in 1:m){
    t[i,j]=j/m;
  }
}

#PDE model 
# Xcurrent=I; Xcured=R; Xdead=D;
#need to estimate I, I', R', D'
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


select_h_one_dimen_one_poly=function(n,m,t,X,h_length,h_candidate,order){
  set.seed(1);
  indexorder=sample(c(1:n),n)
  trainindex=indexorder[1:(n/2)];
  testindex=indexorder[(n/2+1):n]; 
  testerror=array(0,h_length);
  for(hit in 1:h_length){
    h=h_candidate[hit];
    for(i in testindex){
      for(j in 1:m){
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

select_h_one_dimen_sec_poly=function(n,m,t,X,h_length,h_candidate,order){
  set.seed(1);
  indexorder=sample(c(1:n),n)
  trainindex=indexorder[1:(n/2)];
  testindex=indexorder[(n/2+1):n]; 
  testerror=array(0,h_length);
  for(hit in 1:h_length){
    h=h_candidate[hit];
    for(i in testindex){
      for(j in 1:m){
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

##################################################
#no sparsity case
h_length=8;
h_candidate=exp(seq(log(0.3),log(0.6),length.out=h_length))
h_EXcurrent=select_h_one_dimen_one_poly(n,m,t,Xcurrent,h_length,h_candidate,0);
h_EXcurrentp=select_h_one_dimen_one_poly(n,m,t,Xcurrent,h_length,h_candidate,1);
h_EXcuredp=select_h_one_dimen_one_poly(n,m,t,Xcured,h_length,h_candidate,1);
h_EXdeadp=select_h_one_dimen_one_poly(n,m,t,Xdead,h_length,h_candidate,1);
l=20;
H=min(max(h_EXcurrent,h_EXcurrentp,h_EXcuredp,h_EXdeadp),0.45);
betatanchor=seq(H,1-H,length.out=l);

EXcurrent=array(0,l);
EXcurrentp=array(0,l);
for(o in 1:l){
  EXcurrent[o]=one_dimen_one_poly(betatanchor[o],n,m,t,Xcurrent,h_EXcurrent,0);
  EXcurrentp[o]=one_dimen_one_poly(betatanchor[o],n,m,t,Xcurrent,h_EXcurrentp,1);
}
EXcuredp=array(0,l);
for(o in 1:l){
  EXcuredp[o]=one_dimen_one_poly(betatanchor[o],n,m,t,Xcured,h_EXcuredp,1);
}
EXdeadp=array(0,l);
for(o in 1:l){
  EXdeadp[o]=one_dimen_one_poly(betatanchor[o],n,m,t,Xdead,h_EXdeadp,1);
}

hatgamma=sum(EXcuredp*EXcurrent)/sum(EXcurrent*EXcurrent)
hatrho=sum(EXdeadp*EXcurrent)/sum(EXcurrent*EXcurrent)
hatbeta=sum(EXcurrentp*EXcurrent)/sum(EXcurrent*EXcurrent)#+hatgamma+hatrho
cat(c("original data",round(hatbeta,4),round(hatgamma,4),round(hatrho,4)),"\n")


############################################################
#sparse case
parallelobtainhatbeta=function(it){
  index_sparse=array(0,c(n,m_sparse))
  for(i in 1:n){
    temp=sort(sample((1:m),m_sparse));
    for(j in 1:m_sparse){
      index_sparse[i,j]=temp[j];
    }
  }
  t_sparse=array(0,c(n,m_sparse));
  Xcurrent_sparse=array(0,c(n,m_sparse));
  Xcured_sparse=array(0,c(n,m_sparse));
  Xdead_sparse=array(0,c(n,m_sparse));
  for(i in 1:n){
    for(j in 1:m_sparse){
      t_sparse[i,j]=t[i,index_sparse[i,j]];
      Xcurrent_sparse[i,j]=Xcurrent[i,index_sparse[i,j]];
      Xcured_sparse[i,j]=Xcured[i,index_sparse[i,j]];
      Xdead_sparse[i,j]=Xdead[i,index_sparse[i,j]];
  }}
  
  h_EXcurrent_sparse=select_h_one_dimen_one_poly(n,m_sparse,t_sparse,Xcurrent_sparse,h_length,h_candidate,0);
  h_EXcurrentp_sparse=select_h_one_dimen_one_poly(n,m_sparse,t_sparse,Xcurrent_sparse,h_length,h_candidate,1);
  h_EXcuredp_sparse=select_h_one_dimen_one_poly(n,m_sparse,t_sparse,Xcured_sparse,h_length,h_candidate,1);
  h_EXdeadp_sparse=select_h_one_dimen_one_poly(n,m_sparse,t_sparse,Xdead_sparse,h_length,h_candidate,1);
  l=20;
  H=min(max(h_EXcurrent_sparse,h_EXcurrentp_sparse,h_EXcuredp_sparse,h_EXdeadp_sparse),0.45);
  betatanchor=seq(H,1-H,length.out=l);
  
  EXcurrent_sparse=array(0,l);
  EXcurrentp_sparse=array(0,l);
  for(o in 1:l){
    EXcurrent_sparse[o]=one_dimen_one_poly(betatanchor[o],n,m_sparse,t_sparse,Xcurrent_sparse,h_EXcurrent_sparse,0);
    EXcurrentp_sparse[o]=one_dimen_one_poly(betatanchor[o],n,m_sparse,t_sparse,Xcurrent_sparse,h_EXcurrentp_sparse,1);
  }
  
  EXcuredp_sparse=array(0,l);
  for(o in 1:l){
    EXcuredp_sparse[o]=one_dimen_one_poly(betatanchor[o],n,m_sparse,t_sparse,Xcured_sparse,h_EXcuredp_sparse,1);
  }
  
  EXdeadp_sparse=array(0,l);
  for(o in 1:l){
    EXdeadp_sparse[o]=one_dimen_one_poly(betatanchor[o],n,m_sparse,t_sparse,Xdead_sparse,h_EXdeadp_sparse,1);
  }
  
  hatgamma_sparse=sum(EXcuredp_sparse*EXcurrent_sparse)/sum(EXcurrent_sparse*EXcurrent_sparse)
  hatrho_sparse=sum(EXdeadp_sparse*EXcurrent_sparse)/sum(EXcurrent_sparse*EXcurrent_sparse)
  hatbeta_sparse=sum(EXcurrentp_sparse*EXcurrent_sparse)/sum(EXcurrent_sparse*EXcurrent_sparse)#+hatgamma_sparse+hatrho_sparse
  return(as.list(c(hatbeta_sparse,hatgamma_sparse,hatrho_sparse)));
}


relaerror=array(0,c(m,3));
estimates_sparse=array(0,c(m,3));
for(m_sparse in seq(15,37,1)){
  itermax=100;
  hatbeta_sparse_result=array(0,itermax);
  hatgamma_sparse_result=array(0,itermax);
  hatrho_sparse_result=array(0,itermax);
  #parallel
  cl <- makeCluster(100);
  registerDoParallel(cl);
  result <- foreach(it=c(1:itermax), .combine='c') %dopar% parallelobtainhatbeta(it);
  stopCluster(cl);
  for(i in 1:itermax){
    hatbeta_sparse_result[i]=result[[3*i-2]];
    hatgamma_sparse_result[i]=result[[3*i-1]];
    hatrho_sparse_result[i]=result[[3*i]];
  }
  relaerror[m_sparse,1]=round( mean(abs(hatbeta_sparse_result-hatbeta)), 4 )
  relaerror[m_sparse,2]=round( mean(abs(hatgamma_sparse_result-hatgamma)), 4 )  
  relaerror[m_sparse,3]=round( mean(abs(hatrho_sparse_result-hatrho)), 4 )  
  estimates_sparse[m_sparse,1]=round(mean(hatbeta_sparse_result),4)
  estimates_sparse[m_sparse,2]=round(mean(hatgamma_sparse_result),4)
  estimates_sparse[m_sparse,3]=round(mean(hatrho_sparse_result),4)
  cat(c("m'=",m_sparse,
        estimates_sparse[m_sparse,1],
        relaerror[m_sparse,1],
        estimates_sparse[m_sparse,2],
        relaerror[m_sparse,2],
        estimates_sparse[m_sparse,3],
        relaerror[m_sparse,3],
        "\n"))
}

cat(c("original data: ",round(hatbeta,4),round(hatgamma,4),round(hatrho,4)),"\n")


# plot
library("ggplot2")
data=as.data.frame(cbind(relaerror[15:39,],c(15:39)))
ggplot(data)+geom_line(aes(x=V4,y=V1),color="black",size=0.5)+
  geom_line(aes(x=V4,y=V2),color="blue",size=0.5)+
  geom_line(aes(x=V4,y=V3),color="red",size=0.5)+ 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())




