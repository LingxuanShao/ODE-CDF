library("fdapace")
library("expm")
set.seed(97);

#sample size
for(n in c(200)){
  for(scorenumber in c(1,3,5)){
    for(m in c(3,5,7,9,11)){
      set.seed(97);
      itermax=100;
      result_PACEMSE=array(0,itermax);
      it=1;
      while(it<(itermax+1)){
        #generate data
        if(scorenumber==1){zeta=runif(n,1,2);}
        if(scorenumber==3){zeta=rnorm(n,1,0.1);for(i in 1:n){zeta[i]=max(zeta[i],0);}}
        if(scorenumber==5){zeta=rlnorm(n,0.5,0.1);}
        t=array(0,c(n,m));
        X=array(0,c(n,m));
        for(i in 1:n){
          t[i,]=sort(runif(m,0,1));
          for(j in 1:m){
            X[i,j]=t[i,j]*t[i,j]+2*t[i,j]*sqrt(zeta[i]+1)+zeta[i]+3+rnorm(1,0,0.1);
          }
        }
        
        ###PACE method
        err=try(FPCA(Ly=split(X, row(X)),Lt=split(t, row(t)),list(dataType='Sparse')),silent=0);
        if("try-error" %in% class(err)){}
        else{
          temp=FPCA(Ly=split(X, row(X)),Lt=split(t, row(t)),list(dataType='Sparse'));
          temp2=fitted(object=temp);
          PACEMSE=array(0,n);
          for(i in 1:n){
            PACEMSE[i]=sum((temp2[i,]- (temp$workGrid*temp$workGrid+2*temp$workGrid*sqrt(zeta[i]+1)+zeta[i]+3) )**2)/length(temp$workGrid);
          }
          result_PACEMSE[it]=mean(PACEMSE);
          it=it+1;
        }
      }
      result_recoveryMSE=result_PACEMSE[which(result_PACEMSE>0)];
      cat("(2)sf",n,scorenumber,m,round(mean(result_recoveryMSE),4),round(sd(result_recoveryMSE)/sqrt(length(result_recoveryMSE)),4),"\n");
}}}#end for n,m,scorenumber 
