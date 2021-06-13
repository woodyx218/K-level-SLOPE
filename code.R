
library(SLOPE)
library(pixmap)
library(mrbsizeR)
library(bmp)
library(magic)

######### Algorithm 1 & 2 ( function averaging() and PGD_new() respectively )

library(SLOPE)
### state evolution equation
F= function(tau,alpha,prior=x0,iteration=100,sigma=0){
  result=0
  for (i in 1:iteration){
    result=result+mean((prox_sorted_L1(prior+tau*rnorm(p),alpha*tau)-prior)^2)/iteration
  }
  return(sigma^2+result/delta)
}


### alpha to tau calibration
alpha_to_tau=function(alpha_seq,max_iter=100,prior=x0,sigma=0){
  tau=sqrt(sigma^2+mean(prior^2)/delta)
  record_tau=rep(0,max_iter) # initialize
  
  for (t in 1:max_iter){
    tau=sqrt(F(tau,alpha_seq,prior,sigma=sigma))
    record_tau[t]=tau #record each tau
  }
  return(mean(record_tau))
}

### alpha to lambda calibration
alpha_to_lambda=function(alpha_seq,max_iter=100,prior=x0,second_iter=100,sigma=0){
  tau=sqrt(sigma^2+mean(prior^2)/delta)
  record_tau=rep(0,max_iter) # initialize
  
  for (t in 1:max_iter){
    tau=sqrt(F(tau,alpha_seq,prior,sigma=sigma))
    record_tau[t]=tau #record each tau
  }
  tau=mean(record_tau)
  
  E=0
  for (t in 1:second_iter){
    prox_solution=prox_sorted_L1(prior+tau*rnorm(p),alpha_seq*tau)
    E=E+length(unique(abs(prox_solution)))/p/delta/second_iter
  }
  lambda_seq=(1-E)*alpha_seq*tau
  
  return(lambda_seq)
}
### lambda to alpha calibration (implicitly)
lambda_to_alpha=function(lambda_seq,tol=1e-5,prior=x0,sigma=0){
  # standardize lambda sequence
  lambda_seq=sort(lambda_seq,T)
  l=lambda_seq/lambda_seq[1]
  
  
  ### find interval that has larger and smaller at endpoints compared to lambda_seq[1] for bisection
  # starting guess
  alpha1=l/2
  alpha2=l*2
  
  lambda1=alpha_to_lambda(alpha1,sigma=sigma)
  lambda2=alpha_to_lambda(alpha2,sigma=sigma)
  
  while ((lambda1[1]<lambda_seq[1])*(lambda2[1]>lambda_seq[1])==0){
    if (lambda1[1]<lambda_seq[1]){
      alpha1=alpha1*2
      alpha2=alpha2*2
    } else{
      alpha1=alpha1/2
      alpha2=alpha2/2
    }
    lambda1=alpha_to_lambda(alpha1,sigma=sigma)
    lambda2=alpha_to_lambda(alpha2,sigma=sigma)

  }
  
  
  ### bisection to find the alpha_seq which is parallel to lambda_seq
  while ((alpha2[1]-alpha1[1])>tol){
    middle_alpha_seq=(alpha1+alpha2)/2
    middle_lambda=alpha_to_lambda(middle_alpha_seq,sigma=sigma)
    if (middle_lambda[1]>lambda_seq[1]){
      alpha2=middle_alpha_seq
    }else if (middle_lambda[1]<lambda_seq[1]){
      alpha1=middle_alpha_seq
    }else{
      break
    }
  }
  return(middle_alpha_seq)
}

### Algorithm 1
averaging=function(C,iter=100000){
  # aim C decreasing
  l=length(C)
  iii=1
  while(prod((C[2:l]-C[1:(l-1)])<=0)==0){
    # find the first non-decreasing sub-sequence
    start=0
    whether_start=0
    final=0
    for(i in 2:l){
      # start sub-sequence
      if((C[i]>C[i-1])&(whether_start==0)){
        start=i-1
        whether_start=1
      }
      # end sub-sequence
      if((C[i]<C[i-1])&(whether_start==1)){
        final=i-1
        break
      }
      if(i==l){final=i}
    }
    # average the sub-sequence
    C[start:final]=mean(C[start:final])
    
    ### allow fixed iteration break
    iii=iii+1
    if(iii>iter){break}
  }
  return(C)
}


### Algorithm 2
PGD_new=function(alpha_seq,prior=x0,sigma=0,Z_iter=1000,max_iter=50,step=500,noise=0,plotting=FALSE){ 
  
  p=length(prior);n=round(p*delta);
  tau=alpha_to_tau(alpha_seq=alpha_seq,prior=prior,sigma=sigma)
  print(paste0("Initial tau: ",tau))
  tau_record=c()
  

  for(iter in 1:max_iter){
    
    start=Sys.time()
    print(paste0("Iteration ",iter))
    # initialize
    Enumerator=EDenominator=rep(0,p)

    # using mean value over Z_iter repeats to approximate the expectation in Theorem 1
    for (it in 1:Z_iter){
      Z=rnorm(p)
      prox=prox_sorted_L1(prior+tau*Z,alpha_seq*tau)
      abp=abs(prox)
      
      #unique_magnitude= abp[ave(abp, abp, FUN = length) == 1]
      magnitude= ave(abp, abp, FUN = length)
      
      # group all the variables, create a list Group
      refer = unique(abp)
      Group = list()
      for (i in 1:length(refer)){
        Group[[i]] = which(abp == refer[i])
      }
      
      # calculate \sigma(i)
      sigma_i = order(abp,decreasing=TRUE)
      index=c()
      for(i in 1:p){
        # the group which order_i lies in
        index[i] = which(refer == abp[sigma_i[i]])
        Enumerator[i] = Enumerator[i]+mean(((prox-prior)*sign(prox))[Group[[index[i]]]])*tau/Z_iter
      }
      
    }
    
    #### alpha gradient
    agrad = -Enumerator/n
    plot(agrad,ylim=c(-0.003,0.0015));abline(v=sum(abp!=0))
    
    print("gradient done")
    #GD
    new_alpha=pmax(averaging(alpha_seq-step*agrad),0)
    new_tau=alpha_to_tau(new_alpha,prior=prior,sigma=sigma)
    while (new_tau>tau){
      step = step/2
      print("reducing step size...")
      if(step*mean(abs(agrad))<1e-2){break}
      # noise and projection
      new_alpha=pmax(averaging(alpha_seq-step*agrad-rnorm(p,sd=noise)),0)
      new_tau=alpha_to_tau(new_alpha,prior=prior,sigma=sigma)
    }
    
    if(new_tau<tau){tau=new_tau;alpha_seq=new_alpha;print(tau)}
    else{break}
    
    tau_record=c(tau_record,tau)
    print(Sys.time()-start)
    if(plotting==TRUE){plot(alpha_seq)}
  }
  return(list(alpha=alpha_seq,tau=tau_record))
}



######### Experiments using Algorithm 3

# For the sake of simplicity, we only show the code for dependent case and real data. 
# The independent case of synthetic data can be easily modified from the following ARMA dataset
# (by defining independent design A and corresponding response y)

####### Synthetic data

# ARMA (1,1)
library('MASS')
set.seed(918)
# with correlation!
n = 20
p = 50
# using time series to generate the models
A = matrix(rep(0,n*p),nrow = n)
for(i in 1:n){
  A[i,] = arima.sim(list(order = c(1,0,1),'ar' = 0.8,'ma' = 0.8), n = p)
}
beta = rbinom(p,5,0.3)*rnorm(p)
A = scale(A); R = 1
y = A%*%beta + R*rnorm(n); y = scale(y) 
lam1 = 0.01; record=c();p=ncol(A);scale = c();Accu = c(); Length = 100

# Recall that in k-level SLOPE, there are at most 2k-1 parameters to tune. 
# In 2-level SLOPE, there are only three: s1, s2 and p12.
# In 3-level SLOPE, there are only five: s1,s2,s3 and p12, p23
# In the following code, we set s2 = s3 (=0) so that it's a 2-level SLOPE

# using BH SLOPE, using minimax?


record = c()

for(p12 in seq(0.4,0.4,0.1)){
  for(p23 in seq(0.4,0.4,0.1)){
    for(s1 in c(1e-4,1e-3,1e-2,1e-1)){
      for(s2 in seq(0, 0, 0.00005)){
        for(s3 in seq(0, 0, 1e-4)){
          #for(pass in c(1, 100, 200, 300, 500, seq(1000,15000,1000) ) ){
          for(pass in c(10000)){
            
            Accu = 0
            T = 50
            for(t in 1:T){
              
              # using time series to generate the models
              A = matrix(rep(0,n*p),nrow = n)
              for(i in 1:n){
                A[i,] = arima.sim(list(order = c(1,0,1),'ar' = 0.8,'ma' = 0.8), n = p)
              }
              beta = rbinom(p,5,0.3)*rnorm(p)
              A = scale(A); R = 1
              y = A%*%beta + R*rnorm(n); y = scale(y) 
              #lam1 = 0.01; record=c();p=ncol(A);scale = c();Accu = c(); Length = 100
              
              # lambda sequence
              Lda = c(rep(lam1*s1,round(p*p12)), rep(lam1*s2,round(p*p23)-round(p*p12)),rep(lam1*s3,p-round(p*p23)))
            
              # ARIMA(1,0,1) with 'ar' = 0.8,'ma' = 0.8 lasso: lambda = 0.0001 MSE = 0.2985; 2-level-slope: p12 = 0.4 s1 = 0.1 s2 = 0 MSE = 0.2081
              for(k in 1:20){
                A_test = A[(round((k-1)*1)+1):round(k*1),]
                A_train = A[-c((round((k-1)*1)+1):round(k*1)),]
                y_test = y[(round((k-1)*1)+1):round(k*1)]
                y_train = y[-c((round((k-1)*1)+1):round(k*1))]
                S = SLOPE(A_train, y_train, alpha = 1, lambda=Lda, solver = "admm",max_passes = pass) 
                slope=S$coefficients  
                Accu = Accu + mean((A_test%*%slope[-1,1,1]+slope[1,1,1]-y_test)^2)/20
              }
            }
            
            print(c(p12,p23,s1,s2,s3, Accu/T))
            record=c(record, Accu/T)
            
          }
        }
      }
    }
  }
}




####### Real dataset
#### ASCVD
load("cvd.RData")

c_imp = order(-abs(cor(cvd[,1:4216],cvd[,4219])))[1:1000] 
cvd_data_imp = cbind(cvd[,c_imp],cvd[,4219])
A = scale(as.matrix(cvd[,c_imp]))
y = scale(as.vector(cvd[,4219])) 
lam1 = 0.01 ;record=c();p=ncol(A);scale = c();Accu = c(); Length = 100

# Recall that in k-level SLOPE, there are at most 2k-1 parameters to tune. 
# In 2-level SLOPE, there are only three: s1, s2 and p12.
# In 3-level SLOPE, there are only five: s1,s2,s3 and p12, p23
# In the following code, we set s2 = s3 (=0) so that it's a 2-level SLOPE
record = c()
for(p12 in seq(0.3,0.3,0.1)){
  for(p23 in seq(0.3,0.3,0.1)){
    for(s1 in c(1e-7)){
      for(s2 in seq(s1, s1, 0.00005)){
        for(s3 in seq(s1, s1, 0.2)){
          for(pass in c(1,100,200, 300, 500, seq(1000,15000,1000) ) ){
          #for(pass in c(10000 ) ){
            
            # lambda sequence
            Lda = c(rep(lam1*s1,round(p*p12)), rep(lam1*s2,round(p*p23)-round(p*p12)),rep(lam1*s3,p-round(p*p23)))
            ## used for SLOPE-mnimax
            #Lda = c()
            #for(i in 1:p){
            #  Lda[i] = s1*(log(2*p/i)/n)^0.5
            #}
            
            Accu = 0

            # 20-fold CV for slope. To calculate lasso, simply set s1=s2=s3 in the loop. optimal: 0.5276 when lambda is near zero (1e-8)
            # lasso: s1 = 1e-8
            # slope-2: s1=8e-4 s2=0 optimal: 0.4888 p12 = 0.3
            # slope-3: s1=7e-4 s2=3.5e-4 s3=0 optimal: 0.4884 p12 = 0.2, p23 = 0.4
            # slope-minimax: A=1e-5, optimal: 0.559
            for(k in 1:20){
              A_test = A[(round((k-1)*23.6/2)+1):round(k*23.6/2),]
              A_train = A[-c((round((k-1)*23.6/2)+1):round(k*23.6/2)),]
              y_test = y[(round((k-1)*23.6/2)+1):round(k*23.6/2)]
              y_train = y[-c((round((k-1)*23.6/2)+1):round(k*23.6/2))]
              S = SLOPE(A_train, y_train, alpha = 1, lambda= Lda, solver = "admm",max_passes = pass) 
              ## used for SLOPE-BH, alpha = 1e-8, optimal: 0.531
              #S = SLOPE(A_train, y_train, alpha = s1, lambda="bh", max_passes = pass) 
              slope=S$coefficients  
              Accu = Accu + mean((A_test%*%slope[-1,1,1]+slope[1,1,1]-y_test)^2)/20
            }
            
            
            print(c(p12,p23,s1,s2,s3, Accu))
            record=c(record, Accu)
          }
        }
      }
    }
  }
}
 


#### ADNI

load("adni_pred.RData")
adni = sa_pred 
set.seed(925)
index = sample(300) 
A = scale(as.matrix(adni[index,c(1:500)]))
y = 2*as.vector(adni[index,1003]==1)-1

lam1=1;record=c();p=ncol(A);scale = c();Accu = c(); Length = 100
record = c()
for(prop in seq(0.1, 0.1, 0.2)){
  for(s1 in seq( 0.08, 0.08, 0.01)){
    for(s2 in seq(0, 0, 0.01)){ 
      for(pass in c(1,10, 20, 30, 50, seq(100, 1500, 100) ) ){
        
        Accu = 0
        N = 30

        # lasso 0.62 with lambda = 0.03
        # 2-level-SLOPE 0.66 p12 = 0.1, s1 = 0.06/0.08, s2 = 0
        for(k in 1:10){
          A_test = A[(round((k-1)*N)+1):round(k*N),]
          A_train = A[-c((round((k-1)*N)+1):round(k*N)),]
          y_test = y[(round((k-1)*N)+1):round(k*N)]
          y_train = y[-c((round((k-1)*N)+1):round(k*N))]
          S = SLOPE(A_train, y_train, intercept = T, family = "binomial",alpha = 1, lambda=c(rep(lam1*s1,round(p*prop)),rep(lam1*s2,p-round(p*prop))),max_passes = pass) 
          slope=S$coefficients  
          Accu = Accu + mean((A_test%*%slope[-1,1,1] +slope[1,1,1] > 0)==(y_test > 0))/10 
        }
        
        print(c(prop,s1,s2,(Accu)))
        record=c(record,(Accu))
      }
    }
  }
}


