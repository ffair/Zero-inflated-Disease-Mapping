#############################################################################
###         1. Functions to generate the spatial random effect            ###
#############################################################################
phigen=function(tau_in,W,n){
  I=dim(W)[1]
  phi=rep(0,I)
  sn=round(n/3,0)
  res=matrix(0,nrow=sn,ncol=I)
  #using Gibbs for update
  for(t in 1:n)
  {
    for(i in 1:I)
    {
      thesum=sum(W[i,])
      theu=(W[i,]%*%phi)/thesum
      thevar=sqrt(1/(tau_in*thesum))
      phi[i]=rnorm(1,theu,thevar)
    }
    phi=phi-mean(phi)
    if(t>(n-sn))
    {
      res[t-n+sn,]=phi
    }
  }
  phi=apply(res,2,median)
  return(phi)
}

epsiongen=function(delta_in,W,n){
  I=dim(W)[1]
  epsion=rep(0,I)
  sn=round(n/3,0)
  res=matrix(0,nrow=sn,ncol=I)
  #using Gibbs for update
  for(t in 1:n)
  {
    for(i in 1:I)
    {
      thesum=sum(W[i,])
      theu=(W[i,]%*%epsion)/thesum
      thevar=sqrt(1/(delta_in*thesum))
      epsion[i]=rnorm(1,theu,thevar)
    }
    epsion=epsion-mean(epsion)
    if(t>(n-sn))
    {
      res[t-n+sn,]=epsion
    }
  }
  epsion=apply(res,2,median)
  return(epsion)
}



#############################################################################
###                 2. Some fundamental functions for MCMC                ###
#############################################################################
Y0_likelihood=function(i,alpha_in,beta_in,theepsion,thephi,x,z,N){
  ##The liklihood then disease count Y=0
  part1 = log(1+exp(x[i,]%*%alpha_in+theepsion)*exp(-N[i]*exp(z[i,]%*%beta_in+thephi)))
  part2 = -log(1+exp(x[i,]%*%alpha_in+theepsion))
  return (part1+part2)
}

Y1_likelihood =function(i,alpha_in,beta_in,theepsion,thephi,x,z,Y,N){
  ##The liklihood then disease count Y!=0
  part1 = -N[i]*exp(z[i,]%*%beta_in+thephi)
  part2 = Y[i]*(log(N[i])+z[i,]%*%beta_in+thephi)
  part4 = -log(1+exp(-(x[i,]%*%alpha_in+theepsion)))
  return (part1+part2+part4)
}

post_phi1=function(Y,x,z,N,thephi,i,W,epsion_in,phi_in,alpha_in,beta_in,tau_in){
  ##The posterior of spatial random effect phi
  theepsion=epsion_in[i]
  thesum=sum(W[i,])
  phibar=((W[i,]%*%phi_in)/thesum)[1,1]
  part1=-0.5*tau_in*thesum*(thephi-phibar)^2
  part2 = ifelse(Y[i]==0,Y0_likelihood(i,alpha_in,beta_in,theepsion,thephi,x,z,N),
                 Y1_likelihood(i,alpha_in,beta_in,theepsion,thephi,x,z,Y,N))
  res=part1+part2
  return(res)
}

post_epsion1=function(Y,x,z,N,theepsion,i,W,epsion_in,phi_in,alpha_in,beta_in,delta_in){
  ##The posterior of spatial random effect epsion
  thephi=phi_in[i]
  thesum=sum(W[i,])
  epsionbar=((W[i,]%*%epsion_in)/thesum)[1,1]
  part1=-0.5*delta_in*thesum*(theepsion-epsionbar)^2
  part2 = ifelse(Y[i]==0,Y0_likelihood(i,alpha_in,beta_in,theepsion,thephi,x,z,N),
                 Y1_likelihood(i,alpha_in,beta_in,theepsion,thephi,x,z,Y,N))
  res=part1+part2
  return(res)
}

update_phi=function(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in){
  ##Update phi using MH algorithm
  for(i in 1:I)
  {
    thecurr=phi_in[i]
    thestar=rnorm(1,thecurr,sigma_phi_in[i]) 
    logA=post_phi1(Y,x,z,N,thestar,i,W,epsion_in,phi_in,alpha_in,beta_in,tau_in)
    logB=post_phi1(Y,x,z,N,thecurr,i,W,epsion_in,phi_in,alpha_in,beta_in,tau_in)
    logrr=logA-logB
    logu=log(runif(1,0,1))
    if(logrr>logu)
    {
      phi_in[i]=thestar
      Racc_phi_in[i]=Racc_phi_in[i]+1
    }
  }
  phi_in=phi_in-mean(phi_in)
  return(list(phi_in,Racc_phi_in))
}

update_epsion=function(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in){
  ##Update epsion using MH algorithm
  for(i in 1:I)
  {
    thecurr=epsion_in[i]
    thestar=rnorm(1,thecurr,sigma_epsion_in[i]) 
    logA=post_epsion1(Y,x,z,N,thestar,i,W,epsion_in,phi_in,alpha_in,beta_in,delta_in)
    logB=post_epsion1(Y,x,z,N,thecurr,i,W,epsion_in,phi_in,alpha_in,beta_in,delta_in)
    logrr=logA-logB
    logu=log(runif(1,0,1))
    if(logrr>logu)
    {
      epsion_in[i]=thestar
      Racc_epsion_in[i]=Racc_epsion_in[i]+1
    }
  }
  epsion_in=epsion_in-mean(epsion_in)
  return(list(epsion_in,Racc_epsion_in))
}

post_beta1=function(x,z,Y,N,thebeta,i,epsion_in,phi_in,alpha_in,beta_in){
  ##The posterior for coefficient beta
  newbeta=beta_in
  newbeta[i]=thebeta 
  res = vector(length = length(Y))
  for (j in (1:length(Y))){
    theepsion=epsion_in[j]
    thephi=phi_in[j]
    res[j]=ifelse(Y[j]==0,Y0_likelihood(j,alpha_in,newbeta,theepsion,thephi,x,z,N),
                  Y1_likelihood(j,alpha_in,newbeta,theepsion,thephi,x,z,Y,N))
  }
  return (sum(res))
}

post_alpha1=function(x,z,Y,N,thealpha,i,epsion_in,phi_in,alpha_in,beta_in){
  ##The posterior for coefficient alpha
  newalpha=alpha_in
  newalpha[i]=thealpha
  res = vector(length = length(Y))
  for (j in (1:length(Y))){
    theepsion=epsion_in[j]
    thephi=phi_in[j]
    res[j]=ifelse(Y[j]==0,Y0_likelihood(j,newalpha,beta_in,theepsion,thephi,x,z,N),
                  Y1_likelihood(j,newalpha,beta_in,theepsion,thephi,x,z,Y,N))
  }
  return (sum(res))
}

update_beta=function(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in){
  ##Update beta using MH algorithm
  for(i in 1:length(beta_in))
  {
    thecurr=beta_in[i]
    thestar=rnorm(1,thecurr,sigma_beta_in[i])
    logA=post_beta1(x,z,Y,N,thestar,i,epsion_in,phi_in,alpha_in,beta_in)
    logB=post_beta1(x,z,Y,N,thecurr,i,epsion_in,phi_in,alpha_in,beta_in)
    logrr=logA-logB
    logu=log(runif(1,0,1))
    if(logrr>logu)
    {
      beta_in[i]=thestar
      Racc_beta_in[i]=Racc_beta_in[i]+1
    }
  }
  return(list(beta_in,Racc_beta_in))
}

update_alpha=function(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in){
  ##Update alpha using MH algorithm
  for(i in 1:length(alpha_in))
  {
    thecurr=alpha_in[i]
    thestar=rnorm(1,thecurr,sigma_alpha_in[i])
    logA=post_alpha1(x,z,Y,N,thestar,i,epsion_in,phi_in,alpha_in,beta_in)
    logB=post_alpha1(x,z,Y,N,thecurr,i,epsion_in,phi_in,alpha_in,beta_in)
    logrr=logA-logB
    logu=log(runif(1,0,1))
    if(logrr>logu)
    {
      alpha_in[i]=thestar
      Racc_alpha_in[i]=Racc_alpha_in[i]+1
    }
  }
  return(list(alpha_in,Racc_alpha_in))
}

update_tau=function(I,phi_in,a_tau,b_tau,W){
  ##Update tau
  philist=rep(phi_in,I)
  p1=matrix(philist,nrow=I,ncol=I,byrow=T)
  p2=matrix(philist,nrow=I,ncol=I,byrow=F)
  pp=(p1-p2)^2
  ppp=pp[lower.tri(pp)]
  www=W[lower.tri(W)]
  newb=b_tau+sum(ppp*www)/2
  newa=a_tau+I/2
  newtau=rgamma(1,newa,newb)
  return(newtau)
}

update_delta=function(I,epsion_in,a_delta,b_delta,W){
  ##Update delta
  epsionlist=rep(epsion_in,I)
  p1=matrix(epsionlist,nrow=I,ncol=I,byrow=T)
  p2=matrix(epsionlist,nrow=I,ncol=I,byrow=F)
  pp=(p1-p2)^2
  ppp=pp[lower.tri(pp)]
  www=W[lower.tri(W)]
  newb=b_delta+sum(ppp*www)/2
  newa=a_delta+I/2
  newdelta=rgamma(1,newa,newb)
  return(newdelta)
}



#############################################################################
###                       3. MCMC to fit the ZDM model                    ###
#############################################################################
model_new = function(I,Y,z,x,N,W,M,x_M){
  ##MCMC to fit the ZDM model
  ###########################################################################
  ###                      MCMC to fit the ZDM model                      ###
  ###---------------------------------------------------------------------###
  ###     "I": Number of regions in the artificial geographical map       ###
  ###     "Y": Disease count                                              ###
  ###     "z": Covariates in the Poisson count process                    ###
  ###     "x": Covariates in the zero count process                       ###
  ###     "N": Population in risk                                         ###
  ###     "W": the adjacency matrix                                       ###
  ###     "M": Number of covariate in "z"                                 ###
  ###     "x_M": Number of covariate in "x"                               ###
  ###########################################################################
  a_tau=b_tau=1
  a_delta=b_delta=1
  n.T=10000
  n.burn=5000
  n.tune=1000
  times=10
  step = 5
  saveT = (n.T-n.burn)/step
  
  data<-as.data.frame(cbind(Y,x[,2],z[,2:5]))
  names(data)<-c('Y','X','Z2','Z3','Z4','Z5')
  result_zip<-zeroinfl(Y~Z2+Z3+Z4+Z5|X,data = data)
  sum_zip <- summary(result_zip)
  
  sigma_phi=0.5
  sigma_epsion=0.5
  sigma_beta=c(0.06,0.05,0.05,0.05,0.05)
  sigma_alpha=c(0.71,0.72)
  
  Rsigma_phi=rep(0,I)
  Rsigma_epsion=rep(0,I)
  Rsigma_beta=rep(0,M)
  Rsigma_alpha=rep(0,x_M)
  Racc_phi_in=rep(0,I)
  Racc_epsion_in=rep(0,I)
  Racc_beta_in=rep(0,M)
  Racc_alpha_in=rep(0,x_M)
  
  Rphisave=matrix(0,nrow=saveT,ncol=I)
  Repsionsave=matrix(0,nrow=saveT,ncol=I)
  Rbetasave=matrix(0,nrow=saveT,ncol=M)
  Ralphasave=matrix(0,nrow=saveT,ncol=x_M)
  Rtausave=rep(0,saveT)
  Rdeltasave=rep(0,saveT)
  Rprobsave=matrix(0,nrow=saveT,ncol=I) 
  Rquansave=matrix(0,nrow=saveT,ncol=I)
  
  beta_in=sum_zip$coefficients$count[,1]
  alpha_in=sum_zip$coefficients$zero[,1]
  beta_in=as.matrix(beta_in,nrow=M,ncol=1)
  alpha_in = as.matrix(alpha_in,nrow=x_M,ncol=1)
  
  tau_in=rgamma(1,a_tau,b_tau)
  delta_in=rgamma(1,a_delta,b_delta)
  
  phi_in=mvrnorm(I,0,1)[,1]
  epsion_in=mvrnorm(I,0,1)[,1]
  
  sigma_phi_in=rep(sigma_phi,I)
  sigma_epsion_in=rep(sigma_epsion,I)
  sigma_beta_in=sigma_beta
  sigma_alpha_in=sigma_alpha
  
  for(nn in 1:n.burn){
    ######## (1)update each phi using M-H  
    theres=update_phi(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in)
    phi_in=theres[[1]]
    Racc_phi_in=theres[[2]]
    ######## (2)update each epsion using M-H
    theres=update_epsion(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in)
    epsion_in=theres[[1]]
    Racc_epsion_in=theres[[2]]
    ######## (3)update each beta using M-H
    theres=update_beta(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in)
    beta_in=theres[[1]]
    Racc_beta_in=theres[[2]]
    ######## (4)update each alpha using M-H
    theres=update_alpha(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in)
    alpha_in=theres[[1]]
    Racc_alpha_in=theres[[2]]
    ######## (5)update tau
    tau_in=update_tau(I,phi_in,a_tau,b_tau,W)
    ######## (6)update delta
    delta_in=update_delta(I,epsion_in,a_delta,b_delta,W)
  }
  
  for(j in 1:times){
    ##Tune the variance of proposal distribution to control the acceptance rate
    Racc_phi_in=rep(0,I)
    Racc_epsion_in=rep(0,I)
    Racc_beta_in=rep(0,M)
    Racc_alpha_in=rep(0,x_M)
    for(nn in 1:n.tune){
      ######## (1)update each phi using M-H  
      theres=update_phi(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in)
      phi_in=theres[[1]]
      Racc_phi_in=theres[[2]]
      ######## (2)update each epsion using M-H
      theres=update_epsion(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in)
      epsion_in=theres[[1]]
      Racc_epsion_in=theres[[2]]
      ######## (3)update each beta using M-H
      theres=update_beta(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in)
      beta_in=theres[[1]]
      Racc_beta_in=theres[[2]]
      ######## (4)update each alpha using M-H
      theres=update_alpha(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in)
      alpha_in=theres[[1]]
      Racc_alpha_in=theres[[2]]
      ######## (5)update tau
      tau_in=update_tau(I,phi_in,a_tau,b_tau,W)
      ######## (6)update delta
      delta_in=update_delta(I,epsion_in,a_delta,b_delta,W)
    }
    Racc_phi_in=Racc_phi_in/n.tune
    Racc_epsion_in=Racc_epsion_in/n.tune
    Racc_beta_in=Racc_beta_in/n.tune
    Racc_alpha_in=Racc_alpha_in/n.tune
    for(i in 1:I){
      if(Racc_phi_in[i]>0.4){sigma_phi_in[i]=sigma_phi_in[i]+0.12}
      if(Racc_phi_in[i]<0.2){sigma_phi_in[i]=sigma_phi_in[i]-0.08}
      if(Racc_epsion_in[i]>0.4){sigma_epsion_in[i]=sigma_epsion_in[i]+0.28}
      if(Racc_epsion_in[i]<0.2){sigma_epsion_in[i]=sigma_epsion_in[i]-0.06}
    }
    for(i in 1:M){
      if(Racc_beta_in[i]>0.4){sigma_beta_in[i]=sigma_beta_in[i]+0.02}
      if(Racc_beta_in[i]<0.2){sigma_beta_in[i]=sigma_beta_in[i]-0.02}
    }
    for(i in 1:x_M){
      if(Racc_alpha_in[i]>0.4){sigma_alpha_in[i]=sigma_alpha_in[i]+0.08}
      if(Racc_alpha_in[i]<0.2){sigma_alpha_in[i]=sigma_alpha_in[i]-0.08}
    }
  }
  
  Racc_phi_in=rep(0,I)
  Racc_epsion_in=rep(0,I)
  Racc_beta_in=rep(0,M)
  Racc_alpha_in=rep(0,x_M)
  for(nn in 1:n.T){
    ######## (1)update each phi using M-H  
    theres=update_phi(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in)
    phi_in=theres[[1]]
    Racc_phi_in=theres[[2]]
    ######## (2)update each epsion using M-H
    theres=update_epsion(I,W,Y,x,z,N,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in)
    epsion_in=theres[[1]]
    Racc_epsion_in=theres[[2]]
    ######## (3)update each beta using M-H
    theres=update_beta(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in)
    beta_in=theres[[1]]
    Racc_beta_in=theres[[2]]
    ######## (4)update each alpha using M-H
    theres=update_alpha(x,z,Y,N,epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in)
    alpha_in=theres[[1]]
    Racc_alpha_in=theres[[2]]
    ######## (5)update tau
    tau_in=update_tau(I,phi_in,a_tau,b_tau,W)
    ######## (6)update delta
    delta_in=update_delta(I,epsion_in,a_delta,b_delta,W)
    ######## (7)save posteriors after burnin
    if((nn>n.burn)&((nn-n.burn)%%step==0)){
      Rphisave[(nn-n.burn)/step,]=phi_in
      Repsionsave[(nn-n.burn)/step,]=epsion_in
      Rbetasave[(nn-n.burn)/step,]=beta_in
      Ralphasave[(nn-n.burn)/step,]=alpha_in
      Rprobsave[(nn-n.burn)/step,] = exp(x%*%alpha_in+epsion_in)/(1+exp(x%*%alpha_in+epsion_in))
      Rquansave[(nn-n.burn)/step,] = exp(z%*%beta_in+phi_in)
      Rtausave[nn/step] = tau_in
      Rdeltasave[nn/step] = delta_in
    }
  }
  Racc_alpha_in=Racc_alpha_in/n.T
  Racc_beta_in=Racc_beta_in/n.T
  Racc_epsion_in=Racc_epsion_in/n.T
  Racc_phi_in=Racc_phi_in/n.T
  
  beta_each_mean = apply(Rbetasave,2,mean)
  beta_each_sd = apply(Rbetasave,2,sd)
  beta_each_upper = apply(Rbetasave,2,function(x){return(quantile(x,0.95))})
  beta_each_lower = apply(Rbetasave,2,function(x){return(quantile(x,0.05))})
  Rbeta=cbind(beta_each_mean,beta_each_sd,
              beta_each_upper,beta_each_lower)
  
  alpha_each_mean = apply(Ralphasave,2,mean)
  alpha_each_sd = apply(Ralphasave,2,sd)
  alpha_each_upper = apply(Ralphasave,2,function(x){return(quantile(x,0.95))})
  alpha_each_lower = apply(Ralphasave,2,function(x){return(quantile(x,0.05))})
  Ralpha=cbind(alpha_each_mean,alpha_each_sd,
               alpha_each_upper,alpha_each_lower)
  
  prob_each_mean=apply(Rprobsave,2,mean)
  prob_each_sd=apply(Rprobsave,2,sd)
  prob_each_upper=apply(Rprobsave,2,function(x){return(quantile(x,0.95))})
  prob_each_lower=apply(Rprobsave,2,function(x){return(quantile(x,0.05))})
  Rprob=cbind(prob_each_mean,prob_each_sd,prob_each_upper,prob_each_lower)
  
  quan_each_mean=apply(Rquansave,2,mean)
  quan_each_sd=apply(Rquansave,2,sd)
  quan_each_upper=apply(Rquansave,2,function(x){return(quantile(x,0.95))})
  quan_each_lower=apply(Rquansave,2,function(x){return(quantile(x,0.05))})
  Rquan=cbind(quan_each_mean,quan_each_sd,quan_each_upper,quan_each_lower)
  
  tau_each_mean=mean(Rtausave)
  tau_each_sd=sd(Rtausave)
  tau_each_upper=quantile(Rtausave,0.95)
  tau_each_lower=quantile(Rtausave,0.05)
  Rtau=cbind(tau_each_mean,tau_each_sd,tau_each_upper,tau_each_lower)
  
  delta_each_mean=mean(Rdeltasave)
  delta_each_sd=sd(Rdeltasave)
  delta_each_upper=quantile(Rdeltasave,0.95)
  delta_each_lower=quantile(Rdeltasave,0.05)
  Rdelta=cbind(delta_each_mean,delta_each_sd,delta_each_upper,delta_each_lower)
  
  return (list(Ralpha,Rbeta,Rprob,Rquan,Racc_phi_in,Racc_epsion_in,Racc_beta_in,Racc_alpha_in,Rtau,Rdelta))
}



#############################################################################
###                     4. Functions for simulation                       ###
#############################################################################
simul_art_spatial <- function(list_in){
  
  ###########################################################################
  ###                 Simulation function for ZDM, ZIP, POR               ###
  ###---------------------------------------------------------------------###
  ###       "list_in": a list containing simulation design information    ###
  ###                  including correlation "rho", coefficient "alpha"   ###
  ###                  and "beta"                                         ###
  ###########################################################################
  
  
  ###########################################################################
  ###                     4.1. Generate true values                       ###
  ###########################################################################
  ##The adjacency matrix
  W <- read.csv('../../../data/W.csv',header=T)
  W = as.matrix(W)
  I = dim(W)[1]
  ##Generate population
  N_mean = 50000
  N_sd = 10000
  set.seed(5)
  N=rnorm(I,mean=N_mean,sd=N_sd)  
  N=round(N,0)
  ##Generate spatial random effect
  tau_true = 0.1
  delta_true = 0.01
  n=2000
  set.seed(5)
  phi_true = phigen(tau_true,W,n)
  set.seed(5)
  epsion_true = epsiongen(delta_true,W,n)
  ##Generate covariates
  x_M = 2
  M = 5
  x = matrix(NA,I,x_M)
  z = matrix(NA,I,M)
  x[,1] = rep(1,I)
  z[,1] = rep(1,I)
  set.seed(6)
  x[,2] = phigen(1,W,n)
  set.seed(7)
  z[,2] = phigen(2,W,n)
  set.seed(8)
  z[,3] = phigen(1,W,n)
  set.seed(9)
  z[,4] = phigen(0.5,W,n)
  set.seed(10)
  z[,5] = phigen(0.1,W,n)
  
  ##The true value
  beta_true = as.matrix(c(list_in[1],0.1,-0.2,0.3,-0.4),ncol=1)
  alpha_true = as.matrix(c(list_in[2],0.5),ncol=1)
  p_true = exp(x%*%alpha_true+epsion_true)/(1+exp(x%*%alpha_true+epsion_true))
  q_true = exp(z%*%beta_true+phi_true)
  lambda = N*q_true
  q_bar = mean(q_true)
  r_true = q_true/q_bar
  
  ###########################################################################
  ###                4.2. Prepare to start simulation                     ###
  ###########################################################################
  s = 100
  y = matrix(NA,nrow=s,ncol=I)
  
  ##Store the results for ZDM model
  Rbeta_beta_new = list(length = M)
  for (i in 1:M){
    Rbeta_beta_new[[i]] = matrix(NA,nrow=s,ncol=6)
  }
  
  Ralpha_alpha_new = list(length = x_M)
  for (i in 1:x_M){
    Ralpha_alpha_new[[i]] = matrix(NA,nrow=s,ncol=6)
  }
  
  Rprob_prob_new = list(length = I)
  for (i in 1:I){
    Rprob_prob_new[[i]] = matrix(NA,nrow=s,ncol=7)
  }
  
  Rquan_quan_new = list(length = I)
  for (i in 1:I){
    Rquan_quan_new[[i]] = matrix(NA,nrow=s,ncol=7)
  }
  
  Rr_r_new = list(length = I)
  for(i in 1:I){
    Rr_r_new[[i]] = matrix(NA,nrow=s,ncol=4)
  }
  
  Raccbeta_beta_new = matrix(nrow = s,ncol = M)
  Raccalpha_alpha_new = matrix(nrow = s,ncol = x_M)
  Raccphi_phi_new = matrix(nrow = s,ncol = I)
  Raccepsion_epsion_new = matrix(nrow = s,ncol = I)
  
  ##Store the results for ZIP model
  Rbeta_beta_zip = list(length = M)
  for (i in 1:M){
    Rbeta_beta_zip[[i]] = matrix(NA,nrow=s,ncol=6)
  }
  
  Ralpha_alpha_zip = list(length = x_M)
  for (i in 1:x_M){
    Ralpha_alpha_zip[[i]] = matrix(NA,nrow=s,ncol=6)
  }
  
  Rprob_prob_zip = list(length = I)
  for (i in 1:I){
    Rprob_prob_zip[[i]] = matrix(NA,nrow=s,ncol=7)
  }
  
  Rquan_quan_zip = list(length = I)
  for (i in 1:I){
    Rquan_quan_zip[[i]] = matrix(NA,nrow=s,ncol=7)
  }
  
  Rr_r_zip = list(length = I)
  for(i in 1:I){
    Rr_r_zip[[i]] = matrix(NA,nrow=s,ncol=4)
  }
  
  Raccbeta_beta_zip = matrix(nrow = s,ncol = M)
  Raccalpha_alpha_zip = matrix(nrow = s,ncol = x_M)
  
  ##Store the results for POR model
  Rquan_quan_poi = list(length = I)
  for (i in 1:I){
    Rquan_quan_poi[[i]] = matrix(NA,nrow=s,ncol=7)
  }
  
  Rbeta_beta_poi = list(length = M)
  for (i in 1:M){
    Rbeta_beta_poi[[i]] = matrix(NA,nrow=s,ncol=6)
  }
  a = rep(0,s)
  Rr_r_poi = list(length = I)
  for(i in 1:I){
    Rr_r_poi[[i]] = matrix(NA,nrow=s,ncol=4)
  }
  
  ###########################################################################
  ###                     4.3. Start the simulation                       ###
  ###########################################################################
  set.seed(5)
  for(j in 1:s){
    Y=vector(length=I)
    indicator=vector(length=I)
    for(i in c(1:I)){
      indicator[i]=rbinom(1,1,p_true[i])  
    }
    for(i in 1:I){
      Y[i]=ifelse(indicator[i]==0,0,rpois(1,lambda[i]))
    }
    indi = vector(length=I)
    for(i in 1:I){
      indi[i] = ifelse(Y[i]==0,0,1)
    }
    
    ##ZDM
    result_new = model_new(I,Y,z,x,N,W,M,x_M)
    Ralpha = result_new[[1]]
    Rbeta = result_new[[2]]
    Rprob = result_new[[3]]
    Rquan = result_new[[4]]
    Raccphi_phi_new[j,] = result_new[[5]]
    Raccepsion_epsion_new[j,] = result_new[[6]]
    Raccbeta_beta_new[j,] = result_new[[7]]
    Raccalpha_alpha_new[j,] = result_new[[8]]
    
    for (i in 1:M){
      Rbeta_beta_new[[i]][j,1:4] = Rbeta[i,]
      Rbeta_beta_new[[i]][j,5] = (Rbeta[i,1]-beta_true[i,])^2 #MSE
      Rbeta_beta_new[[i]][j,6] = ifelse(beta_true[i]<=Rbeta[i,3]&beta_true[i]>=Rbeta[i,4],1,0)
    }
    
    for (i in 1:x_M){
      Ralpha_alpha_new[[i]][j,1:4] = Ralpha[i,]
      Ralpha_alpha_new[[i]][j,5] = (Ralpha[i,1]-alpha_true[i,])^2 #MSE
      Ralpha_alpha_new[[i]][j,6] = ifelse(alpha_true[i]<=Ralpha[i,3]&alpha_true[i]>=Ralpha[i,4],1,0)
    }
    
    for (i in 1:I){
      Rprob_prob_new[[i]][j,1:4] = Rprob[i,]
      Rprob_prob_new[[i]][j,5] = (Rprob[i,1]-p_true[i,])^2 #MSE
      Rprob_prob_new[[i]][j,6] = abs(Rprob[i,1]-p_true[i,])/p_true[i,] #RSE
      Rprob_prob_new[[i]][j,7] = (log(Rprob[i,1])-log(p_true[i,]))^2
    }
    
    for (i in 1:I){
      Rquan_quan_new[[i]][j,1:4] = Rquan[i,]
      Rquan_quan_new[[i]][j,5] = (Rquan[i,1]-q_true[i,])^2 #MSE
      Rquan_quan_new[[i]][j,6] = abs(Rquan[i,1]-q_true[i,])/q_true[i,] #RSE
      Rquan_quan_new[[i]][j,7] = (log(Rquan[i,1])-log(q_true[i,]))^2
    }
    
    for (i in 1:I){
      Rr_r_new[[i]][j,1] = Rquan[i,1]/mean(Rquan[,1])
      Rr_r_new[[i]][j,2] = (Rr_r_new[[i]][j,1]-r_true[i])^2
      Rr_r_new[[i]][j,3] = (log(Rr_r_new[[i]][j,1])-log(r_true[i,]))^2
    }
    for (i in 1:I){
      if(Rquan_quan_new[[i]][j,1]>1){a[j]=1}
    }
    rm(result_new,Rbeta,Ralpha,Rprob,Rquan)
    
    ##ZIP
    data<-as.data.frame(cbind(Y,x[,2],z[,2:5]))
    names(data)<-c('Y','X','Z2','Z3','Z4','Z5')
    result_zip<-zeroinfl(Y~Z2+Z3+Z4+Z5|X,data = data)
    sum_zip <- summary(result_zip)
    
    for (i in 1:M){
      Rbeta_beta_zip[[i]][j,1] = sum_zip$coefficients$count[i,1]
      Rbeta_beta_zip[[i]][j,2] = sum_zip$coefficients$count[i,2]
      Rbeta_beta_zip[[i]][j,3] = Rbeta_beta_zip[[i]][j,1]-Rbeta_beta_zip[[i]][j,2]*qnorm(0.05)
      Rbeta_beta_zip[[i]][j,4] = Rbeta_beta_zip[[i]][j,1]-Rbeta_beta_zip[[i]][j,2]*qnorm(0.95)
      Rbeta_beta_zip[[i]][j,5] = (Rbeta_beta_zip[[i]][j,1]-beta_true[i,])^2 #MSE
      Rbeta_beta_zip[[i]][j,6] = ifelse(beta_true[i]<=Rbeta_beta_zip[[i]][j,3]&beta_true[i]>=Rbeta_beta_zip[[i]][j,4],1,0)
    }
    
    for (i in 1:x_M){
      Ralpha_alpha_zip[[i]][j,1] = sum_zip$coefficients$zero[i,1]
      Ralpha_alpha_zip[[i]][j,2] = sum_zip$coefficients$zero[i,2]
      Ralpha_alpha_zip[[i]][j,3] = Ralpha_alpha_zip[[i]][j,1]-Ralpha_alpha_zip[[i]][j,2]*qnorm(0.05)
      Ralpha_alpha_zip[[i]][j,4] = Ralpha_alpha_zip[[i]][j,1]-Ralpha_alpha_zip[[i]][j,2]*qnorm(0.95)
      Ralpha_alpha_zip[[i]][j,5] = (Ralpha_alpha_zip[[i]][j,1]-alpha_true[i,])^2 #MSE
      Ralpha_alpha_zip[[i]][j,6] = ifelse(alpha_true[i]<=Ralpha_alpha_zip[[i]][j,3]&alpha_true[i]>=Ralpha_alpha_zip[[i]][j,4],1,0)
    }
    
    alpha_zip = matrix(NA,nrow=x_M,ncol=1)
    alpha_up_zip = matrix(NA,nrow=x_M,ncol=1)
    alpha_low_zip = matrix(NA,nrow=x_M,ncol=1)
    for (i in 1:x_M){
      alpha_zip[i,] = Ralpha_alpha_zip[[i]][j,1]
    }
    for (i in 1:x_M){
      alpha_up_zip[i,] = Ralpha_alpha_zip[[i]][j,3]
    }
    for (i in 1:x_M){
      alpha_low_zip[i,] = Ralpha_alpha_zip[[i]][j,4]
    }
    prob_zip = exp(x%*%alpha_zip)/(1+exp(x%*%alpha_zip))
    prob_up_zip = exp(x%*%alpha_up_zip)/(1+exp(x%*%alpha_up_zip))
    prob_low_zip = exp(x%*%alpha_low_zip)/(1+exp(x%*%alpha_low_zip))
    for (i in 1:I){
      Rprob_prob_zip[[i]][j,1] = prob_zip[i,]
      Rprob_prob_zip[[i]][j,2] = prob_up_zip[i,]
      Rprob_prob_zip[[i]][j,3] = prob_low_zip[i,]
      Rprob_prob_zip[[i]][j,4] = (prob_zip[i,]-p_true[i,])^2 #MSE
      Rprob_prob_zip[[i]][j,5] = abs(prob_zip[i,]-p_true[i,])/p_true[i,]
      Rprob_prob_zip[[i]][j,6] = (log(prob_zip[i,])-log(p_true[i,]))^2
      Rprob_prob_zip[[i]][j,7] = ifelse(p_true[i]<=prob_up_zip[i,]&p_true[i]>=prob_low_zip[i,],1,0)
    }
    
    beta_zip = matrix(NA,nrow=M,ncol=1)
    beta_up_zip = matrix(NA,nrow=M,ncol=1)
    beta_low_zip = matrix(NA,nrow=M,ncol=1)
    for (i in 1:M){
      beta_zip[i,] = Rbeta_beta_zip[[i]][j,1]
    }
    for (i in 1:M){
      beta_up_zip[i,] = Rbeta_beta_zip[[i]][j,3]
    }
    for (i in 1:M){
      beta_low_zip[i,] = Rbeta_beta_zip[[i]][j,4]
    }
    quan_zip = exp(z%*%beta_zip)/N
    quan_up_zip = exp(z%*%beta_up_zip)/N
    quan_low_zip = exp(z%*%beta_low_zip)/N
    
    for (i in 1:I){
      Rquan_quan_zip[[i]][j,1] = quan_zip[i,]
      Rquan_quan_zip[[i]][j,2] = quan_up_zip[i,]
      Rquan_quan_zip[[i]][j,3] = quan_low_zip[i,]
      Rquan_quan_zip[[i]][j,4] = (quan_zip[i,]-q_true[i,])^2 #MSE
      Rquan_quan_zip[[i]][j,5] = abs(quan_zip[i,]-q_true[i,])/q_true[i,]
      Rquan_quan_zip[[i]][j,6] = (log(quan_zip[i,])-log(q_true[i,]))^2
      Rquan_quan_zip[[i]][j,7] = ifelse(q_true[i]<=quan_up_zip[i,]&q_true[i]>=quan_low_zip[i,],1,0)
    }
    
    for (i in 1:I){
      Rr_r_zip[[i]][j,1] = quan_zip[i]/mean(quan_zip)
      Rr_r_zip[[i]][j,2] = (Rr_r_zip[[i]][j,1]-r_true[i])^2
      Rr_r_zip[[i]][j,3] = abs(Rr_r_zip[[i]][j,1]-r_true[i,])/r_true[i,]
      Rr_r_zip[[i]][j,4] = (log(Rr_r_zip[[i]][j,1])-log(r_true[i,]))^2
    }
    for (i in 1:I){
      if(Rquan_quan_zip[[i]][j,1]>1){a[j]=1}
    }
    rm(result_zip,sum_zip)
    
    ##POR
    sum_poi <- summary(glm(Y~z[,2]+z[,3]+z[,4]+z[,5],family=poisson(link="log")))
    for (i in 1:M){
      Rbeta_beta_poi[[i]][j,1] = sum_poi$coefficients[i,1]
      Rbeta_beta_poi[[i]][j,2] = sum_poi$coefficients[i,2]
      Rbeta_beta_poi[[i]][j,3] = Rbeta_beta_poi[[i]][j,1]-Rbeta_beta_poi[[i]][j,2]*qnorm(0.05)
      Rbeta_beta_poi[[i]][j,4] = Rbeta_beta_poi[[i]][j,1]-Rbeta_beta_poi[[i]][j,2]*qnorm(0.95)
      Rbeta_beta_poi[[i]][j,5] = (Rbeta_beta_poi[[i]][j,1]-beta_true[i,])^2
      Rbeta_beta_poi[[i]][j,6] = ifelse(beta_true[i]<=Rbeta_beta_poi[[i]][j,3]&beta_true[i]>=Rbeta_beta_poi[[i]][j,4],1,0)
    }
    
    beta_poi = matrix(NA,nrow=M,ncol=1)
    beta_up_poi = matrix(NA,nrow=M,ncol=1)
    beta_low_poi = matrix(NA,nrow=M,ncol=1)
    for (i in 1:M){
      beta_poi[i,] = Rbeta_beta_poi[[i]][j,1]
    }
    for (i in 1:M){
      beta_up_poi[i,] = Rbeta_beta_poi[[i]][j,3]
    }
    for (i in 1:M){
      beta_low_poi[i,] = Rbeta_beta_poi[[i]][j,4]
    }
    quan_poi = exp(z%*%beta_poi)/N
    quan_up_poi = exp(z%*%beta_up_poi)/N
    quan_low_poi = exp(z%*%beta_low_poi)/N
    for (i in 1:I){
      Rquan_quan_poi[[i]][j,1] = quan_poi[i,]
      Rquan_quan_poi[[i]][j,2] = quan_up_poi[i,]
      Rquan_quan_poi[[i]][j,3] = quan_low_poi[i,]
      Rquan_quan_poi[[i]][j,4] = (quan_poi[i,]-q_true[i,])^2 #MSE
      Rquan_quan_poi[[i]][j,5] = abs(quan_poi[i,]-q_true[i,])/q_true[i,]
      Rquan_quan_poi[[i]][j,6] = (log(quan_poi[i,])-log(q_true[i,]))^2
      Rquan_quan_poi[[i]][j,7] = ifelse(q_true[i]<=quan_up_poi[i,]&q_true[i]>=quan_low_poi[i,],1,0)
    }
    
    for (i in 1:I){
      Rr_r_poi[[i]][j,1] = quan_poi[i]/mean(quan_poi)
      Rr_r_poi[[i]][j,2] = (Rr_r_poi[[i]][j,1]-r_true[i])^2
      Rr_r_poi[[i]][j,3] = abs(Rr_r_poi[[i]][j,1]-r_true[i,])/r_true[i,]
      Rr_r_poi[[i]][j,4] = (log(Rr_r_poi[[i]][j,1])-log(r_true[i,]))^2
    }
    for (i in 1:I){
      if(Rquan_quan_poi[[i]][j,1]>1){a[j]=1}
    }
    rm(sum_poi)
    
    y[j,] = Y
  }
  
  ###########################################################################
  ###                     4.4. Summarize the results                      ###
  ###########################################################################
  r_all<-matrix(nrow=I,ncol=6)
  quan_all<-matrix(nrow = I,ncol = 6)
  prob_all<-matrix(nrow = I,ncol = 4)
  beta_all<-matrix(nrow = M,ncol = 9)
  alpha_all<-matrix(nrow = x_M,ncol = 6)
  for(i in 1:I){
    ##Relative risk
    Rr_r_new[[i]] = Rr_r_new[[i]][which(a==0),]
    Rr_r_zip[[i]] = Rr_r_zip[[i]][which(a==0),]
    Rr_r_poi[[i]] = Rr_r_poi[[i]][which(a==0),]
    ###MSE
    r_all[i,1] = mean(Rr_r_new[[i]][,2])
    r_all[i,2] = mean(Rr_r_zip[[i]][,2])
    r_all[i,3] = mean(Rr_r_poi[[i]][,2])
    ###LMSE
    r_all[i,4] = mean(Rr_r_new[[i]][,3])
    r_all[i,5] = mean(Rr_r_zip[[i]][,4])
    r_all[i,6] = mean(Rr_r_poi[[i]][,4])
    ##Incidence rate
    Rquan_quan_new[[i]] = Rquan_quan_new[[i]][which(a==0),]
    Rquan_quan_zip[[i]] = Rquan_quan_zip[[i]][which(a==0),]
    Rquan_quan_poi[[i]] = Rquan_quan_poi[[i]][which(a==0),]
    ###MSE
    quan_all[i,1] = mean(Rquan_quan_new[[i]][,5])
    quan_all[i,2] = mean(Rquan_quan_zip[[i]][,4])
    quan_all[i,3] = mean(Rquan_quan_poi[[i]][,4])
    ###LMSE
    quan_all[i,4] = mean(Rquan_quan_new[[i]][,7])
    quan_all[i,5] = mean(Rquan_quan_zip[[i]][,6])
    quan_all[i,6] = mean(Rquan_quan_poi[[i]][,6])
    ##Probability of zero count process
    ###MSE
    prob_all[i,1] = mean(Rprob_prob_new[[i]][,5])
    prob_all[i,2] = mean(Rprob_prob_zip[[i]][,4])
    ###LMSE
    prob_all[i,3] = mean(Rprob_prob_new[[i]][,7])
    prob_all[i,4] = mean(Rprob_prob_zip[[i]][,6])
  }
  for(i in 1:M){
    ##Coefficients beta
    Rbeta_beta_new[[i]] = Rbeta_beta_new[[i]][which(a==0),]
    Rbeta_beta_zip[[i]] = Rbeta_beta_zip[[i]][which(a==0),]
    Rbeta_beta_poi[[i]] = Rbeta_beta_poi[[i]][which(a==0),]
    ###MSE
    beta_all[i,1] = mean(Rbeta_beta_new[[i]][,5])
    beta_all[i,2] = mean(Rbeta_beta_zip[[i]][,5])
    beta_all[i,3] = mean(Rbeta_beta_poi[[i]][,5])
    ###Bias^2
    beta_all[i,4] = (mean(Rbeta_beta_new[[i]][,1])-beta_true[i])^2
    beta_all[i,5] = (mean(Rbeta_beta_zip[[i]][,1])-beta_true[i])^2
    beta_all[i,6] = (mean(Rbeta_beta_poi[[i]][,1])-beta_true[i])^2
    ###Var
    beta_all[i,7] = var(Rbeta_beta_new[[i]][,1])*(s-1)/s
    beta_all[i,8] = var(Rbeta_beta_zip[[i]][,1])*(s-1)/s
    beta_all[i,9] = var(Rbeta_beta_poi[[i]][,1])*(s-1)/s
  }
  for(i in 1:x_M){
    ##Coefficients alpha
    ###MSE
    alpha_all[i,1] = mean(Ralpha_alpha_new[[i]][,5])
    alpha_all[i,2] = mean(Ralpha_alpha_zip[[i]][,5])
    ###Bias^2
    alpha_all[i,3] = (mean(Ralpha_alpha_new[[i]][,1])-alpha_true[i])^2
    alpha_all[i,4] = (mean(Ralpha_alpha_zip[[i]][,1])-alpha_true[i])^2
    ###Var
    alpha_all[i,5] = var(Ralpha_alpha_new[[i]][,1])*(s-1)/s
    alpha_all[i,6] = var(Ralpha_alpha_zip[[i]][,1])*(s-1)/s
  }
  colnames(quan_all) = c('new_MSE','zip_MSE','poi_MSE','new_log','zip_log','poi_log')
  colnames(prob_all) = c('new_MSE','zip_MSE','new_log','zip_log')
  colnames(alpha_all) = c('new_MSE','zip_MSE','new_sqBias','zip_sqBias','new_Var','zip_Var')
  colnames(beta_all) = c('new_MSE','zip_MSE','poi_MSE','new_sqBias','zip_sqBias','poi_sqBias','new_Var','zip_Var','poi_Var')
  colnames(r_all) = c('new_MSE','zip_MSE','poi_MSE','new_log','zip_log','poi_log')
  
  beta_integ = apply(beta_all[-1,],2,mean)
  quan_integ = apply(quan_all,2,mean)
  r_integ = apply(r_all,2,mean)
  
  return(list(beta_integ,
              quan_integ,
              r_integ,
              y))
  
}



#############################################################################
###               5. Functions for plot the results                       ###
#############################################################################
plot_beta <- function(n){
  p <- beta_integ_all[n,6:15]
  a = 5
  b = a/5
  data.1<-data.frame(PE=c(-min(p[1],a),min(p[6],a), -min(p[2],a),min(p[7],a),-min(p[5],a),min(p[10],a),-min(p[4],a),min(p[9],a),-min(p[3],a),min(p[8],a)),
                     method=factor(rep(c("ZDM", "ZIP","ZBSJ", "DM","POR"), each = 2), 
                                   levels = c("POR","DM","ZBSJ","ZIP", "ZDM")),
                     type=factor(c(8,7,10,9,4,3,2,1,6,5)))
  t<-expression(paste(rho,"=",0))
  p1<-ggplot(data=data.1, aes(x=method, y=PE, fill=type))+geom_bar(stat = "identity",position = "identity",width = 0.75)+theme(axis.text.x=element_text(angle=90,colour="black"))+coord_flip()
  p1<-p1+scale_fill_brewer(palette = "Paired")+
    theme_bw()+
    scale_y_continuous(limits=c(-a,a), breaks=c(seq(-a,a,b)),
                       labels =c("0.04+", as.character(c(seq(a-b,b,-b),0,seq(b,a-b,b))), "0.04+") )+
    labs(x="",y="")+
    guides(fill=F)+
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(t=10,r=10,b=0,l=0))
  return(p1)
}
