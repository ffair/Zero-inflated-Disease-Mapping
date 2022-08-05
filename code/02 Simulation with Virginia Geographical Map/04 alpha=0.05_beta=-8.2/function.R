#############################################################################
###         1. Functions to generate the spatial random effect            ###
#############################################################################
phigen=function(tau_in,W,n)  
{
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

epsiongen=function(delta_in,W,n)  
{
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
  
  data<-as.data.frame(cbind(Y,x[,2:7],z[,2:7]))
  names(data)<-c('Y','X2','X3','X4','X5','X6','X7','Z2','Z3','Z4','Z5','Z6','Z7')
  result_zip<-zeroinfl(Y~Z2+Z3+Z4+Z5+Z6+Z7|X2+X3+X4+X5+X6+X7,data = data)
  sum_zip <- summary(result_zip)
  
  sigma_phi=1.5 
  sigma_epsion=4.5 
  sigma_beta=rep(0.08,M)
  sigma_alpha=rep(1.5,x_M)
  
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
  beta_in[1]=0
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
      if(Racc_phi_in[i]>0.4){sigma_phi_in[i]=sigma_phi_in[i]+0.2}
      if(Racc_phi_in[i]<0.2){sigma_phi_in[i]=sigma_phi_in[i]-0.08}
      if(Racc_epsion_in[i]>0.4){sigma_epsion_in[i]=sigma_epsion_in[i]+0.28}
      if(Racc_epsion_in[i]<0.2){sigma_epsion_in[i]=sigma_epsion_in[i]-0.12}
    }
    for(i in 1:M){
      if(Racc_beta_in[i]>0.4){sigma_beta_in[i]=sigma_beta_in[i]+0.02}
      if(Racc_beta_in[i]<0.2){sigma_beta_in[i]=sigma_beta_in[i]-0.01}
    }
    for(i in 1:x_M){
      if(Racc_alpha_in[i]>0.4){sigma_alpha_in[i]=sigma_alpha_in[i]+1.5}
      if(Racc_alpha_in[i]<0.2){sigma_alpha_in[i]=sigma_alpha_in[i]-0.2}
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
