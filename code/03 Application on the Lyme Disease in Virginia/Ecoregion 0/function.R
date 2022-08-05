################################################################
###              Some functions for MCMC process             ###
################################################################
Y0_likelihood=function(i,alpha_in,beta_in,theepsion,thephi){
  part1 = ifelse(x[i,]%*%alpha_in+theepsion-N[i]*exp(z[i,]%*%beta_in+thephi)>701,
                 x[i,]%*%alpha_in+theepsion-N[i]*exp(z[i,]%*%beta_in+thephi),
                 log(1+exp(x[i,]%*%alpha_in+theepsion-N[i]*exp(z[i,]%*%beta_in+thephi)))) 
  part2 = ifelse(x[i,]%*%alpha_in+theepsion>701,
                 -(x[i,]%*%alpha_in+theepsion),
                 -log(1+exp(x[i,]%*%alpha_in+theepsion)))
  return (part1+part2)
}

Y1_likelihood =function(i,alpha_in,beta_in,theepsion,thephi){
  part1 = -N[i]*exp(z[i,]%*%beta_in+thephi)
  part2 = Y[i]*(log(N[i])+z[i,]%*%beta_in+thephi)
  part4 = ifelse(-(x[i,]%*%alpha_in+theepsion)>701,
                 (x[i,]%*%alpha_in+theepsion),
                 -log(1+exp(-(x[i,]%*%alpha_in+theepsion))))
  return (part1+part2+part4)               
}

post_phi1=function(thephi,i,W,epsion_in,phi_in,alpha_in,beta_in,tau_in)
{
  theepsion=epsion_in[i]
  thesum=sum(W[i,])
  phibar=((W[i,]%*%phi_in)/thesum)[1,1] 
  part1=-0.5*tau_in*thesum*(thephi-phibar)^2
  part2 = ifelse(Y[i]==0,Y0_likelihood(i,alpha_in,beta_in,theepsion,thephi),
                 Y1_likelihood(i,alpha_in,beta_in,theepsion,thephi))
  res=part1+part2
  return(res)
}

post_epsion1=function(theepsion,i,W,epsion_in,phi_in,alpha_in,beta_in,delta_in)
{
  thephi=phi_in[i]
  thesum=sum(W[i,])
  epsionbar=((W[i,]%*%epsion_in)/thesum)[1,1] 
  part1=-0.5*delta_in*thesum*(theepsion-epsionbar)^2
  part2 = ifelse(Y[i]==0,Y0_likelihood(i,alpha_in,beta_in,theepsion,thephi),
                 Y1_likelihood(i,alpha_in,beta_in,theepsion,thephi))
  res=part1+part2
  return(res)
}

update_phi=function(W,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in)
{
  for(i in 1:I)
  {
    thecurr=phi_in[i]
    thestar=rnorm(1,thecurr,sigma_phi_in[i]) 
    logA=post_phi1(thestar,i,W,epsion_in,phi_in,alpha_in,beta_in,tau_in)
    logB=post_phi1(thecurr,i,W,epsion_in,phi_in,alpha_in,beta_in,tau_in)
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
update_epsion=function(W,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in)
{
  for(i in 1:I)
  {
    thecurr=epsion_in[i]
    thestar=rnorm(1,thecurr,sigma_epsion_in[i]) 
    logA=post_epsion1(thestar,i,W,epsion_in,phi_in,alpha_in,beta_in,delta_in)
    logB=post_epsion1(thecurr,i,W,epsion_in,phi_in,alpha_in,beta_in,delta_in)
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

post_beta1=function(thebeta,i,epsion_in,phi_in,alpha_in,beta_in)
{
  newbeta=beta_in
  newbeta[i]=thebeta 
  res1 = vector(length = length(Y))
  for (j in (1:length(Y))){
    theepsion=epsion_in[j]
    thephi=phi_in[j]
    res1[j]=ifelse(Y[j]==0,Y0_likelihood(j,alpha_in,newbeta,theepsion,thephi),
                   Y1_likelihood(j,alpha_in,newbeta,theepsion,thephi))
  }
  res = sum(res1)+dnorm(thebeta,mean=0,sd=2,log=TRUE)
  return (res)
}

post_alpha1=function(thealpha,i,epsion_in,phi_in,alpha_in,beta_in)
{
  newalpha=alpha_in
  newalpha[i]=thealpha
  res1 = vector(length = length(Y))
  for (j in (1:length(Y))){
    theepsion=epsion_in[j]
    thephi=phi_in[j]
    res1[j]=ifelse(Y[j]==0,Y0_likelihood(j,newalpha,beta_in,theepsion,thephi),
                   Y1_likelihood(j,newalpha,beta_in,theepsion,thephi))
  }
  res2=dnorm(thealpha,mean=0,sd=2,log=TRUE)
  res = sum(res1)+res2
  return (res)
}

update_beta=function(epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in)
{
  for(i in 1:length(beta_in))
  {
    thecurr=beta_in[i]
    thestar=rnorm(1,thecurr,sigma_beta_in[i])
    logA=post_beta1(thestar,i,epsion_in,phi_in,alpha_in,beta_in)
    logB=post_beta1(thecurr,i,epsion_in,phi_in,alpha_in,beta_in)
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

update_alpha=function(epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in)
{
  for(i in 1:length(alpha_in))
  {
    thecurr=alpha_in[i]
    thestar=rnorm(1,thecurr,sigma_alpha_in[i])
    logA=post_alpha1(thestar,i,epsion_in,phi_in,alpha_in,beta_in)
    logB=post_alpha1(thecurr,i,epsion_in,phi_in,alpha_in,beta_in)
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

update_tau=function(phi_in,a_tau,b_tau,W)
{
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

update_delta=function(epsion_in,a_delta,b_delta,W)
{
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