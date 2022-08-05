########################################################################
###   The ZIP model used for comparison needs the R package 'pscl' to 
###   complete the model fitting process, and it can be installed via 'install.packages("pscl")'
###   The DM and ZBSJ models used for comparison needs the R package 'INLA', 
###   and it can be installed via 'install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)'
###   The details of 'INLA' package can be found in 'https://www.r-inla.org/home'
#########################################################################
install.packages("pscl")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
rm(list=ls(all=TRUE))
library(MASS)
library(pscl)
library(INLA)

source("function.R")

################################################################
###                       1. Load data                       ###
################################################################
setwd("../../../data")
lyme <- read.csv("lyme_data.csv",encoding='utf-8')
W <- read.csv("W0.csv")
W = as.matrix(W[,-1])
eco = lyme[,1]
N = lyme$n[which(eco==0)]
Y = lyme$disea[which(eco==0)]
z = cbind(1,lyme[which(eco==0),c(5,6,10,12,15,17)])
z = as.matrix(z)
x = z
I = dim(z)[1]
M = dim(z)[2]
x_M = M



################################################################
###                     2. Model fitting                     ###
###----------------------------------------------------------###
###                         2.1. ZDM                         ###
################################################################
a_tau=b_tau=1 
a_delta=b_delta=1
n.T=10000
n.burn=3000
n.tune=1000
times=10
step = 5
saveT = (n.T-n.burn)/step

data<-as.data.frame(cbind(Y,x[,2:7],z[,2:7]))
names(data)<-c('Y','X2','X3','X4','X5','X6','X7','Z2','Z3','Z4','Z5','Z6','Z7')
set.seed(11)
result_zip<-zeroinfl(Y~Z2+Z3+Z4+Z5+Z6+Z7|X2+X3+X4+X5+X6+X7,data = data)
sum_zip <- summary(result_zip)
sigma_phi=1.5 
sigma_epsion=3 
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

set.seed(8)
tau_in=rgamma(1,a_tau,b_tau)
set.seed(9)
delta_in=rgamma(1,a_delta,b_delta)

set.seed(8)
phi_in=mvrnorm(I,0,1)[,1]
set.seed(9)
epsion_in=mvrnorm(I,0,1)[,1]

sigma_phi_in=rep(sigma_phi,I)
sigma_epsion_in=rep(sigma_epsion,I)
sigma_beta_in=sigma_beta
sigma_alpha_in=sigma_alpha

set.seed(11)
for(nn in 1:n.burn){
  print(nn)
  ######## (1)update each phi using M-H  
  theres=update_phi(W,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in)
  phi_in=theres[[1]]
  Racc_phi_in=theres[[2]]
  ######## (2)update each epsion using M-H
  theres=update_epsion(W,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in)
  epsion_in=theres[[1]]
  Racc_epsion_in=theres[[2]]
  ######## (3)update each beta using M-H
  theres=update_beta(epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in)
  beta_in=theres[[1]]
  Racc_beta_in=theres[[2]]
  ######## (4)update each alpha using M-H
  theres=update_alpha(epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in)
  alpha_in=theres[[1]]
  Racc_alpha_in=theres[[2]]
  ######## (5)update tau
  tau_in=update_tau(phi_in,a_tau,b_tau,W)
  ######## (6)update delta
  delta_in=update_delta(epsion_in,a_delta,b_delta,W)
}

set.seed(11)
for(j in 1:times){
  Racc_phi_in=rep(0,I) 
  Racc_epsion_in=rep(0,I)
  Racc_beta_in=rep(0,M)
  Racc_alpha_in=rep(0,x_M)
  for(nn in 1:n.tune){
    print(nn)
    ######## (1)update each phi using M-H  
    theres=update_phi(W,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in)
    phi_in=theres[[1]]
    Racc_phi_in=theres[[2]]
    ######## (2)update each epsion using M-H
    theres=update_epsion(W,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in)
    epsion_in=theres[[1]]
    Racc_epsion_in=theres[[2]]
    ######## (3)update each beta using M-H
    theres=update_beta(epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in)
    beta_in=theres[[1]]
    Racc_beta_in=theres[[2]]
    ######## (4)update each alpha using M-H
    theres=update_alpha(epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in)
    alpha_in=theres[[1]]
    Racc_alpha_in=theres[[2]]
    ######## (5)update tau
    tau_in=update_tau(phi_in,a_tau,b_tau,W)
    ######## (6)update delta
    delta_in=update_delta(epsion_in,a_delta,b_delta,W)
  }
  Racc_phi_in=Racc_phi_in/n.tune
  Racc_epsion_in=Racc_epsion_in/n.tune
  Racc_beta_in=Racc_beta_in/n.tune
  Racc_alpha_in=Racc_alpha_in/n.tune
  for(i in 1:I){
    if(Racc_phi_in[i]>0.4){sigma_phi_in[i]=sigma_phi_in[i]+0.2}
    if(Racc_phi_in[i]<0.2){sigma_phi_in[i]=sigma_phi_in[i]-0.08}
    if(Racc_epsion_in[i]>0.4){sigma_epsion_in[i]=sigma_epsion_in[i]+0.28}
    if(Racc_epsion_in[i]<0.2){sigma_epsion_in[i]=sigma_epsion_in[i]-0.3}
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
set.seed(11)
for(nn in 1:n.T){
  print(nn)
  ######## (1)update each phi using M-H  
  theres=update_phi(W,epsion_in,phi_in,alpha_in,beta_in,tau_in,sigma_phi_in,Racc_phi_in)
  phi_in=theres[[1]]
  Racc_phi_in=theres[[2]]
  ######## (2)update each epsion using M-H
  theres=update_epsion(W,epsion_in,phi_in,alpha_in,beta_in,delta_in,sigma_epsion_in,Racc_epsion_in)
  epsion_in=theres[[1]]
  Racc_epsion_in=theres[[2]]
  ######## (3)update each beta using M-H
  theres=update_beta(epsion_in,phi_in,alpha_in,beta_in,sigma_beta_in,Racc_beta_in)
  beta_in=theres[[1]]
  Racc_beta_in=theres[[2]]
  ######## (4)update each alpha using M-H
  theres=update_alpha(epsion_in,phi_in,alpha_in,beta_in,sigma_alpha_in,Racc_alpha_in)
  alpha_in=theres[[1]]
  Racc_alpha_in=theres[[2]]
  ######## (5)update tau
  tau_in=update_tau(phi_in,a_tau,b_tau,W)
  ######## (6)update delta
  delta_in=update_delta(epsion_in,a_delta,b_delta,W)
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
beta_each_upper = apply(Rbetasave,2,function(x){return(quantile(x,0.975))})
beta_each_lower = apply(Rbetasave,2,function(x){return(quantile(x,0.025))})
Rbeta=cbind(beta_each_mean,beta_each_sd,
            beta_each_upper,beta_each_lower)

quan_each_mean=apply(Rquansave,2,mean)
quan_each_sd=apply(Rquansave,2,sd)
quan_each_upper=apply(Rquansave,2,function(x){return(quantile(x,0.975))})
quan_each_lower=apply(Rquansave,2,function(x){return(quantile(x,0.025))})
Rquan=cbind(quan_each_mean,quan_each_sd,quan_each_upper,quan_each_lower)

##Log likelihood for ZDM
Y0_likelihood1=function(i,alpha_in,beta_in,theepsion,thephi){
  part1 = ifelse(x[i,]%*%alpha_in+theepsion-N[i]*exp(z[i,]%*%beta_in+thephi)>701,
                 x[i,]%*%alpha_in+theepsion-N[i]*exp(z[i,]%*%beta_in+thephi),
                 log(1+exp(x[i,]%*%alpha_in+theepsion-N[i]*exp(z[i,]%*%beta_in+thephi)))) 
  part2 = ifelse(x[i,]%*%alpha_in+theepsion>701,
                 -(x[i,]%*%alpha_in+theepsion),
                 -log(1+exp(x[i,]%*%alpha_in+theepsion)))
  return (part1+part2)
}
Y1_likelihood1 =function(i,alpha_in,beta_in,theepsion,thephi){
  part1 = -N[i]*exp(z[i,]%*%beta_in+thephi)
  part2 = Y[i]*(log(N[i])+z[i,]%*%beta_in+thephi)
  part4 = ifelse(-(x[i,]%*%alpha_in+theepsion)>701,
                 (x[i,]%*%alpha_in+theepsion),
                 -log(1+exp(-(x[i,]%*%alpha_in+theepsion))))
  part3 = -log(factorial(Y[i]))
  return (part1+part2+part3+part4)               
}
res = vector(length = length(Y))
for (j in (1:length(Y))){
  thephi = phi_in[j]
  theepsion = epsion_in[j]
  res[j]=ifelse(Y[j]==0,Y0_likelihood1(j,alpha_in,beta_in,theepsion,thephi),
                Y1_likelihood1(j,alpha_in,beta_in,theepsion,thephi))
}
loglike_zdm = sum(res)

beta_mean_zdm = beta_each_mean[2:7]
beta_low_zdm = beta_each_lower[2:7]
beta_up_zdm = beta_each_upper[2:7]

################################################################
###                         2.2. ZIP                         ###
################################################################
data<-as.data.frame(cbind(Y,x[,2:7],z[,2:7]))
names(data)<-c('Y','X2','X3','X4','X5','X6','X7','Z2','Z3','Z4','Z5','Z6','Z7')
result_zip<-zeroinfl(Y~Z2+Z3+Z4+Z5+Z6+Z7|X2+X3+X4+X5+X6+X7,data = data)
sum_zip <- summary(result_zip)
beta_zip = sum_zip$coefficients$count
beta_mean_zip = beta_zip[2:7,1]
beta_low_zip = (beta_zip[,1]+qnorm(0.025)*beta_zip[,2])[2:7]
beta_up_zip = (beta_zip[,1]+qnorm(0.975)*beta_zip[,2])[2:7]
loglike_zip = sum_zip$loglik

################################################################
###                         2.3. DM                          ###
################################################################
g <- inla.read.graph(W)
s.index3 = seq(1,I,1)
dt1 = data.frame(Y=Y,z1=z[,2],z2=z[,3],z3=z[,4],z4=z[,5],z5=z[,6],z6=z[,7],N=N)
formula2 <- Y ~ z1 + z2 + z3 + z4 + z5 + z6 + f(s.index3, model="besag", graph=g, hyper=list(prec=list(prior="loggamma", param=c(1, 1))),scale.model = TRUE)
res2 <- inla(formula2, data=dt1, family="poisson", E=N)
sum_dm = summary(res2)
beta_dm = sum_dm$fixed
beta_mean_dm = beta_dm[2:7,1]
beta_low_dm = beta_dm[2:7,3]
beta_up_dm = beta_dm[2:7,5]
loglike_dm = sum_dm$mlik[1]

################################################################
###                        2.4. POR                          ###
################################################################
glmm1 = glm(Y~z[,2]+z[,3]+z[,4]+z[,5]+z[,6]+z[,7],family=poisson(link="log"))
sum_poi <- summary(glmm1)
beta_poi = sum_poi$coefficients
beta_mean_poi = beta_poi[2:7,1]
beta_low_poi = (beta_poi[,1]+qnorm(0.025)*beta_poi[,2])[2:7]
beta_up_poi = (beta_poi[,1]+qnorm(0.975)*beta_poi[,2])[2:7]
loglike_poi = as.numeric(logLik(glmm1))

################################################################
###                       2.5. ZBSJ                          ###
################################################################
s.index4 = c(seq(1,I,1),rep(NA,I))
s.index5 = c(rep(NA,I),seq(1,I,1))
YY = matrix(NA,2*I,2)
indi = vector(length=I)
for(i in 1:I){
  indi[i] = ifelse(Y[i]==0,0,1)
}
YY[1:I,1] = indi
YY[(I+1):dim(YY)[1],2] = Y
dt1 = list(
  YY = YY,
  mu.x = c(rep(1,I),rep(NA,I)),
  mu.z = c(rep(NA,I),rep(1,I)),
  cov.x1 = c(x[,2],rep(NA,I)),
  cov.x2 = c(x[,3],rep(NA,I)),
  cov.x3 = c(x[,4],rep(NA,I)),
  cov.x4 = c(x[,5],rep(NA,I)),
  cov.x5 = c(x[,6],rep(NA,I)),
  cov.x6 = c(x[,7],rep(NA,I)),
  cov.z1 = c(rep(NA,I),z[,2]),
  cov.z2 = c(rep(NA,I),z[,3]),
  cov.z3 = c(rep(NA,I),z[,4]),
  cov.z4 = c(rep(NA,I),z[,5]),
  cov.z5 = c(rep(NA,I),z[,6]),
  cov.z6 = c(rep(NA,I),z[,7]),
  N = c(rep(NA,I),N)
)
pp_bar = sum(Y)/sum(N)
E = N*pp_bar
dt1$E = c(rep(NA,I),E)
formula3 <- YY ~ 0 + mu.x + mu.z + cov.x1 + cov.x2 + cov.x3 + cov.x4 + cov.x5 + cov.x6 + cov.z1 + cov.z2 + cov.z3 + cov.z4 + cov.z5 + cov.z6 +
  f(s.index4, model="bym2", graph=g, scale.model = TRUE, constr=TRUE, 
    hyper=list(phi=list(prior="pc",param=c(0.5,2/3),initial=3),prec=list(prior="pc.prec",param=c(1,0.01),initial=1.5)))+
  f(s.index5,copy='s.index4',fixed=FALSE)
res3 <- inla(formula3, data=dt1, 
             family=c("binomial","poisson"),E=E)
sum_inla = summary(res3)
beta_inla = sum_inla$fixed[c(2,9,10,11,12,13,14),]
beta_mean_inla = beta_inla[2:7,1]
beta_low_inla = beta_inla[2:7,3]
beta_up_inla = beta_inla[2:7,5]
loglike_inla = sum_inla$mlik[1]



################################################################
###          3. The table of results in Ecoregion 0          ###
################################################################
tab = matrix(NA,nrow=6,ncol=15)
tab[,1] = beta_mean_zdm
tab[,2] = beta_low_zdm
tab[,3] = beta_up_zdm
tab[,4] = beta_mean_zip
tab[,5] = beta_low_zip
tab[,6] = beta_up_zip
tab[,7] = beta_mean_inla
tab[,8] = beta_low_inla
tab[,9] = beta_up_inla
tab[,10] = beta_mean_dm
tab[,11] = beta_low_dm
tab[,12] = beta_up_dm
tab[,13] = beta_mean_poi
tab[,14] = beta_low_poi
tab[,15] = beta_up_poi
colnames(tab)=c("ZDM_mean","ZDM_low","ZDM_up","ZIP_mean","ZIP_low","ZIP_up",
         "ZBSJ_mean","ZBSJ_low","ZBSJ_up","DM_mean","DM_low","DM_up",
         "POR_mean","POR_low","POR_up")
tab = round(tab,3)
