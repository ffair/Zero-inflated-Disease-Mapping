########################################################################
###   The ZIP model used for comparison needs the R package 'pscl' to 
###   complete the model fitting process, and it can be installed via 'install.packages("pscl")'
###   The DM and ZBSJ models used for comparison needs the R package 'INLA', 
###   and it can be installed via 'install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)'
###   The details of 'INLA' package can be found in 'https://www.r-inla.org/home'
#########################################################################
install.packages("pscl")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA) #library for fitting ZBSJ and DM
library(pscl) #library for fitting ZIP
library(MASS)
library(parallel) #library for parallel computing
library(ggplot2)

source("function.R")

#############################################################################
###       1. Simulation for ZDM, ZIP, POR using prarllel computing        ###
#############################################################################
nc = detectCores()
clus <- makeCluster(nc)
clusterExport(clus,c("model_new","update_delta","update_tau","update_alpha","update_beta","post_alpha1","post_beta1","update_epsion","update_phi","post_epsion1","post_phi1","Y1_likelihood","Y0_likelihood","phigen","epsiongen"),envir = environment())
clusterEvalQ(clus,{
  library(MASS)
  library(pscl)})
list_in = list(c(0,-6.5,0.95),c(0,-6.5,0.1),c(0,-7.5,0.95),c(0,-7.5,0.1),
               c(0.3,-6.5,0.95),c(0.3,-6.5,0.1),c(0.3,-7.5,0.95),c(0.3,-7.5,0.1),
               c(0.6,-6.5,0.95),c(0.6,-6.5,0.1),c(0.6,-7.5,0.95),c(0.6,-7.5,0.1))
res_art = parLapply(clus, list_in, fun = simul_art)


#############################################################################
###         2. Simulation for DM, ZBSJ using the same true values         ###
###-----------------------------------------------------------------------###
###                      2.1. Generate true values                        ###
#############################################################################
n1 = 12 #Number of simulation settings
W <- read.csv('../../../data/W.csv',header=T)
W = as.matrix(W)
I = dim(W)[1]
##Population
N_mean = 50000
N_sd = 10000
set.seed(5)
N=rnorm(I,mean=N_mean,sd=N_sd)  
N=round(N,0)
##Spatial random effect
tau_true = 0.1
delta_true = 0.01
n=2000
set.seed(5)
phi_true = phigen(tau_true,W,n)
set.seed(5)
epsion_true = epsiongen(delta_true,W,n)
##Covariates
M=5
x_M = 2
x=matrix(NA,I,x_M)  
x[,1]=rep(1,I)
set.seed(5)
x[,2]=rnorm(I,0,1)
##For use of INLA library
g <- inla.read.graph(W)
s.index1 = seq(1,I,1)
s.index2 = c(seq(1,I,1),rep(NA,I))
s.index3 = c(rep(NA,I),seq(1,I,1))

#############################################################################
###                  2.2. Prepare to start the simulation                 ###
#############################################################################
s = 100 

##Store the results of DM model
Rbeta_beta_old = list(length = M)
for (i in 1:M){
  Rbeta_beta_old[[i]] = matrix(NA,nrow=s,ncol=5)
}

Rquan_quan_old = list(length = I)
for (i in 1:I){
  Rquan_quan_old[[i]] = matrix(NA,nrow=s,ncol=3)
}

Rr_r_old = list(length = I)
for(i in 1:I){
  Rr_r_old[[i]] = matrix(NA,nrow=s,ncol=3)
}

##Store the results of the ZBSJ model
Rbeta_beta_inla = list(length = M)
for (i in 1:M){
  Rbeta_beta_inla[[i]] = matrix(NA,nrow=s,ncol=5)
}

Ralpha_alpha_inla = list(length = x_M)
for (i in 1:x_M){
  Ralpha_alpha_inla[[i]] = matrix(NA,nrow=s,ncol=5)
}

Rprob_prob_inla = list(length = I)
for (i in 1:I){
  Rprob_prob_inla[[i]] = matrix(NA,nrow=s,ncol=3)
}

Rquan_quan_inla = list(length = I)
for (i in 1:I){
  Rquan_quan_inla[[i]] = matrix(NA,nrow=s,ncol=3)
}

Rr_r_inla = list(length = I)
for(i in 1:I){
  Rr_r_inla[[i]] = matrix(NA,nrow=s,ncol=3)
}

beta_integ1 = matrix(NA,nrow=n1,ncol=6)
quan_integ1 = matrix(NA,nrow=n1,ncol=4)
r_integ1 = matrix(NA,nrow=n1,ncol=4)

#############################################################################
###                       2.2. Start the simulation                       ###
#############################################################################
set.seed(5)
for(ii in 1:n1){
  z_mean = c(0,0,0,0)
  z_sigma = matrix(list_in[[ii]][1],nrow=4,ncol=4,byrow=TRUE)
  diag(z_sigma) = 1
  set.seed(5)
  z = mvrnorm(I,z_mean,z_sigma)
  z<-cbind(rep(1,I),z)
  beta_true = as.matrix(c(list_in[[ii]][2],0.1,-0.2,0.3,-0.4),ncol=1)
  alpha_true = as.matrix(c(list_in[[ii]][3],0.5),ncol=1)
  tau_true = 0.1
  delta_true = 0.01
  p_true = exp(x%*%alpha_true+epsion_true)/(1+exp(x%*%alpha_true+epsion_true))
  q_true = exp(z%*%beta_true+phi_true)
  lambda = N*q_true
  q_bar = mean(q_true)
  r_true = q_true/q_bar
  
  for(j in 1:s){
    print(j)
    Y=res_art[[ii]][[4]][j,]
    indi = vector(length=I)
    for(i in 1:I){
      indi[i] = ifelse(Y[i]==0,0,1)
    }

    ##DM
    dt1 = data.frame(Y=Y,z1=z[,2],z2=z[,3],z3=z[,4],z4=z[,5],N=N)
    formula2 <- Y ~ z1 + z2 + z3 + z4 + f(s.index1, model="besag", 
                                          graph=g, 
                                          hyper=list(prec=list(prior="loggamma", 
                                                               param=c(1, 1))),
                                          scale.model = TRUE)
    res2 <- inla(formula2, data=dt1, family="poisson", E=N)
    sum_dm = summary(res2)
    Rbeta = sum_dm$fixed
    
    for (i in 1:M){
      Rbeta_beta_old[[i]][j,1:4] = Rbeta[i,c(1,2,3,5)]
      Rbeta_beta_old[[i]][j,5] = (Rbeta[i,1]-beta_true[i,])^2 #MSE
    }
    beta_old = matrix(NA,nrow=M,ncol=1)
    for (i in 1:M){
      beta_old[i,] = Rbeta_beta_old[[i]][j,1]
    }
    quan_old = exp(z%*%beta_old+res2$summary.random$s.index1$mean)
    for (i in 1:I){
      Rquan_quan_old[[i]][j,1] = quan_old[i,]
      Rquan_quan_old[[i]][j,2] = (quan_old[i,]-q_true[i,])^2 #MSE
      Rquan_quan_old[[i]][j,3] = (log(quan_old[i,])-log(q_true[i,]))^2
    }
    for (i in 1:I){
      Rr_r_old[[i]][j,1] = quan_old[i]/mean(quan_old)
      Rr_r_old[[i]][j,2] = (Rr_r_old[[i]][j,1]-r_true[i])^2
      Rr_r_old[[i]][j,3] = (log(Rr_r_old[[i]][j,1])-log(r_true[i,]))^2
    }
    rm(Rbeta,dt1,res2,formula2)
    
    ##ZBSJ
    YY = matrix(NA,2*I,2)
    YY[1:I,1] = indi
    YY[(I+1):dim(YY)[1],2] = Y
    dt1 = list(
      YY = YY,
      mu.x = c(rep(1,I),rep(NA,I)),
      mu.z = c(rep(NA,I),rep(1,I)),
      cov.x = c(x[,2],rep(NA,I)),
      cov.z1 = c(rep(NA,I),z[,2]),
      cov.z2 = c(rep(NA,I),z[,3]),
      cov.z3 = c(rep(NA,I),z[,4]),
      cov.z4 = c(rep(NA,I),z[,5]),
      N = c(rep(NA,I),N)
    )
    pp_bar = sum(Y)/sum(N)
    E = N*pp_bar
    dt1$E = c(rep(NA,I),E)
    formula3 <- YY ~ 0 + mu.x + mu.z + cov.x + cov.z1 + cov.z2 + cov.z3 + cov.z4 +
      f(s.index2, model="bym2", graph=g, scale.model = TRUE, constr=TRUE, 
        hyper=list(phi=list(prior="pc",param=c(0.5,2/3),initial=3),prec=list(prior="pc.prec",param=c(1,0.01),initial=1.5)))+
      f(s.index3,copy='s.index2',fixed=FALSE)
    
    res3 <- inla(formula3, data=dt1, 
                 family=c("binomial","poisson"),E=E,verbose=F,
                 control.predictor=list(compute=TRUE),
                 control.compute=list(dic=TRUE,cpo=TRUE)) #zeroinflatedpoisson1
    sum_inla = summary(res3)
    
    Rbeta = sum_inla$fixed[c(2,seq(4,7,1)),]
    Ralpha = sum_inla$fixed[c(1,3),]
    for (i in 1:M){
      Rbeta_beta_inla[[i]][j,1:4] = Rbeta[i,c(1,2,3,5)]
      Rbeta_beta_inla[[i]][j,5] = (Rbeta[i,1]-beta_true[i,])^2 #MSE
    }
    
    for (i in 1:x_M){
      Ralpha_alpha_inla[[i]][j,1:4] = Ralpha[i,c(1,2,3,5)]
      Ralpha_alpha_inla[[i]][j,5] = (Ralpha[i,1]-alpha_true[i,])^2 #MSE
    }
    
    alpha_inla = matrix(NA,nrow=x_M,ncol=1)
    for (i in 1:x_M){
      alpha_inla[i,] = Ralpha_alpha_inla[[i]][j,1]
    }
    prob_inla = exp(x%*%alpha_inla+res3$summary.random$s.index2$mean[1:100])/(1+exp(x%*%alpha_inla+res3$summary.random$s.index2$mean[1:100]))
    for (i in 1:I){
      Rprob_prob_inla[[i]][j,1] = prob_inla[i,]
      Rprob_prob_inla[[i]][j,2] = (prob_inla[i,]-p_true[i,])^2 #MSE
      Rprob_prob_inla[[i]][j,3] = (log(prob_inla[i,])-log(p_true[i,]))^2
    }
    
    beta_inla = matrix(NA,nrow=M,ncol=1)
    for (i in 1:M){
      beta_inla[i,] = Rbeta_beta_inla[[i]][j,1]
    }
    r_inla = exp(z%*%beta_inla+res3$summary.random$s.index3$mean[101:200])
    for (i in 1:I){
      Rr_r_inla[[i]][j,1] = r_inla[i,]
      Rr_r_inla[[i]][j,2] = (r_inla[i,]-r_true[i,])^2 #MSE
      Rr_r_inla[[i]][j,3] = (log(r_inla[i,])-log(r_true[i,]))^2
    }
    for (i in 1:I){
      Rquan_quan_inla[[i]][j,1] = r_inla[i]*pp_bar
      Rquan_quan_inla[[i]][j,2] = (Rquan_quan_inla[[i]][j,1]-q_true[i])^2
      Rquan_quan_inla[[i]][j,3] = (log(Rquan_quan_inla[[i]][j,1])-log(q_true[i,]))^2
    }
    rm(Rbeta,Ralpha,YY,dt1,res3,formula3)
  }
  
  
  #############################################################################
  ###                       3. Summarize the results                        ###
  ###-----------------------------------------------------------------------###
  ###             3.1. Summarize the results for DM, ZBSJ model             ###
  #############################################################################
  r_all1<-matrix(nrow=I,ncol=4)
  quan_all1<-matrix(nrow = I,ncol = 4)
  prob_all1<-matrix(nrow = I,ncol = 2)
  beta_all1<-matrix(nrow = M,ncol = 6)
  alpha_all1<-matrix(nrow = x_M,ncol = 3)
  for(i in 1:I){
    ##Relative risk
    ###MSE
    r_all1[i,1] = mean(Rr_r_old[[i]][,2])
    r_all1[i,2] = mean(Rr_r_inla[[i]][,2])
    ###LMSE
    r_all1[i,3] = mean(Rr_r_old[[i]][,3])
    r_all1[i,4] = mean(Rr_r_inla[[i]][,3])
    ##Incidence rate
    ###MSE
    quan_all1[i,1] = mean(Rquan_quan_old[[i]][,2])
    quan_all1[i,2] = mean(Rquan_quan_inla[[i]][,2])
    ###LMSE
    quan_all1[i,3] = mean(Rquan_quan_old[[i]][,3])
    quan_all1[i,4] = mean(Rquan_quan_inla[[i]][,3])
    ##Probability of zero count process
    ###MSE
    prob_all1[i,1] = mean(Rprob_prob_inla[[i]][,2])
    ###LMSE
    prob_all1[i,2] = mean(Rprob_prob_inla[[i]][,3])
  }
  for(i in 1:M){
    ###MSE
    beta_all1[i,1] = mean(Rbeta_beta_old[[i]][,5])
    beta_all1[i,2] = mean(Rbeta_beta_inla[[i]][,5])
    ###Bias^2
    beta_all1[i,3] = (mean(Rbeta_beta_old[[i]][,1])-beta_true[i])^2
    beta_all1[i,4] = (mean(Rbeta_beta_inla[[i]][,1])-beta_true[i])^2
    ###Var
    beta_all1[i,5] = var(Rbeta_beta_old[[i]][,1])*(s-1)/s
    beta_all1[i,6] = var(Rbeta_beta_inla[[i]][,1])*(s-1)/s
  }
  for(i in 1:x_M){
    ###MSE
    alpha_all1[i,1] = mean(Ralpha_alpha_inla[[i]][,5])
    ###Bias^2
    alpha_all1[i,2] = (mean(Ralpha_alpha_inla[[i]][,1])-alpha_true[i])^2
    ###Var
    alpha_all1[i,3] = var(Ralpha_alpha_inla[[i]][,1])*(s-1)/s
  }
  colnames(quan_all1) = c('old_MSE','inla_MSE','old_log','inla_log')
  colnames(prob_all1) = c('inla_MSE','inla_log')
  colnames(alpha_all1) = c('inla_MSE','inla_sqBias','inla_Var')
  colnames(beta_all1) = c('old_MSE','inla_MSE','old_sqBias','inla_sqBias','old_Var','inla_Var')
  colnames(r_all1) = c('old_MSE','inla_MSE','old_log','inla_log')
  
  beta_integ1[ii,] = apply(beta_all1[-1,],2,mean)
  quan_integ1[ii,] = apply(quan_all1,2,mean)
  r_integ1[ii,] = apply(r_all1,2,mean)
}

#############################################################################
###                3.2. Summarize the results of 5 models                 ###
#############################################################################
beta_integ_all = matrix(NA,nrow=n1,ncol=15)
quan_integ_all = matrix(NA,nrow=n1,ncol=10)
r_integ_all = matrix(NA,nrow=n1,ncol=10)

for(i in 1:n1){
  beta_integ_all[i,1:3] = res_art[[i]][[1]][1:3]
  beta_integ_all[i,4:5] = beta_integ1[i,1:2]
  beta_integ_all[i,6:8] = res_art[[i]][[1]][4:6]
  beta_integ_all[i,9:10] = beta_integ1[i,3:4]
  beta_integ_all[i,11:13] = res_art[[i]][[1]][7:9]
  beta_integ_all[i,14:15] = beta_integ1[i,5:6]
  
  quan_integ_all[i,1:3] = res_art[[i]][[2]][1:3]
  quan_integ_all[i,4:5] = quan_integ1[i,1:2]
  quan_integ_all[i,6:8] = res_art[[i]][[2]][4:6]
  quan_integ_all[i,9:10] = quan_integ1[i,3:4]
  
  r_integ_all[i,1:3] = res_art[[i]][[3]][1:3]
  r_integ_all[i,4:5] = r_integ1[i,1:2]
  r_integ_all[i,6:8] = res_art[[i]][[3]][4:6]
  r_integ_all[i,9:10] = r_integ1[i,3:4]
}
colnames(quan_integ_all) = c('new_MSE','zip_MSE','poi_MSE','old_MSE','inla_MSE','new_log','zip_log','poi_log','old_log','inla_log')
colnames(beta_integ_all) = c('new_MSE','zip_MSE','poi_MSE','old_MSE','inla_MSE','new_sqBias','zip_sqBias','poi_sqBias','old_sqBias','inla_sqBias','new_Var','zip_Var','poi_Var','old_Var','inla_Var')
colnames(r_integ_all) = c('new_MSE','zip_MSE','poi_MSE','old_MSE','inla_MSE','new_log','zip_log','poi_log','old_log','inla_log')


#############################################################################
###          4. The tables for incidence rate and relative risk           ###
###-----------------------------------------------------------------------###
###  "tab": a list stores the result tables, with each dimension stores a ###
###         kind of simulation setting.                                   ###
#############################################################################
tab = list(length=12)
for(i in 1:12){
  tab[[i]] = matrix(NA,nrow=5,ncol=4)
  tab[[i]][1,1] = quan_integ_all[i,1]
  tab[[i]][1,2] = quan_integ_all[i,6]
  tab[[i]][1,3] = r_integ_all[i,1]
  tab[[i]][1,4] = r_integ_all[i,6]
  tab[[i]][2,1] = quan_integ_all[i,2]
  tab[[i]][2,2] = quan_integ_all[i,7]
  tab[[i]][2,3] = r_integ_all[i,2]
  tab[[i]][2,4] = r_integ_all[i,7]
  tab[[i]][3,1] = quan_integ_all[i,5]
  tab[[i]][3,2] = quan_integ_all[i,10]
  tab[[i]][3,3] = r_integ_all[i,5]
  tab[[i]][3,4] = r_integ_all[i,10]
  tab[[i]][4,1] = quan_integ_all[i,4]
  tab[[i]][4,2] = quan_integ_all[i,9]
  tab[[i]][4,3] = r_integ_all[i,4]
  tab[[i]][4,4] = r_integ_all[i,9]
  tab[[i]][5,1] = quan_integ_all[i,3]
  tab[[i]][5,2] = quan_integ_all[i,8]
  tab[[i]][5,3] = r_integ_all[i,3]
  tab[[i]][5,4] = r_integ_all[i,8]
  rownames(tab[[i]]) = c("ZDM","ZIP","ZBSJ","DM","POR")
  colnames(tab[[i]]) = c("rate_MSE","rate_LMSE","risk_MSE","risk_LMSE")
}


#############################################################################
###          5. The resulting figures for the coefficient beta            ###
#############################################################################
setwd("../../../results/01 Simulation with Artificial Geographical Map/02 Second Situation")
plot_beta(1)
ggsave("rho=0_alpha=0.95_beta=-6.5.png",width=4,height=3)
plot_beta(2)
ggsave("rho=0_alpha=0.1_beta=-6.5.png",width=4,height=3)
plot_beta(3)
ggsave("rho=0_alpha=0.95_beta=-7.5.png",width=4,height=3)
plot_beta(4)
ggsave("rho=0_alpha=0.1_beta=-7.5.png",width=4,height=3)
plot_beta(5)
ggsave("rho=3_alpha=0.95_beta=-6.5.png",width=4,height=3)
plot_beta(6)
ggsave("rho=3_alpha=0.1_beta=-6.5.png",width=4,height=3)
plot_beta(7)
ggsave("rho=3_alpha=0.95_beta=-7.5.png",width=4,height=3)
plot_beta(8)
ggsave("rho=3_alpha=0.1_beta=-7.5.png",width=4,height=3)
plot_beta(9)
ggsave("rho=6_alpha=0.95_beta=-6.5.png",width=4,height=3)
plot_beta(10)
ggsave("rho=6_alpha=0.1_beta=-6.5.png",width=4,height=3)
plot_beta(11)
ggsave("rho=6_alpha=0.95_beta=-7.5.png",width=4,height=3)
plot_beta(12)
ggsave("rho=6_alpha=0.6_beta=-7.5.png",width=4,height=3)







