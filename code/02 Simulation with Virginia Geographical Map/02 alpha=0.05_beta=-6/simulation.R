########################################################################
###   The ZIP model used for comparison needs the R package 'pscl' to 
###   complete the model fitting process, and it can be installed via 'install.packages("pscl")'
###   The DM and ZBSJ models used for comparison needs the R package 'INLA', 
###   and it can be installed via 'install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)'
###   The details of 'INLA' package can be found in 'https://www.r-inla.org/home'
#########################################################################
install.packages("pscl")
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(foreach) #library for parallel computing
library(doParallel) #library for parallel computing
library(parallel) #library for parallel computing
library(pscl) #library for fitting ZIP
library(INLA) #library for fitting ZBSJ
library(MASS)
library(ggplot2)

source("function.R")

#############################################################################
###                      1. Set some true values                          ###
#############################################################################
lyme <- read.csv("../../../data/lyme_data.csv",encoding='utf-8') #Lyme dataset
W <- read.csv("../../../data/W0.csv") #Adajacency matrix in Ecoregion 0 of Virginia
W = as.matrix(W[,-1])
I = dim(W)[1] #Number of regions in Ecoregion 0
eco = lyme[,1]
N = lyme$n[which(eco==0)] #Population
Y = lyme$disea[which(eco==0)] #Disease count
##Covarites
z = cbind(1,lyme[which(eco==0),c(5,6,10,12,15,17)]) 
z = as.matrix(z)
x = z
I = dim(z)[1]
M = dim(z)[2]
x_M = M

##Coefficients alpha and beta
data<-as.data.frame(cbind(Y,x[,2:7],z[,2:7]))
names(data)<-c('Y','X2','X3','X4','X5','X6','X7','Z2','Z3','Z4','Z5','Z6','Z7')
result_zip<-zeroinfl(Y~Z2+Z3+Z4+Z5+Z6+Z7|X2+X3+X4+X5+X6+X7,data = data)
sum_zip <- summary(result_zip)
beta_true = sum_zip$coefficients$count[,1]
alpha_true = sum_zip$coefficients$zero[,1]
list_inin = list(c(-6,5.8),c(-6,-0.05),c(-8.2,0.8),c(-8.2,-0.05))
list_in = list_inin[[2]]
beta_true[1] = list_in[1] #-8.2 -6
alpha_true[1] = list_in[2] #0.8 -0.05
beta_true = as.matrix(beta_true,ncol=1)
alpha_true = as.matrix(alpha_true,ncol=1)

##Spatial random effect
tau_true = 0.1
delta_true = 0.01
n=3000
set.seed(5)
phi_true = phigen(tau_true,W,n)
set.seed(5)
epsion_true = epsiongen(delta_true,W,n)

p_true = exp(x%*%alpha_true+epsion_true)/(1+exp(x%*%alpha_true+epsion_true))
q_true = exp(z%*%beta_true+phi_true)
mean(q_true)
lambda = N*q_true
q_bar = mean(q_true)
r_true = q_true/q_bar

#############################################################################
###                 2. Simulation for ZDM, ZIP and POR                    ###
###-----------------------------------------------------------------------###
###                       2.1. Parallel Computing                         ###
#############################################################################
s=100
nc = detectCores()
cl <- makeCluster(nc)
##Register a cluster for parallel computing
registerDoParallel(cl)
set.seed(5)
res_real <- foreach(j=1:s, .packages = c('pscl','MASS')) %dopar% {
  ##Prepare to strore the results of ZDM
  Rbeta_beta_new = list(length = M)
  for (i in 1:M){
    Rbeta_beta_new[[i]] = matrix(NA,nrow=1,ncol=5)
  }
  
  Ralpha_alpha_new = list(length = x_M)
  for (i in 1:x_M){
    Ralpha_alpha_new[[i]] = matrix(NA,nrow=1,ncol=5)
  }
  
  Rprob_prob_new = list(length = I)
  for (i in 1:I){
    Rprob_prob_new[[i]] = matrix(NA,nrow=1,ncol=6)
  }
  
  Rquan_quan_new = list(length = I)
  for (i in 1:I){
    Rquan_quan_new[[i]] = matrix(NA,nrow=1,ncol=6)
  }
  
  Rr_r_new = list(length = I)
  for(i in 1:I){
    Rr_r_new[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  ##Prepare to strore the results of ZIP
  Rbeta_beta_zip = list(length = M)
  for (i in 1:M){
    Rbeta_beta_zip[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  Ralpha_alpha_zip = list(length = x_M)
  for (i in 1:x_M){
    Ralpha_alpha_zip[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  Rprob_prob_zip = list(length = I)
  for (i in 1:I){
    Rprob_prob_zip[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  Rquan_quan_zip = list(length = I)
  for (i in 1:I){
    Rquan_quan_zip[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  Rr_r_zip = list(length = I)
  for(i in 1:I){
    Rr_r_zip[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  ##Prepare to strore the results of POR
  Rquan_quan_poi = list(length = I)
  for (i in 1:I){
    Rquan_quan_poi[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  Rbeta_beta_poi = list(length = M)
  for (i in 1:M){
    Rbeta_beta_poi[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
  Rr_r_poi = list(length = I)
  for(i in 1:I){
    Rr_r_poi[[i]] = matrix(NA,nrow=1,ncol=3)
  }
  
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
  Y_zero_rate = sum(ifelse(Y==0,1,0))
  
  ##Model fitting of ZDM
  result_new = model_new(I,Y,z,x,N,W,M,x_M)
  Ralpha = result_new[[1]]
  Rbeta = result_new[[2]]
  Rprob = result_new[[3]]
  Rquan = result_new[[4]]
  
  for (i in 1:M){
    Rbeta_beta_new[[i]][1:4] = Rbeta[i,]
    Rbeta_beta_new[[i]][5] = (Rbeta[i,1]-beta_true[i,])^2 #MSE
  }
  
  for (i in 1:x_M){
    Ralpha_alpha_new[[i]][1:4] = Ralpha[i,]
    Ralpha_alpha_new[[i]][5] = (Ralpha[i,1]-alpha_true[i,])^2 #MSE
  }
  
  for (i in 1:I){
    Rprob_prob_new[[i]][1:4] = Rprob[i,]
    Rprob_prob_new[[i]][5] = (Rprob[i,1]-p_true[i,])^2 #MSE
    Rprob_prob_new[[i]][6] = (log(Rprob[i,1])-log(p_true[i,]))^2
  }
  
  for (i in 1:I){
    Rquan_quan_new[[i]][1:4] = Rquan[i,]
    Rquan_quan_new[[i]][5] = (Rquan[i,1]-q_true[i,])^2 #MSE
    Rquan_quan_new[[i]][6] = (log(Rquan[i,1])-log(q_true[i,]))^2
  }
  
  for (i in 1:I){
    Rr_r_new[[i]][1] = Rquan[i,1]/mean(Rquan[,1])
    Rr_r_new[[i]][2] = (Rr_r_new[[i]][1]-r_true[i])^2
    Rr_r_new[[i]][3] = (log(Rr_r_new[[i]][1])-log(r_true[i,]))^2
  }
  rm(result_new,Rbeta,Ralpha,Rprob,Rquan)
  
  ##Model fitting of ZIP
  data<-as.data.frame(cbind(Y,x[,2:7],z[,2:7]))
  names(data)<-c('Y','X2','X3','X4','X5','X6','X7','Z2','Z3','Z4','Z5','Z6','Z7')
  result_zip<-zeroinfl(Y~Z2+Z3+Z4+Z5+Z6+Z7|X2+X3+X4+X5+X6+X7,data = data)
  sum_zip <- summary(result_zip)
  
  for (i in 1:M){
    Rbeta_beta_zip[[i]][1] = sum_zip$coefficients$count[i,1]
    Rbeta_beta_zip[[i]][2] = sum_zip$coefficients$count[i,2]
    Rbeta_beta_zip[[i]][3] = (Rbeta_beta_zip[[i]][1]-beta_true[i,])^2 #MSE
  }
  
  for (i in 1:x_M){
    Ralpha_alpha_zip[[i]][1] = sum_zip$coefficients$zero[i,1]
    Ralpha_alpha_zip[[i]][2] = sum_zip$coefficients$zero[i,2]
    Ralpha_alpha_zip[[i]][3] = (Ralpha_alpha_zip[[i]][1]-alpha_true[i,])^2 #MSE
  }
  
  alpha_zip = matrix(NA,nrow=x_M,ncol=1)
  for (i in 1:x_M){
    alpha_zip[i,] = Ralpha_alpha_zip[[i]][1]
  }
  prob_zip = exp(x%*%alpha_zip)/(1+exp(x%*%alpha_zip))
  for (i in 1:I){
    Rprob_prob_zip[[i]][1] = prob_zip[i,]
    Rprob_prob_zip[[i]][2] = (prob_zip[i,]-p_true[i,])^2 #MSE
    Rprob_prob_zip[[i]][3] = (log(prob_zip[i,])-log(p_true[i,]))^2
  }
  
  beta_zip = matrix(NA,nrow=M,ncol=1)
  for (i in 1:M){
    beta_zip[i,] = Rbeta_beta_zip[[i]][1]
  }
  quan_zip = exp(z%*%beta_zip)/N
  for (i in 1:I){
    Rquan_quan_zip[[i]][1] = quan_zip[i,]
    Rquan_quan_zip[[i]][2] = (quan_zip[i,]-q_true[i,])^2 #MSE
    Rquan_quan_zip[[i]][3] = (log(quan_zip[i,])-log(q_true[i,]))^2
  }
  
  for (i in 1:I){
    Rr_r_zip[[i]][1] = quan_zip[i]/mean(quan_zip)
    Rr_r_zip[[i]][2] = (Rr_r_zip[[i]][1]-r_true[i])^2
    Rr_r_zip[[i]][3] = (log(Rr_r_zip[[i]][1])-log(r_true[i,]))^2
  }
  rm(result_zip,sum_zip)
  
  ##Model fitting of POR
  sum_poi <- summary(glm(Y~z[,2]+z[,3]+z[,4]+z[,5]+z[,6]+z[,7],family=poisson(link="log")))
  for (i in 1:M){
    Rbeta_beta_poi[[i]][1] = sum_poi$coefficients[i,1]
    Rbeta_beta_poi[[i]][2] = sum_poi$coefficients[i,2]
    Rbeta_beta_poi[[i]][3] = (Rbeta_beta_poi[[i]][1]-beta_true[i,])^2
  }
  
  beta_poi = matrix(NA,nrow=M,ncol=1)
  for (i in 1:M){
    beta_poi[i,] = Rbeta_beta_poi[[i]][1]
  }
  quan_poi = exp(z%*%beta_poi)/N
  for (i in 1:I){
    Rquan_quan_poi[[i]][1] = quan_poi[i,]
    Rquan_quan_poi[[i]][2] = (quan_poi[i,]-q_true[i,])^2 #MSE
    Rquan_quan_poi[[i]][3] = (log(quan_poi[i,])-log(q_true[i,]))^2
  }
  
  for (i in 1:I){
    Rr_r_poi[[i]][1] = quan_poi[i]/mean(quan_poi)
    Rr_r_poi[[i]][2] = (Rr_r_poi[[i]][1]-r_true[i])^2
    Rr_r_poi[[i]][3] = (log(Rr_r_poi[[i]][1])-log(r_true[i,]))^2
  }
  rm(sum_poi)
  
  list(Rbeta_beta_new,
       Ralpha_alpha_new,
       Rquan_quan_new,
       Rprob_prob_new,
       Rr_r_new,
       Rbeta_beta_zip,
       Ralpha_alpha_zip,
       Rquan_quan_zip,
       Rprob_prob_zip,
       Rr_r_zip,
       Rbeta_beta_poi,
       Rquan_quan_poi,
       Rr_r_poi,
       Y)
}

##Close the cluster for parallel computing
stopImplicitCluster()
stopCluster(cl)

#############################################################################
###                      2.2. Summarize the results                       ###
#############################################################################
y = matrix(NA,nrow=s,ncol=I)

Rbeta_beta_new = list(length = M)
for (i in 1:M){
  Rbeta_beta_new[[i]] = matrix(NA,nrow=s,ncol=5)
}
for (i in 1:M){
  for (j in 1:s){
    Rbeta_beta_new[[i]][j,] = res_real[[j]][[1]][[i]]
  }
}

Ralpha_alpha_new = list(length = x_M)
for (i in 1:x_M){
  Ralpha_alpha_new[[i]] = matrix(NA,nrow=s,ncol=5)
}
for (i in 1:x_M){
  for (j in 1:s){
    Ralpha_alpha_new[[i]][j,] = res_real[[j]][[2]][[i]]
  }
}

Rprob_prob_new = list(length = I)
for (i in 1:I){
  Rprob_prob_new[[i]] = matrix(NA,nrow=s,ncol=6)
}
for (i in 1:I){
  for (j in 1:s){
    Rprob_prob_new[[i]][j,] = res_real[[j]][[4]][[i]]
  }
}

Rquan_quan_new = list(length = I)
for (i in 1:I){
  Rquan_quan_new[[i]] = matrix(NA,nrow=s,ncol=6)
}
for (i in 1:I){
  for (j in 1:s){
    Rquan_quan_new[[i]][j,] = res_real[[j]][[3]][[i]]
  }
}

Rr_r_new = list(length = I)
for(i in 1:I){
  Rr_r_new[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:I){
  for (j in 1:s){
    Rr_r_new[[i]][j,] = res_real[[j]][[5]][[i]]
  }
}

Rbeta_beta_zip = list(length = M)
for (i in 1:M){
  Rbeta_beta_zip[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:M){
  for (j in 1:s){
    Rbeta_beta_zip[[i]][j,] = res_real[[j]][[6]][[i]]
  }
}

Ralpha_alpha_zip = list(length = x_M)
for (i in 1:x_M){
  Ralpha_alpha_zip[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:M){
  for (j in 1:s){
    Ralpha_alpha_zip[[i]][j,] = res_real[[j]][[7]][[i]]
  }
}

Rprob_prob_zip = list(length = I)
for (i in 1:I){
  Rprob_prob_zip[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:I){
  for (j in 1:s){
    Rprob_prob_zip[[i]][j,] = res_real[[j]][[9]][[i]]
  }
}

Rquan_quan_zip = list(length = I)
for (i in 1:I){
  Rquan_quan_zip[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:I){
  for (j in 1:s){
    Rquan_quan_zip[[i]][j,] = res_real[[j]][[8]][[i]]
  }
}

Rr_r_zip = list(length = I)
for(i in 1:I){
  Rr_r_zip[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:I){
  for (j in 1:s){
    Rr_r_zip[[i]][j,] = res_real[[j]][[10]][[i]]
  }
}

Rquan_quan_poi = list(length = I)
for (i in 1:I){
  Rquan_quan_poi[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:I){
  for (j in 1:s){
    Rquan_quan_poi[[i]][j,] = res_real[[j]][[12]][[i]]
  }
}

Rbeta_beta_poi = list(length = M)
for (i in 1:M){
  Rbeta_beta_poi[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:M){
  for (j in 1:s){
    Rbeta_beta_poi[[i]][j,] = res_real[[j]][[11]][[i]]
  }
}

Rr_r_poi = list(length = I)
for(i in 1:I){
  Rr_r_poi[[i]] = matrix(NA,nrow=s,ncol=3)
}
for (i in 1:I){
  for (j in 1:s){
    Rr_r_poi[[i]][j,] = res_real[[j]][[13]][[i]]
  }
}

y = matrix(NA,nrow=s,ncol=I)
for (i in 1:s){
  y[i,] = res_real[[i]][[14]]
}

r_all<-matrix(nrow=I,ncol=6)
quan_all<-matrix(nrow = I,ncol = 6)
prob_all<-matrix(nrow = I,ncol = 4)
beta_all<-matrix(nrow = M,ncol = 9)
alpha_all<-matrix(nrow = x_M,ncol = 6)
for(i in 1:I){
  ##Relative risk
  ###MSE
  r_all[i,1] = mean(Rr_r_new[[i]][,2])
  r_all[i,2] = mean(Rr_r_zip[[i]][,2])
  r_all[i,3] = mean(Rr_r_poi[[i]][,2])
  ###LMSE
  r_all[i,4] = mean(Rr_r_new[[i]][,3])
  r_all[i,5] = mean(Rr_r_zip[[i]][,3])
  r_all[i,6] = mean(Rr_r_poi[[i]][,3])
  ##Disease incidence
  ###MSE
  quan_all[i,1] = mean(Rquan_quan_new[[i]][,5])
  quan_all[i,2] = mean(Rquan_quan_zip[[i]][,2])
  quan_all[i,3] = mean(Rquan_quan_poi[[i]][,2])
  ###LMSE
  quan_all[i,4] = mean(Rquan_quan_new[[i]][,6])
  quan_all[i,5] = mean(Rquan_quan_zip[[i]][,3])
  quan_all[i,6] = mean(Rquan_quan_poi[[i]][,3])
  ###MSE
  prob_all[i,1] = mean(Rprob_prob_new[[i]][,5])
  prob_all[i,2] = mean(Rprob_prob_zip[[i]][,2])
  ###LMSE
  prob_all[i,3] = mean(Rprob_prob_new[[i]][,6])
  prob_all[i,4] = mean(Rprob_prob_zip[[i]][,3])
}
for(i in 1:M){
  #MSE
  beta_all[i,1] = mean(Rbeta_beta_new[[i]][,5])
  beta_all[i,2] = mean(Rbeta_beta_zip[[i]][,3])
  beta_all[i,3] = mean(Rbeta_beta_poi[[i]][,3])
  #Bias^2
  beta_all[i,4] = (mean(Rbeta_beta_new[[i]][,1])-beta_true[i])^2
  beta_all[i,5] = (mean(Rbeta_beta_zip[[i]][,1])-beta_true[i])^2
  beta_all[i,6] = (mean(Rbeta_beta_poi[[i]][,1])-beta_true[i])^2
  #Var
  beta_all[i,7] = var(Rbeta_beta_new[[i]][,1])*(s-1)/s
  beta_all[i,8] = var(Rbeta_beta_zip[[i]][,1])*(s-1)/s
  beta_all[i,9] = var(Rbeta_beta_poi[[i]][,1])*(s-1)/s
}
for(i in 1:x_M){
  #MSE
  alpha_all[i,1] = mean(Ralpha_alpha_new[[i]][,5])
  alpha_all[i,2] = mean(Ralpha_alpha_zip[[i]][,3])
  #Bias^2
  alpha_all[i,3] = (mean(Ralpha_alpha_new[[i]][,1])-alpha_true[i])^2
  alpha_all[i,4] = (mean(Ralpha_alpha_zip[[i]][,1])-alpha_true[i])^2
  #Var
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



#############################################################################
###                  3. Simulation for ZBSJ and POR                       ###
#############################################################################
n1 = 1
##For use of INLA
g <- inla.read.graph(W)
s.index1 = seq(1,I,1)
s.index2 = c(seq(1,I,1),rep(NA,I))
s.index3 = c(rep(NA,I),seq(1,I,1))

s = 100

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
###                         3.1. Start simulation                         ###
#############################################################################
set.seed(7)
for(ii in 1:n1){
  for(j in 1:s){
    print(j)
    Y = y[j,]
    indi = vector(length=I)
    for(i in 1:I){
      indi[i] = ifelse(Y[i]==0,0,1)
    }
    Y_zero_rate = sum(ifelse(Y==0,1,0))
    
    ##Model fitting of DM
    dt1 = data.frame(Y=Y,z1=z[,2],z2=z[,3],z3=z[,4],z4=z[,5],z5=z[,6],z6=z[,7],N=N)
    formula2 <- Y ~ z1 + z2 + z3 + z4 + z5 + z6 + f(s.index1, model="besag", graph=g, hyper=list(prec=list(prior="loggamma", param=c(1, 1))),scale.model = TRUE)
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
    
    ##Model fitting of ZBSJ
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
      f(s.index2, model="bym2", graph=g, scale.model = TRUE, constr=TRUE, 
        hyper=list(phi=list(prior="pc",param=c(0.5,2/3),initial=3),prec=list(prior="pc.prec",param=c(1,0.01),initial=1.5)))+
      f(s.index3,copy='s.index2',fixed=FALSE)
    res3 <- inla(formula3, data=dt1, 
                 family=c("binomial","poisson"),E=E) #zeroinflatedpoisson1
    sum_inla = summary(res3)
    
    Rbeta = sum_inla$fixed[c(2,seq(9,14,1)),]
    Ralpha = sum_inla$fixed[c(1,seq(3,8,1)),]
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
    prob_inla = exp(x%*%alpha_inla+res3$summary.random$s.index2$mean[1:583])/(1+exp(x%*%alpha_inla+res3$summary.random$s.index2$mean[1:583]))
    for (i in 1:I){
      Rprob_prob_inla[[i]][j,1] = prob_inla[i,]
      Rprob_prob_inla[[i]][j,2] = (prob_inla[i,]-p_true[i,])^2 #MSE
      Rprob_prob_inla[[i]][j,3] = (log(prob_inla[i,])-log(p_true[i,]))^2
    }
    
    beta_inla = matrix(NA,nrow=M,ncol=1)
    for (i in 1:M){
      beta_inla[i,] = Rbeta_beta_inla[[i]][j,1]
    }
    r_inla = exp(z%*%beta_inla+res3$summary.random$s.index3$mean[584:1166])
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
    ##Disease incidence
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
    #MSE
    beta_all1[i,1] = mean(Rbeta_beta_old[[i]][,5])
    beta_all1[i,2] = mean(Rbeta_beta_inla[[i]][,5])
    #Bias^2
    beta_all1[i,3] = (mean(Rbeta_beta_old[[i]][,1])-beta_true[i])^2
    beta_all1[i,4] = (mean(Rbeta_beta_inla[[i]][,1])-beta_true[i])^2
    #Var
    beta_all1[i,5] = var(Rbeta_beta_old[[i]][,1])*(s-1)/s
    beta_all1[i,6] = var(Rbeta_beta_inla[[i]][,1])*(s-1)/s
  }
  for(i in 1:x_M){
    #MSE
    alpha_all1[i,1] = mean(Ralpha_alpha_inla[[i]][,5])
    #Bias^2
    alpha_all1[i,2] = (mean(Ralpha_alpha_inla[[i]][,1])-alpha_true[i])^2
    #Var
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
###                    3.2. Summarize the results                         ###
#############################################################################
beta_integ_all = matrix(NA,nrow=n1,ncol=15)
quan_integ_all = matrix(NA,nrow=n1,ncol=10)
r_integ_all = matrix(NA,nrow=n1,ncol=10)

for(i in 1:n1){
  beta_integ_all[i,1:3] = beta_integ[1:3]
  beta_integ_all[i,4:5] = beta_integ1[i,1:2]
  beta_integ_all[i,6:8] = beta_integ[4:6]
  beta_integ_all[i,9:10] = beta_integ1[i,3:4]
  beta_integ_all[i,11:13] = beta_integ[7:9]
  beta_integ_all[i,14:15] = beta_integ1[i,5:6]
  
  quan_integ_all[i,1:3] = quan_integ[1:3]
  quan_integ_all[i,4:5] = quan_integ1[i,1:2]
  quan_integ_all[i,6:8] = quan_integ[4:6]
  quan_integ_all[i,9:10] = quan_integ1[i,3:4]
  
  r_integ_all[i,1:3] = r_integ[1:3]
  r_integ_all[i,4:5] = r_integ1[i,1:2]
  r_integ_all[i,6:8] = r_integ[4:6]
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
tab = matrix(NA,nrow=5,ncol=4)
tab[1,1] = quan_integ_all[i,1]
tab[1,2] = quan_integ_all[i,6]
tab[1,3] = r_integ_all[i,1]
tab[1,4] = r_integ_all[i,6]
tab[2,1] = quan_integ_all[i,2]
tab[2,2] = quan_integ_all[i,7]
tab[2,3] = r_integ_all[i,2]
tab[2,4] = r_integ_all[i,7]
tab[3,1] = quan_integ_all[i,5]
tab[3,2] = quan_integ_all[i,10]
tab[3,3] = r_integ_all[i,5]
tab[3,4] = r_integ_all[i,10]
tab[4,1] = quan_integ_all[i,4]
tab[4,2] = quan_integ_all[i,9]
tab[4,3] = r_integ_all[i,4]
tab[4,4] = r_integ_all[i,9]
tab[5,1] = quan_integ_all[i,3]
tab[5,2] = quan_integ_all[i,8]
tab[5,3] = r_integ_all[i,3]
tab[5,4] = r_integ_all[i,8]
rownames(tab) = c("ZDM","ZIP","ZBSJ","DM","POR")
colnames(tab) = c("rate_MSE","rate_LMSE","risk_MSE","risk_LMSE")



#############################################################################
###          5. The resulting figures for the coefficient beta            ###
#############################################################################
setwd("../../../results/02 Simulation with Virginia Geographical Map")
p <- beta_integ_all[6:15]
a = 0.1
b = a/2
data.1<-data.frame(PE=c(-min(p[1],a),min(p[6],a), -min(p[2],a),min(p[7],a),-min(p[5],a),min(p[10],a),-min(p[4],a),min(p[9],a),-min(p[3],a),min(p[8],a)),
                   method=factor(rep(c("ZDM", "ZIP","ZBSJ", "DM","POR"), each = 2), 
                                 levels = c("POR","DM","ZBSJ","ZIP", "ZDM")),
                   type=factor(c(8,7,10,9,4,3,2,1,6,5)))
t<-expression(paste(rho,"=",0))
p1<-ggplot(data=data.1, aes(x=method, y=PE, fill=type))+geom_bar(stat = "identity",position = "identity",width = 0.75)+theme(axis.text.x=element_text(angle=90,colour="black"))+coord_flip()
p1<-p1+scale_fill_brewer(palette = "Paired")+
  theme_bw()+
  scale_y_continuous(limits=c(-a,a), breaks=c(seq(-a,a,b)),
                     labels =c("0.1+", as.character(c(seq(a-b,b,-b),0,seq(b,a-b,b))), "0.1+") )+
  labs(x="",y="")+
  guides(fill=F)+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t=10,r=10,b=0,l=0))
ggsave("alpha=0.05_beta=-6.png",width=4,height=3)

