setwd("/Users/wendy/Dropbox/Dynamic expectile factor model/Code")
rm(list = ls())
graphics.off()
libraries = c("svd","fda","vars", "stats", "tseries",  "dygraphs", "CCA", "ggplot2","quantreg","mboost", "BayesX")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})

#install.packages("~/Dropbox/My Mac (王冰菱’s MacBook Pro (2))/Downloads/expectreg_0.39.tar", repos = NULL,
                 #type = "source")
#"qgraph"
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

library(parallel)
library(expectreg)
source("fpca.R", encoding="UTF-8")
source("PrincipalDirection_modified.R")

N_seq=c(100,100,100,100,1000,1000,1000,1000)
T_seq=c(100,200,500,1000,1000,500,200,100)
#tau_seq = c(0.05,0.5,0.95)
re_list =c()
#for(j in 1:length(tau_seq)){
for (q in 1:8){
nsimu = 10
MSE_model = rep(0,nsimu)
MSE_f = rep(0,nsimu)
MSE_l = rep(0,nsimu)
MSE_fo = rep(0,nsimu)
MSE_lo = rep(0,nsimu)
Rsquared_f1 = rep(0,nsimu)
Rsquared_f2 = rep(0,nsimu)
Rsquared_l1 = rep(0,nsimu)
Rsquared_l2 = rep(0,nsimu)

tau = 0.05
alpha = tau-0.5
#tau = tau_seq[j]
#N=100
#T=200
N = N_seq[2]
T = T_seq[2]
### start the loop of the simulation ###
for (simu_i in 1: nsimu)
{
  cat("Simu:",simu_i)
  ### simulate data matrix 
  ### simulate error term
  e = array(0,c(N, T))
  sigma = 0.1
  ### where lth column of e is an i.i.d. sample from N(0, sigma), for l from 1 to N.
  for (i in  1:T)
  {
    one_col_e = rnorm(N, 0, sigma)
    #one_col_e = rt(N, 5)
    e[,i] = one_col_e
  }
  e = as.matrix(e)
  
  #set.seed(123)
  #Just generate random AR(1) time series; based on this, I want to estimate the parameters
  t_seq = seq(0,1,length.out = T)
  mu = 1+ t_seq + exp(-(t_seq-0.6)^2/0.05)
  f_1 = arima.sim(n=T, list(ar=c(0.5)))
  f_2 = arima.sim(n=T, list(ar=c(0.3,0.2))) 
  l_1 = rnorm(N, 0, 3)
  l_2 = rnorm(N, 0, 2)
 
  lt_1 = as.matrix(sign(l_1[1])*l_1)
  lt_2 = as.matrix(sign(l_2[1])*l_2)
  signf1 = ifelse(sign(lt_1[1])==sign(f_1[1]),1,-1)
  signf2 = ifelse(sign(lt_2[1])==sign(f_2[1]),1,-1)
  ft_1 = signf1*f_1
  ft_2 = signf2*f_2
  
  e_et = rep(0, T)
  for (j in 1:T){
    e_et[j] = expectile(e[,j],alpha)
  }
  
  return_in = lt_1%*%ft_1 + lt_2%*%ft_2 +e
  
  e_returnt =  rep(0, T)
  for (j in 1:T){
    e_returnt[j]=expectile(return_in[,j],alpha)
  }
  
  ### calculate the expectile of errors and returns
  e_returnn =  rep(0, N)
  e_en = rep(0, N)
  for (i in 1:N){
    e_en[i] = expectile(e[i,],alpha)
    e_returnn[i] = expectile(return_in[i,],alpha)
  }
 
  return_s = as.matrix(return_in-e_returnt)
  
  
  #### pec case
  return=return_s
  
  set.seed(598)
  output.pec0 =	pec.k(return,nk=2,alpha=-0.45)
  output.pec = lapply(alpha,pec.k,Y=return, nk=2,reset.tol=50,lab.ini=-output.pec0[[4]])
  basis.pec = lapply(output.pec,listBasis,basis0=output.pec[[1]][[2]],num=2)
  lambdas = output.pec[[1]][[5]]
  PCs = matrix(unlist(basis.pec), N, 2)
  loadings_n= matrix(0,N,2)
  loadings_n[,1] = PCs[,1]/sqrt(lambdas[1])
  loadings_n[,2] = PCs[,2]/sqrt(lambdas[2])
  scores_pec = rbind(basis.pec[[1]][,1]%*%return, basis.pec[[1]][,2]%*%(return-apply(basis.pec[[1]][,1]%*%return, 2, "*",basis.pec[[1]][,1])))
  
  F1 = as.matrix(sign(scores_pec[1,])*scores_pec[1,])
  F2 = as.matrix(sign(scores_pec[2,])*scores_pec[2,])
  factors = cbind(F1/sqrt(lambdas[1]), F2/sqrt(lambdas[2]))
  
  loadings = matrix(0, N, 2)
  inter = rep(0, N)
  #pred = matrix(0, nrow(return_in), ncol(return_in))
  for (i in 1:N){
    fit = lm(as.numeric(return[i,])~factors[,1]+factors[,2])
    loadings[i,] = fit$coefficients[-1]
    inter[i] = as.numeric(fit$coefficients[1])
    #pred[i,] = as.numeric(unlist(predict(fit)[1]))
    #newfactor[i,] = pre[1]
  }
  
  signl1 = ifelse(sign(loadings[1,1])==sign(factors[1,1]),1,-1)
  signl2 = ifelse(sign(loadings[1,2])==sign(factors[1,2]),1,-1)
  loadings[,1] = signl1*loadings[,1]
  loadings[,2] = signl2*loadings[,2]
  
  plot(factors[,1],type='l')
  lines(ft_1,col='blue')
  plot(factors[,2],type ='l')
  lines(ft_2+1,col='blue')
  plot(loadings[,1],type='l')
  lines(lt_1,col='red')
  plot(loadings[,2],type='l')
  lines(lt_2,col='red')

  
  ### identification, make sure the same loading and corresponding factor have the same sign
  #est_model = loadings_n%*%scores_pec
  est_model = loadings%*%t(factors)+inter
  F1 = factors[,1]
  F2 = factors[,2]
  #ft_1 = ft_1+e_return
  #ft_2 = ft_2+e_return
  lto = cbind(lt_1, lt_2)
  fto = cbind(ft_1, ft_2)
  factorso = cbind(F1, F2)
  loadingso= loadings
  MSE_lo[simu_i] = sqrt(mean((lto-loadingso)^2))
  MSE_fo[simu_i] = sqrt(mean((fto-factorso)^2))
  
  ft_1 = (ft_1 - min(ft_1))/(max(ft_1) - min(ft_1))
  F1 = (F1 - min(F1))/(max(F1) - min(F1))
  ft_2 = (ft_2 - min(ft_2))/(max(ft_2) - min(ft_2))
  F2 = (F2 - min(F2))/(max(F2) - min(F2))
  lt_1 = (lt_1 - min(lt_1))/(max(lt_1) - min(lt_1))
  loadings[,1] = (loadings[,1] - min(loadings[,1]))/(max(loadings[,1]) - min(loadings[,1]))
  lt_2 = (lt_2 - min(lt_2))/(max(lt_2) - min(lt_2))
  loadings[,2] = (loadings[,2] - min(loadings[,2]))/(max(loadings[,2]) - min(loadings[,2]))
  lt = cbind(lt_1, lt_2)
  ft = cbind(ft_1, ft_2)
  factors = cbind(F1, F2)
  
  MSE_l[simu_i] = sqrt(mean((lt-loadings)^2))
  MSE_f[simu_i] = sqrt(mean((ft-factors)^2))
  MSE_model[simu_i] = sqrt(mean((return-est_model)^2))
  
  #MSE_loading[simu_i]= sqrt(mean((lt_1-loadings[1,])^2+(lt_2-loadings[2,])^2))
  #MSE_factor[simu_i] = sqrt(mean(cbind((ft_1-F1)^2,(ft_2-F2)^2)))
  
  Rsquared_f1[simu_i] = summary(lm(ft_1~F1))$adj.r.squared
  Rsquared_f2[simu_i] = summary(lm(ft_2~F2))$adj.r.squared
  Rsquared_l1[simu_i] = summary(lm(lt_1~loadings[,1]))$adj.r.squared
  Rsquared_l2[simu_i] = summary(lm(lt_2~loadings[,2]))$adj.r.squared
}
Result = cbind(Rsquared_f1,Rsquared_f2,Rsquared_l1,Rsquared_l2,MSE_l, MSE_f, MSE_model,MSE_lo, MSE_fo)
colnames(Result) = c("Rsquared_f1","Rsquared_f2","Rsquared_l1","Rsquared_l2","MSE_l","MSE_f", "MSE_model","MSE_lo","MSE_fo")

re = colMeans(Result)
reT= t(re)
print(tau)
print(N)
print(T)
re_list = append(re_list,reT)
}
#}

new_re = t(array(re_list,c(9,8)))
write.csv(new_re, 're.csv')
######################
ts_mean_return = colMeans(return_in)
ts_mean_return=  (ts_mean_return  - min(ts_mean_return))/(max(ts_mean_return) - min(ts_mean_return))


plot.vec(return_in, type = "l", lwd = 2, ylab = "", xlab = "", color = expgr)
plot(factors[,1], type = "l", lwd = 2, col = "green")
lines(ft_1,type = "l",col = "blue")
plot(factors[,2], type = "l", lwd = 2, col = "green")
lines(ft_2,type = "l",col = "red")


apply(resid.pec(tempsm, output.pec[[2]], fitted = T, tau), 1, mean) + 
  apply(resid.pec(tempsm, output.pec[[2]], tau), 1, expectile, alpha )

plot(lt_1,type = "l",col = "red",lwd=2,ylab ="", xlab="")
lines(loadings[,1],type = "l",col = "blue")
######################
