##### simulation for many times ##############
setwd("/Users/wendy/Dropbox/Dynamic expectile factor model/")
rm(list = ls())
graphics.off()
libraries = c("svd","fda","vars", "stats", "tseries",  "dygraphs", "CCA", "ggplot2","quantreg","mboost", "BayesX")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})

install.packages("~/Dropbox/My Mac (王冰菱’s MacBook Pro (2))/Downloads/expectreg_0.39.tar", repos = NULL,
                 type = "source")
#"qgraph"
lapply(libraries, library, quietly = TRUE, character.only = TRUE)



library(parallel)
library(expectreg)

nsimu = 50
MSE_model = rep(0,nsimu)
MSE_f = rep(0,nsimu)
MSE_l = rep(0,nsimu)
Rsquared_f1 = rep(0,nsimu)
Rsquared_f2 = rep(0,nsimu)
Rsquared_l1 = rep(0,nsimu)
Rsquared_l2 = rep(0,nsimu)

  tau = 0.05
  #tau = 0.5
  #tau = 0.95
### start the loop of the simulation ###
for (simu_i in 1: nsimu)
{
cat("Simu:",simu_i)
 ### simulate data matrix 
  N = 100
  T = 1000
  ### simulate error term
  e = matrix(0, T, 1)
  sigma = 0.1
  #sigma = 0.01
  ### where lth column of e is an i.i.d. sample from N(0, sigma), for l from 1 to N.
  for (l in  1:N)
  {
    #one_col_e = rnorm(T, 0, sigma)
    one_col_e = rt(T, 5)
    e = cbind(e, one_col_e)
  }
  e = e[,-1]
#set.seed(123)
#Just generate random AR(1) time series; based on this, I want to estimate the parameters
f_1 = arima.sim(n=T, list(ar=c(0.5)))
f_2 = arima.sim(n=T, list(ar=c(0.3,.2))) 
l_1 = rnorm(N, 0, 3)
l_2 = rnorm(N, 0, 2)

lt_1 = sign(l_1[1])*l_1
lt_2 = sign(l_2[1])*l_2
signf1 = ifelse(sign(lt_1[1])==sign(f_1[1]),1,-1)
signf2 = ifelse(sign(lt_2[1])==sign(f_2[1]),1,-1)
ft_1 = signf1*f_1
ft_2 = signf2*f_2

return_in = ft_1%*%t(lt_1)+ft_2%*%t(lt_2)+e
n = ncol(return_in)

#### expectile regression case
exp_seq = rep(0,N)
for (j in 1:ncol(return_in)){
exp_seq[j] = expectile(return_in[,j],tau)
}

diff = matrix(0, nrow(return_in), ncol(return_in))
for (j in 1:ncol(return_in)){
diff[,j] = return_in[,j]-exp_seq[j]
}
k = 2
tau_matrix = ifelse(diff>0, -tau, tau-1)
#### first derivative of loss function ####
expt = 2*diff*tau_matrix
cov_tau_matrix = cov(expt)
##### Spectral Decomposition of a Matrix ###########
e_val = eigen(cov_tau_matrix, symmetric = TRUE)$value
e_vec = eigen(cov_tau_matrix, symmetric = TRUE)$vector

loadings = t(sqrt(n)*e_vec[,1:k])
#plot(e_val)
sum(e_val[1:k])/sum(e_val)
#pro = rep(0,N)
#for (j in 1: N){
#pro[j] = sum(e_val[1:j])/sum(e_val)
#}
#plot(pro)
###loadings = rbind(rnorm(N, 0, 6), rnorm(N, 0, 3))

loadings[1,] = sign(loadings[1,1])*loadings[1,]
loadings[2,] = sign(loadings[2,1])*loadings[2,]

factors = matrix(0, nrow(return_in), k)
inter = rep(0, nrow(return_in))
#pred = matrix(0, nrow(return_in), ncol(return_in))
for (i in 1:nrow(return_in)){
fit = expectreg.ls(as.numeric(return_in[i,])~loadings[1,]+loadings[2,], estimate="laws",expectiles=tau)
factors[i,] = as.numeric(fit$coefficients)
inter[i] = as.numeric(fit$intercepts)
#pred[i,] = as.numeric(unlist(predict(fit)[1]))
#newfactor[i,] = pre[1]
}
### identification, make sure the same loading and corresponding factor have the same sign
signf1 = ifelse(sign(loadings[1,1])==sign(factors[1,1]),1,-1)
signf2 = ifelse(sign(loadings[2,1])==sign(factors[1,2]),1,-1)
factors[,1] = signf1*factors[,1]
factors[,2] = signf2*factors[,2]

est_model = factors%*%loadings+inter
F1 = factors[,1]
F2 = factors[,2]
ft_1 = (ft_1 - min(ft_1))/(max(ft_1) - min(ft_1))
F1 = (F1 - min(F1))/(max(F1) - min(F1))
ft_2 = (ft_2 - min(ft_2))/(max(ft_2) - min(ft_2))
F2 = (F2 - min(F2))/(max(F2) - min(F2))
lt_1 = (lt_1 - min(lt_1))/(max(lt_1) - min(lt_1))
loadings[1,] = (loadings[1,] - min(loadings[1,]))/(max(loadings[1,]) - min(loadings[1,]))
lt_2 = (lt_2 - min(lt_2))/(max(lt_2) - min(lt_2))
loadings[2,] = (loadings[2,] - min(loadings[2,]))/(max(loadings[2,]) - min(loadings[2,]))
lt = rbind(lt_1, lt_2)
ft = cbind(ft_1, ft_2)
factors = cbind(F1, F2)
MSE_l[simu_i] =sqrt(mean((lt-loadings)^2))
MSE_f[simu_i] = sqrt(mean((ft-factors)^2))
MSE_model[simu_i] = sqrt(mean((return_in-est_model)^2))
#MSE_loading[simu_i]= sqrt(mean((lt_1-loadings[1,])^2+(lt_2-loadings[2,])^2))
#MSE_factor[simu_i] = sqrt(mean(cbind((ft_1-F1)^2,(ft_2-F2)^2)))

Rsquared_f1[simu_i] = summary(lm(ft_1~F1))$adj.r.squared
Rsquared_f2[simu_i] = summary(lm(ft_2~F2))$adj.r.squared
Rsquared_l1[simu_i] = summary(lm(lt_1~loadings[1,]))$adj.r.squared
Rsquared_l2[simu_i] = summary(lm(lt_2~loadings[2,]))$adj.r.squared
}
Result = cbind(Rsquared_f1,Rsquared_f2,Rsquared_l1,Rsquared_l2,MSE_l, MSE_f, MSE_model)
colnames(Result) = c("Rsquared_f1","Rsquared_f2","Rsquared_l1","Rsquared_l2","MSE_l","MSE_f", "MSE_model")

re = colMeans(Result)
write.csv(re, 're.csv')
re

write.csv(Rsquared, file = paste("Rsquared_Normal","_tau",tau,"_N",N,"_T",T,".csv", sep =""))


i = 1
plot(return_in[,i])
lines(est_model[,i])
plot(return_in[i,])
lines(est_model[i,])
lt_1 = (lt_1 - min(lt_1))/(max(lt_1) - min(lt_1))
loadings[1,] = (loadings[1,] - min(loadings[1,]))/(max(loadings[1,]) - min(loadings[1,]))
plot(lt_1, type = "l",col = "red",lwd=5,ylab ="", xlab="",cex.axis = 2.5)#
lines(loadings[1,], type = "l",col = "blue",lwd=6,ylab ="", xlab="")
#print(MSE_factor[simu_i])
MSE_loading[simu_i] =sqrt(sum((lt_1-loadings[1,])^2+(lt_2-loadings[2,])^2))
#print(MSE_loading[simu_i])
expectile(a, tau)
