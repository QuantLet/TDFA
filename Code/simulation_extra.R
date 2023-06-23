setwd("/Users/wendy/Dropbox/Dynamic expectile factor model/Code")
rm(list = ls())
graphics.off()
libraries = c("svd","fda","vars", "stats", "tseries",  "dygraphs", "CCA", "ggplot2","quantreg","mboost", "BayesX","
Trig")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})


#"qgraph"
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

install.packages('expectreg')
library(parallel)
library(expectreg)
source("fpca.R", encoding="UTF-8")
source("PrincipalDirection.R")

#tau = 0.05
#tau = 0.5
tau = 0.95

N = 100
T = 1000
### simulate error term
e = array(0,c(N, T))
sigma = 0.5
#sigma = 0.01
### where lth column of e is an i.i.d. sample from N(0, sigma), for l from 1 to N.
for (i in  1:T)
{
  one_col_e = rnorm(N, 0, sigma)
  #one_row_e = rt(T, 5)
  e[,i] = one_col_e
}
e = as.matrix(e)


t_seq = seq(0,1,length.out = T)
mu = 1+ t_seq + exp(-(t_seq-0.6)^2/0.05)
phi1 = sqrt(2)*sin(2*pi*t_seq)
phi2 = sqrt(2)*cos(2*pi*t_seq)
#phi = cbind(phi1,phi2)
#matplot(phi)
alpha1 = rnorm(N, 0, 6)
alpha2 = rnorm(N, 0, 3)
return_in = mu + alpha1%*%t(phi1) + alpha2%*%t(phi2) + e


l_1 = alpha1
l_2 = alpha2
lt_1 = sign(l_1[1])*l_1
lt_2 = sign(l_2[1])*l_2
loadings = cbind(lt_1,lt_2)
 

f_1 = phi1
f_2 = phi2
signf1 = ifelse(sign(lt_1[1])==sign(f_1[1]),1,-1)
signf2 = ifelse(sign(lt_2[1])==sign(f_2[1]),1,-1)
ft_1 = signf1*f_1
ft_2 = signf2*f_2
factors = rbind(ft_1,ft_2)

return = mu + loadings%*%factors + e

alpha = tau-0.5

### calculate the expectile of errors and returns
e_e = rep(0, T)
e_return =  rep(0, T)
for (t in 1:T){
  e_e[t] = expectile(e[,t],alpha)
  e_return[t] = expectile(return[,t],alpha)
}

return_s = as.matrix(return-e_return)


#### pec case
set.seed(1234)
output.pec0 =	pec.k(return_s,nk=2,alpha=0.45)
output.pec = lapply(alpha,pec.k,Y=return_s, nk=2,reset.tol=50,lab.ini=-output.pec0[[4]])
basis.pec = lapply(output.pec,listBasis,basis0=output.pec[[1]][[2]],num=2)
PCs = matrix(unlist(basis.pec), nrow(return),2)
scores_pec = rbind(basis.pec[[1]][,1]%*%return_s, basis.pec[[1]][,2]%*%(return_s-apply(basis.pec[[1]][,1]%*%return_s, 2, "*",basis.pec[[1]][,1])))


F1 = as.matrix(sign(scores_pec[1,1])*scores_pec[1,])
F2 = as.matrix(sign(scores_pec[2,1])*scores_pec[2,])
factors = cbind(F1, F2)
#signf1 = ifelse(sign(loadings[1,1])==sign(factors[1,1]),1,-1)
#signf2 = ifelse(sign(loadings[2,1])==sign(factors[1,2]),1,-1)
#factors[,1] = signf1*factors[,1]
#factors[,2] = signf2*factors[,2]


#factors = t(scores_pec)
loadings = matrix(0, nrow(return), 2)
inter = rep(0, nrow(return))
#pred = matrix(0, nrow(return_in), ncol(return_in))
for (i in 1:nrow(return)){
  fit = lm(as.numeric(return_s[i,])~factors[,1]+factors[,2])
  #loadings[i,] = fit$coefficients
  loadings[i,] = fit$coefficients[-1]
  inter[i] = as.numeric(fit$coefficients[1])
  #pred[i,] = as.numeric(unlist(predict(fit)[1]))
  #newfactor[i,] = pre[1]
}


### identification, make sure the same loading and corresponding factor have the same sign
est_model = loadings%*%t(factors)+inter
F1 = factors[,1]
F2 = factors[,2]
ft_1 = ft_1+e_e
ft_1 = (ft_1 - min(ft_1))/(max(ft_1) - min(ft_1))
F1 = (F1 - min(F1))/(max(F1) - min(F1))
ft_2 = ft_2+e_e
ft_2 = (ft_2 - min(ft_2))/(max(ft_2) - min(ft_2))
F2 = (F2 - min(F2))/(max(F2) - min(F2))
lt_1 = (lt_1 - min(lt_1))/(max(lt_1) - min(lt_1))
loadings[,1] = (loadings[,1] - min(loadings[,1]))/(max(loadings[,1]) - min(loadings[,1]))
lt_2 = (lt_2 - min(lt_2))/(max(lt_2) - min(lt_2))
loadings[,2] = (loadings[,2] - min(loadings[,2]))/(max(loadings[,2]) - min(loadings[,2]))
lt = cbind(lt_1, lt_2)
ft = cbind(ft_1, ft_2)
factors = cbind(F1,F2)
MSE_l[simu_i] =sqrt(mean((lt-loadings)^2))
MSE_f[simu_i] = sqrt(mean((ft-factors)^2))
MSE_model[simu_i] = sqrt(mean((return_s-est_model)^2))

######################
plot(ft_2,type = "l",col = "red",lwd=2,ylab ="", xlab="")
lines(factors[,2],type = "l",col = "blue")

plot(lt_1,type = "l",col = "red",lwd=2,ylab ="", xlab="")
lines(loadings[,1],type = "l",col = "blue")
######################

  