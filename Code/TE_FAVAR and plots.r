
###################### Stress test macro variables ###########################
rm(list = ls())
graphics.off()
libraries = c("fda","vars", "stats", "tseries",  "dygraphs", "CCA", "expectreg", "qgraph", "ggplot2")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)
setwd("/Users/wendy/Dropbox/Dynamic expectile factor model/Code")
#in mac computer "qgraph","ggplot2", "GGally",
#setwd("....")
source("fpca.R", encoding="UTF-8")
source("PrincipalDirection.R")
#source("gfevd.r")
#sort(apply(return_in[,2:26], 1, sd))
## for  period 1991--1996
#st = 1
#ed = 71
## for  period 1997--2000 crisis
#st = 72
#ed = 119
## for  period 2001--2006 
#st = 120
#ed = 191
## for crisis period 2008
#st = 192
#ed = 239
## for stable period
#st = 240
#ed = 287
## other period ###
st = 1
ed = 287
#ed = 310
#[2:26]
#return_in = read.csv("49_Industry_Portfolios_returns.csv")[,-1]
#rownames(return_in) = return_in
return_in=as.matrix(read.csv(("100_Financial_firms_returns_abb.csv"))[2:101])
##return_in=read.csv(("100_Financial_firms_returns_abb.csv"))[77:101]
rownames(return_in)= read.csv(("100_Financial_firms_returns_abb.csv"))[,1]
dt = as.Date(rownames(return_in),format = "%d/%m/%Y")#[25:length(rownames(return_in))][st:ed]
return_in = return_in[st:ed,]
#return_in = t(return_in)
  tau = 0.05
  #tau = 0.5
  #tau = 0.95
  alpha = -0.45 ### alpha + 0.5 is the expectile level tau
  #alpha	= 0
  #alpha	= 0.45
    return =  t(as.matrix(return_in))
#for (j in 1:ncol(return)) {
# return[, j] = (return[, j] - min(return[, j]))/(max(return[, j]) - min(return[, j]))
#}
    set.seed(1234)
    output.pec0 =	pec.k(return,nk=2,alpha=-0.44)
    output.pec = lapply(alpha,pec.k,Y=return, nk=2,reset.tol=50,lab.ini=-output.pec0[[4]])
    #PCs = output.pec[[1]][[2]] #The first two estimated PCs
    basis.pec = lapply(output.pec,listBasis,basis0=output.pec[[1]][[2]],num=2)
    PCs = matrix(unlist(basis.pec), nrow(return),2)
    scores_pec = rbind(basis.pec[[1]][,1]%*%return, basis.pec[[1]][,2]%*%(return-apply(basis.pec[[1]][,1]%*%return, 2, "*",basis.pec[[1]][,1])))

#F1 = abs(F1)
#F2 = abs(F2)
#F1 = as.matrix(scores_pec[1,])
#F2 = as.matrix(scores_pec[2,])
 
######### expectile regression ##########
#scores_pec = matrix(0,ncol(PCs),nrow(return_in))
#for (i in 1:nrow(return_in)){
#scores_pec[,i] = as.numeric(expectreg.ls(as.numeric(return_in[i,])~PCs[,1]+PCs[,2], estimate="laws",expectiles=tau)$coefficients)
##scores_pec[,i] = as.numeric(expectreg.ls(temp[,i]~PCs[,1]+PCs[,2], estimate="laws",expectiles=0.5)$coefficients)
##scores_pec[,i] = as.numeric(expectreg.ls(temp[,i]~PCs[,1]+PCs[,2], estimate="laws",expectiles=0.95)$coefficients)
#}
##}

    F1 = as.matrix(sign(scores_pec[1,1])*scores_pec[1,])
    F2 = as.matrix(sign(scores_pec[2,1])*scores_pec[2,])
    scores_pec = t(cbind(F1, F2))
    PCs[,1] = sign(PCs[1,1])*PCs[,1]
    PCs[,2] = sign(PCs[1,2])*PCs[,2]
    PCs = cbind(PCs[,1], PCs[,2])

############## plot mean and score 1  ######################
ts_mean_return = colMeans(return)
scores_pec[1,] =  (scores_pec[1,]  - min(scores_pec[1,]))/(max(scores_pec[1,]) - min(scores_pec[1,]))
ts_mean_return=  (ts_mean_return  - min(ts_mean_return))/(max(ts_mean_return) - min(ts_mean_return))
#dev.new()
#pdf(file = paste("FC1scores_Meanreturns_new111", tau, ".pdf", sep =""),onefile = F, width = 16, height = 8) 
pdf(file = paste("FC1scores_Meanreturns_new111", tau, ".pdf", sep =""),onefile = F, width = 8, height = 8) 

#plot(dt, ts_mean_return, type ="l", col ="cornflowerblue", lwd = 3,xlab ="", ylab ="", cex.lab = 1.5,cex.axis = 2.5, font.axis = 1.5)#, ylim = range(-0.5,0.3)
#lines(dt, scores_pec[1,], type ="l", col ="grey", lwd = 5)
plot(dt, ts_mean_return, pch = 20, col ="black", ylim = range(-0,1), lwd = 5,xlab ="", ylab ="", cex.lab = 1.5,cex.axis = 2.5, font.axis = 1.5)
lines(dt, scores_pec[1,], type ="l", col ="grey", lwd = 5)#
#lines(dt, scores_pec[1,], type ="l", col ="red", lwd = 6)#
#lines(dt, scores_pec[1,], type ="l", col ="blue", lwd = 3)#

##lines(dt, scores_pec[2,], type ="l", col ="blue", lwd = 3)#, ylim = range(-0.5,0.3)
dev.off()
ts_sd_return = apply(return,2, sd)
cor((ts_mean_return),scores_pec[1,], method = "pearson")
#cor((ts_sd_return),scores_pec[1,], method = "pearson")

#F1 =log(abs(scores_pec[1,]))
#F2 =log(abs(scores_pec[2,]))
#F1 = log((scores_pec[1,])^2)
#F2 = log((scores_pec[2,])^2)
scores_pec[2,] =  (scores_pec[2,]  - min(scores_pec[2,]))/(max(scores_pec[2,]) - min(scores_pec[2,]))

#pdf(file = paste("FCs2_05_005.pdf", sep =""),onefile = F, width = 8, height = 8) 
plot(dt,F1, type = "l",col = "red",lwd=5,ylab ="", xlab="",  ylim =range(-3,3),cex.axis = 2.5)#
#plot(dt,scores_pec[2,], type = "l",col = "blue",lwd=3,ylab ="", xlab="", cex.lab = 1.5,cex.axis = 2.5, font.axis = 1.5)#
lines(dt,scores_pec[2,], type = "l",col = "red",lwd=6,ylab ="", xlab="")
#lines(dt,F1, type = "l",col = "blue",lwd=2,ylab ="", xlab="")
#dev.off()
##################################
   ess	= matrix(NA,length(tau),1)
for (k in 1:length(tau)){
tau0		= tau[[k]]
sstot		<-  sum(return^2*(sign(return)*(tau0 - 0.5) + 0.5))
ess[k,1] 	<- 1-resid.pec(return,output.pec[[k]],tau=tau0,rss=T)/sstot
}
ess	<- data.frame(round(ess,2))
rownames(ess) = "tau=0.05"
colnames(ess)	= "PEC"
ess



#############################################################
macros = read.csv(file="Macro economic variables_new.csv") 
UNE = as.vector(macros[,2]) # nonstationary
CPI = as.vector(macros[,3]) # nonstationary
THR = as.vector(macros[,4]) # nonstationary
GDP = as.vector(read.csv("GDP_cubic_growth_new.csv")[,2][-1]) #  stationary, GDP_cubic_growth rate
NGD = as.vector(read.csv("nominal_GDP_cubic_growth_new.csv")[,2][-1])
macros1 = read.csv(file="Macro economic variables.csv") 
SP  = as.vector(macros1[,4]) # nonstationary
macros2 = read.csv(file="9_macro_variables.csv") 
INC =  as.vector(macros2[,2]) # nonstationary
FIV =  as.vector(macros2[,8])[-1] # stationary
TEN =  as.vector(macros2[,3])[-1] # stationary
HPI =  as.vector(macros2[,4]) # nonstationary 
VIX =  as.vector(macros2[,6])
PRI =  as.vector(macros2[,5])
NIN =  as.vector(macros2[,7])
adf.test(F1)
adf.test(F2)
adf.test(UNE)
adf.test(GDP)
adf.test(UNE)
adf.test(THR)
adf.test(TEN)
adf.test(CPI)
adf.test(INC)
adf.test(SP)
adf.test(HPI)
adf.test(diff(CPI))
adf.test(diff(THR))
adf.test(diff(SP))
UNE = UNE[-1]
#UNE = diff(UNE)
#UNE =  diff(UNE)/UNE[-length(UNE)]*100
#CPI = CPI[-1]
#CPI = diff(CPI)
CPI = diff(CPI)/CPI[-length(CPI)]*100
adf.test(CPI)
THR = (THR)[-1]
#THR  = diff(THR)
#THR = diff(THR)/THR[-length(THR)]*100
SP =  diff(log(SP))
#SP =  diff(SP)
#SP =  SP[-1]
#SP = diff(SP)/SP[-length(SP)]*100
adf.test(SP)
#INC = INC[-1]
#INC =  diff(INC)
INC =  diff(INC)/INC[-length(INC)]*100
adf.test(INC)
NIN =  diff(NIN)/NIN[-length(NIN)]*100
#HPI = HPI[-1]
HPI = diff(HPI)
#HPI = diff(HPI)/HPI[-length(HPI)]*100
adf.test(HPI)
PRI = PRI[-1]

VARmodel = cbind(F1, F2, GDP[st:ed], CPI[st:ed], INC[st:ed], UNE[st:ed], THR[st:ed], TEN[st:ed], SP[st:ed], HPI[st:ed] )#
colnames(VARmodel) = c("F1","F2", "GDP", "CPI", "INC", "UNE", "THR", "TEN", "SP","HPI")

VARselect(VARmodel, lag.max =20)#
varest1 = VAR(VARmodel, p = 15) #, type = "both"
#varest1 = VAR(VARmodel, p = 7) 
#varest1 = VAR(VARmodel, p = 10) ### expectile 005
resid <- residuals(varest1)
root = roots(varest1) ##the var process is stationary
summary(root)

### we found that if order 4 is optimal. If we choose large order, the result of the impulse response analysis fluctuated and the trend is unclear.

serial.test(varest1,lags.pt = 20, type = "PT.asymptotic")
serial.test(varest1,lags.pt = 20,  type = "PT.adjusted")#
serial.test(varest1,  lags.pt = 16,type = "BG")#
serial.test(varest1,lags.pt = 16,  type = "ES")#

#IRF_ALL = irf(varest1, impulse = NULL, response = NULL, n.ahead = 27,ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95)
#plot(IRF_ALL,plot.type = "multiple")

causality(varest1,cause = "F1")$Granger # yes  0.007358
causality(varest1,cause = "F2")$Granger # yes 3.822e-05
causality(varest1,cause = "GDP")$Granger # yes 4.645e-05
causality(varest1,cause = "CPI")$Granger # no  0.3228
causality(varest1,cause = "INC")$Granger # no 0.3188
causality(varest1,cause = "UNE")$Granger # no 0.262
causality(varest1,cause = "THR")$Granger # yes 0.00087
causality(varest1,cause = "TEN")$Granger # no 0.1246
causality(varest1,cause = "SP")$Granger # no 0.2189
causality(varest1,cause = "HPI")$Granger # yes 0.007765

######################################################
n  = nrow(return_in)
dt_in = read.csv("dt_forecast.csv")[,1]
dt_f = as.Date(dt_in,format = "%d/%m/%Y")
fw = 12 ### one year forecasting
#fw = 6 ### one year forecasting
nv = 10 ### number of variables
ws = 48 ### two years moving window
total.connect = matrix(0, (n - ws), 1)
OUT_link = matrix(0, (n - ws), nv)
IN_link  = matrix(0, (n - ws), nv)
total_OUT= matrix(0, (n - ws), 1)
total_IN= matrix(0, (n - ws), 1)
total.direct = matrix(0, nv,nv)

dt_fc = dt_f[(1+ws+fw):(n+fw)]

for (i in 1:(n - ws)) {
    print(i)
    wF1 = F1[i:(i + ws)]
    wF2 = F2[i:(i + ws)]
    wGDP = GDP[i:(i + ws)]
    wCPI = CPI[i:(i + ws)]
    wUNE = UNE[i:(i + ws)]
    wINC = INC[i:(i + ws)]
    wTHR = THR[i:(i + ws)]
    wTEN = TEN[i:(i + ws)]
    wSP = SP[i:(i + ws)]
    wHPI = HPI[i:(i + ws)]
wVAR = cbind(wF1, wF2, wGDP, wCPI, wINC, wUNE, wTHR, wTEN, wSP, wHPI)#
#wVAR = cbind(wF1, wF2, wGDP, wCPI, wINC, wUNE, wTHR, wHPI)#
#colnames(wVAR) = c("F1", "F2", "GDP", "CPI", "INC", "UNE", "THR", "HPI")

# standardize macro variables
#    for (j in 1:ncol(wVAR)) {
#      wVAR[, j] = (wVAR[, j] - min(wVAR[, j]))/(max(wVAR[, j]) - min(wVAR[, j]))
#    }
colnames(wVAR) = c("F1","F2", "GDP", "CPI", "INC", "UNE", "THR", "TEN", "SP", "HPI")

 #### the order selection ################
Ord_Select =as.numeric(unlist(VARselect(wVAR, lag.max = 20)[1]))
#for (l in Ord_Select){
for (l in 1:20){
varest2 = VAR(wVAR, p = l) ## tau = 0.05, type = "both"
resid <- residuals(varest2)
PT_as = as.numeric(unlist(serial.test(varest2, lags.pt = 16, type = "PT.asymptotic")[2])[3])
PT_ad = as.numeric(unlist(serial.test(varest2, lags.pt = 16, type = "PT.adjusted")[2])[3])
BG = as.numeric(unlist(serial.test(varest2, lags.pt = 16, type = "BG")[2])[3])
ES = as.numeric(unlist(serial.test(varest2, lags.pt = 16, type = "ES")[2])[4])
if( PT_as >=0.05 | PT_ad>=0.05 | BG>=0.05 | ES>=0.05)
break
l = l+1
}
#print(l)
#cat("PT_as:",PT_as)
#cat("PT_ad:", PT_ad)
#cat("BG:", BG)
#cat("ES:",ES)
ord = l
  varest = VAR(wVAR, p = ord) ## tau = 0.05, type = "both"
  resid <- residuals(varest)
  vd = fevd(varest, n.ahead = fw)
  #plot(vd, plot.type = "multiple",cex.lab=1.5,cex.axis=2, cex.main = 2, ylab = "", xlab = "")
  vd_matrix <- matrix(cbind(unlist(print(vd[1])), unlist(print(vd[2])), unlist(print(vd[3])), unlist(print(vd[4])), unlist(print(vd[5])), unlist(print(vd[6])), unlist(print(vd[7])), unlist(print(vd[8])), unlist(print(vd[9])), unlist(print(vd[10]))), ncol = fw, byrow = TRUE)#
  #vd_matrix <- matrix(cbind(unlist(print(vd[1])), unlist(print(vd[2])), unlist(print(vd[3])), unlist(print(vd[4])), unlist(print(vd[5])), unlist(print(vd[6])), unlist(print(vd[7])), unlist(print(vd[8]))), ncol = fw, byrow = TRUE)
  
  adj_vd = matrix(vd_matrix[, fw], nv,nv)
  diag(adj_vd)=0
  print(adj_vd)
  OUT_link[i,] = rowSums(adj_vd)
  total_OUT[i] = sum(rowSums(adj_vd))
  IN_link[i,]  = colSums(adj_vd)
  total_IN[i] = sum(colSums(adj_vd))  
  total.connect[i] = (sum(adj_vd))/nv
  total.direct = total.direct + adj_vd
  ##### add more variables, first use these four to do net work
  #col = c(rep("red", 2), rep("blue", 2), rep("green4", 2), rep("mediumorchid4", 
  #2))
  #col = c("red", "orange", "yellow", "greenyellow", "limegreen", "cyan", "dodgerblue", "mediumblue", "purple",  "violet" )
  #names =  colnames(wVAR)
  #pdf(file = paste("FA_Network",i,".pdf", sep =""),onefile = F, width = 8, height = 8) 
  #plot_g = qgraph(adj_vd, layout = "circle", layoutScale = c(1.2, 1.2), label.font = 2, label.cex = 2, shape = "circle", labels = names, esize = 5,   maximum = max(adj_vd), color = "white", node.width = 0.8, label.cex = 1.8, label.color = col,  edge.color = col, curve = 1, border.width = 3, border.color = col, asize = 5)
  ##plot_g = qgraph(adj_vd, layout = "spring", layoutScale = c(1.2, 1.2), label.font = 2, label.cex = 2, shape = "circle", labels = names, esize = 5,   maximum = max(adj_vd), color = "white", node.width = 0.8, label.cex = 1.8, label.color = col,  edge.color = col, curve = 1, border.width = 1.2, border.color = col, asize = 5)
  #text(x = -0.9, y = 1.2, labels = substitute(paste(d), list(d = as.character(dt_fc[i]))), ,xpd = NA, cex = 1.5)
  ##text(x = -0.9, y = 1, labels = substitute(paste(d), list(d = as.character(dt[i+ws+fw]))),xpd = NA, cex = 1.5)
  #dev.off()
}

##################### write plots in a folder ###################

############### plot estimated factors #########################
pdf(file = paste("FCs_005.pdf", sep =""),onefile = F, width = 8, height = 8) 
plot(dt,F1, type = "l",col = "red",lwd=5,ylab ="", xlab="",  ylim =range(-3,3),cex.axis = 2.5)#
#plot(dt,F2, type = "l",col = "blue",lwd=5,ylab ="", xlab="",  ylim =range(-3,3),cex.axis = 2.5)#
lines(dt,F2, type = "l",col = "red",lwd=2,ylab ="", xlab="")
#lines(dt,F1, type = "l",col = "blue",lwd=2,ylab ="", xlab="")
dev.off()


##### plot factor loadings #####################################
pdf(file = paste("FC_loadings_new_005.pdf", sep =""),onefile = F, width = 8, height = 8) 
lc = ncol(return_in)/4
color = c(rep("red",lc),rep("blue",lc),rep("green4",lc),rep("mediumorchid4",lc))
plot(PCs,type="n",  ylim =range(-0.6,0.6),xlim=range(-0.6,0.6),cex.axis = 1.5, xlab ="F 1", ylab ="F 2", cex.lab = 1.5,cex.axis = 1.8, font.axis = 1.5) 
text(PCs,labels=colnames(return_in),col = color,cex=2)
abline(h=0.0,v=0.0)
dev.off()

######################## plot correlation ########################
#
correlation = cor(t(scores_pec),t(return))
#correlation = cor((PCs),(return))
ucircle = cbind(cos((0:360)/180*pi),sin((0:360)/180*pi))
pdf(file = paste("cor_","tau", ".pdf", sep =""),onefile = F, width = 8, height = 8) 
plot(ucircle,type="l",lty="solid",col="black",xlab="F 1",ylab="F 2",,cex.axis = 2, cex.lab = 1.5, font.axis = 1.5)
abline(h=0.0,v=0.0)
lc = ncol(return_in)/4
#lc = nrow(return_in)/4
color = c(rep("red",lc),rep("blue",lc),rep("green4",lc),rep("mediumorchid4",lc))
text(t(correlation),labels=rownames(return),col = color,cex=2)
dev.off()

##### plot factor loadings ###########################
#pdf(file = paste("FC1_loadings_crisis.pdf", sep =""),onefile = F, width = 8, height = 8) 
lc = ncol(return_in)/4
color = c(rep("red",lc),rep("blue",lc),rep("green4",lc),rep("mediumorchid4",lc))
pdf(file = paste("FC1_loadings", tau,".pdf", sep =""),onefile = F, width = 8, height = 8) 
#plot(PCs,type="n",  ylim =range(-0.6,0.6),xlim=range(-0.6,0.6),cex.axis = 1.5, xlab ="F 1", ylab ="F 2", cex.lab = 1.5,cex.axis = 1.8, font.axis = 1.5) # set up plot 
plot(PCs,type="n",  ylim =range(-0.65,0.6),xlim=range(-0.6,0.6),cex.axis = 1.5, xlab ="F 1", ylab ="F 2", cex.lab = 1.5,cex.axis = 1.8, font.axis = 1.5) # set up plot 
text(PCs,labels=rownames(return),col = color,cex=2)
abline(h=0.0,v=0.0)
dev.off()
#write.csv(PCs, file = "PC_changedata_005.csv")

############## plot mean and score 1  ######################
ts_mean_return = colMeans(return)
scores_pec[1,] =  (scores_pec[1,]  - min(scores_pec[1,]))/(max(scores_pec[1,]) - min(scores_pec[1,]))
ts_mean_return=  (ts_mean_return  - min(ts_mean_return))/(max(ts_mean_return) - min(ts_mean_return))
dev.new()
pdf(file = paste("FC1scores_Meanreturns_new", tau, ".pdf", sep =""),onefile = F, width = 16, height = 8) 
#plot(dt, ts_mean_return, type ="l", col ="cornflowerblue", lwd = 3,xlab ="", ylab ="", cex.lab = 1.5,cex.axis = 2.5, font.axis = 1.5)#, ylim = range(-0.5,0.3)
#lines(dt, scores_pec[1,], type ="l", col ="grey", lwd = 5)
plot(dt, ts_mean_return, pch = 20, col ="black", ylim = range(-0,1), lwd = 5,xlab ="", ylab ="", cex.lab = 1.5,cex.axis = 2.5, font.axis = 1.5)
lines(dt, scores_pec[1,], type ="l", col ="grey", lwd = 5)#
#lines(dt, scores_pec[2,], type ="l", col ="blue", lwd = 3)#, ylim = range(-0.5,0.3)
dev.off()
ts_sd_return = apply(return,2, sd)
cor((ts_mean_return),scores_pec[1,], method = "pearson")
#cor((ts_sd_return),scores_pec[1,], method = "pearson")

#plot(dt, ts_mean_return, pch = 20, col ="black", ylim = range(-1,1), lwd = 5,xlab ="", ylab ="", cex.lab = 1.5,cex.axis = 2.5, font.axis = 1.5)

####################### plot against expectile #################
########
expec_ts_all = as.matrix(read.csv("expectile005_all_firms.csv")[,2])
#expec_ts_all = as.matrix(read.csv("expectile005_all_firms.csv")[,2])
cor(expec_ts_all,scores_pec[1,], method = "pearson")
#cor(expec_ts_all, PCs[,1])
ts_mean_return = rowMeans(return)
cor(ts_mean_return, PCs[,1])
dev.new()
scores_pec[1,]=  (scores_pec[1,]  - min(scores_pec[1,]))/(max(scores_pec[1,]) - min(scores_pec[1,]))
expec_ts_all=  (expec_ts_all  - min(expec_ts_all))/(max(expec_ts_all) - min(expec_ts_all))
pdf(file = paste("FC1scores_Expectreturns", tau,".pdf", sep =""),onefile = F, width = 16, height = 8) 
plot(dt,scores_pec[1,],col="grey",type ="l", lwd = 3, ylim = range(0, 1),cex.axis=2.5, ylab="",xlab="")
lines(dt,expec_ts_all, col = "red",lwd = 5)
dev.off()

########################## acf plots #############
#pdf(file = "acf_05_%03d.pdf", width = 14, height = 8, onefile = FALSE)
par ( mfrow =c(5 ,8))
acf(resid, ylab = "",  cex.axis =  2.0, lwd = 5, xlab="", cex.main = 1.8)[1]
#dev.off()


######################### IRF plots ##############
nf = 27
nv = 10
for (i in colnames(VARmodel)){
IRF = irf(varest1, impulse = i, response = c("F1", "F2", "GDP", "CPI", "INC", "UNE", "THR", "TEN", "SP","HPI"), n.ahead = nf,ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95)
pdf(file = paste("IRF_new", tau, i,".pdf", sep =""),onefile = F, width = 8, height = 10) 
plot(IRF, lwd=3, cex.axis = 3, font.axis = 2,cex.lab=1.4,sub ="", cex.main=2)
dev.off()
}

########################### prediction #############
######### Forecast ################################
require(forecast)
library(forecast)
pred = predict(varest1, n.ahead=27, ci=0.95,)
pdf(file = "forecast1_005_%03d_12.pdf", width = 16, height = 8, onefile = FALSE)
#par ( mfrow =c(5 ,8))
plot(pred,cex.lab=3,cex.axis=3, cex.main = 2, ylab = "", xlab = "", lwd = 5, lty =1, plot.type = "single", ylim = range(-3,3)) # "multiple"
dev.off()

######################### total connection ###############
  write.csv(total.connect,file = paste("DEFM_tc_","ws",ws,"fw",fw,".csv", sep =""))
  write.csv(OUT_link,file = paste("OUT_link_new","ws",ws,"fw",fw,".csv", sep =""))
  write.csv(IN_link,file = paste("IN_link_new","ws",ws,"fw",fw,".csv", sep =""))
  write.csv(total_OUT,file = paste("total_OUT_new","ws",ws,"fw",fw,".csv", sep =""))
  write.csv(total_IN,file = paste("total_IN_new","ws",ws,"fw",fw,".csv", sep =""))
  write.csv(total.direct,file = paste("total_d_c_new","ws",ws,"fw",fw,".csv", sep =""))

  total.direct = t(total.direct)
  new_direct = rbind(cbind(total.direct,rowSums(total.direct)), c(colSums(total.direct),sum(colSums(total.direct))))
  colnames(new_direct) = c(colnames(wVAR), "Total IN")
  rownames(new_direct) = c(colnames(wVAR), "Total OUT")
  library(xtable)
  xtable(new_direct, digits = 0)
  dt_fc = dt_f[(ws+1):287]

  pdf(file = paste("total_new",tau,"ws",ws,"fw",fw,".pdf", sep =""),onefile = F, width = 16, height = 8) 
  plot(dt_fc,total.connect, type ="l", lwd = 5,cex.axis=2.5, ylab="",xlab="", ylim = range(0.35,0.9))
  dev.off()

############# total direct connectedness plot #######################################
total_direct = t(total.direct)
dat1 = data.frame(
         EMIT=factor(rep(colnames(wVAR), 10), levels = colnames(wVAR)), 
         RECEIVE=factor(rep(colnames(wVAR), c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10)), levels = colnames(wVAR)),
         value=c(as.numeric(round(total_direct[,1], digits = 0)), as.numeric(round(total_direct[,2], digits = 0)),as.numeric(round(total_direct[,3], digits = 0)),as.numeric(round(total_direct[,4], digits = 0)),as.numeric(round(total_direct[,5], digits = 0)),as.numeric(round(total_direct[,6], digits = 0)),as.numeric(round(total_direct[,7], digits = 0)),as.numeric(round(total_direct[,8], digits = 0)),as.numeric(round(total_direct[,9], digits = 0)),as.numeric(round(total_direct[,10], digits = 0)))
)

p1 = ggplot(dat1, aes(x=EMIT, y=RECEIVE, fill=value)) +
     theme_bw() +
     geom_tile() + 
     theme(axis.text=element_text(size=15,face="bold"),
     axis.title=element_text(size=16,face="bold")) +
     geom_text(aes(label=paste(value)),size=6) +
     scale_fill_gradient2(midpoint=0, low="#B2182B", high="#2166AC")
    

pdf(file = paste("total_direct",tau,"ws",ws,"fw",fw,".pdf", sep =""),onefile = F, width = 10, height = 8) 
p1
dev.off()



####### IN link overtime dynamic graph ################
#IN_link = read.csv("IN_link_005_newws48fw12.csv")[,-1]
IN_link = as.data.frame(IN_link)
  col = c("red", "orange", "yellow", "greenyellow", "limegreen", "cyan", "dodgerblue", "mediumblue",   "violet","purple" ) #""mediumorchid", "purple","deeppink"
rownames(IN_link) = dt_fc
  dygraph(IN_link) %>%
  #dyAxis("x", label = dt_fc)
  dySeries("V1", label = colnames(VARmodel)[1], color = col[1], strokeWidth = 2) %>%
  dySeries("V2", label = colnames(VARmodel)[2], color = col[2], strokeWidth = 2) %>%
  dySeries("V3", label = colnames(VARmodel)[3], color = col[3], strokeWidth = 2) %>%
  dySeries("V4", label = colnames(VARmodel)[4], color = col[4], strokeWidth = 2) %>%
  dySeries("V5", label = colnames(VARmodel)[5], color = col[5], strokeWidth = 2) %>% 
  dySeries("V6", label = colnames(VARmodel)[6], color = col[6], strokeWidth = 2) %>%  
  dySeries("V7", label = colnames(VARmodel)[7], color = col[7], strokeWidth = 2) %>%
  dySeries("V8", label = colnames(VARmodel)[8], color = col[8], strokeWidth = 2) %>%
  dySeries("V9", label = colnames(VARmodel)[9], color = col[9], strokeWidth = 2) %>%
  dySeries("V10", label = colnames(VARmodel)[10], color = col[10], strokeWidth = 2) %>%
  dyOptions(stackedGraph = TRUE) %>%
  dyAxis("x",axisLineWidth = 2,labelWidth =2, pixelsPerLabel = 50) %>%
  dyAxis("y",axisLineWidth = 2,labelWidth =2) %>%
  dyRangeSelector(height = 10)
#, colors = co


####### OUT link overtime dynamic graph ################
#OUT_link = read.csv("OUT_link_05_newws48fw12.csv")[,-1]
OUT_link = as.data.frame(OUT_link)
  col = c("red", "orange", "yellow", "greenyellow", "limegreen", "cyan", "dodgerblue", "mediumblue", "purple",  "violet" ) #""mediumorchid", "purple","deeppink"
rownames(OUT_link) = dt_fc
  dygraph(OUT_link) %>%
  #dyAxis("x", label = dt_fc)
  dySeries("V1", label = colnames(VARmodel)[1], color = col[1], strokeWidth = 2) %>%
  dySeries("V2", label = colnames(VARmodel)[2], color = col[2], strokeWidth = 2) %>%
  dySeries("V3", label = colnames(VARmodel)[3], color = col[3], strokeWidth = 2) %>%
  dySeries("V4", label = colnames(VARmodel)[4], color = col[4], strokeWidth = 2) %>%
  dySeries("V5", label = colnames(VARmodel)[5], color = col[5], strokeWidth = 2) %>% 
  dySeries("V6", label = colnames(VARmodel)[6], color = col[6], strokeWidth = 2) %>%  
  dySeries("V7", label = colnames(VARmodel)[7], color = col[7], strokeWidth = 2) %>%
  dySeries("V8", label = colnames(VARmodel)[8], color = col[8], strokeWidth = 2) %>%
  dySeries("V9", label = colnames(VARmodel)[9], color = col[9], strokeWidth = 2) %>%
  dySeries("V10", label = colnames(VARmodel)[10], color = col[10], strokeWidth = 2) %>%
  dyOptions(stackedGraph = TRUE) %>%
  dyAxis("x",axisLineWidth = 2,labelWidth =2, pixelsPerLabel = 50) %>%
  dyAxis("y",axisLineWidth = 2,labelWidth =2) %>%
  dyRangeSelector(height = 10)

print(new_direct)

library(moments)
## tau = 0.5
mean(F1)      ## 0.0551875
sd(F1)        ## 0.4913917
skewness(F1)  ## -1.261906
kurtosis(F1)  ## 7.854962

mean(F2)      ##  0.01215134
sd(F2)        ## 0.2292202
skewness(F2)  ## 1.375039
kurtosis(F2)  ## 25.53147

## tau = 0.05
mean(F1)      ## 0.04670503
sd(F1)        ## 0.4759937
skewness(F1)  ## -1.080036
kurtosis(F1)  ## 8.084804

mean(F2)      ##  0.02496107
sd(F2)        ## 0.2224448
skewness(F2)  ## 3.921703
kurtosis(F2)  ## 45.95963

## mean of returns 

mean(ts_mean_return)   ## 0.0551875
sd(ts_mean_return)     ## 0.4913917
skewness(ts_mean_return)  ##-1.261906
kurtosis(ts_mean_return)  ## 7.854962

