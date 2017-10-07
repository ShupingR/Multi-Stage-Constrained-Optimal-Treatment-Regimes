setwd('/Users/shuping.ruan/GitHub/research-github/Multi-Stage-Constrained-Optimal-Treatment-Regimes/test/')
source('myPredict.R')

library('FlexCoDE') # conditional density estimation
library('stats') # univariate normal
library('MASS') # multivariate noraml
library('iqLearn')# iqLearn
library('nloptr')


########################
# treatment assignment #
########################


#values <- function(tau, filename){

  #tau  <-read.table("tau.txt",header = FALSE, sep = ",", dec = ".")
  tau <- c(10,2,4,10,5,6)
  tau1 <- tau[1:2]
  tau2 <- tau[3:6]
  load("LTrain.Rdata")
  load("fit.Rdata")
  LTrain <- as.data.frame(LTrain)
  # browser()
  intTrain = rep(1,length(LTrain$X1))
  T1.pred = 
    1*( as.matrix(data.frame(intTrain, LTrain$X1)) %*% t(t(tau1)) > 0) 
  T1.pred[T1.pred == 0] <- -1
  L1.pred <- as.matrix(data.frame(T1.pred, LTrain$X1))
  #plot(fitX2,L1.pred, NaN)
  X2.pred <- predict.FlexCoDE(fitX2, L1.pred, B=length(LTrain$X1))$z
  # H2.pred <- data.frame(intTrain, X2.pred, T1.pred, LTrain$X1)
  # T2.pred <- 1*(as.matrix(H2.pred) %*% t(t(tau2)) > 0)
  # T2.pred[T2.pred == 0] <- -1
  # L2.pred <- data.frame(T2.pred, X2.pred, T1.pred, LTrain$X1)
  # Y.pred <- predict.FlexCoDE(fitY, as.matrix(L2.pred))$z
  # Z.pred <- predict.FlexCoDE(fitZ, as.matrix(L2.pred))$z
  # aveY <- mean(Y.pred)
  # aveZ <- mean(Z.pred)
#  write(aveY, '~/GitHub/research-github/Multi-Stage-Constrained-Optimal-Treatment-Regimes/test/aveY_1.txt', append=FALSE)
#  write(aveZ, '~/GitHub/research-github/Multi-Stage-Constrained-Optimal-Treatment-Regimes/test/aveZ_1.txt', append=FALSE)
#  print(filename)
#  save(T1.pred, X2.pred, T2.pred, Y.pred, Z.pred, file = filename)
#  return(aveY)
#}
#values(c(-10,-2,-4,-10,-5,-6), 'out1.Rdata')
#values(c(10,2,4,10,5,6), 'out2.Rdata')
