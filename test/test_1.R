library('R.matlab') # call matlab from R
library('FlexCoDE') # conditional density estimation
library('stats') # univariate normal
library('MASS') # multivariate noraml
library('iqLearn')# iqLearn
library('ggplot2')
library('nloptr')
set.seed(2017)

# generate dataset
# estimate optimal regime/index parm on ing dataset and
# apply to testset

n <- 1000
int <- rep(1, n) # intercept
muX1 <- 1 # mean
sigmaX1 <- 1 # covariance 
X1 <- rnorm(n, mean = muX1, sd = sigmaX1)
H1 <- data.frame(int, X1)
A1 <- sample(0:1, n, replace=TRUE)
A1[A1 == 0] <- -1
X2Beta0 <- c(0.5, 0.75)
X2Beta1 <- c(0.25, 0.5)
dim(X2Beta0) <- c(2,1)
dim(X2Beta1) <- c(2,1)
muEpX2 = 0
SigmaEpX2 = 1
EpX2 <- rnorm(n, muEpX2, SigmaEpX2)
X2 = as.matrix(H1) %*% X2Beta0 + A1 * (as.matrix(H1) %*% X2Beta1) + EpX2 
int2 <- int
H2 = data.frame(int2, X2)
Ybeta0 = c(30, 3);
Ybeta1 = c(5, -1.5);
Zbeta0 = c(15, 1);
Zbeta1 = c(3, -0.5);
dim(Ybeta0) <- c(2,1)
dim(Ybeta1) <- c(2,1)
dim(Zbeta0) <- c(2,1)
dim(Zbeta1) <- c(2,1)
muEpYZ = c(0, 0);
SigmaEpYZ = matrix(c(1 , 0.7, 0.7, 1), 2, 2)
A2 <- sample(0:1, n, replace=TRUE)
A2[A2 == 0] <- -1
EpYZ =  mvrnorm(n, muEpYZ, SigmaEpYZ)
EpY = EpYZ[,1]
EpZ = EpYZ[,2]
Y = as.matrix(H2) %*% Ybeta0 + A2 * (as.matrix(H2) %*% Ybeta1) + EpY
Z = as.matrix(H2) %*% Zbeta0 + A2 * (as.matrix(H2) %*% Zbeta1) + EpZ
yzData <- data.frame(Y,Z)
firstStageCovars <- data.frame(H1,A1)
secondStageCovars <- data.frame(H2, A2)
fullHistoryCovars <- data.frame(H2, A2, H1, A1)
#ggplot(data, aes(x=Y, y=Z)) + 
#  geom_point() + 
#
#geom_rug(col="steelblue",alpha=0.1, size=1.5)

# generate data
# determine sample sizes
nTrain=round(0.7*n)
nValidation=round(0.25*n)
nTest=n-nTrain-nValidation

# split data
randomIndex=sample(1:n)

yzTrain=yzData[randomIndex[1:nTrain],]
yzValidation=yzData[randomIndex[(nTrain+1):(nTrain+nValidation)],]
yzTest=yzData[randomIndex[(nTrain+nValidation+1):n],]

firstStageTrain=firstStageCovars[randomIndex[1:nTrain],]
firstStageValidation=firstStageCovars[randomIndex[(nTrain+1):(nTrain+nValidation)],]
firstStageTest=firstStageCovars[randomIndex[(nTrain+nValidation+1):n],]

secondStageTrain=secondStageCovars[randomIndex[1:nTrain],]
secondStageValidation=secondStageCovars[randomIndex[(nTrain+1):(nTrain+nValidation)],]
secondStageTest=secondStageCovars[randomIndex[(nTrain+nValidation+1):n],]

fullHistoryTrain=fullHistoryCovars[randomIndex[1:nTrain],]
fullHistoryValidation=fullHistoryCovars[randomIndex[(nTrain+1):(nTrain+nValidation)],]
fullHistoryTest=fullHistoryCovars[randomIndex[(nTrain+nValidation+1):n],]

# Fit sparse additive FlexCoDE, Y on full history 
fit_y_history=fitFlexCoDE(fullHistoryTrain[,c(2:3,5:6)],yzTrain[,1],
                          fullHistoryValidation[,c(2:3,5:6)],yzValidation[,1],
                          fullHistoryTest[,c(2:3,5:6)],yzTest[,1],nIMax = 30,
                          regressionFunction = regressionFunction.SpAM)
fit_y_history$estimatedRisk
print(fit_y_history)
plot(fit_y_history,as.matrix(fullHistoryTest[,c(2:3,5:6)]),yzTest[,1])

# Fit sparse additive FlexCoDE, Z on full history
fit_z_history=fitFlexCoDE(fullHistoryTrain[,c(2:3,5:6)],yzTrain[,2],
                          fullHistoryValidation[,c(2:3,5:6)],yzValidation[,2],
                          fullHistoryTest[,c(2:3,5:6)],yzTest[,2],nIMax = 30,
                          regressionFunction = regressionFunction.SpAM)
fit_z_history$estimatedRisk
print(fit_z_history)
plot(fit_z_history, as.matrix(fullHistoryTest[,c(2:3,5:6)]), yzTest[,2])

# Fit X2 on X1, A1
fit_x2_first=fitFlexCoDE(firstStageTrain[,2:3], secondStageTrain[,2],
                         firstStageValidation[,2:3], secondStageValidation[,2],
                         firstStageTest[,2:3], secondStageTest[,2], nIMax = 30,
                         regressionFunction = regressionFunction.SpAM)
fit_x2_first$estimatedRisk
print(fit_x2_first)
plot(fit_x2_first,firstStageTest[,2:3], secondStageTest[,2])

########################
# treatment assignment #
########################
tau1 = c(1,2)
tau2 = c(2,1)

T1 = 1*( as.matrix(firstStageTrain[,1:2]) %*% tau1 > 0) 
T1[T1 == 0] <- -1
H1T = data.frame(firstStageTrain[,2], T1)
X2.predict = predict.FlexCoDE(fit_x2_first, as.matrix(firstStageCovars[,2:3]))$z
H2.predict= data.frame(int, X2.predict)
T2 = 1*(as.matrix(H2.predict) %*% tau2 > 0)
T2[T2 == 0] <- -1
Y.predict = predict.FlexCoDE(fit_y_history, 
            as.matrix(data.frame(X1[1:700], T1, X2.predict[1:700], T2[1:700])))$z

###################### 
## New test dataset ##
######################

# to sample from estimated CDF
n <- 1000
int <- rep(1, n) # intercept
muX1 <- 1 # mean
sigmaX1 <- 1 # covariance 
X1 <- rnorm(n, mean = muX1, sd = sigmaX1)
H1 <- data.frame(int, X1)
A1 <- sample(0:1, n, replace=TRUE)
A1[A1 == 0] <- -1
X2Beta0 <- c(0.5, 0.75)
X2Beta1 <- c(0.25, 0.5)
dim(X2Beta0) <- c(2,1)
dim(X2Beta1) <- c(2,1)
muEpX2 = 0
SigmaEpX2 = 1
EpX2 <- rnorm(n, muEpX2, SigmaEpX2)
X2 = as.matrix(H1) %*% X2Beta0 + A1 * (as.matrix(H1) %*% X2Beta1) + EpX2 
int2 <- int
H2 = data.frame(int2, X2)
Ybeta0 = c(30, 3);
Ybeta1 = c(5, -1.5);
Zbeta0 = c(15, 1);
Zbeta1 = c(3, -0.5);
dim(Ybeta0) <- c(2,1)
dim(Ybeta1) <- c(2,1)
dim(Zbeta0) <- c(2,1)
dim(Zbeta1) <- c(2,1)
muEpYZ = c(0, 0);
SigmaEpYZ = matrix(c(1 , 0.7, 0.7, 1), 2, 2)
A2 <- sample(0:1, n, replace=TRUE)
A2[A2 == 0] <- -1
EpYZ =  mvrnorm(n, muEpYZ, SigmaEpYZ  )
EpY = EpYZ[,1]
EpZ = EpYZ[,2]
Y = as.matrix(H2) %*% Ybeta0 + A2 * (as.matrix(H2) %*% Ybeta1) + EpY
Z = as.matrix(H2) %*% Zbeta0 + A2 * (as.matrix(H2) %*% Zbeta1) + EpZ
yzData <- data.frame(Y,Z)
firstStageCovars <- data.frame(H1, A1)
secondStageCovars <- data.frame(H2, A2)
fullHistoryCovars <- data.frame(H2, A2, H1, A1) 
# based on X1, A1 to estimate the density of X2
X2.predict = predict.FlexCoDE(fit_x2_first, as.matrix(firstStageCovars[,2:3]))
# based on X2, A2, X1, A1 to estimate the density of Y or Z
Y.predict = predict.FlexCoDE(fit_y_history, as.matrix(fullHistoryTrain[,c(2:3,5:6)]))$z
# based on X2, A2, X1, A1 to estimate the density of X1
Z.predict = predict.FlexCoDE(fit_z_history, as.matrix(fullHistoryTrain[,c(2:3,5:6)]))$z