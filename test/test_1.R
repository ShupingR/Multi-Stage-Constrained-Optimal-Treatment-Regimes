library('FlexCoDE') # conditional density estimation
library('stats') # univariate normal
library('MASS') # multivariate noraml
library('iqLearn')# iqLearn
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
A1 <- sample(0:1, n, replace=TRUE)
A1[A1 == 0] <- -1
muEpX2 = 0
SigmaEpX2 = 1
EpX2 <- rnorm(n, muEpX2, SigmaEpX2)
H1 <- data.frame(int, X1)
X2Beta0 <- c(0.5, 0.75)
X2Beta1 <- c(0.25, 0.5)
dim(X2Beta0) <- c(2,1)
dim(X2Beta1) <- c(2,1)
X2 = as.matrix(H1) %*% X2Beta0 + A1 * (as.matrix(H1) %*% X2Beta1) + EpX2 
H2 <- data.frame(int, X1, A1, X2)
A2 <- sample(0:1, n, replace=TRUE)
A2[A2 == 0] <- -1
muEpYZ = c(0, 0);
SigmaEpYZ = matrix(c(1, 0.7, 0.7, 1), 2, 2)
EpYZ =  mvrnorm(n, muEpYZ, SigmaEpYZ)
EpY = EpYZ[,1]
EpZ = EpYZ[,2]
Ybeta0 = c(5, 3, 1, 10);
Ybeta1 = c(5, 6, 2, 15);
Zbeta0 = c(5, 1.5, 1, 5);
Zbeta1 = c(5, 3, 2, 7.5);
dim(Ybeta0) <- c(2,1,2,1)
dim(Ybeta1) <- c(2,1,2,1)
dim(Zbeta0) <- c(2,1,2,1)
dim(Zbeta1) <- c(2,1,2,1)
Y = as.matrix(H2) %*% Ybeta0 + A2 * (as.matrix(H2) %*% Ybeta1) + EpY
Z = as.matrix(H2) %*% Zbeta0 + A2 * (as.matrix(H2) %*% Zbeta1) + EpZ
yzData <- data.frame(Y,Z)
L <- data.frame(Y, Z, A2, X2, A1, X1)

nTrain=round(0.7*n)
nValidat=round(0.25*n)
nTest=n-nTrain-nValidat

# split data
randomIndex=sample(1:n)
LTrain = L[randomIndex[1:nTrain],]
LValidat = L[randomIndex[(nTrain+1):(nTrain+nValidat)],]
LTest = L[randomIndex[(nTrain+nValidat+1):n],]

L2Train = as.matrix(LTrain[,c("X2","A2", "X1", "A1")])
L2Validat = as.matrix(LValidat[,c("X2","A2", "X1", "A1")])
L2Test = as.matrix(LTest[,c("X2","A2", "X1", "A1")])

# Fit sparse additive FlexCoDE, Y on full history 
fitY = fitFlexCoDE(L2Train, LTrain$Y,
                   L2Validat, LValidat$Y,
                   L2Test, LTest$Y,
                   nIMax = 30, regressionFunction = regressionFunction.SpAM)
fitY$estimatedRisk
print(fitY)
plot(fitY,L2Test, LTest$Y)
 
# Fit sparse additive FlexCoDE, Z on full history
fitZ = fitFlexCoDE(L2Train, LTrain$Z,
                    L2Validat, LValidat$Z,
                    L2Test, LTest$Z,
                    nIMax = 30, regressionFunction = regressionFunction.SpAM)
fitZ$estimatedRisk
print(fitZ)
plot(fitZ, L2Test, LTest$Z)

# Fit X2 on X1, A1
L1Train = as.matrix(LTrain[,c("X1", "A1")])
L1Validat = as.matrix(LValidat[,c("X1", "A1")])
L1Test = as.matrix(LTest[,c("X1", "A1")])

fitX2=fitFlexCoDE(L1Train, LTrain$Y,
                  L1Validat, LValidat$Y,
                  L1Test, LTest$Y,
                  nIMax = 30, regressionFunction = regressionFunction.SpAM)
fitX2$estimatedRisk
print(fitX2)
plot(fitX2, L1Test, LTest$X2)

## save 
save(fitY, fitZ, fitX2, file = "fit.Rdata")
save(L1Train, file="L1Train.Rdata")
