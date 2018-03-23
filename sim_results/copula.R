library (mvtnorm);
n = 500;
sig = diag(2);
sig[1,2] = 0.9;
sig[2,1] = sig[1,2];
xy = rmvnorm(n, sigma=sig);
x = xy[,1]
y = xy[,2]
estCopula = function(x, y, M=1000){
  xTrans = qnorm(rank(x)/(length(x)+1)); 
  yTrans = qnorm(rank(y)/(length(y)+1)); 
  rho = cor(xTrans, yTrans);
  
  copulaCor = diag(2);
  copulaCor[1,2] = rho;
  copulaCor[2,1] = rho;
  
  xySample = rmvnorm(M, sigma=copulaCor);
  retSample = xySample;
  for (m in 1:M){
    #browser();
    retSample[m,1] = x[which.min(abs(xySample[m,1] - xTrans))];
    retSample[m,2] = y[which.min(abs(xySample[m,2] - yTrans))];
  }   
  return(retSample);
}

retSample = estCopula(x, y, 1000)
