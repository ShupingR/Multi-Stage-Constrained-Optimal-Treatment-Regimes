	## @file cqRun.R

	library (Rcpp);
	library (parallel);
	dyn.load ("test.so");

	kappa = seq (12, 20, by=.25);
	lk = length (kappa);

	myFn = function (i){
		kappa = seq (12, 20, by=.25);
		lk = length (kappa);
		n = as.integer (250);
		m = as.integer (10000);
		J = as.integer (500);
		set.seed (10);
		seeds = abs (floor (10000*rnorm (lk)));
		mu = 1; 
		beta10 = c (.5, .75);
		beta11 = c (.25, .5);
		gamma0 = c (0.5, -0.1); 
		gamma1 = c (0.2, -0.1);
		betaY20 = c (30, 2); 
		betaY21 = c (5, -1.5);
		betaZ20 = c (15, 1); 
		betaZ21 = c (3, -.5);  
		sigma = matrix (c (1, .7, .7, 1), 2, 2);
		es = eigen (sigma);
		sqrtSigma = es$vectors %*% sqrt (diag (es$values)) %*% t (es$vectors); 
		N = as.integer (2500);
		alpha = .25;

		op = .Call ("genEx", n, m, mu, beta10, beta11, gamma0, 
			gamma1, 
			betaY20, betaY21, betaZ20, betaZ21, sqrtSigma, J, N, 
			kappa[i], alpha, seeds[i])

		return (list ("etaHat1"=op$etaHat1, "etaHat2"=op$etaHat2, "feasible"=op$feasible, "meanYopt"=op$meanYopt, "meanZopt"=op$meanZopt, "meanY"=op$meanY, "meanZ"=op$meanZ, "uncstY"=op$uncstY, "uncstZ"=op$uncstZ, "mxYmnZ"=op$mxYmnZ))
	}

results = mclapply (1:lk, myFn, mc.cores=lk);

feasVec = rep (0, lk);
EYopt = rep (0, lk);
EZopt = rep (0, lk);
EY = rep (0, lk);
EZ = rep (0, lk);
uncEY = rep (0, lk);
uncEZ = rep (0, lk);
maxYmeanZ = rep (0, lk);
for (k in 1:lk){
	feasVec[k] = results[[k]]$feasible;
	EYopt[k] = results[[k]]$meanYopt;
	EZopt[k] = results[[k]]$meanZopt;
	EY[k] = results[[k]]$meanY;
	EZ[k] = results[[k]]$meanZ;
	uncEY[k] = results[[k]]$uncstY;
	uncEZ[k] = results[[k]]$uncstZ;
	maxYmeanZ[k] = results[[k]]$mxYmnZ;
}

rMat = cbind (feasVec, EYopt, EZopt, EY, EZ, uncEY, uncEZ, 
	maxYmeanZ);
rMat = as.data.frame (rMat);

write.table (rMat, file="simResults.txt");



## plot results
if (2 < 1){

	setwd ("/Users/kalinnMB/Desktop")
	sim = read.table ("simResults.txt", header=T)
	kappa = seq (12, 20, by=.25)
	nonzero = which (sim[,2]!=0)
	par (mar=c(4.5, 4.5, 1, 2))
	plot (kappa[nonzero], sim[nonzero,2], ylab=expression (paste ("E(Y) under ", hat(pi)[opt])), ylim=c(0, 37), xlim=c(12,20), xlab="kappa")
	points (kappa[-nonzero], sim[-nonzero,2], col="red")
	lines (kappa, sim[,2])
	legend ("bottomright", c("Feasible regime exists", "No feasible regime exists"), col=c("black", "red"), pch=c(1,1), lty=c(-1, -1))
}

# system.time (.Call ("genEx", n, mu, beta10, beta11, gamma0, 
	# gamma1, betaY20, 
# betaY21, betaZ20, betaZ21, sqrtSigma, J, N, kappa))



#myn = new ("f2Model", testEst, f1test, f2test, prSeq, J, N, kappa);
    
#runCQ =
#    system.time (cqEval (mu=mu, beta10=beta10, beta11=beta11, gamma0=gamma0,
#    gamma1=gamma1, betaY20=betaY20, betaY21=betaY21, betaZ20=betaZ20,
#    betaZ21=betaZ21, xiDistn=xiDistn, epsDistn=epsDistn,
#    testSize=testSize, trSize=trSize, J=J, N=N, prSeq=prSeq,
#    kappa=kappa, methodY=methodY, methodZ=methodZ, methodYZ=methodYZ,
#    estGyz=estGyz, M=M, ncores=ncores, sparseSize=sparseSize));

