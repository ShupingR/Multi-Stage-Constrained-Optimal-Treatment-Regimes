function [X1test, EpX2test, EpYtest, EpZtest] = ...
         genTestDataset(testseed, ntest, ...
                        muX1, SigmaX1, ...
                        muEpX2,SigmaEpX2, ...
                        muEpYZ, SigmaEpYZ)
    %================================================
    % @Fun: to generate all the random parts in the testdataset
    % @Input: 
    %     testseed: seed for random generator
    %     ntest: sample size of the test dataset                                   
    %                                                          
    %     n : resapme size                                                              
    %================================================
    % training dataset and apply estimated optimal regime/tauSol to testset
    rng(testseed,'twister');
    % ntest = 10000;
    % muX1 = 1;
    % SigmaX1 = 1; % covariance matrix
    X1test = mvnrnd(muX1, SigmaX1, ntest);
    % H1test = [ones(ntest, 1), X1test];

    % X2Beta0 = [0.5; 0.75];
    % X2Beta1 = [0.25; 0.5];
    % muEpX2 = 0;
    % SigmaEpX2 = 1;
    EpX2test =  mvnrnd(muEpX2, SigmaEpX2, ntest);
    % d1 =  ( H1test * tauSol1 > 0 ) -  ( H1test * tauSol1 <= 0 );
    % X2test = H1test * X2Beta0 + d1 .* (H1test * X2Beta1) + EpX2test; 
    % H2test = [ones(ntest, 1), X2test];

    % Ybeta0 = [30 ; 3];
    % Ybeta1 = [5 ; -1.5];
    % Zbeta0 = [15;  1];
    % Zbeta1 = [3;  -0.5];
    % muEpYZ = [ 0, 0 ];
    % SigmaEpYZ = [1 , 0.7 ; 0.7, 1];
    % d2 =  ( H2test * tauSol2 > 0 ) -  ( H2test * tauSol2 <= 0 );
    EpYZtest =  mvnrnd( muEpYZ, SigmaEpYZ, ntest );
    EpYtest = EpYZtest( : , 1);
    EpZtest = EpYZtest( : , 2);
    % Ytest = H2test * Ybeta0 + d2 .* (H2test * Ybeta1) + EpYtest;
    % Ztest = H2test * Zbeta0 + d2 .* (H2test * Zbeta1) + EpZtest;
    
end