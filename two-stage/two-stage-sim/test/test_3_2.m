%% multiple start fimincon only
% check simulation results of different tau, if the tauSol is the best
% pY = 0? pY sum to 1?
% nonparameteric conditional

%% output
% tauSolZ is the min Z
% tauSolY is the max Y
% tauSol is the constrained problem solution
% NaN mean infeasible
%profile on;
%%
clc
clear

parfor rep = 1:150
    testseed = 350 + rep;
    disp(rep);
    rng(testseed, 'twister');
    ns = 5; % num of randon starts
    nk = 30; % num of the kappas for constraint
    %% declare vectors for plotting
    meanYtestdList = zeros(1 , nk); % list of mean Y under regime d
    meanZtestdList = zeros(1 , nk); % list of means Z under regime d

    %% Generate Training Dataset
    n = 1000; % training data sample size
    muX1 = 1; 
    SigmaX1 = 1; % covariance matrix of X1
    X1 = mvnrnd( muX1, SigmaX1, n ); 
    H1 = [ones(n, 1), X1];

    A1 = randi( 0:1, [n,1]); % A1 compeletely random assignment 
    A1(A1 == 0)  = -1;

    X2Beta0 = [0.5; 0.75];
    X2Beta1 = [0.25; 0.5];
    muEpX2 = 0;
    SigmaEpX2 = 1;
    EpX2 =  mvnrnd( muEpX2, SigmaEpX2, n );
    X2 = H1 * X2Beta0 + A1 .* (H1 * X2Beta1) + EpX2; 
    H2 = [ones(n, 1), X2];

    A2 = randi( 0:1, [n,1] ); % A2 compeletely random assignment 
    A2( A2 == 0 )  = -1;

    Ybeta0 = [30 ; 3];
    Ybeta1 = [5 ; -1.5];
    Zbeta0 = [15;  1];
    Zbeta1 = [3;  -0.5];
    muEpYZ = [ 0, 0 ];
    SigmaEpYZ = [1 , 0.7 ; 0.7, 1];
    EpYZ =  mvnrnd( muEpYZ, SigmaEpYZ, n );
    EpY = EpYZ( : , 1);
    EpZ = EpYZ( : , 2);
    Y = H2 * Ybeta0 + A2 .* (H2 * Ybeta1) + EpY;
    Z = H2 * Zbeta0 + A2 .* (H2 * Zbeta1) + EpZ;

    %%
    % remember to set seed for simulation to avoid stochastic 
    % fimsearch use nealder-mean simplex method, gradient free

    %%
    % Multiple Start ininital point
    tau0 = [-0.5, 0.5, 0.5, -0.5 ];

    % display('Z')
    % find min Z for later check if it satisfies the constraint, aka <= kappa
    %%tic;
    objZ = @(tau) preEstPrY (tau, Z, H2, A2, H1, A1, n, 1) ;
    optsZ = optimoptions(@fminunc,'Algorithm', 'quasi-newton', 'Display','off' , 'FinDiffRelStep', 1e-2);
    problemZ = createOptimProblem('fminunc', 'objective', objZ, 'x0', tau0,  'options', optsZ);
    msZ = MultiStart('StartPointsToRun', 'all', 'Display','off');
    [tauSolZ, fvalZ, exitflagZ] = run(msZ, problemZ, ns);
    %toc;

    % display('Y')
    % find max Y for later check if its corresponding Z satisfies the constraint, aka <= kappa
    %tic;
    objY = @(tau) preEstPrY (tau, Y, H2, A2, H1, A1, n, -1) ;
    optsY = optimoptions(@fminunc,'Algorithm', 'quasi-newton', 'Display','off' , 'FinDiffRelStep', 1e-2);
    problemY = createOptimProblem('fminunc', 'objective', objY, 'x0', tau0,  'options', optsY);
    msY = MultiStart('StartPointsToRun', 'all', 'Display','off');
    [tauSolY, fvalY, exitflagY] = run(msY, problemY, ns);
    fvalY = -fvalY;
    fvalZmaxY = preEstPrY (tauSolY, Z, H2, A2, H1, A1, n, 1) ;
    %toc;

    %%
    % display('C')
    % Create the kappaList for the constraint
    stdZ = std(Z);
    kappaList = linspace (min(Z) - 0.5 * stdZ , max(Z) + 0.5 * stdZ, nk);
    kappaListFile =  'test7_kappaList.txt';
    dlmwrite(kappaListFile, kappaList, '-append');
    %reach = 0;
    %profile on
    %write out solutions
    %parpool(2);
     %tic;

    for k = 1:nk 
        kappa = kappaList(k);
        %----------------------------------------------------------------------
        % output file names
        optTauHatFile =  strcat('test7_optTauHat_kappa_', num2str(kappa), '.txt');
            % estimated optimal regime indexing parameters
        optObjValFile =  strcat('test7_optObjVal_kappa_', num2str(kappa), '.txt');
            % objective function value under the solution above
        optConValFile =  strcat('test7_optContVal_kappa_', num2str(kappa), '.txt');
            % constraint function value under the solution above
        exitFlagFile =  strcat('test7_exitFlagHat_kappa_', num2str(kappa), '.txt');
            % exit flag of solving the problem above
        optMeanYTestFile =  strcat('test7_optMeanYTest_kappa_', num2str(kappa), '.txt');
            % estimated meanY under the estimated constrained opt regime
        optMeanZTestFile =  strcat('test7_optMeanZTest_kappa_', num2str(kappa), '.txt');
            % estimated meanZ under the estimated constrained opt regime
        %----------------------------------------------------------------------   
        % check if minMeanZ satisfies the constraint, aka <= kappa
        % if minMeanZ > kappa, problem infeasible. We skip to the next kappa
        if (fvalZ > kappa) 
            tauSol = [ NaN, NaN, NaN, NaN ]; 
            fval = NaN;
            con = NaN;
            exitflag = NaN;
            meanYtestdList(k) = NaN;
            meanZtestdList(k) = NaN;
            dlmwrite(optTauHatFile, tauSol, '-append');
            dlmwrite(optObjValFile, fval, '-append');
            dlmwrite(optConValFile, con, '-append');
            dlmwrite(exitFlagFile, exitflag, '-append');
            continue;
        end

        %----------------------------------------------------------------------   
        % check if the correspond meanZ of maxMeanY satisfies the constraint,
        % aka <= kappa, if it satisfies, then this is the solution
        if (fvalZmaxY <= kappa) 
            % reach = reach + 1;
            tauSol = tauSolY;
            con = fvalZmaxY;
            exitflag = exitflagY;
            tauSol1 = tauSol(1:2)';
            tauSol2 = tauSol(3:4)';
            meanYZ = testset(testseed, tauSol1, tauSol2);
            meanYtestdList(k) = meanYZ(1);
            meanZtestdList(k) = meanYZ(2);

            dlmwrite(optTauHatFile, tauSol, '-append');
            dlmwrite(optObjValFile, fvalY, '-append');
            dlmwrite(optConValFile, con, '-append');
            dlmwrite(exitFlagFile, exitflag, '-append');

            % if maxY is achieved for three times, we stop calculating further,
            % not suitable for parfor
            % if ( reach > 3) 
            %    break;
            % end
            %
            continue;
        end

        %----------------------------------------------------------------------
        % if neither the situation above is satisfied, we solve the problem 
        % using constrained optimization nonlinear objective function 
        % (negative mean Y to be minimized)
        obj = @(tau) preEstPrY (tau, Y, H2, A2, H1, A1, n, -1) ;
        %Nonlinear Inequality and Equality Constraints
        constraint = @(tau) preEstPrZ(tau, Z, H2, A2, H1, A1, kappa, n);
        % fmincon options
        % options = optimoptions(@fmincon,'Display','iter-detailed','Algorithm','interior-point' , 'FiniteDifferenceStepSize', 1e-2);  
        opts = optimoptions(@fmincon,'Display','iter-detailed','Algorithm','interior-point' , 'FinDiffRelStep', 1e-2);  
        % opts = optimset('Algorithm', 'interior-point', 'FinDiffRelStep',1e-2); 
        problem = createOptimProblem('fmincon', 'objective', obj, 'x0', tau0,  'nonlcon', constraint, 'options', opts);
        ms = MultiStart('StartPointsToRun', 'all', 'Display','off');
        [tauSol, fval, exitflag] = run(ms, problem, ns);
        con = preEstPrY (tauSol, Z, H2, A2, H1, A1, n, 1) ;
        fval = -fval;
        dlmwrite(optTauHatFile, tauSol, '-append');
        dlmwrite(optObjValFile, fval, '-append');
        dlmwrite(optConValFile, con, '-append');
        dlmwrite(exitFlagFile, exitflag, '-append');

       %----------------------------------------------------------------------
       % apple the estimated optimal regime / tauSol above to test dataset
       % tauSol1 = tauSol(1:2)';
       % tauSol2 = tauSol(3:4)';
       % meanYZ = testset(testseed, tauSol1, tauSol2);
       % meanYtestdList(k) = meanYZ(1);
       % meanZtestdList(k) = meanYZ(2);
    end
    %toc;
    
end
% 
% % plot 
% %% calculate the meanYMax and meanZYmax on the test dataset
% meanYZtestY = testset(testseed, tauSolY(1:2)', tauSolY(3:4)');
% meanYtestY = meanYZtestY(1);
% meanZtestY = meanYZtestY(2);
% 
% figure% new figure window
% stem(kappaList, meanYtestdList, 'color', 'r');
% hold on 
% stem(kappaList, meanZtestdList, 'color', 'b');
% hlineY = refline([0 meanYtestY]);
% set(hlineY,'Color','r', 'LineStyle', '-.');
% hlineZ = refline([0 meanZtestY]);
% set(hlineZ,'Color','b', 'LineStyle', '-.');
% hold off % reset hold state
% 
% saveas(gcf,'test_new_7.png');
% %profile off;
% %profsave(profile('info'),'test7');
% delete(gcp)
