
% Create the kappaList for the constraint
nk = 10;
kappaList = linspace (10:50, nk);
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
    optTauHatFile =  strcat('test7_optTauHat_', num2str(kappa), '.txt');
        % estimated optimal regime indexing parameters
    optObjValFile =  strcat('test7_optObjVal_', num2str(kappa), '.txt');
        % objective function value under the solution above
    optConValFile =  strcat('test7_optContVal_', num2str(kappa), '.txt');
        % constraint function value under the solution above
    exitFlagFile =  strcat('test7_exitFlagHat_', num2str(kappa), '.txt');
        % exit flag of solving the problem above
    optMeanYTestFile =  strcat('test7_optMeanYTest_', num2str(kappa), '.txt');
        % estimated meanY under the estimated constrained opt regime
    optMeanZTestFile =  strcat('test7_optMeanZTest_', num2str(kappa), '.txt');
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
    tauSol1 = tauSol(1:2)';
    tauSol2 = tauSol(3:4)';
    meanYZ = testset(testseed, tauSol1, tauSol2);
    meanYtestdList(k) = meanYZ(1);
    meanZtestdList(k) = meanYZ(2);
end
%toc;
% plot 


%% calculate the meanYMax and meanZYmax on the test dataset
meanYZtestY = testset(testseed, tauSolY(1:2)', tauSolY(3:4)');
meanYtestY = meanYZtestY(1);
meanZtestY = meanYZtestY(2);

figure% new figure window
stem(kappaList, meanYtestdList, 'color', 'r');
hold on 
stem(kappaList, meanZtestdList, 'color', 'b');
hlineY = refline([0 meanYtestY]);
set(hlineY,'Color','r', 'LineStyle', '-.');
hlineZ = refline([0 meanZtestY]);
set(hlineZ,'Color','b', 'LineStyle', '-.');
hold off % reset hold state

saveas(gcf,'test7.png');
%profile off;
%profsave(profile('info'),'test7');